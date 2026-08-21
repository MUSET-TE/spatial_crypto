"""
Preprocessing for "Spatial Effects and Uncertainty in Cryptocurrencies:
the case of Network ARCH Models" (Finance Research Letters).

Builds the balanced daily panel used in the paper and writes the log-return
matrices that modelling.R consumes.

Steps, in the order the paper describes them:
  1. read the three raw series (price, volume, market capitalisation)
  2. drop the six coins excluded by the retention rule (Supplementary Table S1)
  3. unify missingness: a coin-day missing in any one series is missing in all three
  4. choose the analysis window that maximises the number of fully observed
     coin-day observations, n * T
  5. keep only the coins fully observed inside that window
  6. take log returns and drop the first row

Output: data/processed/returns_prices.csv, returns_volumes.csv, panel_info.json

Usage:  python preprocessing.py
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

RAW = Path("data/raw")
OUT = Path("data/processed")
SEED = 123

# Excluded under the retention rule; see Supplementary Table S1 for the reasons.
DROPPED = [
    "jasmycoin",
    "bitget.token",
    "xdce.crowd.sale",
    "tokenize.xchange",
    "whitebit",
    "weth",
]

NON_COIN = ["X", "Date"]


def r_make_names(cols) -> list:
    """Reproduce R read.csv column-name sanitising.

    The raw headers are CoinGecko ids containing hyphens, for example
    "bitget-token" and "xdce-crowd-sale". R's read.csv rewrites any character
    that is not alphanumeric, underscore or dot into a dot, and names the
    unnamed index column "X". The published pipeline therefore works with dotted
    names, and the exclusion list above is written that way. Normalising here
    keeps Python and R on exactly the same coin identifiers; without it three of
    the six excluded coins are silently retained and the panel comes out with 53
    coins instead of 50.
    """
    out = []
    for c in cols:
        s = str(c)
        if s == "" or s.startswith("Unnamed:"):
            out.append("X")
        else:
            out.append(re.sub(r"[^A-Za-z0-9_.]", ".", s))
    return out


def load(name: str) -> pd.DataFrame:
    df = pd.read_csv(RAW / f"{name}_top_100.csv")
    df.columns = r_make_names(df.columns)
    return df.drop(columns=[c for c in DROPPED if c in df.columns])


def choose_window(mcap: pd.DataFrame) -> tuple[int, int]:
    """Return (n_coins, T_days) maximising n * T.

    Coins are ordered by their number of valid days, descending. For each k the
    number of rows complete across the first k coins is counted, and k * T_k is
    maximised. Ties go to the smaller k, matching R's which.max.
    """
    valid_days = mcap.notna().sum().sort_values(ascending=False, kind="mergesort")
    ordered = mcap[valid_days.index]
    complete = ordered.notna().cumprod(axis=1)      # 1 while every coin so far is present
    t_by_k = complete.sum(axis=0).to_numpy()        # complete rows for each k
    k_grid = np.arange(1, ordered.shape[1] + 1)
    product = k_grid * t_by_k
    best = int(np.argmax(product))                  # first maximum, as in R
    return int(k_grid[best]), int(t_by_k[best])


def log_returns(prices: pd.DataFrame, coins: list[str]) -> pd.DataFrame:
    """log(x_t / x_{t-1}), first row dropped."""
    return np.log(prices[coins] / prices[coins].shift(1)).iloc[1:].reset_index(drop=True)


def impute(mat: pd.DataFrame, rng: np.random.Generator, label: str) -> pd.DataFrame:
    """Defensive replacement of non-finite values and exact zeros by N(0, 0.001).

    Kept for fidelity with the published pipeline. On this data it never fires:
    the counters below are zero for both series.
    """
    arr = mat.to_numpy(dtype=float)
    bad = ~np.isfinite(arr)
    zero = arr == 0
    n_bad, n_zero = int(bad.sum()), int(zero.sum())
    if n_bad:
        arr[bad] = rng.normal(0.0, 0.001, n_bad)
    if n_zero:
        arr[zero] = rng.normal(0.0, 0.001, n_zero)
    print(
        f"  {label}: {arr.size} cells, non-finite {n_bad} "
        f"({100 * n_bad / arr.size:.3f}%), exact zeros {n_zero} "
        f"({100 * n_zero / arr.size:.3f}%)"
    )
    return pd.DataFrame(arr, columns=mat.columns)


def main() -> None:
    rng = np.random.default_rng(SEED)
    OUT.mkdir(parents=True, exist_ok=True)

    prices, volumes, mcap = load("prices"), load("volumes"), load("market_cap")

    # coin columns, in the order they appear in the market-cap file, which is
    # descending market capitalisation; the intersection guards against a coin
    # being present in one file only
    coins = [
        c
        for c in mcap.columns
        if c not in NON_COIN and c in prices.columns and c in volumes.columns
    ]
    print(f"coins after the exclusion rule: {len(coins)}")

    # a coin-day missing anywhere is missing everywhere
    missing = mcap[coins].isna() | prices[coins].isna() | volumes[coins].isna()
    for df in (mcap, prices, volumes):
        df[coins] = df[coins].mask(missing)

    n_sel, t_sel = choose_window(mcap[coins])
    total_rows = len(mcap)
    window = slice(total_rows - t_sel, total_rows)          # the trailing t_sel rows
    print(f"window: n*T maximised at n={n_sel}, T={t_sel} (rows {window.start + 1}-{window.stop})")

    mcap_w = mcap.iloc[window]
    retained = [c for c in coins if mcap_w[c].notna().all()]
    print(f"coins fully observed inside the window: {len(retained)}")

    prices_w = prices.iloc[window].reset_index(drop=True)
    volumes_w = volumes.iloc[window].reset_index(drop=True)

    print("imputation counters (expected to be zero):")
    rp = impute(log_returns(prices_w, retained), rng, "returns_prices")
    rv = impute(log_returns(volumes_w, retained), rng, "returns_volumes")

    rp.to_csv(OUT / "returns_prices.csv", index=False)
    rv.to_csv(OUT / "returns_volumes.csv", index=False)

    dates = None
    if "Date" in prices.columns:
        d = prices.iloc[window]["Date"].reset_index(drop=True)
        dates = [str(d.iloc[1]), str(d.iloc[-1])]     # first return day, last day

    info = {
        "n_coins": len(retained),
        "n_days_levels": int(t_sel),
        "n_days_returns": int(rp.shape[0]),
        "coins": retained,
        "dropped": DROPPED,
        "return_window": dates,
    }
    (OUT / "panel_info.json").write_text(json.dumps(info, indent=2))

    print(
        f"\nwritten: returns_prices.csv and returns_volumes.csv, "
        f"{rp.shape[0]} days x {rp.shape[1]} coins"
    )
    if dates:
        print(f"return window: {dates[0]} to {dates[1]}")


if __name__ == "__main__":
    main()
