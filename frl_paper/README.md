# Spatial Effects and Uncertainty in Cryptocurrencies: the case of Network ARCH Models

Replication code and data for the paper published in **Finance Research Letters**.

Repository: https://github.com/MUSET-TE/spatial_crypto

## Abstract

Cryptocurrency markets form tightly interconnected networks in which volatility spills
immediately across coins, yet conventional GARCH models treat each asset in isolation.
Building on the Spatial ARCH approach, we model these spillovers with a dynamic network
ARCH framework, defining neighbours through return correlation adjacency matrices. Using a
balanced daily panel of 50 cryptocurrencies over 1517 trading days (2021–2025), we estimate
volatility dependence across network, time and network diffusion effects by
quasi-maximum likelihood after a multivariate ARMA mean filter. Price volatility exhibits
strong contemporaneous network spillover (ρ ≈ 0.43), whereas volume volatility is largely
idiosyncratic (ρ ≈ 0.12). Diagnostics confirm the model properly addresses the
cryptocurrencies features, as it removes remaining ARCH dependence from residuals.

**Keywords:** Volatility networks, Cryptocurrency, Spatial ARCH
**JEL:** C58, C32, C33, C55, E42

## Repository structure

```
preprocessing.py      builds the balanced panel and the log-return matrices
modelling.R           estimates the two-stage model and reproduces the paper tables
requirements.txt      Python dependencies
data/raw/             the three raw daily series, top 100 coins by market capitalisation
data/processed/       output of preprocessing.py, input to modelling.R
```

Two code files, deliberately: preprocessing in Python, modelling in R. Run them in that
order.

## How to run

```bash
pip install -r requirements.txt
python preprocessing.py            # about 10 seconds
Rscript modelling.R                # about 30 minutes, or 10 with QUICK <- TRUE
```

`modelling.R` needs the R package `Rsolnp`:

```r
install.packages("Rsolnp")
```

## What preprocessing.py does

1. Reads the three raw series: price, trading volume and market capitalisation.
2. Drops the six coins excluded by the retention rule (Supplementary Table S1):
   `whitebit`, `bitget.token`, `jasmycoin`, `tokenize.xchange`, `xdce.crowd.sale`, `weth`.
3. Unifies missingness: a coin-day missing in any one series is treated as missing in all
   three, so the three panels stay aligned.
4. Chooses the analysis window maximising the number of fully observed coin-day
   observations, n × T, rather than imputing gaps.
5. Keeps only the coins fully observed inside that window.
6. Takes log returns.

Result: **50 coins × 1516 return days**, 2021-01-07 to 2025-03-02, with no missing values.

One implementation note that matters for reproduction. The raw headers are CoinGecko ids
containing hyphens, such as `bitget-token`. R rewrites those to dots on import, so the
published pipeline works with dotted names and the exclusion list is written that way.
`preprocessing.py` reproduces that renaming before applying the exclusions; without it
three of the six excluded coins are silently retained and the panel comes out with 53
coins instead of 50. The Python output has been checked against the R pipeline
element by element and agrees to 1e-15.

## What modelling.R does

**Stage 1, the mean filter.** For each coin, an ARMA model of its log returns with the
contemporaneous and once-lagged neighbour returns, `W r_t` and `W r_{t-1}`, as external
regressors. This removes own-asset and cross-asset dependence in the conditional mean, so
that what remains is the unpredictable part. Prices select AR(1), volumes AR(2).

**Stage 2, the volatility network.** The residuals are squared, logged and globally
demeaned, and then

```
(I - rho W) y*_t = gamma y*_{t-1} + delta W y*_{t-1} + u_t
```

is estimated by quasi-maximum likelihood, including the Jacobian term `log|I - rho W|`.
Here `rho` is the contemporaneous spillover between neighbouring coins, `gamma` the
coin's own temporal persistence and `delta` the transmission of neighbours' past
volatility.

**Weight matrices.** Three definitions of neighbourhood, all built from the return
correlations and row-standardised: significance-thresholded at p < 0.01, absolute
correlation thresholded at 0.15, and k nearest neighbours with k = 5. The script fits all
three for each series, which is the comparison the reported specification was selected on:
the matrix minimising the negative log-likelihood, subject to the residual temporal
dependence diagnostics.

**Dynamic stability.** The script also reports the spectral radius of the reduced-form
transition matrix `(I - rho W)^-1 (gamma I + delta W)`, which must be strictly below one.
Because `W` is row-standardised it has 1 as an eigenvalue, so this equals
`(gamma + delta) / (1 - rho)`, equivalently `rho + gamma + delta < 1`. Both selected
specifications satisfy it; the script flags any specification that does not.

## Expected results

| Series | Weight matrix | rho | gamma | delta | −logL |
|---|---|---|---|---|---|
| Prices, AR(1) | threshold, selected | 0.428 | 0.175 | 0.337 | 174832 |
| Volumes, AR(2) | k-NN k=5, selected | 0.116 | 0.183 | 0.015 | 174131 |

Spectral radius 0.895 for prices and 0.224 for volumes, so both are dynamically stable.

## What is reproduced, and what is not

`modelling.R` reproduces exactly, to the digits printed in the paper: the parameter
estimates and their standard errors, the log-likelihood, AIC and BIC, the Ljung-Box shares
for residuals and squared residuals, and the spectral radius used for the stability check.

The ARMA orders are taken as given, not re-derived. `modelling.R` fits AR(1) for prices
and AR(2) for volumes, which are the orders selected in the paper. The selection itself
searched a grid of AR 1 to 7 and MA 0 to 4 across all three weight matrices and stopped at
the first specification whose residual diagnostics were satisfied; that search is described
in the paper but is not reproduced here, because it involves 210 full estimations and would
take many hours. Anyone wanting to verify the choice of orders rather than the estimates
given those orders has to re-run the grid.

Two entries of the published Table 2 are not re-derived here, Moran's I and the ARCH-LM
share. The standard centred Moran statistic gives values close to zero, as the paper
reports, but not the same second decimal, and the normalisation used in the original
estimation notebook could not be reconstructed from the paper alone. Rather than print a
number that disagrees with the published table, the script omits both. Both versions
support the same reading, namely that no material spatial dependence remains in the
residuals.

## Notes on the numbers

The optimiser is `Rsolnp::solnp`, which minimises, so the reported `-logL` is minus the
log-likelihood and smaller is better. Standard errors come from the inverse Hessian at the
optimum. Small differences in the last digits across machines and R versions are expected,
since stage 1 fits 50 ARMA models per specification by maximum likelihood.

## Data

`data/raw/` holds three daily series for the top 100 cryptocurrencies by market
capitalisation: closing price, trading volume and market capitalisation. The first column
is a row index, the second is the date, and every remaining column is one coin identified
by its CoinGecko id.

The series run from 2015-03-06 to 2025-03-02, 3648 daily observations, and were retrieved
from the public CoinGecko API on 31 October 2025. The data are publicly available, so the
panel can be re-retrieved from the same source and checked against the copy included here.

Note that `extract_data/` at the repository root is a different extraction, built with the
R package `crypto2` and therefore sourced from CoinMarketCap. It is not the origin of the
files used for this paper.

## Environment

The published results were produced with:

| | version |
|---|---|
| R | 4.2.2 (2022-10-31 ucrt), x86_64-w64-mingw32 |
| Rsolnp | 1.16 |
| Python | 3.12.7 |
| pandas | 2.2.2 |
| numpy | 1.26.4 |

## Citation

Please cite the published article. Add the DOI here once it is assigned.

```
Mare, C., Belbe, S., Catineanu, M., Otto, P. (2026).
Spatial Effects and Uncertainty in Cryptocurrencies: the case of Network ARCH Models.
Finance Research Letters. DOI to be added.
Code and data: https://github.com/MUSET-TE/spatial_crypto
```
