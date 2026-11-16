# ------------------------------------------------------------------------------
# NARCH_data
#
# Builds the balanced panel behind the network ARCH study of cryptocurrency
# volatility, and the two descriptive figures that motivate it.
#
# The raw files hold the top 100 cryptocurrencies by market capitalisation, with
# daily price, trading volume and market capitalisation. They are not a rectangle:
# coins are listed at different dates, some are delisted, some have gaps. Two bad
# options present themselves and both are refused here.
#
#   imputing the gaps   inserts smooth, predictable movement. The object under
#                       study is the co-movement of volatility, so imputation
#                       would manufacture the very pattern being tested.
#   shortening to zero  a window short enough to be complete for all 100 coins
#                       leaves too little data to estimate anything.
#
# Instead a coin is retained only if it is fully observed on all three series
# over a common window, and the window itself is chosen to maximise the number of
# fully observed coin-day observations, that is coins x days.
# ------------------------------------------------------------------------------

# ---- raw input ---------------------------------------------------------------
market_cap <- read.csv("market_cap_top_100.csv")
prices     <- read.csv("prices_top_100.csv")
volumes    <- read.csv("volumes_top_100.csv")

# Six coins are excluded under a single rule: too short a history to span any
# usable window, a discontinuous listing, or values that are non-finite once the
# log transform is applied.
#
#   whitebit           907 valid days
#   bitget.token      1293 valid days
#   jasmycoin         1377 valid days
#   tokenize.xchange  discontinuous listing
#   xdce.crowd.sale   discontinuous listing, a rebranded XinFin asset
#   weth              long runs of zeros, non-finite under the log transform
#
# The names carry dots, not hyphens, and that is not cosmetic. The raw coin ids
# use hyphens ("usd-coin"), and read.csv rewrites every hyphen as a dot when it
# makes syntactic column names. An exclusion list written with hyphens therefore
# matches nothing, the six coins stay in, and the panel silently comes out with
# 53 coins instead of 50. No error is raised: the code runs, the tables look
# plausible, and they are wrong. Count your columns after every filtering step.
drop_coins <- c("jasmycoin", "bitget.token", "xdce.crowd.sale",
                "tokenize.xchange", "whitebit", "weth")

market_cap <- market_cap[, !(names(market_cap) %in% drop_coins)]
prices     <- prices[,     !(names(prices)     %in% drop_coins)]
volumes    <- volumes[,    !(names(volumes)    %in% drop_coins)]

non_coin <- c("X", "Date")
coins <- Reduce(intersect, list(setdiff(names(market_cap), non_coin),
                                setdiff(names(prices),     non_coin),
                                setdiff(names(volumes),    non_coin)))
cat(sprintf("coins after the exclusion rule: %d\n", length(coins)))

# ---- a cell counts as observed only if all three series are observed ----------
missing_any <- is.na(market_cap[, coins]) | is.na(prices[, coins]) | is.na(volumes[, coins])
market_cap[, coins][missing_any] <- NA

M <- apply(as.matrix(market_cap[, coins]), 2, as.numeric)

# ---- window selection --------------------------------------------------------
# Order coins by how much history they have, then ask, for each k, how many
# complete days the best k coins share. The objective is k x days, so the choice
# trades breadth against length explicitly instead of by eye.
valid_days <- colSums(!is.na(M))
M_sorted   <- M[, order(valid_days, decreasing = TRUE)]
objective  <- sapply(seq_len(ncol(M)),
                     function(k) k * sum(rowSums(is.na(M_sorted[, 1:k, drop = FALSE])) == 0))

k_best <- which.max(objective)
n_days <- sum(rowSums(is.na(M_sorted[, 1:k_best, drop = FALSE])) == 0)
window <- (nrow(M) - n_days + 1):nrow(M)

# Every coin fully observed on that window is retained, which can exceed k_best:
# k_best only locates the window, it does not fix the membership.
retained <- coins[colSums(is.na(M[window, , drop = FALSE])) == 0]

cat(sprintf("window: %d consecutive days, objective maximised at k = %d\n", n_days, k_best))
cat(sprintf("coins retained: %d\n", length(retained)))
if (all(c("Date") %in% names(prices)))
  cat(sprintf("period: %s to %s\n",
              prices$Date[window[1]], prices$Date[window[length(window)]]))

# ---- log-returns --------------------------------------------------------------
# Prices and volumes are not analysed in levels. Levels are not comparable across
# coins, since one trades in tens of thousands of dollars and another in
# thousandths of a cent, and they carry trends. Logarithmic changes are additive
# over time and symmetric in gains and losses: +50% followed by -50% returns to
# the starting point in logs, but loses a quarter in simple percentages.
#
# For volume the same transform is a growth rate of trading activity, not a
# return. It is treated in parallel with prices throughout the study.
log_returns <- function(df) {
  X <- as.matrix(df[window, retained, drop = FALSE])
  storage.mode(X) <- "numeric"
  log(X[-1, , drop = FALSE] / X[-nrow(X), , drop = FALSE])
}

returns_prices  <- log_returns(prices)
returns_volumes <- log_returns(volumes)

# The original pipeline guards against non-finite values and exact zeros left by
# the log transform. On the retained panel the guard never fires, which is worth
# reporting rather than leaving implicit: the balance of the panel is genuine and
# no value anywhere in it has been replaced.
cat(sprintf("guard, non-finite cells: prices %d, volumes %d\n",
            sum(!is.finite(returns_prices)), sum(!is.finite(returns_volumes))))
cat(sprintf("guard, exact zeros:      prices %d, volumes %d\n",
            sum(returns_prices == 0), sum(returns_volumes == 0)))
cat(sprintf("panel: %d days x %d coins, %d observations per series\n",
            nrow(returns_prices), ncol(returns_prices), length(returns_prices)))
cat(sprintf("average pairwise return correlation: prices %.2f, volumes %.2f\n",
            mean(cor(returns_prices)[upper.tri(cor(returns_prices))]),
            mean(cor(returns_volumes)[upper.tri(cor(returns_volumes))])))

write.csv(returns_prices,  "returns_prices.csv",  row.names = FALSE)
write.csv(returns_volumes, "returns_volumes.csv", row.names = FALSE)

# ---- figure --------------------------------------------------------------------
# Left: price levels of the ten largest coins, on a log scale, so series spanning
# five orders of magnitude can share an axis. Right: daily log-returns of bitcoin.
# The point of the right panel is visible without any model, namely that turbulent
# days arrive in clusters rather than uniformly. That single fact is the reason the
# whole GARCH family exists, and the reason this study has anything to measure.
top10 <- retained[1:10]
dates <- as.Date(prices$Date[window])
P10   <- as.matrix(prices[window, top10, drop = FALSE])
storage.mode(P10) <- "numeric"

palette10 <- c("#E4457A", "#2E9E4F", "#2E7FD4", "#E8872B", "#C64BD0",
               "#39B7C4", "#8FA01F", "#1E2761", "#A0522D", "#7A7A7A")

png("NARCH_data.png", width = 2600, height = 1150, res = 170)
par(mfrow = c(1, 2), mar = c(4.0, 4.4, 2.6, 1.0), cex.main = 1.0)

matplot(dates, P10, type = "l", lty = 1, lwd = 1.3, col = palette10, log = "y",
        xaxt = "n", xlab = "", ylab = "price (USD, log scale)",
        main = "(a) Price levels, ten largest coins")
axis.Date(1, dates, format = "%Y")
legend("bottomright", legend = top10, col = palette10, lty = 1, lwd = 1.6,
       cex = 0.62, ncol = 2, bty = "n")

plot(dates[-1], returns_prices[, "bitcoin"], type = "l", lwd = 0.5, col = "#2E7FD4",
     xaxt = "n", xlab = "", ylab = "daily log-return",
     main = "(b) Bitcoin daily log-returns")
axis.Date(1, dates[-1], format = "%Y")
abline(h = 0, col = "grey60", lwd = 0.8)

dev.off()
cat("wrote NARCH_data.png, returns_prices.csv, returns_volumes.csv\n")
