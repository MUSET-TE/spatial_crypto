# ------------------------------------------------------------------------------
# NARCH_meanfilter
#
# Stage 1 of the two-stage pipeline: the spatiotemporal ARMA mean filter.
#
# Why this stage exists at all is the point worth understanding, because it is the
# methodological contribution of the paper and it is easy to skip when reading.
#
# The neighbourhood matrix W is built from the correlation of returns, in
# NARCH_weights. Suppose we now fitted a spatial volatility model directly to raw
# returns. The spatial parameter would pick up the fact that the coins move
# together, which is precisely what we inserted into W by hand. The estimate would
# be large and would mean nothing. The argument is circular.
#
# Stage 1 breaks the circle. For each coin separately, it removes everything in the
# return that was predictable: the coin's own ARMA dynamics, and the contemporaneous
# and lagged neighbour averages. What survives is the innovation, the part of the
# return explained neither by the coin's own past nor by what its neighbours did
# today or yesterday.
#
#   phi(L) r_it = c_i + theta(L) eps_it + lambda_1i (W r_t)_i + lambda_2i (W r_t-1)_i
#
# Only then does stage 2 ask whether the magnitude of those innovations travels
# between neighbours. Return co-movement has already been taken out, so a spatial
# parameter estimated downstream can no longer be measuring it.
#
# The contemporaneous neighbour term deserves a note. Including a same-period
# regressor in a mean equation is normally hard to defend, because an exchange
# closes and transmission shows up the next day. Cryptocurrency markets trade
# continuously, with no closing auction, no weekend and no time zone, so same-day
# transmission in the mean is a real possibility rather than an artefact.
# ------------------------------------------------------------------------------

returns_prices  <- as.matrix(read.csv("returns_prices.csv",  check.names = FALSE))
returns_volumes <- as.matrix(read.csv("returns_volumes.csv", check.names = FALSE))
coins <- colnames(returns_prices)
cat(sprintf("panel: %d days x %d coins\n", nrow(returns_prices), length(coins)))

# ---- neighbourhood matrices (see NARCH_weights for the full treatment) ---------
row_standardise <- function(W) {
  W[!is.finite(W)] <- 0
  diag(W) <- 0
  rs <- rowSums(W)
  rs[rs == 0] <- 1          # an isolated coin keeps an all-zero row
  W / rs
}

W_threshold <- function(R, tau = 0.15) {
  C <- cor(R); C[!is.finite(C)] <- 0; diag(C) <- 0
  row_standardise((abs(C) > tau) * abs(C))
}

W_knn <- function(R, k = 5) {
  C <- cor(R); C[!is.finite(C)] <- 0; diag(C) <- 0
  n <- ncol(R); W <- matrix(0, n, n)
  for (i in seq_len(n)) {
    ord <- order(abs(C[i, ]), decreasing = TRUE)
    ord <- ord[ord != i]
    W[i, head(ord, k)] <- abs(C[i, head(ord, k)])
  }
  row_standardise(W)
}

# ---- the filter ----------------------------------------------------------------
# Estimated coin by coin by maximum likelihood. The two external regressors are the
# row-standardised neighbour average of returns, contemporaneously and at one lag.
lag1 <- function(x) c(NA, head(x, -1))

mean_filter <- function(R, W, order) {
  WR <- R %*% t(W)
  E  <- matrix(NA_real_, nrow(R), ncol(R), dimnames = dimnames(R))
  isolated <- character(0)
  failed   <- character(0)
  for (i in seq_len(ncol(R))) {
    y  <- R[, i]
    xr <- cbind(neighbour_t = WR[, i], neighbour_t1 = lag1(WR[, i]))
    informative <- apply(xr, 2, function(z) sd(z, na.rm = TRUE) > 0)
    if (!any(informative)) {
      xr <- NULL
      isolated <- c(isolated, colnames(R)[i])
    } else {
      xr <- xr[, informative, drop = FALSE]
    }
    ok  <- if (is.null(xr)) !is.na(y) else complete.cases(y, xr)
    fit <- tryCatch(arima(y[ok], order = order, include.mean = TRUE, method = "ML",
                          xreg = if (is.null(xr)) NULL else xr[ok, , drop = FALSE]),
                    error = function(e) NULL)
    if (is.null(fit)) failed <- c(failed, colnames(R)[i])
    else E[which(ok), i] <- as.numeric(residuals(fit))
  }
  attr(E, "isolated") <- isolated
  attr(E, "failed")   <- failed
  E
}

# The selected orders come from a grid over AR 1 to 7 and MA 0 to 4, searched
# jointly with the three weight matrices and scored by the stage 2 objective
# subject to the residual diagnostics. See NARCH_estimation for that search.
jobs <- list(
  prices  = list(R = returns_prices,  W = W_threshold(returns_prices),  order = c(1, 0, 0)),
  volumes = list(R = returns_volumes, W = W_knn(returns_volumes),       order = c(2, 0, 0))
)

for (series in names(jobs)) {
  j <- jobs[[series]]
  cat(sprintf("\n%s: ARMA(%d,%d) with contemporaneous and lagged neighbour returns\n",
              series, j$order[1], j$order[3]))
  E <- mean_filter(j$R, j$W, j$order)
  if (length(attr(E, "failed")))
    cat(sprintf("  filter did not converge for: %s\n", paste(attr(E, "failed"), collapse = ", ")))

  keep <- apply(E, 1, function(z) all(is.finite(z)))
  cat(sprintf("  innovations: %d usable days of %d, one lost to the lagged regressor\n",
              sum(keep), nrow(E)))
  cat(sprintf("  mean innovation %.2e, sd %.4f\n", mean(E[keep, ]), sd(E[keep, ])))

  # How much of the return the filter actually removed, coin by coin.
  share <- 1 - apply(E[keep, ], 2, var) / apply(j$R[keep, ], 2, var)
  cat(sprintf("  variance removed by the filter: median %.1f%%, range %.1f%% to %.1f%%\n",
              100 * median(share), 100 * min(share), 100 * max(share)))

  write.csv(E[keep, ], sprintf("innovations_%s.csv", series), row.names = FALSE)
}

# ---- figure ---------------------------------------------------------------------
# What the filter achieves here, and what it does not.
#
# It is tempting to show the autocorrelation of returns before and after filtering,
# but on this data that figure would mislead. Cryptocurrency returns are already
# close to serially uncorrelated, so there is very little mean dependence to strip
# out and the two autocorrelation functions look much alike. What the filter
# actually removes is the contemporaneous and lagged neighbour component of the
# level, and that shows up in variance rather than in autocorrelation.
#
# Left: the share of each coin's return variance removed by the filter, sorted.
# Right: the autocorrelation of the squared innovations, left almost untouched.
# That is the whole point. Dependence in the mean has gone, dependence in the
# magnitude has not, and the surviving dependence is what stage 2 models.
E_p <- as.matrix(read.csv("innovations_prices.csv", check.names = FALSE))
R_p <- returns_prices[(nrow(returns_prices) - nrow(E_p) + 1):nrow(returns_prices), ]
share_p <- sort(1 - apply(E_p, 2, var) / apply(R_p, 2, var))

png("NARCH_meanfilter.png", width = 2400, height = 1050, res = 170)
par(mfrow = c(1, 2), mar = c(6.8, 4.4, 3.0, 1.0), cex.main = 1.0)

bp <- barplot(100 * share_p, col = "#2E7FD4", border = NA, las = 2,
              cex.names = 0.42, names.arg = names(share_p),
              ylab = "share of return variance removed (%)", ylim = c(0, 100),
              main = "(a) What the mean filter removes, coin by coin")
abline(h = 100 * median(share_p), lty = 2, lwd = 1.6, col = "#E8872B")
text(bp[length(bp)], 100 * median(share_p) + 4,
     sprintf("median %.0f%%", 100 * median(share_p)), adj = c(1, 0),
     col = "#E8872B", cex = 0.85, font = 2)

par(mar = c(4.2, 4.4, 3.0, 1.0))
acf_sq <- acf(E_p[, "bitcoin"]^2, lag.max = 20, plot = FALSE)
plot(acf_sq$lag, acf_sq$acf, type = "h", lwd = 6, col = "#1E2761",
     xlab = "lag", ylab = "autocorrelation",
     ylim = c(-0.05, max(0.25, acf_sq$acf[-1])),
     main = "(b) Squared innovations: dependence survives")
abline(h = c(-1.96, 1.96) / sqrt(nrow(E_p)), lty = 2, col = "#E8872B")
abline(h = 0, col = "grey50")

dev.off()
cat("\nwrote NARCH_meanfilter.png, innovations_prices.csv, innovations_volumes.csv\n")
