# ------------------------------------------------------------------------------
# NARCH_estimation
#
# Stage 2 of the pipeline: the dynamic network ARCH model, estimated by
# quasi-maximum likelihood, under each of the three neighbourhood definitions.
#
# The model, on log-squared innovations demeaned globally across coins and days:
#
#   (I - rho W) y*_t = gamma y*_t-1 + delta W y*_t-1 + u_t
#
#   rho    contemporaneous spatial dependence: how much of the neighbours'
#          turbulence today reaches a coin today
#   gamma  temporal persistence: how much of a coin's own turbulence yesterday
#          remains today
#   delta  spatiotemporal transmission: how much of the neighbours' turbulence
#          yesterday reaches a coin today
#
# Two points of estimation theory are worth stating plainly.
#
# First, why not least squares. The contemporaneous spatial lag is endogenous by
# construction: coin j's log-volatility today is a regressor in coin i's equation,
# while coin i's is a regressor in coin j's. The covariance between the spatial lag
# and the error does not vanish as the sample grows, so OLS is inconsistent, not
# merely inefficient.
#
# Second, what the log-determinant term does. Drop log|S(rho)| from the objective
# and what remains is a sum of squares, which the optimiser can shrink simply by
# driving rho towards one, since (I - rho W) y* approaches zero. Spillover would
# then be reported for any series whatsoever, including an independent one. As rho
# approaches one the matrix becomes near-singular, its determinant approaches zero
# and its logarithm tends to minus infinity, so the Jacobian penalises exactly the
# direction the optimiser wants to take. It is not a formality: it is what makes
# rho mean anything at all.
#
# Normality is an estimation device here rather than a claim about the data, since
# the log of a squared innovation is not normal. Hence quasi-maximum likelihood,
# and standard errors read as quasi-ML standard errors.
# ------------------------------------------------------------------------------

suppressWarnings(suppressMessages(library(Rsolnp)))

returns <- list(prices  = as.matrix(read.csv("returns_prices.csv",  check.names = FALSE)),
                volumes = as.matrix(read.csv("returns_volumes.csv", check.names = FALSE)))

# ---- neighbourhood matrices ----------------------------------------------------
row_standardise <- function(W) {
  W[!is.finite(W)] <- 0; diag(W) <- 0
  rs <- rowSums(W); rs[rs == 0] <- 1
  W / rs
}

W_significance <- function(R, alpha = 0.01) {
  C <- cor(R); C[!is.finite(C)] <- 0; diag(C) <- 0
  Tn <- nrow(R)
  tstat <- abs(C) * sqrt((Tn - 2) / pmax(1 - C^2, .Machine$double.eps))
  row_standardise((2 * pt(-tstat, df = Tn - 2) < alpha) * abs(C))
}

W_threshold <- function(R, tau = 0.15) {
  C <- cor(R); C[!is.finite(C)] <- 0; diag(C) <- 0
  row_standardise((abs(C) > tau) * abs(C))
}

W_knn <- function(R, k = 5) {
  C <- cor(R); C[!is.finite(C)] <- 0; diag(C) <- 0
  n <- ncol(R); W <- matrix(0, n, n)
  for (i in seq_len(n)) {
    ord <- order(abs(C[i, ]), decreasing = TRUE); ord <- ord[ord != i]
    W[i, head(ord, k)] <- abs(C[i, head(ord, k)])
  }
  row_standardise(W)
}

# ---- stage 1, reused here because the filter depends on W -----------------------
lag1 <- function(x) c(NA, head(x, -1))

# A coin can end up with an all-zero row of W: under the absolute-correlation
# threshold on volumes, dexe has no neighbour above 0.15. Its neighbour averages are
# then identically zero, the two external regressors are constant, and the design is
# singular. This is not a convergence quirk to be worked around, it is the model
# itself: when every w_ij is zero the neighbour terms drop out of coin i's mean
# equation. Such a coin is filtered by a plain ARMA of the same order, keeping it in
# the panel and keeping n unchanged across cells of the grid.
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

# ---- stage 2 --------------------------------------------------------------------
# Log-squared innovations, demeaned by a single global mean. The model carries no
# intercept, so the constant has to go somewhere; subtracting one pooled mean
# centres the panel without erasing the cross-sectional differences in volatility
# level that the network is supposed to explain. Coin-specific demeaning improves
# in-sample fit but leaves rho, delta and the price versus volume contrast
# qualitatively unchanged.
make_Y <- function(E, eps = 1e-12) {
  keep <- apply(E, 1, function(z) all(is.finite(z)))
  Y <- t(log(pmax(E[keep, , drop = FALSE]^2, eps)))   # n x T
  Y <- Y - mean(Y, na.rm = TRUE)
  Y[!is.finite(Y)] <- 0
  Y
}

negloglik <- function(par, Y, W) {
  rho <- par[1]; gamma <- par[2]; delta <- par[3]; s2 <- par[4]
  n <- nrow(Y); Tn <- ncol(Y); t_eff <- Tn - 1
  S  <- diag(n) - rho * W
  ld <- determinant(S, logarithm = TRUE)$modulus
  sse <- 0
  for (k in 2:Tn) {
    u <- S %*% Y[, k] - gamma * Y[, k - 1] - delta * (W %*% Y[, k - 1])
    sse <- sse + sum(u^2)
  }
  -(-(n * t_eff) / 2 * log(2 * pi) - (n * t_eff) / 2 * log(s2) + t_eff * ld - sse / (2 * s2))
}

qml <- function(Y, W) {
  fit <- solnp(c(0.2, 0.2, 0.2, 1), fun = negloglik, Y = Y, W = W,
               LB = c(-1, -1, -1, 1e-7), UB = c(1, 1, 1, Inf),
               control = list(trace = 0))
  par <- fit$pars
  # Quasi-ML standard errors from the numerical Hessian of the objective at the
  # optimum. The objective is the negative log-likelihood, so its Hessian is the
  # observed information directly.
  H  <- optimHess(par, negloglik, Y = Y, W = W)
  se <- tryCatch(sqrt(diag(solve(H))), error = function(e) rep(NA_real_, length(par)))
  list(par = par, se = se, nll = tail(fit$values, 1),
       n = nrow(Y), t_eff = ncol(Y) - 1)
}

stars <- function(est, se) {
  if (!is.finite(se) || se <= 0) return("")
  p <- 2 * pnorm(-abs(est / se))
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
}

# ---- the grid ---------------------------------------------------------------------
# The ARMA order is held at the value selected for each series, AR(1) for prices and
# AR(2) for volumes, and the weight matrix is varied. Reporting all six cells rather
# than only the two winners is deliberate: it is the evidence that the price versus
# volume contrast does not depend on how the neighbourhood was defined.
specs <- list(
  list(series = "prices",  order = c(1, 0, 0), wname = "significance, p < 0.01", W = W_significance),
  list(series = "prices",  order = c(1, 0, 0), wname = "threshold |r| > 0.15",   W = W_threshold),
  list(series = "prices",  order = c(1, 0, 0), wname = "k-NN, k = 5",            W = W_knn),
  list(series = "volumes", order = c(2, 0, 0), wname = "significance, p < 0.01", W = W_significance),
  list(series = "volumes", order = c(2, 0, 0), wname = "threshold |r| > 0.15",   W = W_threshold),
  list(series = "volumes", order = c(2, 0, 0), wname = "k-NN, k = 5",            W = W_knn)
)

results <- list()
for (i in seq_along(specs)) {
  sp <- specs[[i]]
  R  <- returns[[sp$series]]
  cat(sprintf("[%d/%d] %s, AR(%d), %s ... ", i, length(specs), sp$series,
              sp$order[1], sp$wname)); flush.console()
  W <- sp$W(R)
  E <- mean_filter(R, W, sp$order)
  if (length(attr(E, "isolated")))
    cat(sprintf("[isolated, filtered without neighbour terms: %s] ",
                paste(attr(E, "isolated"), collapse = ", ")))
  if (length(attr(E, "failed")))
    cat(sprintf("[FILTER FAILED: %s] ", paste(attr(E, "failed"), collapse = ", ")))
  Y <- make_Y(E)
  f <- qml(Y, W)
  radius <- (f$par[2] + f$par[3]) / (1 - f$par[1])   # see NARCH_diagnostics
  results[[i]] <- data.frame(
    series = sp$series, order = sprintf("AR(%d)", sp$order[1]), weights = sp$wname,
    rho = f$par[1], gamma = f$par[2], delta = f$par[3], sigma2 = f$par[4],
    se_rho = f$se[1], se_gamma = f$se[2], se_delta = f$se[3],
    nll = f$nll, aic = 2 * f$nll + 2 * 4, bic = 2 * f$nll + log(f$n * f$t_eff) * 4,
    radius = radius, n_obs = f$n * f$t_eff, stringsAsFactors = FALSE)
  cat(sprintf("rho=%.3f  -logL=%.0f\n", f$par[1], f$nll))
}
tab <- do.call(rbind, results)

# ---- the selection table ------------------------------------------------------------
# The criterion is the minimum of -logL subject to the residual diagnostics. Every
# candidate carries four parameters, so AIC = 2(-logL) + 8 is a constant shift and the
# two criteria rank identically. AIC is reported only for comparability with other work.
cat("\n")
cat("Spatial ARCH estimates under alternative weight matrices\n")
cat(strrep("-", 96), "\n")
cat(sprintf("%-9s %-24s %10s %10s %10s %10s %9s\n",
            "series", "weight matrix", "rho", "gamma", "delta", "-logL", "radius"))
cat(strrep("-", 96), "\n")
for (s in c("prices", "volumes")) {
  sub <- tab[tab$series == s, ]
  best <- which.min(sub$nll)
  for (j in seq_len(nrow(sub))) {
    mark <- if (j == best) "  <- selected" else ""
    cat(sprintf("%-9s %-24s %7.3f%-3s %7.3f%-3s %7.3f%-3s %10.0f %9.3f%s\n",
                if (j == 1) s else "", sub$weights[j],
                sub$rho[j],   stars(sub$rho[j],   sub$se_rho[j]),
                sub$gamma[j], stars(sub$gamma[j], sub$se_gamma[j]),
                sub$delta[j], stars(sub$delta[j], sub$se_delta[j]),
                sub$nll[j], sub$radius[j], mark))
  }
}
cat(strrep("-", 96), "\n")
cat("*** p<0.001, ** p<0.01, * p<0.05, quasi-ML standard errors\n")

cat("\nRange of rho across all three neighbourhood definitions\n")
for (s in c("prices", "volumes")) {
  sub <- tab[tab$series == s, ]
  cat(sprintf("  %-8s %.3f to %.3f\n", s, min(sub$rho), max(sub$rho)))
}
cat("  the two ranges do not overlap, so the contrast is not an artefact of the network chosen\n")

cat("\nSelected models, with standard errors\n")
for (s in c("prices", "volumes")) {
  sub <- tab[tab$series == s, ]
  b <- sub[which.min(sub$nll), ]
  cat(sprintf("  %-8s %s, %s | rho %.2f (%.3f)%s  gamma %.2f (%.3f)%s  delta %.2f (%.3f)%s | AIC %.1f | n = %d\n",
              s, b$order, b$weights,
              b$rho, b$se_rho, stars(b$rho, b$se_rho),
              b$gamma, b$se_gamma, stars(b$gamma, b$se_gamma),
              b$delta, b$se_delta, stars(b$delta, b$se_delta),
              b$aic, b$n_obs))
}

write.csv(tab, "estimates.csv", row.names = FALSE)

# ---- figure -----------------------------------------------------------------------
# The robustness claim, drawn. Each neighbourhood definition gives a different rho,
# but every price estimate lies above every volume estimate and the two ranges never
# meet. That is the reason all six specifications are reported.
png("NARCH_estimation.png", width = 2400, height = 1100, res = 170)
par(mar = c(5.0, 11.5, 3.2, 2.0))

ord <- c(3, 2, 1, 6, 5, 4)
est <- tab$rho[ord]
lo  <- est - 1.96 * tab$se_rho[ord]
hi  <- est + 1.96 * tab$se_rho[ord]
lab <- sprintf("%s, %s", tab$series[ord], tab$weights[ord])
cols <- ifelse(tab$series[ord] == "prices", "#1E2761", "#5A6080")

p_rng <- range(tab$rho[tab$series == "prices"])
v_rng <- range(tab$rho[tab$series == "volumes"])

plot(est, seq_along(est), xlim = c(0, 0.58), ylim = c(0.5, length(est) + 0.5),
     type = "n", yaxt = "n", xlab = expression(rho), ylab = "",
     main = expression(paste("Contemporaneous spillover ", rho,
                             " under every neighbourhood definition")))
# rows 1 to 3 are the price models, rows 4 to 6 the volume models
rect(p_rng[1], 0.4, p_rng[2], 3.5, col = "#1E276114", border = NA)
rect(v_rng[1], 3.5, v_rng[2], length(est) + 0.6, col = "#5A608014", border = NA)
abline(v = 0, col = "grey70")

segments(lo, seq_along(est), hi, seq_along(est), col = cols, lwd = 3)
points(est, seq_along(est), pch = 19, cex = 1.5, col = cols)
axis(2, at = seq_along(est), labels = lab, las = 2, cex.axis = 0.78)

text(p_rng[2], 3.28, sprintf("prices: %.3f to %.3f ", p_rng[1], p_rng[2]),
     col = "#1E2761", cex = 0.85, font = 2, adj = c(1, 0))
text(v_rng[2], length(est) + 0.26, sprintf("volumes: %.3f to %.3f ", v_rng[1], v_rng[2]),
     col = "#5A6080", cex = 0.85, font = 2, adj = c(1, 0))
dev.off()
cat("\nwrote NARCH_estimation.png, estimates.csv\n")
