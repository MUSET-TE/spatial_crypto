## =============================================================================
## Modelling for "Spatial Effects and Uncertainty in Cryptocurrencies:
## the case of Network ARCH Models" (Finance Research Letters)
##
## Reads the balanced panel produced by preprocessing.py and estimates the
## two-stage dynamic Spatial ARCH model.
##
##   Stage 1  per-coin ARMA mean filter with the contemporaneous and once-lagged
##            neighbour returns, W r_t and W r_{t-1}, as external regressors
##   Stage 2  log-squared residuals, globally demeaned, then
##              (I - rho W) y*_t = gamma y*_{t-1} + delta W y*_{t-1} + u_t
##            estimated by QML, including the Jacobian term log|I - rho W|
##
## Outputs
##   Table 1    estimates under all three weight matrices (main text)
##   Supp.      with RUN_SUPPLEMENTARY, the rolling-window network (Table S2)
##              and the coin-specific demeaning (Table S3)
##   Table 2    the two selected models: estimates, fit statistics and the
##              Ljung-Box diagnostics. Moran's I and the ARCH-LM share printed in
##              the published Table 2 come from the original estimation notebook
##              and are not re-derived here, because the normalisation used there
##              could not be reconstructed from the paper.
##   Stability  spectral radius of the reduced-form transition matrix
##
## Requires: Rsolnp. Run with the working directory set to the repository root.
## Set QUICK to TRUE to fit only the two selected models (about 10 minutes)
## instead of all six candidates (about 30 minutes).
## =============================================================================

QUICK <- FALSE

## Set to TRUE to also reproduce the Supplementary tables: the rolling-window
## network and the coin-specific demeaning. Adds roughly 40 minutes.
RUN_SUPPLEMENTARY <- FALSE

suppressWarnings(suppressMessages({library(Rsolnp)}))
set.seed(123)

Rp <- as.matrix(read.csv("data/processed/returns_prices.csv",  check.names = FALSE))
Rv <- as.matrix(read.csv("data/processed/returns_volumes.csv", check.names = FALSE))
stopifnot(all(dim(Rp) == dim(Rv)), all(is.finite(Rp)), all(is.finite(Rv)))
cat(sprintf("panel: %d days x %d coins\n\n", nrow(Rp), ncol(Rp)))

## ------------------------------------------------------------- weight matrices
rowstd <- function(W) {
  W[!is.finite(W)] <- 0
  diag(W) <- 0
  rs <- rowSums(W)
  rs[rs == 0] <- 1          # a coin with no neighbour keeps an all-zero row
  W / rs
}

W_thr <- function(R, tau = 0.15) {
  C <- cor(R); C[!is.finite(C)] <- 0; diag(C) <- 0
  rowstd((abs(C) > tau) * abs(C))
}

W_knn <- function(R, k = 5) {
  C <- cor(R); C[!is.finite(C)] <- 0; diag(C) <- 0
  n <- ncol(R); W <- matrix(0, n, n)
  for (i in seq_len(n)) {
    o <- order(abs(C[i, ]), decreasing = TRUE)
    o <- o[o != i]
    W[i, head(o, k)] <- abs(C[i, head(o, k)])
  }
  rowstd(W)
}

W_sig <- function(R, p_thr = 0.01) {
  C <- cor(R); C[!is.finite(C)] <- 0; diag(C) <- 0
  P <- outer(colnames(R), colnames(R),
             Vectorize(function(i, j) cor.test(R[, i], R[, j])$p.value))
  rowstd((P < p_thr) * abs(C))
}

## --------------------------------------------------------- stage 1, the filter
lag1 <- function(x) c(NA, head(x, -1))

mean_filter <- function(R, W, order) {
  WR <- R %*% t(W)
  E  <- matrix(NA_real_, nrow(R), ncol(R))
  fallback <- character(0)
  for (i in seq_len(ncol(R))) {
    y  <- R[, i]
    x1 <- WR[, i]
    ## an all-zero W row makes the neighbour regressor constant, so its
    ## coefficient is not identified; fit a plain ARMA for that coin instead
    flat <- !any(is.finite(x1)) || isTRUE(all.equal(sd(x1, na.rm = TRUE), 0))
    if (flat) {
      ok  <- is.finite(y)
      fallback <- c(fallback, colnames(R)[i])
      fit <- tryCatch(arima(y[ok], order = order, include.mean = TRUE, method = "ML"),
                      error = function(e) NULL)
    } else {
      X   <- cbind(Wr_t = x1, Wr_lag1 = lag1(x1))
      ok  <- complete.cases(y, X)
      fit <- tryCatch(arima(y[ok], order = order, xreg = X[ok, ],
                            include.mean = TRUE, method = "ML"),
                      error = function(e) NULL)
    }
    if (!is.null(fit)) E[which(ok), i] <- as.numeric(residuals(fit))
  }
  list(E = E, fallback = fallback)
}

## ------------------------------------------------------------ stage 2, the QML
make_Y <- function(E, eps = 1e-12) {
  keep <- apply(E, 1, function(z) all(is.finite(z)))
  Z <- log(pmax(E[keep, , drop = FALSE]^2, eps))
  Y <- t(Z)
  Y <- Y - mean(Y, na.rm = TRUE)      # global demeaning
  Y[!is.finite(Y)] <- 0
  Y
}

qml <- function(Y, W) {
  n <- nrow(Y); Tn <- ncol(Y); t <- Tn - 1
  nll <- function(p) {
    rho <- p[1]; g <- p[2]; d <- p[3]; s <- p[4]
    S  <- diag(n) - rho * W
    ld <- determinant(S, logarithm = TRUE)$modulus
    sse <- 0
    for (k in 2:Tn) {
      u   <- S %*% Y[, k] - g * Y[, k - 1] - d * (W %*% Y[, k - 1])
      sse <- sse + sum(u^2)
    }
    ## returns MINUS the log-likelihood, so solnp minimises it
    -(-(n * t) / 2 * log(2 * pi) - (n * t) / 2 * log(s) + t * ld - sse / (2 * s))
  }
  o  <- solnp(c(0.2, 0.2, 0.2, 1), fun = nll,
              LB = c(-1, -1, -1, 1e-7), UB = c(1, 1, 1, Inf),
              control = list(trace = 0))
  se <- tryCatch(sqrt(diag(solve(o$hessian))), error = function(e) rep(NA, 4))
  list(theta = o$pars, se = se, nll = o$values[length(o$values)])
}

stars <- function(e, s) {
  t <- abs(e / s)
  if (is.na(t)) "" else if (t > 3.291) "***" else if (t > 2.576) "**" else if (t > 1.960) "*" else ""
}

## Spectral radius of (I - rho W)^-1 (gamma I + delta W). Dynamic stability needs
## it strictly below one. Because W is row-standardised it has 1 as an eigenvalue,
## so this quantity always equals (gamma + delta) / (1 - rho).
spectral_radius <- function(W, rho, g, d) {
  n <- nrow(W)
  A <- solve(diag(n) - rho * W) %*% (g * diag(n) + d * W)
  max(Mod(eigen(A, only.values = TRUE)$values))
}

## ----------------------------------------------------------------- diagnostics
diagnostics <- function(Y, W, theta, lag = 20) {
  n <- nrow(Y); Tn <- ncol(Y)
  S <- diag(n) - theta[1] * W
  U <- matrix(NA_real_, n, Tn - 1)
  for (k in 2:Tn) {
    U[, k - 1] <- S %*% Y[, k] - theta[2] * Y[, k - 1] - theta[3] * (W %*% Y[, k - 1])
  }
  lb <- function(x) tryCatch(Box.test(x, lag = lag, type = "Ljung-Box")$p.value,
                             error = function(e) NA)
  p_lv <- apply(U,   1, lb)
  p_sq <- apply(U^2, 1, lb)
  list(share_lb_lv = mean(p_lv < 0.05, na.rm = TRUE),
       share_lb_sq = mean(p_sq < 0.05, na.rm = TRUE))
}

## ------------------------------------------------------------------ estimation
fit_one <- function(R, order, W, label) {
  mf <- mean_filter(R, W, order)
  Y  <- make_Y(mf$E)
  r  <- qml(Y, W)
  sr <- spectral_radius(W, r$theta[1], r$theta[2], r$theta[3])
  cat(sprintf("  %-34s rho=%.3f%-3s gamma=%.3f%-3s delta=%.3f%-3s  -logL=%7.0f  sr=%.3f%s\n",
              label,
              r$theta[1], stars(r$theta[1], r$se[1]),
              r$theta[2], stars(r$theta[2], r$se[2]),
              r$theta[3], stars(r$theta[3], r$se[3]),
              r$nll, sr, if (sr >= 1) "  << UNSTABLE" else ""))
  list(theta = r$theta, se = r$se, nll = r$nll, sr = sr, Y = Y, W = W,
       fallback = mf$fallback)
}

series <- list(
  list(name = "prices",  R = Rp, order = c(1, 0, 0), sel = "threshold"),
  list(name = "volumes", R = Rv, order = c(2, 0, 0), sel = "k-NN (k=5)")
)

cat("=== Table 1: estimates under alternative weight matrices ===\n")
cat("Each candidate uses the SAME matrix in both stages, which is the comparison\n")
cat("the specification was selected on. sr is the spectral radius.\n\n")

results <- list()
for (s in series) {
  cat(sprintf("%s, AR(%d)\n", s$name, s$order[1]))
  Ws <- if (QUICK) {
    setNames(list(if (s$sel == "threshold") W_thr(s$R) else W_knn(s$R)), s$sel)
  } else {
    list("significance (p<0.01)" = W_sig(s$R),
         "threshold"             = W_thr(s$R),
         "k-NN (k=5)"            = W_knn(s$R))
  }
  for (nm in names(Ws)) {
    tag <- if (nm == s$sel) " (selected)" else ""
    results[[paste(s$name, nm)]] <- fit_one(s$R, s$order, Ws[[nm]], paste0(nm, tag))
  }
  cat("\n")
}

cat("=== Table 2: the two selected models, with diagnostics ===\n")
for (s in series) {
  f <- results[[paste(s$name, s$sel)]]
  d <- diagnostics(f$Y, f$W, f$theta)
  cat(sprintf("\n%s (%s)\n", s$name, s$sel))
  cat(sprintf("  rho      %.4f (%.4f)%s\n", f$theta[1], f$se[1], stars(f$theta[1], f$se[1])))
  cat(sprintf("  gamma    %.4f (%.4f)%s\n", f$theta[2], f$se[2], stars(f$theta[2], f$se[2])))
  cat(sprintf("  delta    %.4f (%.4f)%s\n", f$theta[3], f$se[3], stars(f$theta[3], f$se[3])))
  cat(sprintf("  sigma2_u %.4f\n", f$theta[4]))
  cat(sprintf("  -logL %.1f   AIC %.1f   BIC %.1f\n", f$nll, 2 * f$nll + 8,
              2 * f$nll + 4 * log(nrow(f$Y) * (ncol(f$Y) - 1))))
  cat(sprintf("  Ljung-Box share significant, residuals         %.2f\n", d$share_lb_lv))
  cat(sprintf("  Ljung-Box share significant, squared residuals %.2f\n", d$share_lb_sq))
  cat(sprintf("  spectral radius %.3f -> %s\n", f$sr,
              if (f$sr < 1) "dynamically stable" else "NOT stable"))
  if (length(f$fallback)) {
    cat(sprintf("  plain-ARMA fallback used for: %s\n", paste(f$fallback, collapse = ", ")))
  }
}

cat("\n=== Dynamic stability across all fitted specifications ===\n")
for (k in names(results)) {
  f <- results[[k]]
  cat(sprintf("  %-34s spectral radius %.3f  %s\n", k, f$sr,
              if (f$sr < 1) "stable" else "UNSTABLE"))
}

## =============================================================================
## Supplementary analyses
## =============================================================================
if (RUN_SUPPLEMENTARY) {

  ## --- Table S2: static versus rolling-window (time-varying) network ---------
  ## The weight matrix is rebuilt from trailing return correlations and the model
  ## is re-estimated with time-varying weights. The mean filter is held fixed at
  ## the static specification, so the filtered series is identical across rows
  ## and the comparison is a clean one of in-sample fit.

  build_rolling_W <- function(R, kept, L, rule) {
    static <- if (rule == "thr") W_thr(R) else W_knn(R)
    lapply(seq_along(kept), function(k) {
      tt <- kept[k]
      lo <- tt - L
      if (lo < 1) static else {
        sub <- R[lo:(tt - 1), , drop = FALSE]
        tryCatch(if (rule == "thr") W_thr(sub) else W_knn(sub),
                 error = function(e) static)
      }
    })
  }

  qml_tv <- function(Y, Wlist) {
    n <- nrow(Y); Tn <- ncol(Y); t <- Tn - 1
    nll <- function(p) {
      rho <- p[1]; g <- p[2]; d <- p[3]; s <- p[4]
      sse <- 0; ld <- 0
      for (k in 2:Tn) {
        Wt <- Wlist[[k]]
        S  <- diag(n) - rho * Wt
        ld  <- ld + determinant(S, logarithm = TRUE)$modulus
        u   <- S %*% Y[, k] - g * Y[, k - 1] - d * (Wt %*% Y[, k - 1])
        sse <- sse + sum(u^2)
      }
      -(-(n * t) / 2 * log(2 * pi) - (n * t) / 2 * log(s) + ld - sse / (2 * s))
    }
    o  <- solnp(c(0.2, 0.2, 0.2, 1), fun = nll,
                LB = c(-1, -1, -1, 1e-7), UB = c(1, 1, 1, Inf),
                control = list(trace = 0))
    se <- tryCatch(sqrt(diag(solve(o$hessian))), error = function(e) rep(NA, 4))
    list(theta = o$pars, se = se, nll = o$values[length(o$values)])
  }

  make_Y_keep <- function(E, eps = 1e-12) {
    keep <- apply(E, 1, function(z) all(is.finite(z)))
    Z <- log(pmax(E[keep, , drop = FALSE]^2, eps))
    Y <- t(Z); Y <- Y - mean(Y, na.rm = TRUE); Y[!is.finite(Y)] <- 0
    list(Y = Y, kept = which(keep))
  }

  cat("\n=== Table S2: static versus rolling-window network ===\n")
  cat("For a row-standardised W the spectral radius equals (gamma+delta)/(1-rho),\n")
  cat("which is reported so the stability condition can be checked on every row.\n\n")

  for (s in series) {
    rule <- if (s$sel == "threshold") "thr" else "knn"
    W0   <- if (rule == "thr") W_thr(s$R) else W_knn(s$R)
    mf   <- mean_filter(s$R, W0, s$order)
    YY   <- make_Y_keep(mf$E)
    cat(sprintf("%s\n", s$name))
    r  <- qml(YY$Y, W0)
    sr <- (r$theta[2] + r$theta[3]) / (1 - r$theta[1])
    cat(sprintf("  %-18s rho=%.3f gamma=%.3f delta=%.3f  -logL=%7.0f AIC=%8.0f  sr=%.3f%s\n",
                "static", r$theta[1], r$theta[2], r$theta[3], r$nll, 2 * r$nll + 8, sr,
                if (sr >= 1) "  << UNSTABLE" else ""))
    for (L in c(14, 21, 28)) {
      Wl <- build_rolling_W(s$R, YY$kept, L, rule)
      rr <- qml_tv(YY$Y, Wl)
      sr <- (rr$theta[2] + rr$theta[3]) / (1 - rr$theta[1])
      cat(sprintf("  %-18s rho=%.3f gamma=%.3f delta=%.3f  -logL=%7.0f AIC=%8.0f  sr=%.3f%s\n",
                  sprintf("rolling %d days", L),
                  rr$theta[1], rr$theta[2], rr$theta[3], rr$nll, 2 * rr$nll + 8, sr,
                  if (sr >= 1) "  << UNSTABLE" else ""))
    }
    cat("\n")
  }

  ## --- Table S3: global versus coin-specific demeaning -----------------------
  ## Coin-specific demeaning subtracts each coin's own mean log-squared residual,
  ## that is, it introduces spatial fixed effects.

  make_Y_coin <- function(E, eps = 1e-12) {
    keep <- apply(E, 1, function(z) all(is.finite(z)))
    Z <- log(pmax(E[keep, , drop = FALSE]^2, eps))
    Y <- t(Z)
    Y <- Y - rowMeans(Y, na.rm = TRUE)
    Y[!is.finite(Y)] <- 0
    Y
  }

  cat("=== Table S3: global versus coin-specific demeaning ===\n\n")
  for (s in series) {
    W0 <- if (s$sel == "threshold") W_thr(s$R) else W_knn(s$R)
    mf <- mean_filter(s$R, W0, s$order)
    cat(sprintf("%s\n", s$name))
    for (nm in c("global", "coin-specific")) {
      Y  <- if (nm == "global") make_Y(mf$E) else make_Y_coin(mf$E)
      r  <- qml(Y, W0)
      sr <- spectral_radius(W0, r$theta[1], r$theta[2], r$theta[3])
      cat(sprintf("  %-15s rho=%.3f gamma=%.3f delta=%.3f  -logL=%7.0f  sr=%.3f%s\n",
                  nm, r$theta[1], r$theta[2], r$theta[3], r$nll, sr,
                  if (sr >= 1) "  << UNSTABLE" else ""))
    }
    cat("\n")
  }
}

cat("\nDONE\n")
