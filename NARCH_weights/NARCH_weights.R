# ------------------------------------------------------------------------------
# NARCH_weights
#
# Spatial weight matrices for the cryptocurrency volatility network.
#
# Cryptocurrencies have no geographical location, so the neighbourhood structure
# has to be built from the data. Here it comes from the correlation matrix of
# daily log returns, under two competing definitions:
#
#   threshold   coins i and j are neighbours if |cor(r_i, r_j)| > 0.15
#   k-NN        each coin is linked to its k = 5 most correlated coins
#
# Both matrices are row-standardised: every row sums to one, so a weight w_ij is
# the share of coin i's neighbourhood accounted for by coin j. That also means
# the size of a weight is governed by how many neighbours a coin has, roughly
# 1/45 under the threshold rule and 1/5 under k-NN, which is why the two panels
# below carry separate colour bars.
# ------------------------------------------------------------------------------

# ---- data: daily log returns, 50 coins x 1516 days ---------------------------
returns_prices  <- as.matrix(read.csv("returns_prices.csv",  check.names = FALSE))
returns_volumes <- as.matrix(read.csv("returns_volumes.csv", check.names = FALSE))
coins <- colnames(returns_prices)
n     <- length(coins)

# ---- weight matrix builders --------------------------------------------------

# Row-standardise: zero diagonal, every row sums to one. A coin left without any
# neighbour keeps an all-zero row rather than producing a division by zero.
row_standardise <- function(W) {
    W[!is.finite(W)] <- 0
    diag(W)          <- 0
    row_sums         <- rowSums(W)
    row_sums[row_sums == 0] <- 1
    W / row_sums
}

# Absolute correlation threshold
weights_threshold <- function(returns, tau = 0.15) {
    C <- cor(returns)
    C[!is.finite(C)] <- 0
    diag(C)          <- 0
    row_standardise((abs(C) > tau) * abs(C))
}

# k nearest neighbours by absolute correlation
weights_knn <- function(returns, k = 5) {
    C <- cor(returns)
    C[!is.finite(C)] <- 0
    diag(C)          <- 0
    W <- matrix(0, ncol(returns), ncol(returns))
    for (i in seq_len(ncol(returns))) {
        ranked    <- order(abs(C[i, ]), decreasing = TRUE)
        neighbours <- head(ranked[ranked != i], k)
        W[i, neighbours] <- abs(C[i, neighbours])
    }
    row_standardise(W)
}

W_prices  <- weights_threshold(returns_prices,  tau = 0.15)   # price model
W_volumes <- weights_knn(returns_volumes,       k   = 5)      # volume model

# ---- what the two networks look like -----------------------------------------
describe <- function(W, label) {
    links <- W[W > 0]
    cat(sprintf("%-10s density %4.1f%%   neighbours per coin (median) %2d   median weight %.4f   max %.4f\n",
                label, 100 * length(links) / length(W),
                median(rowSums(W > 0)), median(links), max(links)))
}
describe(W_prices,  "prices")
describe(W_volumes, "volumes")

# ---- figure ------------------------------------------------------------------
# Colour is scaled by quantile within each panel. A linear scale would wash the
# price panel out: 99 per cent of its weights sit below 0.032 while a handful of
# stablecoin pairs reach 0.73, because those coins have as few as two neighbours.

palette_blues <- function(k) {
    colorRampPalette(c("#deebf7", "#9ecae1", "#4292c6", "#2171b5", "#08306b"))(k)
}

quantile_breaks <- function(W, n_bins = 100) {
    links  <- W[W > 0]
    breaks <- unique(as.numeric(quantile(links, seq(0, 1, length.out = n_bins + 1))))
    if (length(breaks) < 3) breaks <- seq(min(links), max(links), length.out = n_bins + 1)
    breaks[1] <- breaks[1] - 1e-12
    breaks
}

label_every <- 5
ticks       <- seq(1, n, by = label_every)

draw_panel <- function(W, title) {
    breaks <- quantile_breaks(W)
    Wp     <- W
    Wp[Wp <= 0] <- NA                       # no link stays blank
    image(1:n, 1:n, t(Wp[n:1, ]), col = palette_blues(length(breaks) - 1),
          breaks = breaks, axes = FALSE, xlab = "", ylab = "",
          main = title, font.main = 1, cex.main = 1.95)
    axis(1, at = ticks,         labels = coins[ticks], las = 2, cex.axis = 1.5, tcl = -0.4, mgp = c(3, 0.7, 0))
    axis(2, at = n + 1 - ticks, labels = coins[ticks], las = 2, cex.axis = 1.5, tcl = -0.4, mgp = c(3, 0.7, 0))
    mtext("coin j (neighbour)", side = 1, line = 13.4, cex = 1.72)
    mtext("coin i (receiving)", side = 2, line = 13.8, cex = 1.72)
    box(col = "grey55")
}

draw_colour_bar <- function(W) {
    breaks <- quantile_breaks(W)
    n_bins <- length(breaks) - 1
    image(1, 1:n_bins, matrix(1:n_bins, 1), col = palette_blues(n_bins),
          axes = FALSE, xlab = "", ylab = "")
    probs  <- c(0, 0.25, 0.50, 0.75, 0.90, 1)
    labels <- formatC(as.numeric(quantile(W[W > 0], probs)), format = "f", digits = 3)
    labels[length(labels)] <- paste0(labels[length(labels)], " (max)")
    axis(4, at = 1 + probs * (n_bins - 1), labels = labels, las = 1, cex.axis = 1.4,
         tcl = -0.3, mgp = c(3, 0.6, 0))
    mtext(expression(w[ij]), side = 3, line = 0.7, cex = 1.72)
    box(col = "grey55")
}

png("NARCH_weights.png", width = 4900, height = 2450, res = 300)
layout(matrix(1:5, 1, 5), widths = c(6, 1.15, 0.5, 6, 1.15))
par(mar = c(15.2, 15.6, 3.6, 0.4)); draw_panel(W_prices,  "(a) Prices: correlation-thresholded, |r| > 0.15")
par(mar = c(15.2,  0.4, 3.6, 5.4)); draw_colour_bar(W_prices)
par(mar = c(0, 0, 0, 0));           plot.new()
par(mar = c(15.2, 15.6, 3.6, 0.4)); draw_panel(W_volumes, "(b) Volumes: k-nearest-neighbour, k = 5")
par(mar = c(15.2,  0.4, 3.6, 7.2)); draw_colour_bar(W_volumes)
dev.off()

cat("written NARCH_weights.png\n")
