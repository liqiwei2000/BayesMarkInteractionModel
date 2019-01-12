# ============================================================================================================
# Set working directory
setwd("~/Dropbox/research_qbrc/bayesian_cin/shared/paper/aoas/github");

# Load libraries
# (empty)

# Load functions
Rcpp::sourceCpp('functions_x2.cpp');

# Load settings
shape = c(1, 3, 2, 4, 0, 8);
# ============================================================================================================


# Figure S2
n <- 100;
c_max <- 0.4;
Q <- 3;
set.seed(2);
x <- runif(n, 0, 1);
y <- runif(n, 0, 1);
z <- sample(1:Q, n, replace = TRUE);
build <- dist_list(x, y, c_max);
edge_all <- build$edge;
distance_all <- build$distance;
edge_03 <- edge_all[distance_all <= 0.3,];
distance_03 <- distance_all[distance_all <= 0.3];
edge_02 <- edge_all[distance_all <= 0.2,];
distance_02 <- distance_all[distance_all <= 0.2];
edge_01 <- edge_all[distance_all <= 0.1,];
distance_01 <- distance_all[distance_all <= 0.1];
# Figure S2 (a)
par(pty = "s", mar =c(2.5, 3, 1, 1), las = 1);
plot(x, y, xlab = "", ylab = "", pch = shape[z], col = z, main = "", cex.lab = 1.5, cex.axis = 1.5, cex = 1);
for (i in 1:dim(edge_03)[1]) {
  lines(c(x[edge_03[i, 1]], x[edge_03[i, 2]]), c(y[edge_03[i, 1]], y[edge_03[i, 2]]), col = 8);
}
points(x, y, pch = shape[z], col = z, cex = 1);
# Figure S2 (b)
par(pty = "s", mar =c(2.5, 3, 1, 1), las = 1);
plot(x, y, xlab = "", ylab = "", pch = shape[z], col = z, main = "", cex.lab = 1.5, cex.axis = 1.5, cex = 1);
for (i in 1:dim(edge_02)[1]) {
  lines(c(x[edge_02[i, 1]], x[edge_02[i, 2]]), c(y[edge_02[i, 1]], y[edge_02[i, 2]]), col = 8);
}
points(x, y, pch = shape[z], col = z, cex = 1);
# Figure S2 (c)
par(pty = "s", mar =c(2.5, 3, 1, 1), las = 1);
plot(x, y, xlab = "", ylab = "", pch = shape[z], col = z, main = "", cex.lab = 1.5, cex.axis = 1.5, cex = 1);
for (i in 1:dim(edge_01)[1]) {
  lines(c(x[edge_01[i, 1]], x[edge_01[i, 2]]), c(y[edge_01[i, 1]], y[edge_01[i, 2]]), col = 8);
}
points(x, y, pch = shape[z], col = z, cex = 1);