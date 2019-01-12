# ============================================================================================================
# Set working directory
setwd("~/Dropbox/research_qbrc/bayesian_cin/shared/paper/aoas/github");

# Load libraries
library(spatstat);
library(lattice);

# Load functions
Rcpp::sourceCpp('functions_x2.cpp');
source('functions.R');

# Load settings
shape = c(1, 3, 2, 4, 0, 8);
# ============================================================================================================



# ============================================================================================================
# Figure S1
c <- 0.05;
d <- seq(0, 1, by = 0.001);
e1000 <- exp(-500*d);
e100 <- exp(-200*d);
e50 <- exp(-50*d);
e20 <- exp(-20*d);
par(pty = "m", mar = c(5, 5, 1, 1), las = 1);
plot(NULL, xlim = c(0, 0.3), ylim = c(0, 1), xlab = expression(d), ylab = expression(exp(-lambda*d)), main = "", cex.axis = 1.5, cex.lab = 1.5, cex = 1.5);
lines(d, e20, type = "l", lty = 1, lwd = 3, col = 1)
lines(d, e50, type = "l", lty = 6, lwd = 3, col = 2)
lines(d, e100, type = "l", lty = 5, lwd = 3, col = 3)
lines(d, e1000, type = "l", lty = 3, lwd = 3, col = 4)
legend("topright", inset = c(0.05, 0), bty = "n", cex = 1.5, c(expression(lambda==20), expression(lambda==50), expression(lambda==200), expression(lambda==500)), lty = c(1, 6, 5, 3), col = 1:4, lwd = 1.5);
# ============================================================================================================



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
# ============================================================================================================



# ============================================================================================================
# Figure S3
load("data_simulated_example_high_repulsion.Rdata")
x <- data_x[[1]];
y <- data_y[[1]];
z <- data_z[[1]];
Q <- 2
# Figure S3(a)
par(pty = "s", mar=c(2.5, 5, 1, 1), xpd = FALSE, las = 1);
plot(x, y, xlim = c(0, 1), ylim = c(0, 1), col = z, pch = shape[z], main = "", xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 0.5);
# Figure S3(e)
par(pty = "s", mar=c(2.5, 5, 1, 1), xpd = FALSE, las = 1);
plot(x, y, xlim = c(0.4, 0.6), ylim = c(0.4, 0.6), col = z, pch = shape[z], main = "", xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 1.5);
# Figure S3(e)
temp <- ppp(x, y, marks = as.factor(z));
r <- seq(0, 1, by = 0.01);
mcf_11 <- markconnect(temp, "1", "1", r);
mcf_12 <- markconnect(temp, "1", "2", r);
mcf_21 <- markconnect(temp, "2", "1", r);
mcf_22 <- markconnect(temp, "2", "2", r);
par(pty = "s", mar=c(2.5, 5, 1, 1), xpd = FALSE, las = 1);
plot(r, mcf_11$iso, type = "l", ylab = "", xlab = "", ylim = c(0, 1), cex = 1.5, cex.axis = 1.5, cex.lab = 1.5, lwd = 1.5);
lines(r, mcf_12$iso + mcf_21$iso, lty = 6, col = 3, lwd = 1.5);
lines(r, mcf_22$iso, lty = 5, col = 2, lwd = 1.5);
# Figure S3(m)
temp <- ppp(x, y, marks = as.factor(z));
r <- seq(0, 0.1, by = 0.001);
mk_11 <- Kcross(temp, "1", "1", r);
mk_12 <- Kcross(temp, "1", "2", r);
mk_21 <- Kcross(temp, "2", "1", r);
mk_22 <- Kcross(temp, "2", "2", r);
min <- min(mk_11$iso - pi*r*r, mk_12$iso - pi*r*r, mk_22$iso - pi*r*r);
max <- max(mk_11$iso - pi*r*r, mk_12$iso - pi*r*r, mk_22$iso - pi*r*r);
min <- -0.035;
max <- 0.035;
par(pty = "s", mar=c(2.5, 5, 1, 1), xpd = FALSE, las = 1);
plot(NULL, ylab = "", xlab = "", xlim = c(min(r), max(r)), ylim = c(min, max), cex = 1.5, cex.axis = 1.5, cex.lab = 1.5);
lines(r, mk_11$iso - pi*r*r, lty = 1, col = 1, lwd = 1.5);
lines(r, mk_12$iso - pi*r*r, lty = 6, col = 3, lwd = 1.5);
lines(r, mk_22$iso - pi*r*r, lty = 5, col = 2, lwd = 1.5);
lines(r[which(r < sqrt(-min/pi))], - pi*r[which(r < sqrt(-min/pi))]*r[which(r < sqrt(-min/pi))], lty = 7, col = 8, lwd = 1.5);

load("data_simulated_example_low_repulsion.Rdata")
x <- data_x[[1]];
y <- data_y[[1]];
z <- data_z[[1]];
Q <- 2
# Figure S3(b)
par(pty = "s", mar=c(2.5, 5, 1, 1), xpd = FALSE, las = 1);
plot(x, y, xlim = c(0, 1), ylim = c(0, 1), col = z, pch = shape[z], main = "", xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 0.5);
# Figure S3(f)
par(pty = "s", mar=c(2.5, 5, 1, 1), xpd = FALSE, las = 1);
plot(x, y, xlim = c(0.4, 0.6), ylim = c(0.4, 0.6), col = z, pch = shape[z], main = "", xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 1.5);
# Figure S3(j)
temp <- ppp(x, y, marks = as.factor(z));
r <- seq(0, 1, by = 0.01);
mcf_11 <- markconnect(temp, "1", "1", r);
mcf_12 <- markconnect(temp, "1", "2", r);
mcf_21 <- markconnect(temp, "2", "1", r);
mcf_22 <- markconnect(temp, "2", "2", r);
par(pty = "s", mar=c(2.5, 5, 1, 1), xpd = FALSE, las = 1);
plot(r, mcf_11$iso, type = "l", ylab = "", xlab = "", ylim = c(0, 1), cex = 1.5, cex.axis = 1.5, cex.lab = 1.5, lwd = 1.5);
lines(r, mcf_12$iso + mcf_21$iso, lty = 6, col = 3, lwd = 1.5);
lines(r, mcf_22$iso, lty = 5, col = 2, lwd = 1.5);
# Figure S3(n)
temp <- ppp(x, y, marks = as.factor(z));
r <- seq(0, 0.1, by = 0.001);
mk_11 <- Kcross(temp, "1", "1", r);
mk_12 <- Kcross(temp, "1", "2", r);
mk_21 <- Kcross(temp, "2", "1", r);
mk_22 <- Kcross(temp, "2", "2", r);
min <- min(mk_11$iso - pi*r*r, mk_12$iso - pi*r*r, mk_22$iso - pi*r*r);
max <- max(mk_11$iso - pi*r*r, mk_12$iso - pi*r*r, mk_22$iso - pi*r*r);
min <- -0.035;
max <- 0.035;
par(pty = "s", mar=c(2.5, 5, 1, 1), xpd = FALSE, las = 1);
plot(NULL, ylab = "", xlab = "", xlim = c(min(r), max(r)), ylim = c(min, max), cex = 1.5, cex.axis = 1.5, cex.lab = 1.5);
lines(r, mk_11$iso - pi*r*r, lty = 1, col = 1, lwd = 1.5);
lines(r, mk_12$iso - pi*r*r, lty = 6, col = 3, lwd = 1.5);
lines(r, mk_22$iso - pi*r*r, lty = 5, col = 2, lwd = 1.5);
lines(r[which(r < sqrt(-min/pi))], - pi*r[which(r < sqrt(-min/pi))]*r[which(r < sqrt(-min/pi))], lty = 7, col = 8, lwd = 1.5);

load("data_simulated_example_low_attraction.Rdata")
x <- data_x[[1]];
y <- data_y[[1]];
z <- data_z[[1]];
Q <- 2
# Figure S3(c)
par(pty = "s", mar=c(2.5, 5, 1, 1), xpd = FALSE, las = 1);
plot(x, y, xlim = c(0, 1), ylim = c(0, 1), col = z, pch = shape[z], main = "", xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 0.5);
# Figure S3(g)
par(pty = "s", mar=c(2.5, 5, 1, 1), xpd = FALSE, las = 1);
plot(x, y, xlim = c(0.4, 0.6), ylim = c(0.4, 0.6), col = z, pch = shape[z], main = "", xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 1.5);
# Figure S3(k)
temp <- ppp(x, y, marks = as.factor(z));
r <- seq(0, 1, by = 0.01);
mcf_11 <- markconnect(temp, "1", "1", r);
mcf_12 <- markconnect(temp, "1", "2", r);
mcf_21 <- markconnect(temp, "2", "1", r);
mcf_22 <- markconnect(temp, "2", "2", r);
par(pty = "s", mar=c(2.5, 5, 1, 1), xpd = FALSE, las = 1);
plot(r, mcf_11$iso, type = "l", ylab = "", xlab = "", ylim = c(0, 1), cex = 1.5, cex.axis = 1.5, cex.lab = 1.5, lwd = 1.5);
lines(r, mcf_12$iso + mcf_21$iso, lty = 6, col = 3, lwd = 1.5);
lines(r, mcf_22$iso, lty = 5, col = 2, lwd = 1.5);
# Figure S3(o)
temp <- ppp(x, y, marks = as.factor(z));
r <- seq(0, 0.1, by = 0.001);
mk_11 <- Kcross(temp, "1", "1", r);
mk_12 <- Kcross(temp, "1", "2", r);
mk_21 <- Kcross(temp, "2", "1", r);
mk_22 <- Kcross(temp, "2", "2", r);
min <- -0.00045;
max <- 0.00045;
par(pty = "s", mar=c(2.5, 5, 1, 1), xpd = FALSE, las = 1);
plot(NULL, ylab = "", xlab = "", xlim = c(min(r), max(r)), ylim = c(min, max), cex = 1.5, cex.axis = 1.5, cex.lab = 1.5);
lines(r, mk_11$iso - pi*r*r, lty = 1, col = 1, lwd = 1.5);
lines(r, mk_12$iso - pi*r*r, lty = 6, col = 3, lwd = 1.5);
lines(r, mk_22$iso - pi*r*r, lty = 5, col = 2, lwd = 1.5);
lines(r[which(r < sqrt(-min/pi))], - pi*r[which(r < sqrt(-min/pi))]*r[which(r < sqrt(-min/pi))], lty = 7, col = 8, lwd = 1.5);

load("data_simulated_example_high_attraction.Rdata")
x <- data_x[[1]];
y <- data_y[[1]];
z <- data_z[[1]];
Q <- 2
# Figure S3(d)
par(pty = "s", mar=c(2.5, 5, 1, 1), xpd = FALSE, las = 1);
plot(x, y, xlim = c(0, 1), ylim = c(0, 1), col = z, pch = shape[z], main = "", xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 0.5);
# Figure S3(h)
par(pty = "s", mar=c(2.5, 5, 1, 1), xpd = FALSE, las = 1);
plot(x, y, xlim = c(0.4, 0.6), ylim = c(0.4, 0.6), col = z, pch = shape[z], main = "", xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 1.5);
# Figure S3(l)
temp <- ppp(x, y, marks = as.factor(z));
r <- seq(0, 1, by = 0.01);
mcf_11 <- markconnect(temp, "1", "1", r);
mcf_12 <- markconnect(temp, "1", "2", r);
mcf_21 <- markconnect(temp, "2", "1", r);
mcf_22 <- markconnect(temp, "2", "2", r);
par(pty = "s", mar=c(2.5, 5, 1, 1), xpd = FALSE, las = 1);
plot(r, mcf_11$iso, type = "l", ylab = "", xlab = "", ylim = c(0, 1), cex = 1.5, cex.axis = 1.5, cex.lab = 1.5, lwd = 1.5);
lines(r, mcf_12$iso + mcf_21$iso, lty = 6, col = 3, lwd = 1.5);
lines(r, mcf_22$iso, lty = 5, col = 2, lwd = 1.5);
legend("topright", bty = "n", cex = 1.5, c(expression(MCF["1,1"](d)), expression(MCF["1,2"](d)), expression(MCF["2,2"](d))), lty = c(1, 6, 5), col = c(1, 3, 2), lwd = 1.5);
# Figure S3(p)
temp <- ppp(x, y, marks = as.factor(z));
r <- seq(0, 0.1, by = 0.001);
mk_11 <- Kcross(temp, "1", "1", r);
mk_12 <- Kcross(temp, "1", "2", r);
mk_21 <- Kcross(temp, "2", "1", r);
mk_22 <- Kcross(temp, "2", "2", r);
min <- -0.00045;
max <- 0.00045;
par(pty = "s", mar=c(2.5, 5, 1, 1), xpd = FALSE, las = 1);
plot(NULL, ylab = "", xlab = "", xlim = c(min(r), max(r)), ylim = c(min, max), cex = 1.5, cex.axis = 1.5, cex.lab = 1.5);
lines(r, mk_11$iso - pi*r*r, lty = 1, col = 1, lwd = 1.5);
lines(r, mk_12$iso - pi*r*r, lty = 6, col = 3, lwd = 1.5);
lines(r, mk_22$iso - pi*r*r, lty = 5, col = 2, lwd = 1.5);
lines(r[which(r < sqrt(-min/pi))], - pi*r[which(r < sqrt(-min/pi))]*r[which(r < sqrt(-min/pi))], lty = 7, col = 8, lwd = 1.5);
legend("topright", inset = c(0, 0.3), bty = "n", cex = 1.2, c(expression(K["off,off"](d)-pi*d^2), expression(K["off,on"](d)-pi*d^2), expression(K["on,on"](d)-pi*d^2), expression(-pi*d^2)), lty = c(1, 6, 5, 7), col = c(1, 3, 2, 8), lwd = 1.5);
# ============================================================================================================



# ============================================================================================================
# Figure S4
iter <- 50000;
burn <- iter/2;
input_path <- "~/Dropbox/research_qbrc/bayesian_cin/local/results_homogeneous_poisson_lambda=60_c=0.03/Rdata/";
interaction <- "li";
files <- list.files(input_path, pattern = paste0("interaction=", interaction))
omega_1 <- NULL;
theta_11 <- NULL;
theta_12 <- NULL;
lambda <- NULL;
for (f in files){
  load(paste0(input_path, f));
  omega_1 <- c(omega_1, mean(Y$omega[burn:iter, 1]));
  theta_11 <- c(theta_11, mean(Y$theta[burn:iter, 1]));
  theta_12 <- c(theta_12, mean(Y$theta[burn:iter, 2]));
  lambda <- c(lambda, mean(Y$lambda[burn:iter]));
}
results <- data.frame(cbind(omega_1, theta_11, theta_12));
par(mar=c(2.5, 2.5, 1, 1));
boxplot(results, boxwex = 0.15, at = (1:3) - 0.2, names = c("", "", ""), border = 2, cex = 0.5, ylim = c(-1.5, 3.5), xlim = c(0.5, 3.5), xaxt = "n", yaxt = "n");

input_path <- "~/Dropbox/research_qbrc/bayesian_cin/local/results_homogeneous_poisson_lambda=60_c=0.05/Rdata/";
files <- list.files(input_path, pattern = paste0("interaction=", interaction))
omega_1 <- NULL;
theta_11 <- NULL;
theta_12 <- NULL;
lambda <- NULL;
for (f in files){
  load(paste0(input_path, f));
  omega_1 <- c(omega_1, mean(Y$omega[burn:iter, 1]));
  theta_11 <- c(theta_11, mean(Y$theta[burn:iter, 1]));
  theta_12 <- c(theta_12, mean(Y$theta[burn:iter, 2]));
  lambda <- c(lambda, mean(Y$lambda[burn:iter]));
}
results <- data.frame(cbind(omega_1, theta_11, theta_12));
boxplot(results, add = TRUE, boxwex = 0.15, at = 1:3, border = 3, cex = 0.5, names = c(expression(omega[1]), expression(theta[11]), expression(theta[12])), cex.lab = 1.5, cex.axis = 1.5);

input_path <- "~/Dropbox/research_qbrc/bayesian_cin/local/results_homogeneous_poisson_lambda=60_c=0.1/Rdata/";
files <- list.files(input_path, pattern = paste0("interaction=", interaction))
omega_1 <- NULL;
theta_11 <- NULL;
theta_12 <- NULL;
lambda <- NULL;
for (f in files){
  load(paste0(input_path, f));
  omega_1 <- c(omega_1, mean(Y$omega[burn:iter, 1]));
  theta_11 <- c(theta_11, mean(Y$theta[burn:iter, 1]));
  theta_12 <- c(theta_12, mean(Y$theta[burn:iter, 2]));
  lambda <- c(lambda, mean(Y$lambda[burn:iter]));
}
results <- data.frame(cbind(omega_1, theta_11, theta_12));
boxplot(results, add = TRUE, boxwex = 0.15, at = (1:3) + 0.2, border = 4, cex = 0.5, xaxt = "n", yaxt = "n");

lines(c(0.6, 1.4), c(1, 1), lwd = 1.5, lty = 3)
if (interaction == "ha") {
  Thetat <- interaction2Theta(matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2)); # High attraction
} else if (interaction == "la") {
  Thetat <- interaction2Theta(matrix(c(0.7, 0.3, 0.3, 0.7), nrow = 2)); # Low attraction
} else if (interaction == "cr") {
  Thetat <- interaction2Theta(matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2)); # Complete random
} else if (interaction == "li") {
  Thetat <- interaction2Theta(matrix(c(0.3, 0.7, 0.7, 0.3), nrow = 2)); # Low inhibition
} else if (interaction == "hi") {
  Thetat <- interaction2Theta(matrix(c(0.1, 0.9, 0.9, 0.1), nrow = 2)); # High inhibition
}
lines(c(1.6, 2.4), c(Thetat[1, 1], Thetat[1, 1]), lwd = 1.5, lty = 3)
lines(c(2.6, 3.4), c(Thetat[1, 2], Thetat[1, 2]), lwd = 1.5, lty = 3)
# legend("topright", bty = "n", inset = c(0.05, 0), cex = 1.5, c(expression(c==0.03), expression(c==0.05), expression(c==0.10)), pch = rep(0, 3), col = c(2, 3, 4));
# ============================================================================================================



# ============================================================================================================
# Figure S5
load(paste0("data_spatstat_amacrine.Rdata"));
Q <- max(data_list$class);
temp <- ppp(data_list$x, data_list$y, marks = as.factor(data_list$class));

# Figure S5(a)
r <- seq(0, 0.5, by = 0.01);
mcf_11 <- markconnect(temp, "1", "1", r);
mcf_12 <- markconnect(temp, "1", "2", r);
mcf_21 <- markconnect(temp, "2", "1", r);
mcf_22 <- markconnect(temp, "2", "2", r);
par(pty = "s", xpd = TRUE, mar=c(2.5, 5, 1, 1), las = 1)
plot(NULL, ylab = "", xlab = "", xlim = c(min(r), max(r)), ylim = c(0, 1), cex = 1.5, cex.axis = 1.5, cex.lab = 1.5);
lines(r, mcf_11$iso, lty = 1, col = 1, lwd = 1.5);
lines(r, mcf_12$iso + mcf_21$iso, lty = 6, col = 3, lwd = 1.5);
lines(r, mcf_22$iso, lty = 5, col = 2, lwd = 1.5);
legend("topright", bty = "n", cex = 1.2, c(expression(MCF["off,off"](d)), expression(MCF["off,on"](d)), expression(MCF["on,on"](d))), lty = c(1, 6, 5), col = c(1, 3, 2), lwd = 1.5);
# Figure S5(d)
r <- seq(0, 0.06, by = 0.0001);
mk_11 <- Kcross(temp, "1", "1", r);
mk_12 <- Kcross(temp, "1", "2", r);
mk_21 <- Kcross(temp, "2", "1", r);
mk_22 <- Kcross(temp, "2", "2", r);
par(pty = "s", xpd = TRUE, mar=c(2.5, 5, 1, 1), las = 1)
min <- min(mk_11$iso - pi*r*r, mk_12$iso - pi*r*r, mk_22$iso - pi*r*r);
max <- max(mk_11$iso - pi*r*r, mk_12$iso - pi*r*r, mk_22$iso - pi*r*r);
plot(NULL, ylab = "", xlab = "", xlim = c(min(r), max(r)), ylim = c(min, max), cex = 1.5, cex.axis = 1.5, cex.lab = 1.5);
lines(r, mk_11$iso - pi*r*r, lty = 1, col = 1, lwd = 1.5);
lines(r, mk_12$iso - pi*r*r, lty = 6, col = 3, lwd = 1.5);
lines(r, mk_22$iso - pi*r*r, lty = 5, col = 2, lwd = 1.5);
lines(r[which(r < sqrt(-min/pi))], - pi*r[which(r < sqrt(-min/pi))]*r[which(r < sqrt(-min/pi))], lty = 7, col = 8, lwd = 1.5);
legend("topleft", bty = "n", cex = 1.2, c(expression(K["off,off"](d)-pi*d^2), expression(K["off,on"](d)-pi*d^2), expression(K["on,on"](d)-pi*d^2), expression(-pi*d^2)), lty = c(1, 6, 5, 7), col = c(1, 3, 2, 8), lwd = 1.5);

load(paste0("data_spatstat_betacells.Rdata"));
Q <- max(data_list$class);
temp <- ppp(data_list$x, data_list$y, marks = as.factor(data_list$class));

# Figure S5(b)
r <- seq(0, 0.5, by = 0.01);
mcf_11 <- markconnect(temp, "1", "1", r);
mcf_12 <- markconnect(temp, "1", "2", r);
mcf_21 <- markconnect(temp, "2", "1", r);
mcf_22 <- markconnect(temp, "2", "2", r);
par(pty = "s", xpd = TRUE, mar=c(2.5, 5, 1, 1), las = 1)
plot(NULL, ylab = "", xlab = "", xlim = c(min(r), max(r)), ylim = c(0, 1), cex = 1.5, cex.axis = 1.5, cex.lab = 1.5);
lines(r, mcf_11$iso, lty = 1, col = 1, lwd = 1.5);
lines(r, mcf_12$iso + mcf_21$iso, lty = 6, col = 3, lwd = 1.5);
lines(r, mcf_22$iso, lty = 5, col = 2, lwd = 1.5);
legend("topright", bty = "n", cex = 1.2, c(expression(MCF["off,off"](d)), expression(MCF["off,on"](d)), expression(MCF["on,on"](d))), lty = c(1, 6, 5), col = c(1, 3, 2), lwd = 1.5);

# Figure S5(e)
r <- seq(0, 0.06, by = 0.0001);
mk_11 <- Kcross(temp, "1", "1", r);
mk_12 <- Kcross(temp, "1", "2", r);
mk_21 <- Kcross(temp, "2", "1", r);
mk_22 <- Kcross(temp, "2", "2", r);
par(pty = "s", xpd = TRUE, mar=c(2.5, 5, 1, 1), las = 1)
min <- min(mk_11$iso - pi*r*r, mk_12$iso - pi*r*r, mk_22$iso - pi*r*r);
max <- max(mk_11$iso - pi*r*r, mk_12$iso - pi*r*r, mk_22$iso - pi*r*r);
plot(NULL, ylab = "", xlab = "", xlim = c(min(r), max(r)), ylim = c(min, max), cex = 1.5, cex.axis = 1.5, cex.lab = 1.5);
lines(r, mk_11$iso - pi*r*r, lty = 1, col = 1, lwd = 1.5);
lines(r, mk_21$iso - pi*r*r, lty = 6, col = 3, lwd = 1.5);
lines(r, mk_22$iso - pi*r*r, lty = 5, col = 2, lwd = 1.5);
lines(r[which(r < sqrt(-min/pi))], - pi*r[which(r < sqrt(-min/pi))]*r[which(r < sqrt(-min/pi))], lty = 7, col = 8, lwd = 1.5);
legend("bottomleft", bty = "n", cex = 1.2, c(expression(K["off,off"](d)-pi*d^2), expression(K["off,on"](d)-pi*r^2), expression(K["on,on"](d)-pi*r^2), expression(-pi*d^2)), lty = c(1, 6, 5, 7), col = c(1, 3, 2, 8), lwd = 1.5);

load(paste0("data_spatstat_lansing.Rdata"));
Q <- max(data_list$class);
temp <- ppp(data_list$x, data_list$y, marks = as.factor(data_list$class));

# Figure S5(c)
r <- seq(0, 0.5, by = 0.01);
mcf_11 <- markconnect(temp, "1", "1", r);
mcf_22 <- markconnect(temp, "2", "2", r);
mcf_33 <- markconnect(temp, "3", "3", r);
mcf_44 <- markconnect(temp, "4", "4", r);
mcf_55 <- markconnect(temp, "5", "5", r);
mcf_66 <- markconnect(temp, "6", "6", r);
par(pty = "s", xpd = TRUE, mar=c(2.5, 5, 1, 1), las = 1)
plot(r, mcf_11$iso, type = "l", ylab = "", xlab = "", ylim = c(0, 1), cex = 1.5, cex.axis = 1.5, cex.lab = 1.5, lwd = 1.5);
lines(r, mcf_22$iso, lty = 6, col = 2, lwd = 1.5);
lines(r, mcf_33$iso, lty = 5, col = 3, lwd = 1.5);
lines(r, mcf_44$iso, lty = 2, col = 4, lwd = 1.5);
lines(r, mcf_55$iso, lty = 3, col = 5, lwd = 1.5);
lines(r, mcf_66$iso, lty = 4, col = 6, lwd = 1.5);
legend("topright", bty = "n", cex = 1.2, c(expression(MCF["black oak,black oak"](d)), expression(MCF["hickory,hickory"](d)), expression(MCF["maple,maple"](d)), expression(MCF["misc,misc"](d)), expression(MCF["red oak,red oak"](d)), expression(MCF["white oak,white oak"](d))), lty = c(1, 6, 5, 2, 3, 4), col = 1:6, lwd = 1.5);

# Figure S5(f)
r <- seq(0, 0.06, by = 0.0001);
mk_11 <- Kcross(temp, "1", "1", r);
mk_22 <- Kcross(temp, "2", "2", r);
mk_33 <- Kcross(temp, "3", "3", r);
mk_44 <- Kcross(temp, "4", "4", r);
mk_55 <- Kcross(temp, "5", "5", r);
mk_66 <- Kcross(temp, "6", "6", r);
min <- min(mk_11$iso, mk_22$iso, mk_33$iso, mk_44$iso, mk_55$iso, mk_66$iso);
max <- max(mk_11$iso, mk_22$iso, mk_33$iso, mk_44$iso, mk_55$iso, mk_66$iso);
par(pty = "s", xpd = TRUE, mar=c(2.5, 5, 1, 1), las = 1)
plot(NULL, ylab = "", xlab = "", xlim = c(min(r), max(r)), ylim = c(min, max), cex = 1.5, cex.axis = 1.5, cex.lab = 1.5);
lines(r, mk_11$iso, lty = 1, col = 1, lwd = 1.5);
lines(r, mk_22$iso, lty = 6, col = 2, lwd = 1.5);
lines(r, mk_33$iso, lty = 5, col = 3, lwd = 1.5);
lines(r, mk_44$iso, lty = 2, col = 4, lwd = 1.5);
lines(r, mk_55$iso, lty = 3, col = 5, lwd = 1.5);
lines(r, mk_66$iso, lty = 4, col = 6, lwd = 1.5);
legend("topleft", bty = "n", cex = 1.2, c(expression(K["black oak,black oak"](d)), expression(K["hickory,hickory"](d)), expression(K["maple,maple"](d)), expression(K["misc,misc"](d)), expression(K["red oak,red oak"](d)), expression(K["white oak,white oak"](d))), lty = c(1, 6, 5, 2, 3, 4), col = c(1:6), lwd = 1.5);
# ============================================================================================================



# ============================================================================================================
# Figure S6
load("result_amacrine_chain_1.Rdata");
M_1 <- Y;
load("result_amacrine_chain_2.Rdata");
M_2 <- Y;
load("result_amacrine_chain_3.Rdata");
M_3 <- Y;
load("result_amacrine_chain_4.Rdata");
M_4 <- Y;
theta<- rbind(M_1$theta[burn:iter,], M_2$theta[burn:iter,], M_3$theta[burn:iter,], M_4$theta[burn:iter,]);
print(array2matrix_r(round(colMeans(theta), 3), Q)); # 0.350, -4.024, -4.024, 1.000
phi <- matrix(NA, nrow = dim(theta)[1], ncol = Q*Q);
for (i in 1:dim(theta)[1]) {
  phi[i,] <- c(Theta2interaction(array2matrix_r(theta[i,], Q)));
}
Phi <- Theta2interaction(array2matrix_r(round(colMeans(theta), 3), Q));
colnames(Phi) <- c("off", "on");
rownames(Phi) <- c("off", "on");
print(round(Phi, 3)) # 0.012 0.999 0.993 0.007
CR <- matrix(NA, nrow = Q, ncol = Q);
count <- 1;
for (j in 1:Q) {
  for (i in 1:Q) {
    CR[j, i] <- paste0("[", round(quantile(phi[, count], 0.025), 3), ",", round(quantile(phi[, count], 0.975), 3), "]");
    count <- count + 1;
  }
}

NC <- matrix(c("black", "white", "white", "black"), nrow = 2);
myPanel <- function(x, y, z, ...) {
  panel.levelplot(x, y, z, ...);
  panel.text(x, y, t(CR)[cbind(x, y)], cex = 1.5, col = NC[cbind(x, y)]);
}
levelplot(t(Phi), col.regions = gray(0:100/100)[100:1], xlab = "", ylab = "", scales = list(x = list(cex = 1.5), y = list(cex = 1.5)), colorkey = list(labels = list(cex = 1.5)), panel = myPanel);
# ============================================================================================================



# ============================================================================================================
# Figure S7
load("result_betacells_chain_1.Rdata");
M_1 <- Y;
load("result_betacells_chain_2.Rdata");
M_2 <- Y;
load("result_betacells_chain_3.Rdata");
M_3 <- Y;
load("result_betacells_chain_4.Rdata");
M_4 <- Y;
theta<- rbind(M_1$theta[burn:iter,], M_2$theta[burn:iter,], M_3$theta[burn:iter,], M_4$theta[burn:iter,]);
print(array2matrix_r(round(colMeans(theta), 3), Q)); # 0.650, -3.104, -3.104, 1.000
phi <- matrix(NA, nrow = dim(theta)[1], ncol = Q*Q);
for (i in 1:dim(theta)[1]) {
  phi[i,] <- c(Theta2interaction(array2matrix_r(theta[i,], Q)));
}
Phi <- Theta2interaction(array2matrix_r(round(colMeans(theta), 3), Q));
colnames(Phi) <- c("off", "on");
rownames(Phi) <- c("off", "on");
print(round(Phi, 3)) # 0.023 0.977 0.984 0.016
CR <- matrix(NA, nrow = Q, ncol = Q);
count <- 1;
for (j in 1:Q) {
  for (i in 1:Q) {
    CR[j, i] <- paste0("[", round(quantile(phi[, count], 0.025), 3), ",", round(quantile(phi[, count], 0.975), 3), "]");
    count <- count + 1;
  }
}

NC <- matrix(c("black", "white", "white", "black"), nrow = 2);
myPanel <- function(x, y, z, ...) {
  panel.levelplot(x, y, z, ...);
  panel.text(x, y, t(CR)[cbind(x, y)], cex = 1.5, col = NC[cbind(x, y)]);
}
levelplot(t(Phi), col.regions = gray(0:100/100)[100:1], xlab = "", ylab = "", scales = list(x = list(cex = 1.5), y = list(cex = 1.5)), colorkey = list(labels = list(cex = 1.5)), panel = myPanel);
# ============================================================================================================



# ============================================================================================================
# Figure S8
load("result_lansing_chain_1.Rdata");
M_1 <- Y;
load("result_lansing_chain_2.Rdata");
M_2 <- Y;
load("result_lansing_chain_3.Rdata");
M_3 <- Y;
load("result_lansing_chain_4.Rdata");
M_4 <- Y;
theta<- rbind(M_1$theta[burn:iter,], M_2$theta[burn:iter,], M_3$theta[burn:iter,], M_4$theta[burn:iter,]);
print(array2matrix_r(round(colMeans(theta), 3), Q)); 
# -0.066 0.978 1.449  3.836 0.970 0.997
#  0.978 0.570 1.332  1.166 1.003 1.200
#  1.449 1.332 0.495  0.955 1.108 1.202
#  3.836 1.166 0.955 -0.092 1.044 1.187
#  0.970 1.003 1.108  1.044 0.535 1.266
#  0.997 1.200 1.202  1.187 1.266 1.000
phi <- matrix(NA, nrow = dim(theta)[1], ncol = Q*Q);
for (i in 1:dim(theta)[1]) {
  phi[i,] <- c(Theta2interaction(array2matrix_r(theta[i,], Q)));
}
Phi <- Theta2interaction(array2matrix_r(round(colMeans(theta), 3), Q));
colnames(Phi) <- c("black oak", "hickory", "maple", "misc", "red oak", "white oak");
rownames(Phi) <- c("black oak", "hickory", "maple", "misc", "red oak", "white oak");
print(round(Phi, 3));
#     0.436   0.172 0.111 0.009   0.165     0.192
#       0.154   0.259 0.124 0.126   0.160     0.156
#        0.096   0.121 0.287 0.156   0.144     0.156
#          0.009   0.143 0.181 0.444   0.153     0.158
#       0.155   0.168 0.155 0.142   0.255     0.146
#     0.151   0.138 0.142 0.123   0.123     0.191
CR <- matrix(NA, nrow = Q, ncol = Q);
count <- 1;
for (j in 1:Q) {
  for (i in 1:Q) {
    CR[j, i] <- paste0("[", round(quantile(phi[, count], 0.025), 3), ",", round(quantile(phi[, count], 0.975), 3), "]");
    count <- count + 1;
  }
}
NC <- matrix("black", nrow = Q, ncol = Q);
diag(NC) <- "white";
NC[Q, Q] <- "black";
myPanel <- function(x, y, z, ...) {
  panel.levelplot(x, y, z, ...);
  panel.text(x, y, t(CR)[cbind(x, y)], cex = 1, col = t(NC)[cbind(x, y)]);
}
levelplot(t(Phi), at = seq(0, 1, by = 0.2), col.regions = gray(0:100/100)[100:1], xlab = "", ylab = "", scales = list(x = list(cex = 1.5), y = list(cex = 1.5)), colorkey = list(at = seq(0, 1, by = 0.2), labels = list(cex = 1.5)), panel = myPanel);
# ============================================================================================================



# ============================================================================================================
# Figrue S10
load("data_nlst_examples.Rdata");
id <- 1;
Q <- max(data_list[[id]]$class);

# Figure S10(a)
par(pty = "s", xpd = TRUE, mar=c(2.5, 4, 1, 1), las = 1)
plot(data_list[[id]]$x, data_list[[id]]$y, xlim = c(0, 1), ylim = c(0, 1), col = data_list[[id]]$class, pch = shape[data_list[[id]]$class], main = "", xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 0.2);

# Figure S10(b)
temp <- ppp(data_list[[id]]$x, data_list[[id]]$y, marks = as.factor(data_list[[id]]$class));
r <- seq(0, 0.05, by = 0.001);
mcf_11 <- markconnect(temp, "1", "1", r);
mcf_12 <- markconnect(temp, "1", "2", r);
mcf_13 <- markconnect(temp, "1", "3", r);
mcf_21 <- markconnect(temp, "2", "1", r);
mcf_22 <- markconnect(temp, "2", "2", r);
mcf_23 <- markconnect(temp, "2", "3", r);
mcf_31 <- markconnect(temp, "3", "1", r);
mcf_32 <- markconnect(temp, "3", "2", r);
mcf_33 <- markconnect(temp, "3", "3", r);
par(pty = "s", xpd = TRUE, mar=c(2.5, 4, 1, 1), las = 1)
plot(r, mcf_11$iso, type = "l", ylab = "", xlab = "", ylim = c(0, 1), cex = 1.5, cex.axis = 1.5, cex.lab = 1.5, lwd = 1.5);
lines(r, mcf_12$iso + mcf_21$iso, lty = 6, col = 4, lwd = 1.5);
lines(r, mcf_13$iso + mcf_31$iso, lty = 5, col = 5, lwd = 1.5);
lines(r, mcf_22$iso, lty = 2, col = 2, lwd = 1.5);
lines(r, mcf_23$iso + mcf_32$iso, lty = 3, col = 6, lwd = 1.5);
lines(r, mcf_33$iso, lty = 4, col = 3, lwd = 1.5);

# Figure S10(c)
temp <- ppp(data_list[[id]]$x, data_list[[id]]$y, marks = as.factor(data_list[[id]]$class));
r <- seq(0, 0.05, by = 0.001);
mk_11 <- Kcross(temp, "1", "1", r);
mk_12 <- Kcross(temp, "1", "2", r);
mk_13 <- Kcross(temp, "1", "3", r);
mk_21 <- Kcross(temp, "2", "1", r);
mk_22 <- Kcross(temp, "2", "2", r);
mk_23 <- Kcross(temp, "2", "3", r);
mk_31 <- Kcross(temp, "3", "1", r);
mk_32 <- Kcross(temp, "3", "2", r);
mk_33 <- Kcross(temp, "3", "3", r);
min <- -0.0025;
max <- 0.062;
par(pty = "s", xpd = TRUE, mar=c(2.5, 4, 1, 1), las = 1)
plot(NULL, ylab = "", xlab = "", xlim = c(min(r), max(r)), ylim = c(min, max), cex = 1.5, cex.axis = 1.5, cex.lab = 1.5);
lines(r, mk_11$iso - pi*r*r, lty = 1, col = 1, lwd = 1.5);
lines(r, mk_12$iso - pi*r*r, lty = 6, col = 4, lwd = 1.5);
lines(r, mk_13$iso - pi*r*r, lty = 5, col = 5, lwd = 1.5);
lines(r, mk_22$iso- pi*r*r, lty = 2, col = 2, lwd = 1.5);
lines(r, mk_23$iso - pi*r*r, lty = 3, col = 6, lwd = 1.5);
lines(r, mk_33$border- pi*r*r, lty = 4, col = 3, lwd = 1.5);
lines(r[which(r < sqrt(-min/pi))], - pi*r[which(r < sqrt(-min/pi))]*r[which(r < sqrt(-min/pi))], lty = 7, col = 8, lwd = 1.5);

# Figure S10(d)
load("result_nlst_examples.Rdata")
iter <- 50000;
burn <- iter/2;
Y <- result_list[[id]];
omega <- Y$omega[burn:iter,];
theta <- Y$theta[burn:iter,];
lambda <- Y$lambda[burn:iter]
d <- seq(0, 0.05, by = 0.001);
prob_mean <- matrix(NA, nrow = length(d), ncol = Q*Q);
count <- 1;
for (dd in d) {
  print(dd);
  temp <- matrix(NA, nrow = iter - burn + 1, ncol = Q*Q)
  for (ii in 1:(iter - burn + 1)) {
    temp[ii,] <- c(t(Theta2interaction_2(omega[ii,], array2matrix_r(theta[ii,], Q), lambda[ii], dd)));
  }
  prob_mean[count,] <- colMeans(temp);
  count <- count + 1;
}
par(pty = "s", xpd = TRUE, mar=c(2.5, 4, 1, 1), las = 1)
plot(NULL, xlim = c(min(d), max(d)), ylim = c(0, 1), xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 1.5)
lines(d, prob_mean[, 2], lwd = 1.5, col = 4, lty = 2);
lines(d, prob_mean[, 4], lwd = 1.5, col = 4, lty = 3);
lines(d, prob_mean[, 3], lwd = 1.5, col = 5, lty = 2);
lines(d, prob_mean[, 7], lwd = 1.5, col = 5, lty = 3);
lines(d, prob_mean[, 6], lwd = 1.5, col = 6, lty = 2);
lines(d, prob_mean[, 8], lwd = 1.5, col = 6, lty = 3);
lines(d, prob_mean[, 1], lwd = 1.5, col = 1, lty = 1);
lines(d, prob_mean[, 5], lwd = 1.5, col = 2, lty = 5);
lines(d, prob_mean[, 9], lwd = 1.5, col = 3, lty = 6);

id <- 2;
Q <- max(data_list[[id]]$class);

# Figure S10(e)
par(pty = "s", xpd = TRUE, mar=c(2.5, 4, 1, 1), las = 1)
plot(data_list[[id]]$x, data_list[[id]]$y, xlim = c(0, 1), ylim = c(0, 1), col = data_list[[id]]$class, pch = shape[data_list[[id]]$class], main = "", xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 0.2);

# Figure S10(f)
temp <- ppp(data_list[[id]]$x, data_list[[id]]$y, marks = as.factor(data_list[[id]]$class));
r <- seq(0, 0.05, by = 0.001);
mcf_11 <- markconnect(temp, "1", "1", r);
mcf_12 <- markconnect(temp, "1", "2", r);
mcf_13 <- markconnect(temp, "1", "3", r);
mcf_21 <- markconnect(temp, "2", "1", r);
mcf_22 <- markconnect(temp, "2", "2", r);
mcf_23 <- markconnect(temp, "2", "3", r);
mcf_31 <- markconnect(temp, "3", "1", r);
mcf_32 <- markconnect(temp, "3", "2", r);
mcf_33 <- markconnect(temp, "3", "3", r);
par(pty = "s", xpd = TRUE, mar=c(2.5, 4, 1, 1), las = 1)
plot(r, mcf_11$iso, type = "l", ylab = "", xlab = "", ylim = c(0, 1), cex = 1.5, cex.axis = 1.5, cex.lab = 1.5, lwd = 1.5);
lines(r, mcf_12$iso + mcf_21$iso, lty = 6, col = 4, lwd = 1.5);
lines(r, mcf_13$iso + mcf_31$iso, lty = 5, col = 5, lwd = 1.5);
lines(r, mcf_22$iso, lty = 2, col = 2, lwd = 1.5);
lines(r, mcf_23$iso + mcf_32$iso, lty = 3, col = 6, lwd = 1.5);
lines(r, mcf_33$iso, lty = 4, col = 3, lwd = 1.5);
legend("topright", bty = "n", cex = 1.2, c(expression(MCF["lym,lym"](d)), expression(MCF["lym,str"](d)), expression(MCF["lym,tum"](d)), expression(MCF["str,str"](d)), expression(MCF["str,tum"](d)), expression(MCF["tum,tum"](d))), lty = c(1, 6, 5, 2, 3, 4), col = c(1, 4, 5, 2, 6, 3), lwd = 1.5);

# Figure S10(c)
temp <- ppp(data_list[[id]]$x, data_list[[id]]$y, marks = as.factor(data_list[[id]]$class));
r <- seq(0, 0.05, by = 0.001);
mk_11 <- Kcross(temp, "1", "1", r);
mk_12 <- Kcross(temp, "1", "2", r);
mk_13 <- Kcross(temp, "1", "3", r);
mk_21 <- Kcross(temp, "2", "1", r);
mk_22 <- Kcross(temp, "2", "2", r);
mk_23 <- Kcross(temp, "2", "3", r);
mk_31 <- Kcross(temp, "3", "1", r);
mk_32 <- Kcross(temp, "3", "2", r);
mk_33 <- Kcross(temp, "3", "3", r);
min <- -0.0025;
max <- 0.062;
par(pty = "s", xpd = TRUE, mar=c(2.5, 4, 1, 1), las = 1)
plot(NULL, ylab = "", xlab = "", xlim = c(min(r), max(r)), ylim = c(min, max), cex = 1.5, cex.axis = 1.5, cex.lab = 1.5);
lines(r, mk_11$iso - pi*r*r, lty = 1, col = 1, lwd = 1.5);
lines(r, mk_12$iso - pi*r*r, lty = 6, col = 4, lwd = 1.5);
lines(r, mk_13$iso - pi*r*r, lty = 5, col = 5, lwd = 1.5);
lines(r, mk_22$iso- pi*r*r, lty = 2, col = 2, lwd = 1.5);
lines(r, mk_23$iso - pi*r*r, lty = 3, col = 6, lwd = 1.5);
lines(r, mk_33$border- pi*r*r, lty = 4, col = 3, lwd = 1.5);
lines(r[which(r < sqrt(-min/pi))], - pi*r[which(r < sqrt(-min/pi))]*r[which(r < sqrt(-min/pi))], lty = 7, col = 8, lwd = 1.5);
legend("topleft", bty = "n", cex = 1.2, c(expression(K["lym,lym"](d)-pi*d^2), expression(K["lym,str"](d)-pi*d^2), expression(K["lym,tum"](d)-pi*d^2), expression(K["str,str"](d)-pi*d^2), expression(K["str,tum"](d)-pi*d^2), expression(K["tum,tum"](d)-pi*d^2), expression(-pi*d^2)), lty = c(1, 6, 5, 2, 3, 4, 7), col = c(1, 4, 5, 2, 6, 3, 8), lwd = 1.5);

# Figure S10(d)
load("result_nlst_examples.Rdata")
iter <- 50000;
burn <- iter/2;
Y <- result_list[[id]];
omega <- Y$omega[burn:iter,];
theta <- Y$theta[burn:iter,];
lambda <- Y$lambda[burn:iter]
d <- seq(0, 0.05, by = 0.001);
prob_mean <- matrix(NA, nrow = length(d), ncol = Q*Q);
count <- 1;
for (dd in d) {
  print(dd);
  temp <- matrix(NA, nrow = iter - burn + 1, ncol = Q*Q)
  for (ii in 1:(iter - burn + 1)) {
    temp[ii,] <- c(t(Theta2interaction_2(omega[ii,], array2matrix_r(theta[ii,], Q), lambda[ii], dd)));
  }
  prob_mean[count,] <- colMeans(temp);
  count <- count + 1;
}
par(pty = "s", xpd = TRUE, mar=c(2.5, 4, 1, 1), las = 1)
plot(NULL, xlim = c(min(d), max(d)), ylim = c(0, 1), xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 1.5)
lines(d, prob_mean[, 2], lwd = 1.5, col = 4, lty = 2);
lines(d, prob_mean[, 4], lwd = 1.5, col = 4, lty = 3);
lines(d, prob_mean[, 3], lwd = 1.5, col = 5, lty = 2);
lines(d, prob_mean[, 7], lwd = 1.5, col = 5, lty = 3);
lines(d, prob_mean[, 6], lwd = 1.5, col = 6, lty = 2);
lines(d, prob_mean[, 8], lwd = 1.5, col = 6, lty = 3);
lines(d, prob_mean[, 1], lwd = 1.5, col = 1, lty = 1);
lines(d, prob_mean[, 5], lwd = 1.5, col = 2, lty = 5);
lines(d, prob_mean[, 9], lwd = 1.5, col = 3, lty = 6);
legend("topright", inset = c(0, -0.04), bty = "n", cex = 1.2, c(expression(MIF["lym | lym"](d)), expression(MIF["str | lym"](d)), expression(MIF["tum | lym"](d)), expression(MIF["lym | str"](d)), expression(MIF["str | str"](d)), expression(MIF["tum | str"](d))), lty = c(1, 2, 2, 3, 5, 2), col = c(1, 4, 5, 4, 2, 6), lwd = 1.5);
legend("bottomright", inset = c(0, 0.01), bty = "n", cex = 1.2, c(expression(MIF["lym | tum"](d)), expression(MIF["str | tum"](d)), expression(MIF["tum | tum"](d))), lty = c(3, 3, 6), col = c(5, 6, 3), lwd = 1.5);
# ============================================================================================================



# ============================================================================================================
# Figure S11
load("result_nlst_summary.Rdata");
input_path <- "~/Dropbox/research_qbrc/bayesian_cin/local/results_nlst_mcf/";
mode <- "iso";
r <- seq(0, 0.2, by = 0.01);

# Initialize results
files <- list.files(input_path);
lym_lym_results <- matrix(NA, nrow = dim(data)[1], ncol = length(r));
colnames(lym_lym_results) <- r;
str_str_results <- matrix(NA, nrow = dim(data)[1], ncol = length(r));
colnames(str_str_results) <- r;
tum_tum_results <- matrix(NA, nrow = dim(data)[1], ncol = length(r));
colnames(tum_tum_results) <- r;
lym_str_results <- matrix(NA, nrow = dim(data)[1], ncol = length(r));
colnames(lym_str_results) <- r;
lym_tum_results <- matrix(NA, nrow = dim(data)[1], ncol = length(r));
colnames(lym_tum_results) <- r;
str_tum_results <- matrix(NA, nrow = dim(data)[1], ncol = length(r));
colnames(str_tum_results) <- r;
# Read information from the data
for (f in files) {
  # print(f);
  i <- as.numeric(substring(f, unlist(gregexpr("=", f)) + 1, nchar(f) - 6));
  
  load(paste0(input_path, f));
  
  if (mode == "iso") {
    results <- results_iso;
  } else if (mode == "trans") {
    results <- results_trans;
  }
  lym_lym_results[i,] <- results["lym-lym",];
  str_str_results[i,] <- results["str-str",];
  tum_tum_results[i,] <- results["tum-tum",];
  lym_str_results[i,] <- results["lym-str",];
  lym_tum_results[i,] <- results["lym-tum",];
  str_tum_results[i,] <- results["str-tum",];
}
# Visualize the result
results <- str_tum_results;
mcf_lwr <- NULL;
mcf_med <- NULL;
mcf_upp <- NULL;
mcf_ave <- NULL;
for (i in 1:length(r)) {
  mcf_lwr <- c(mcf_lwr, quantile(results[, i], 0.25, na.rm = TRUE));
  # mcf_med <- c(mcf_med, quantile(results[, i], 0.5, na.rm = TRUE));
  mcf_upp <- c(mcf_upp, quantile(results[, i], 0.75, na.rm = TRUE));
  mcf_ave <- c(mcf_ave, mean(results[, i], na.rm = TRUE));
}
par(pty = "m", mar = c(5, 5, 1, 1), las = 1);
plot(NULL, xlim = c(0, max(r)), ylim = c(0, 1), main = "", xlab = "d", ylab = "", cex.axis = 1.5, cex.lab = 1.5);
for (i in 1:dim(results)[1]) {
  lines(r, results[i,], col = 8)
}
lines(r, mcf_upp, lty = 2, lwd = 1.5);
lines(r, mcf_lwr, lty = 2, lwd = 1.5);
lines(r, mcf_ave, lwd = 1.5, col = 2);
# ============================================================================================================



# ============================================================================================================
# Figure S11(a)
load("result_nlst_summary.Rdata");
data <- data[order(as.numeric(data$image_id)),] 
r <- seq(0, 0.1, by = 0.01);
P <- matrix(NA, nrow = 6, ncol = length(r));
for (j in 1:length(r)) {
  results <- cbind(lym_lym_results[, j], str_str_results[, j], tum_tum_results[, j], lym_str_results[, j], lym_tum_results[, j], str_tum_results[, j]);
  colnames(results) <- c("mcf_lym_lym", "mcf_str_str", "mcf_tum_tum", "mcf_lym_str", "mcf_lym_tum", "mcf_str_tum");
  data_2 <- cbind(data, results);
  fit <- coxph(Surv(survival_time_new, dead) ~ mcf_lym_lym + mcf_lym_str + mcf_lym_tum + mcf_str_str + mcf_str_tum + freq.lym + freq.str + age + female + tobacco + cluster(patient_id), data = data_2);
  P[,j] <- summary(fit)$coefficient[1:6, 6]
}
par(pty = "m", mar = c(2.5, 5, 1, 1), las = 1);
plot(NULL, xlim = c(0, max(r)), ylim = c(0, 1), xlab = "", ylab = "P-value", cex.axis = 1.5, cex.lab = 1.5)
lines(r, P[1,], col = 1, lwd = 1.5, lty = 1);
lines(r, P[2,], col = 4, lwd = 1.5, lty = 6);
lines(r, P[3,], col = 5, lwd = 1.5, lty = 5);
lines(r, P[4,], col = 2, lwd = 1.5, lty = 2);
lines(r, P[5,], col = 6, lwd = 1.5, lty = 3);
# lines(r, P[6,], col = 3, lwd = 1.5, lty = 4);
# legend("bottomright", inset = c(0.05, 0), cex = 1.1, bty = "n", c(expression(MCF["lym,str"]), expression(MCF["lym,tum"]), expression(MCF["str,tum"])), lty = c(1, 6, 5), col = 1:3, lwd = rep(1.5, 3));
legend("bottomleft", bty = "n", cex = 1.2, c(expression(MCF["lym,lym"]), expression(MCF["lym,str"]), expression(MCF["lym,tum"])), lty = c(1, 6, 5), col = c(1, 4, 5), lwd = 1.5);
legend("bottomleft", inset = c(0.3, 0), bty = "n", cex = 1.2, c(expression(MCF["str,str"]), expression(MCF["str,tum"])), lty = c(2, 3), col = c(2, 6), lwd = 1.5);

# Figure S11(b)
input_path <- "~/Dropbox/research_qbrc/bayesian_cin/local/results_nlst_mk/";
mode <- "iso";
r <- seq(0, 0.2, by = 0.01);
# Initialize results
files <- list.files(input_path);
lym_lym_results <- matrix(NA, nrow = dim(data_info), ncol = length(r));
colnames(lym_lym_results) <- r;
str_str_results <- matrix(NA, nrow = dim(data_info), ncol = length(r));
colnames(str_str_results) <- r;
tum_tum_results <- matrix(NA, nrow = dim(data_info), ncol = length(r));
colnames(tum_tum_results) <- r;
lym_str_results <- matrix(NA, nrow = dim(data_info), ncol = length(r));
colnames(lym_str_results) <- r;
lym_tum_results <- matrix(NA, nrow = dim(data_info), ncol = length(r));
colnames(lym_tum_results) <- r;
str_tum_results <- matrix(NA, nrow = dim(data_info), ncol = length(r));
colnames(str_tum_results) <- r;
for (f in files) {
  print(f);
  i <- as.numeric(substring(f, unlist(gregexpr("=", f)) + 1, nchar(f) - 6));
  
  load(paste0(input_path, f));
  
  if (mode == "iso") {
    results <- results_iso;
  } else if (mode == "trans") {
    results <- results_trans;
  }
  index <- which(is.na(rowSums(results)));
  if (length(index) > 0) {
    results[index,] <- results_border[index,];
  }
  
  lym_lym_results[i,] <- results["lym-lym",];
  str_str_results[i,] <- results["str-str",];
  tum_tum_results[i,] <- results["tum-tum",];
  lym_str_results[i,] <- results["lym-str",];
  lym_tum_results[i,] <- results["lym-tum",];
  str_tum_results[i,] <- results["str-tum",];
}
load("result_nlst_summary.Rdata");
data <- data[order(as.numeric(data$image_id)),] 
r <- seq(0, 0.1, by = 0.01);
P <- matrix(NA, nrow = 6, ncol = length(r));
for (j in 1:length(r)) {
  results <- cbind(lym_lym_results[, j], str_str_results[, j], tum_tum_results[, j], lym_str_results[, j], lym_tum_results[, j], str_tum_results[, j]);
  colnames(results) <- c("mk_lym_lym", "mk_str_str", "mk_tum_tum", "mk_lym_str", "mk_lym_tum", "mk_str_tum");
  data_2 <- cbind(data, results);
  fit <- coxph(Surv(survival_time_new, dead) ~ mk_lym_lym + mk_lym_str + mk_lym_tum + mk_str_str + mk_str_tum + mk_tum_tum + freq.lym + freq.str + age + female + tobacco + cluster(patient_id), data = data_2);
  P[,j] <- summary(fit)$coefficient[1:6, 6]
}
par(pty = "m", mar = c(2.5, 5, 1, 1), las = 1);
plot(NULL, xlim = c(0, max(r)), ylim = c(0, 1), xlab = "", ylab = "P-value", cex.axis = 1.5, cex.lab = 1.5)
lines(r, P[1,], col = 1, lwd = 1.5, lty = 1);
lines(r, P[2,], col = 4, lwd = 1.5, lty = 6);
lines(r, P[3,], col = 5, lwd = 1.5, lty = 5);
lines(r, P[4,], col = 2, lwd = 1.5, lty = 2);
lines(r, P[5,], col = 6, lwd = 1.5, lty = 3);
lines(r, P[6,], col = 3, lwd = 1.5, lty = 4);
legend("bottomleft", inset = c(0.15, 0), bty = "n", cex = 1.2, c(expression(K["lym,lym"]), expression(K["lym,str"]), expression(K["lym,tum"])), lty = c(1, 6, 5), col = c(1, 4, 5), lwd = 1.5);
legend("bottomleft", inset = c(0.40, 0), bty = "n", cex = 1.2, c(expression(K["str,str"]), expression(K["str,tum"]), expression(K["tum,tum"])), lty = c(2, 3, 4), col = c(2, 6, 3), lwd = 1.5);
# ============================================================================================================



# ============================================================================================================
# Table S6
load("result_nlst_summary.Rdata");
data <- data[order(as.numeric(data$image_id)),] 
input_path <- "~/Dropbox/research_qbrc/bayesian_cin/local/results_nlst_mcf/";
mode <- "iso";
r <- seq(0, 0.2, by = 0.01);

# Initialize results
files <- list.files(input_path);
lym_lym_results <- matrix(NA, nrow = dim(data)[1], ncol = length(r));
colnames(lym_lym_results) <- r;
str_str_results <- matrix(NA, nrow = dim(data)[1], ncol = length(r));
colnames(str_str_results) <- r;
tum_tum_results <- matrix(NA, nrow = dim(data)[1], ncol = length(r));
colnames(tum_tum_results) <- r;
lym_str_results <- matrix(NA, nrow = dim(data)[1], ncol = length(r));
colnames(lym_str_results) <- r;
lym_tum_results <- matrix(NA, nrow = dim(data)[1], ncol = length(r));
colnames(lym_tum_results) <- r;
str_tum_results <- matrix(NA, nrow = dim(data)[1], ncol = length(r));
colnames(str_tum_results) <- r;
# Read information from the data
for (f in files) {
  # print(f);
  i <- as.numeric(substring(f, unlist(gregexpr("=", f)) + 1, nchar(f) - 6));
  
  load(paste0(input_path, f));
  
  if (mode == "iso") {
    results <- results_iso;
  } else if (mode == "trans") {
    results <- results_trans;
  }
  lym_lym_results[i,] <- results["lym-lym",];
  str_str_results[i,] <- results["str-str",];
  tum_tum_results[i,] <- results["tum-tum",];
  lym_str_results[i,] <- results["lym-str",];
  lym_tum_results[i,] <- results["lym-tum",];
  str_tum_results[i,] <- results["str-tum",];
}
j<-11;
results <- cbind(lym_lym_results[, j], str_str_results[, j], tum_tum_results[, j], lym_str_results[, j], lym_tum_results[, j], str_tum_results[, j]);
colnames(results) <- c("mcf_lym_lym", "mcf_str_str", "mcf_tum_tum", "mcf_lym_str", "mcf_lym_tum", "mcf_str_tum");
data_2 <- cbind(data, results);
fit <- coxph(Surv(survival_time_new, dead) ~ mcf_lym_lym + mcf_lym_str + mcf_lym_tum + mcf_str_str + mcf_tum_tum + freq.lym + freq.str + age + female + tobacco + cluster(patient_id), data = data_2);
print(summary(fit));
print(round(summary(fit)$coef[, c(1, 2, 3, 6)],3))
# ============================================================================================================



# ============================================================================================================
# Table S7
input_path <- "~/Dropbox/research_qbrc/bayesian_cin/local/results_nlst_mk/";
mode <- "iso";
r <- seq(0, 0.2, by = 0.01);
# Initialize results
files <- list.files(input_path);
lym_lym_results <- matrix(NA, nrow = dim(data_info), ncol = length(r));
colnames(lym_lym_results) <- r;
str_str_results <- matrix(NA, nrow = dim(data_info), ncol = length(r));
colnames(str_str_results) <- r;
tum_tum_results <- matrix(NA, nrow = dim(data_info), ncol = length(r));
colnames(tum_tum_results) <- r;
lym_str_results <- matrix(NA, nrow = dim(data_info), ncol = length(r));
colnames(lym_str_results) <- r;
lym_tum_results <- matrix(NA, nrow = dim(data_info), ncol = length(r));
colnames(lym_tum_results) <- r;
str_tum_results <- matrix(NA, nrow = dim(data_info), ncol = length(r));
colnames(str_tum_results) <- r;
for (f in files) {
  print(f);
  i <- as.numeric(substring(f, unlist(gregexpr("=", f)) + 1, nchar(f) - 6));
  
  load(paste0(input_path, f));
  
  if (mode == "iso") {
    results <- results_iso;
  } else if (mode == "trans") {
    results <- results_trans;
  }
  index <- which(is.na(rowSums(results)));
  if (length(index) > 0) {
    results[index,] <- results_border[index,];
  }
  
  lym_lym_results[i,] <- results["lym-lym",];
  str_str_results[i,] <- results["str-str",];
  tum_tum_results[i,] <- results["tum-tum",];
  lym_str_results[i,] <- results["lym-str",];
  lym_tum_results[i,] <- results["lym-tum",];
  str_tum_results[i,] <- results["str-tum",];
}
load("result_nlst_summary.Rdata");
data <- data[order(as.numeric(data$image_id)),] 
j<-11;
results <- cbind(lym_lym_results[, j], str_str_results[, j], tum_tum_results[, j], lym_str_results[, j], lym_tum_results[, j], str_tum_results[, j]);
colnames(results) <- c("mk_lym_lym", "mk_str_str", "mk_tum_tum", "mk_lym_str", "mk_lym_tum", "mk_str_tum");
data_2 <- cbind(data, results);
fit <- coxph(Surv(survival_time_new, dead) ~ mk_lym_lym + mk_lym_str + mk_lym_tum + mk_str_str + mk_str_tum + mk_tum_tum + freq.lym + freq.str + age + female + tobacco + cluster(patient_id), data = data_2);
print(summary(fit));
print(round(summary(fit)$coef[, c(1, 2, 3, 6)],3))
# ============================================================================================================



# ============================================================================================================
# Table S8

# ============================================================================================================
