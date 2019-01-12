# ============================================================================================================
# Set working directory
setwd("~/Dropbox/research_qbrc/bayesian_cin/shared/paper/aoas/github");

# Load libraries
library(spatstat);

# Load functions
source('functions.R');

# Load settings
shape = c(1, 3, 2, 4, 0, 8);
# ============================================================================================================


# Figure S10
load("data/data_nlst_examples.Rdata");
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
load("result/result_nlst_examples.Rdata")
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
load("result/result_nlst_examples.Rdata")
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
