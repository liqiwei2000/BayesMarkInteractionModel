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


# FIG 2 [Left]
load("result/result_amacrine_chain_1.Rdata");
M_1 <- Y;
load("result/result_amacrine_chain_2.Rdata");
M_2 <- Y;
load("result/result_amacrine_chain_3.Rdata");
M_3 <- Y;
load("result/result_amacrine_chain_4.Rdata");
M_4 <- Y;
omega <- rbind(M_1$omega[burn:iter,], M_2$omega[burn:iter,], M_3$omega[burn:iter,], M_4$omega[burn:iter,]);
theta <- rbind(M_1$theta[burn:iter,], M_2$theta[burn:iter,], M_3$theta[burn:iter,], M_4$theta[burn:iter,]);
lambda <- c(M_1$lambda[burn:iter], M_2$lambda[burn:iter], M_3$lambda[burn:iter], M_4$lambda[burn:iter]);
d <- seq(0, 0.5, by = 0.01);
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
par(pty = "s", xpd = TRUE, mar=c(2.5, 3, 1, 1), las = 1);
plot(NULL, xlim = c(min(d), max(d)), ylim = c(0, 1), xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 1.5, lwd = 1.5);
lines(d, prob_mean[, 2], lwd = 1.5, col = 3, lty = 6);
lines(d, prob_mean[, 3], lwd = 1.5, col = 3, lty = 3);
lines(d, prob_mean[, 1], lwd = 1.5, col = 1, lty = 1);
lines(d, prob_mean[, 4], lwd = 1.5, col = 2, lty = 5);
legend("topright", bty = "n", cex = 1.2, c(expression(MIF["off | off"]), expression(MIF["on | off"]), expression(MIF["off | on"]), expression(MIF["on | on"])), lty = c(1, 6, 3, 5), col = c(1, 3, 3, 2), lwd = 1.5);


# FIG 2 [Middle]
load("result/result_betacells_chain_1.Rdata");
M_1 <- Y;
load("result/result_betacells_chain_2.Rdata");
M_2 <- Y;
load("result/result_betacells_chain_3.Rdata");
M_3 <- Y;
load("result/result_betacells_chain_4.Rdata");
M_4 <- Y;
omega <- rbind(M_1$omega[burn:iter,], M_2$omega[burn:iter,], M_3$omega[burn:iter,], M_4$omega[burn:iter,]);
theta <- rbind(M_1$theta[burn:iter,], M_2$theta[burn:iter,], M_3$theta[burn:iter,], M_4$theta[burn:iter,]);
lambda <- c(M_1$lambda[burn:iter], M_2$lambda[burn:iter], M_3$lambda[burn:iter], M_4$lambda[burn:iter]);
d <- seq(0, 0.5, by = 0.01);
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
par(pty = "s", xpd = TRUE, mar=c(2.5, 3, 1, 1), las = 1);
plot(NULL, xlim = c(min(d), max(d)), ylim = c(0, 1), xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 1.5, lwd = 1.5);
lines(d, prob_mean[, 2], lwd = 1.5, col = 3, lty = 6);
lines(d, prob_mean[, 3], lwd = 1.5, col = 3, lty = 3);
lines(d, prob_mean[, 1], lwd = 1.5, col = 1, lty = 1);
lines(d, prob_mean[, 4], lwd = 1.5, col = 2, lty = 5);
legend("topright", bty = "n", cex = 1.2, c(expression(MIF["off | off"]), expression(MIF["on | off"]), expression(MIF["off | on"]), expression(MIF["on | on"])), lty = c(1, 6, 3, 5), col = c(1, 3, 3, 2), lwd = 1.5);

# FIG 2 [Right]
load("result/result_lansing_chain_1.Rdata");
M_1 <- Y;
load("result/result_lansing_chain_2.Rdata");
M_2 <- Y;
load("result/result_lansing_chain_3.Rdata");
M_3 <- Y;
load("result/result_lansing_chain_4.Rdata");
M_4 <- Y;
omega <- rbind(M_1$omega[burn:iter,], M_2$omega[burn:iter,], M_3$omega[burn:iter,], M_4$omega[burn:iter,]);
theta <- rbind(M_1$theta[burn:iter,], M_2$theta[burn:iter,], M_3$theta[burn:iter,], M_4$theta[burn:iter,]);
lambda <- c(M_1$lambda[burn:iter], M_2$lambda[burn:iter], M_3$lambda[burn:iter], M_4$lambda[burn:iter]);
d <- seq(0, 0.5, by = 0.01);
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
par(pty = "s", xpd = TRUE, mar=c(2.5, 3, 1, 1), las = 1);
plot(d, prob_mean[, 1], type = "l", ylim = c(0, 1), xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 1.5, lwd = 1.5)
lines(d, prob_mean[, 8], lwd = 1.5, col = 2, lty = 6);
lines(d, prob_mean[, 15], lwd = 1.5, col = 3, lty = 5);
lines(d, prob_mean[, 22], lwd = 1.5, col = 4, lty = 2);
lines(d, prob_mean[, 29], lwd = 1.5, col = 5, lty = 3);
lines(d, prob_mean[, 36], lwd = 1.5, col = 6, lty = 4);
legend("topright", bty = "n", cex = 1.1, c(expression(MIF["black oak | black oak"]), expression(MIF["hickory | hickory"]), expression(MIF["maple | maple"]), expression(MIF["misc | misc"]), expression(MIF["red oak | red oak"]), expression(MIF["white oak | white oak"])), lty = c(1, 6, 5, 2, 3, 4), col = 1:6, lwd = 1.5);
