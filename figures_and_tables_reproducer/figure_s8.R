# ============================================================================================================
# Set working directory
setwd("~/Dropbox/research_qbrc/bayesian_cin/shared/paper/aoas/github");

# Load libraries
library(lattice);

# Load functions
source('functions.R');

# Load settings
# (empty)
# ============================================================================================================


# Figure S8
load("result/result_lansing_chain_1.Rdata");
M_1 <- Y;
load("result/result_lansing_chain_2.Rdata");
M_2 <- Y;
load("result/result_lansing_chain_3.Rdata");
M_3 <- Y;
load("result/result_lansing_chain_4.Rdata");
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
plot(levelplot(t(Phi), at = seq(0, 1, by = 0.2), col.regions = gray(0:100/100)[100:1], xlab = "", ylab = "", scales = list(x = list(cex = 1.5), y = list(cex = 1.5)), colorkey = list(at = seq(0, 1, by = 0.2), labels = list(cex = 1.5)), panel = myPanel));
