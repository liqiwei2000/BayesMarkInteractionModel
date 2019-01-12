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


# Figure S7
load("result/result_betacells_chain_1.Rdata");
M_1 <- Y;
load("result/result_betacells_chain_2.Rdata");
M_2 <- Y;
load("result/result_betacells_chain_3.Rdata");
M_3 <- Y;
load("result/result_betacells_chain_4.Rdata");
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
plot(levelplot(t(Phi), col.regions = gray(0:100/100)[100:1], xlab = "", ylab = "", scales = list(x = list(cex = 1.5), y = list(cex = 1.5)), colorkey = list(labels = list(cex = 1.5)), panel = myPanel));
