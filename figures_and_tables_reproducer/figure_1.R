# ============================================================================================================
# Set working directory
setwd("~/Dropbox/research_qbrc/bayesian_cin/shared/paper/aoas/github");

# Load libraries
# (empty)

# Load functions
# (empty)

# Load settings
shape = c(1, 3, 2, 4, 0, 8);
# ============================================================================================================


# FIG 1 [Left]
load("data/data_spatstat_amacrine.Rdata");
Q <- max(data_list$class);
par(pty = "s", xpd = TRUE, mar=c(2.5, 3, 1, 1), las = 1);
plot(data_list$x, data_list$y, xlim = c(0, 1), ylim = c(0, 1), col = data_list$class, pch = shape[data_list$class], main = "", xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 1.2);
legend("bottomright", inset = c(-0.03, 0.7), legend = data_list$read.me[2:(Q + 1), 1], col = 1:Q, pch = shape[1:Q], bty = "n", cex = 1.5);

# FIG 1 [Middle]
load("data/data_spatstat_betacells.Rdata");
Q <- max(data_list$class);
par(pty = "s", xpd = TRUE, mar=c(2.5, 3, 1, 1), las = 1)
plot(data_list$x, data_list$y, xlim = c(0, 1), ylim = c(0, 1), col = data_list$class, pch = shape[data_list$class], main = "", xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 1.2);
legend("bottomright", inset = c(-0.03, 0.7), legend = data_list$read.me[2:(Q + 1), 1], col = 1:Q, pch = shape[1:Q], bty = "n", cex = 1.5);

# FIG 1 [Right]
load("data/data_spatstat_lansing.Rdata");
Q <- max(data_list$class);
par(pty = "s", xpd = TRUE, mar=c(2.5, 3, 1, 1), las = 1)
plot(data_list$x, data_list$y, xlim = c(0, 1), ylim = c(0, 1), col = data_list$class, pch = shape[data_list$class], main = "", xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 0.5);
# legend("bottomright", inset = c(-0.23, 0), legend = data_list$read.me[2:(Q + 1), 1], col = 1:Q, pch = shape[1:Q], bty = "n", cex = 0.75);