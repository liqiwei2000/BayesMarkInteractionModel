# ============================================================================================================
# Set working directory
setwd("~/Dropbox/research_qbrc/bayesian_cin/shared/paper/aoas/github");

# Load libraries
# (empty)

# Load functions
# (empty)

# Load settings
# (empty)
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
