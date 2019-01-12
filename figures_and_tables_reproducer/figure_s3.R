# ============================================================================================================
# Set working directory
setwd("~/Dropbox/research_qbrc/bayesian_cin/shared/paper/aoas/github");

# Load libraries
library(spatstat);

# Load functions
# (empty)

# Load settings
shape = c(1, 3, 2, 4, 0, 8);
# ============================================================================================================


# Figure S3
load("data/data_simulated_example_high_repulsion.Rdata")
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

load("data/data_simulated_example_low_repulsion.Rdata")
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

load("data/data_simulated_example_low_attraction.Rdata")
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

load("data/data_simulated_example_high_attraction.Rdata")
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
