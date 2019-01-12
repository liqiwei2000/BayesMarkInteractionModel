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


# Figure S5
load(paste0("data/data_spatstat_amacrine.Rdata"));
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

load(paste0("data/data_spatstat_betacells.Rdata"));
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

load(paste0("data/data_spatstat_lansing.Rdata"));
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
