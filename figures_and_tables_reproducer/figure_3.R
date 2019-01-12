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


# FIG 3 [Left]
load("data/data_nlst_examples.Rdata");
id <- 1;
Q <- max(data_list[[id]]$class);
par(pty = "s", xpd = TRUE, mar=c(2.5, 3, 1, 1), las = 1)
plot(data_list[[id]]$x, data_list[[id]]$y, xlim = c(0, 1), ylim = c(0, 1), col = data_list[[id]]$class, pch = shape[data_list[[id]]$class], main = "", xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 0.2);
# legend("bottomright", inset = c(-0.19, 0), legend = c("lym", "str", "tum"), col = 1:Q, pch = shape[1:Q], bty = "n", cex = 1);

# FIG 3 [Right]
load("data/data_nlst_examples.Rdata");
id <- 2;
Q <- max(data_list[[id]]$class);
par(pty = "s", xpd = TRUE, mar=c(2.5, 3, 1, 1), las = 1)
plot(data_list[[id]]$x, data_list[[id]]$y, xlim = c(0, 1), ylim = c(0, 1), col = data_list[[id]]$class, pch = shape[data_list[[id]]$class], main = "", xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 0.2);
# legend("bottomright", inset = c(-0.19, 0), legend = c("lym", "str", "tum"), col = 1:Q, pch = shape[1:Q], bty = "n", cex = 1);
