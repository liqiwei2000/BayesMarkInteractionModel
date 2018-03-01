# ***README***
# The following script is used to analyze the marked spatial point pattern data by 
# plotting its mark connection functions powered by R package spatstat.

# Before running the following code, please first load the data using data_loader.R. The
# necessary inputs should be x, y, and z, which denote the x, y coordinates, and marks of
# each point.
# ***END***

# Load libraries
require(spatstat);



# ========================================================================================
# ========================================================================================
# Plot the data
# ========================================================================================
# ========================================================================================
par(mfcol = c(1, 1));
par(pty = "s", mar=c(2.5, 2.5, 1, 1));
plot(x, y, xlim = c(0, 1), col = z, pch = z, ylim = c(0, 1), main = "", xlab = "", 
     ylab = "", cex.lab = 1.5, cex.axis = 1.5, cex = 0.5);



# ========================================================================================
# ========================================================================================
# Plot mark connection functions
# ========================================================================================
# ========================================================================================
data <- ppp(x, y, marks = as.factor(z));
par(mfcol = c(1, 1));
par(pty = "s", mar=c(2.5, 2.5, 1, 1));
plot(alltypes(data, markconnect));



