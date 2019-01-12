# ============================================================================================================
# Set working directory
setwd("~/Dropbox/research_qbrc/bayesian_cin/shared/paper/aoas/github");

# Load libraries
library(survival);
library(plyr);
library(mclust);
library(fmsb);

# Load functions
source('functions.R');

# Load settings
# (empty)
# ============================================================================================================


# FIG 5 [Left]
load("result/result_nlst_summary.Rdata");
Q <- 3;
interaction <- matrix(NA, ncol = Q*Q, nrow = dim(data)[1])
for (i in 1:dim(data)[1]) {
  interaction[i,] <- c(Theta2interaction(array2matrix_r(c(data$lym_lym[i], data$lym_str[i], data$lym_tum[i], data$str_str[i], data$str_tum[i], data$tum_tum[i]), Q = 3)))
}
colnames(interaction) <- c("lym.given.lym", "str.given.lym", "tum.given.lym", "lym.given.str", "str.given.str", "tum.given.str", "lym.given.tum", "str.given.tum", "tum.given.tum")
proportion <- matrix(NA, ncol = Q, nrow = dim(data)[1])
for (i in 1:dim(data)[1]) {
  proportion[i,] <- omega2proportion(c(data$lym[i], data$str[i], data$tum[i]))
}
colnames(proportion) <- c("pi_lym", "pi_str", "pi_tum");
data <- cbind(data, proportion, interaction);
data[, c("pi_lym", "pi_str", "str.given.lym", "tum.given.lym", "lym.given.str", "tum.given.str", "lym.given.tum", "str.given.tum")] <- data[, c("pi_lym", "pi_str", "str.given.lym", "tum.given.lym", "lym.given.str", "tum.given.str", "lym.given.tum", "str.given.tum")]*100;

## survival_time_new is defined as the time between biopsy and death or the end of study, which comes first
## survival_time_new is defined as the time between the beginning of enroll and death or the end of study, which comes first
data_patient <- ddply(data, c("patient_id", "survival_time_new", "dead"), summarise, str.given.lym = mean(str.given.lym), tum.given.lym = mean(tum.given.lym), lym.given.str = mean(lym.given.str), tum.given.str = mean(tum.given.str), lym.given.tum = mean(lym.given.tum), str.given.tum = mean(str.given.tum), pi_lym = mean(pi_lym), pi_str = mean(pi_str));
m = Mclust(data_patient[,4:11], modelNames = "VVV"); # ellipsoidal, varying volume, shape, and orientation
par(pty = "s", xpd = F, mar = c(4.5, 4.5, 1, 1));
plot(1:length(m$BIC[,1]), m$BIC[, 1], xlab = "Number of groups", ylab = "BIC", main = "", xaxt="n", yaxt="n", type = "b", cex.lab = 1.5, cex.axis = 1.5, lwd = 1.5, col = 1, cex = 1.5);
axis(2, at = seq(-8000, -7000, by = 200), labels = seq(-8000, -7000, by = 200), cex.axis = 1);
axis(1, at = 1:length(m$BIC[,1]), labels = 1:length(m$BIC[,1]), cex.axis = 1.5);
abline(v = 3, col = 2, lty = 2, lwd = 1.5)
points(3, m$BIC[3, 1], cex = 1.5, col = 2)

# FIG 5 [Middle]
group <- summary(m)$classification;
data_patient$group <- group;
data_group <- data_patient[, 4:12]
data_group <- ddply(data_group, c("group"), summarise, str.given.lym = mean(str.given.lym), tum.given.lym = mean(tum.given.lym), lym.given.str = mean(lym.given.str), tum.given.str = mean(tum.given.str), lym.given.tum = mean(lym.given.tum), str.given.tum = mean(str.given.tum), pi_lym = mean(pi_lym), pi_str = mean(pi_str));
data_group <- rbind(rep(100, 9), data_group);
data_group <- rbind(rep(0, 9), data_group)
par(pty = "s", xpd = TRUE, mar = c(1, 1, 1, 1));
radarchart(data_group[, -1], pcol = c(1, 2, 3), plwd = 1.5, plty = c(1, 5, 6), axislabcol = 1, vlabels = c(expression(phi["str,lym"]), expression(phi["tum,lym"]), expression(phi["lym,str"]), expression(phi["tum,str"]), expression(phi["lym,tum"]), expression(phi["str,tum"]), expression(pi["lym"]), expression(pi["str"])), vlcex = 1.5)

# FIG 5 [Right]
group <- summary(m)$classification;
surv <- survfit(Surv(survival_time_new, dead) ~ group, data = data_patient);
par(pty = "s", mar = c(4.5, 4.5, 1, 1), las = 1);
plot(surv, lwd = 1.5, col = 1:3, lty = c(1, 5, 6), xlab = "Days", ylab = "Survival probability", cex.axis = 1.5, cex.lab = 1.5, cex = 1.5)
legend("bottomright", inset = c(0, 0), bty = "n", c("Group 1", "Group 2", "Group 3"), cex = 1.75, lwd = 1.5, col = 1:3, lty = c(1, 5, 6));
logrank <- survdiff(Surv(survival_time_new, dead) ~ group, data = data_patient);
print(pchisq(logrank$chisq,1, lower.tail=F))


