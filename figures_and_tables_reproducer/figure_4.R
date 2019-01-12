# ============================================================================================================
# Set working directory
setwd("~/Dropbox/research_qbrc/bayesian_cin/shared/paper/aoas/github");

# Load libraries
library(survival);
library(plyr);

# Load functions
source('functions.R');

# Load settings
# (empty)
# ============================================================================================================


# FIG 4
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

n <- 1585;
risk <- rep(NA, n);
for (i in 1:n) {
  print(i);
  ## survival_time_new is defined as the time between biopsy and death or the end of study, which comes first
  ##survival_time_new is defined as the time between the beginning of enroll and death or the end of study, which comes first
  fit <- coxph(Surv(survival_time_new, dead) ~ str.given.lym + tum.given.lym + lym.given.str + tum.given.str + lym.given.tum + str.given.tum + lambda + pi_lym + pi_str + age + female + tobacco, data = data[-i,])
  risk[i] <- predict(fit, data[i, c("str.given.lym", "tum.given.lym", "lym.given.str", "tum.given.str", "lym.given.tum", "str.given.tum", "lambda", "pi_lym", "pi_str", "age", "female", "tobacco")], type = "risk")
}
data$risk <- risk;
data_patient <- ddply(data, c("patient_id", "survival_time_new", "dead"), summarise, risk = mean(risk));
data_patient$group <- data_patient$risk >= median(data_patient$risk);
surv <- survfit(Surv(survival_time_new, dead) ~ group, data = data_patient);
par(pty = "s", xpd = TRUE, mar=c(4.5, 4.5, 1, 1), las = 1);
plot(surv, lwd = 1.5, col = 1:3, lty = c(1, 5, 6), xlab = "Days", ylab = "Survival probability", cex.axis = 1.5, cex.lab = 1.5, cex = 1.5)
legend("bottomright", inset = c(0.0, 0), bty = "n", c("Low risk", "High risk"), cex = 1.75, lwd = 1.5, col = 1:2, lty = c(1, 5));
logrank <- survdiff(Surv(survival_time_new, dead) ~ group, data = data_patient);
print(pchisq(logrank$chisq,1, lower.tail=F))
