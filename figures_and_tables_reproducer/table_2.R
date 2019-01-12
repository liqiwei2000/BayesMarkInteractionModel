# ============================================================================================================
# Set working directory
setwd("~/Dropbox/research_qbrc/bayesian_cin/shared/paper/aoas/github");

# Load libraries
library(survival);

# Load functions
source('functions.R');

# Load settings
# (empty)
# ============================================================================================================


# TABLE 2
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
##survival_time_new is defined as the time between the beginning of enroll and death or the end of study, which comes first
fit <- coxph(Surv(survival_time_new, dead) ~ str.given.lym + tum.given.lym + lym.given.str + tum.given.str + lym.given.tum + str.given.tum + lambda + pi_lym + pi_str + age + female + tobacco + cluster(patient_id), data = data)
print(summary(fit));
# print(round(summary(fit)$coef[, c(1, 2, 3, 6)],3))
