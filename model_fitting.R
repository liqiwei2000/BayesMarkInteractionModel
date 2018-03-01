# ***README***
# The following script is used to fit the marked spatial point pattern data to the model 
# proposed in the submitted manuscript titled "A Bayesian Mark Interaction Model for 
# Analysis of Tumor Pathology Images".

# Before running the following code, please first load the data using data_loader.R. The
# necessary inputs should be x, y, and z, which denote the x, y coordinates, and marks of
# each point.
# ***END***

# Load libraries
require(spatstat);
require(Rcpp);

# Load functions
source('functions.R');
Rcpp::sourceCpp('functions_x2.cpp');



# ========================================================================================
# ========================================================================================
# Load algorithm setting
# ========================================================================================
# ========================================================================================
Q <- 2;          # Number of marks
c <- 0.1;        # Tuning parameter c that defines the neighborhood for each point
iter <- 10000;   # Number of iterations
burn <- iter/2;  # Number of burn-in 



# ========================================================================================
# ========================================================================================
# Load hyperparameters
# ========================================================================================
# ========================================================================================
# Prior for theta
mu_theta <- 0;
sigma_theta <- 1;
# Prior for omega
mu_omega <- 1;
sigma_omega <- 1;
# Prior for lambda
a <- 0.001;
b <- 0.001;



# ========================================================================================
# ========================================================================================
# Load initial configuration
# ========================================================================================
# ========================================================================================
theta_start <- rnorm(Q*(Q + 1)/2, mu_theta, sigma_theta);
theta_start[length(theta_start)] <- 1;
Theta_start <- array2matrix_r(theta_start, Q);
omega_start <- c(1, 1);
lambda_start <- rgamma(1, 1, 1);



# ========================================================================================
# ========================================================================================
# Preprocess data (DO NOT EDIT THE CODE IN THIS BLOCK)
# ========================================================================================
# ========================================================================================
id_start <- 1;
id_end <- length(z);
build <- dist_list(x, y, c);
edge <- build$edge;
distance <- build$distance;
duplicate <- build$duplicate;
flag_start <- build$flag_start;
flag_end <- build$flag_end;



# ========================================================================================
# ========================================================================================
# Implement MCMC algorithm (DO NOT EDIT THE CODE IN THIS BLOCK)
# ========================================================================================
# ========================================================================================
start_time <- proc.time();
Y <- model_estimator(z, edge, distance, duplicate, id_start, id_end, flag_start, flag_end, Theta_start, omega_start, lambda_start, mu_theta, sigma_theta, mu_omega, sigma_omega, a, b, iter, burn);
end_time <- proc.time();
time <- end_time - start_time;



# ========================================================================================
# ========================================================================================
# Summarize result
# ========================================================================================
# ========================================================================================
print(paste0("Runtime = ", round(time[3], 1), "s"));
# print(paste0("Acceptance rate (omega): ", round(Y$accept_omega, 3)));
# print(paste0("Acceptance rate (theta): ", round(Y$accept_theta, 3)));
# print(paste0("Acceptance rate (lambda): ", round(Y$accept_lambda, 3)));
print(paste0("Estimated parameters (omega) = ", paste0(round(colMeans(Y$omega[burn:iter,]), 3), collapse = "", sep = ", ")));
print(paste0("Estimated interpretable parameters (pi) = ", paste0(round(omega2proportion(colMeans(Y$omega[burn:iter,])), 3), collapse = "", sep = ", ")));
print(paste0("Estimated parameters (theta) = ", paste0(round(colMeans(Y$theta[burn:iter,]), 3), collapse = "", sep = ", ")));
print("Estimated interpretable parameters (Phi) =");
print(round(Theta2interaction(array2matrix_r(colMeans(Y$theta[burn:iter,]), Q)), 3));
print(paste0("Estimated parameters (lambda) = ", paste0(round(mean(Y$lambda[burn:iter]), 3), collapse = "", sep = ", ")));



# ========================================================================================
# ========================================================================================
# Plot mark interaction functions
# ========================================================================================
# ========================================================================================
omega <- Y$omega[burn:iter,];
theta <- Y$theta[burn:iter,];
lambda <- Y$lambda[burn:iter];
d <- seq(0, 0.5, by = 0.01);
prob_mean <- matrix(NA, nrow = length(d), ncol = Q*Q);
count <- 1;
for (dd in d) {
  temp <- matrix(NA, nrow = iter - burn + 1, ncol = Q*Q)
  for (ii in 1:(iter - burn + 1)) {
    temp[ii,] <- c(t(Theta2interaction_2(omega[ii,], array2matrix_r(theta[ii,], Q), lambda[ii], dd)));
  }
  prob_mean[count,] <- colMeans(temp);
  count <- count + 1;
}

par(mfcol = c(Q, Q));
par(pty = "s", mar=c(2.5, 2.5, 1, 1));
count <- 1;
for (q in 1:Q) {
  for (qq in 1:Q) {
    plot(d, prob_mean[, count], type = "l", xlim = c(min(d), max(d)), ylim = c(0, 1), xlab = "", ylab = "", main = paste0("MIF ", q, "|", qq), cex.lab = 1.5, cex.axis = 1.5, cex = 1.5, lwd = 1.5);
    count <- count + 1;
  }
}


