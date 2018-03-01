# ***README***
# The following script is used to generate the marked spatial point pattern data analyzed 
# in the submitted manuscript titled "A Bayesian Mark Interaction Model for Analysis of 
# Tumor Pathology Images".

# After loading the necessary libraries and functions, please only run the code of one 
# block each time. The outputs are x, y, and z, which denote the x, y coordinates, and 
# marks of each point. Marks should be labelled as natural numbers 1, 2, 3...
# ***END***

# Load libraries
require(spatstat);
require(Rcpp);

# Load functions
source('functions.R');
Rcpp::sourceCpp('functions_x2.cpp');



# ========================================================================================
# ========================================================================================
# Generate simulated data
# ========================================================================================
# ========================================================================================
# Load settings
Q <- 2;
homogenous <- TRUE; # If TRUE, then from the homogeneous Poisson process; 
                    # otherwise from the log Gaussian Cox process
iter <- 10000; # Iterations for Monte Carlo simulation
seed <- 1;
lambda_true <- 60;                                  # True value of lambda
c_true <- 0.05;                                      # True value of c
Phi_true <- matrix(c(0.1, 0.9, 0.9, 0.1), nrow = 2); # True value of Phi (high inhibition)
pi_true <- c(0.5, 0.5);                              # True value of pi
Theta_true <- interaction2Theta(Phi_true);           # Obtain the true value of Theta
omega_true <- proportion2omega(pi_true);             # Obtain the true value of omega
# Generate points from a certain process
if (homogenous) {
  temp <- rpoispp(2000, win = owin(c(0, 1),c(0, 1)));
} else {
  temp <- rLGCP("exp", function(x, y){6 + abs(x - 0.3) + abs(y - 0.3)}, var = 1, 
                scale = 1, win = owin(c(0, 1), c(0, 1)));
}
x <- temp$x;
y <- temp$y;
# Simulate marks
z <- z_generator(x, y, Q, Theta_true, lambda_true, c_true, omega_true, iter, seed);



# ========================================================================================
# ========================================================================================
# Load amacrine data, as shown in Figure 1(a)
# ========================================================================================
# ========================================================================================
# Load data
data("amacrine"); # Q = 2, n = 294
data <- amacrine;
# Rescale data to the unit square
R <- max(max(data$x) - min(data$x), max(data$y) - min(data$y)); 
x <- (data$x - min(data$x))/R;
y <- (data$y - min(data$y))/R;
z <- as.integer(data$marks);



# ========================================================================================
# ========================================================================================
# Load betacells data, as shown in Figure 2(a)
# ========================================================================================
# ========================================================================================
# Load data
data("betacells"); # Q = 2, n = 135
data <- betacells;
# Rescale data to the unit square
R <- max(max(data$x) - min(data$x), max(data$y) - min(data$y)); 
x <- (data$x - min(data$x))/R;
y <- (data$y - min(data$y))/R;
z <- as.integer(data$marks[, 1]);



# ========================================================================================
# ========================================================================================
# Load lansing data, as shown in Figure 3(a)
# ========================================================================================
# ========================================================================================
# Load data
data("lansing"); # Q = 6, n = 2251
data <- lansing;
# Rescale data to the unit square
R <- max(max(data$x) - min(data$x), max(data$y) - min(data$y)); 
x <- (data$x - min(data$x))/R;
y <- (data$y - min(data$y))/R;
z <- as.integer(data$marks);



# ========================================================================================
# ========================================================================================
# Load NLST data, as shown in Figure 4(a)
# ========================================================================================
# ========================================================================================
data <- read.table("real_data_fig_4a.txt", header = TRUE); # Q = 3, n = 7558
x <- data$x;
y <- data$y;
z <- data$z;



# ========================================================================================
# ========================================================================================
# Load NLST data, as shown in Figure 4(b)
# ========================================================================================
# ========================================================================================
data <- read.table("real_data_fig_4b.txt", header = TRUE); # Q = 3, n = 8835
x <- data$x;
y <- data$y;
z <- data$z;


