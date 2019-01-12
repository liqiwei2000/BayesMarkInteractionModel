proportion2omega = function(pi) {
  return (-log(pi) + log(pi)[length(pi)] + 1);
}

omega2proportion = function(omega) {
  return (exp(-omega)/sum(exp(-omega)));
}

interaction2Theta = function(interaction) {
  temp <- -log(interaction);
  Q <- dim(temp)[1];
  temp[, Q] <- temp[, Q] - temp[Q, Q] + 1
  for (q in 1:(Q - 1)) {
    temp[, q] <- temp[, q] - temp[Q, q] + temp[q, Q];
  }
  return (temp)
}

Theta2interaction = function(Theta) {
  return (t(exp(-Theta)/(colSums(exp(-Theta)))));
}

Theta2interaction_2 = function(Theta, lambda, d) {
  return (t(exp(-Theta*exp(-lambda*d))/(colSums(exp(-Theta*exp(-lambda*d))))));
}

z2interaction = function(x, y, z, Q, c) {
  build <- dist_list(x, y, c);
  z_edge <- cbind(z[build$edge[, 1]], z[build$edge[, 2]]);
  interaction <- matrix(0, nrow = Q, ncol = Q);
  for (q in 1:Q) {
    temp <- rep(0, Q);
    for (i in 1:dim(z_edge)[1]) {
      if (z_edge[i, 1] == q) {
        temp[z_edge[i, 2]] <- temp[z_edge[i, 2]] + 1;
      } else if (z_edge[i, 2] == q) {
        temp[z_edge[i, 1]] <- temp[z_edge[i, 1]] + 1;
      }
    }
    interaction[, q] <- temp/sum(temp)
  }
  return (interaction);

}

z_generator = function(x, y, Q, Theta, lambda, c, omega, iter, seed) {
  set.seed(seed);
  
  # Generate imaging data
  build <- dist_list(x, y, c);
  edge <- build$edge;
  distance <- build$distance;
  duplicate <- build$duplicate;
  flag_start <- build$flag_start;
  flag_end <- build$flag_end;
  
  # Generate cell types
  n <- length(x);
  z_s <- sample(1:Q, n, replace = TRUE, prob = exp(-omega)/sum(exp(-omega)));
  z <- z_s;
  # energy_s <- logenergy(z_s, edge, distance, duplicate, Thetat, lambdat);
  # energy <- rep(NA, iter);
  count <- 10;
  for (i in 1:iter) {
    if (i/iter*100 == count) {
      print(paste0("Simulating the data ------ ", count, "% has been done"));
      count <- count + 10;
    }
    z <- model(z, edge, distance, 1, n, flag_start, flag_end, Theta, omega, lambda);
    # energy[i] <- logenergy(z, edge, distance, duplicate, Thetat, lambdat);
  }
  return (z); 
}

array2matrix_r = function(theta, Q) {
  Theta <- matrix(0, nrow = Q, ncol = Q);
  count <- 1;
  for (q in 1:Q) {
    for (qq in q:Q) {
      Theta[q, qq] <- theta[count];
      Theta[qq, q] <- theta[count];
      count <- count + 1;
    }
  }
  return (Theta);
}

color = function(i, j) {
  col <- NA;
  if (i == j) {
    col = i;
  } else {
    if ((i = 1 && j == 2) || (i = 2 && j == 1)) {
      col = 8;
    }
    if ((i = 1 && j == 3) || (i = 3 && j == 1)) {
      col = 4;
    }
    if ((i = 3 && j == 2) || (i = 2 && j == 3)) {
      col = 6;
    }
  }
  return (col);
}