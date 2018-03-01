#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double logenergy(arma::vec z, NumericVector zi, IntegerMatrix E, NumericVector d, LogicalVector dd, NumericMatrix Theta, NumericVector omega, double lambda);
arma::vec model(arma::vec z, IntegerMatrix E, NumericVector d, int id_start, int id_end, IntegerVector flag_start, IntegerVector flag_end, NumericMatrix Theta, NumericVector omega, double lambda);

// [[Rcpp::export]]
Rcpp::List model_estimator(arma::vec z, IntegerMatrix E, NumericVector d, LogicalVector dd, IntegerVector id_start, IntegerVector id_end, IntegerVector flag_start, IntegerVector flag_end, NumericMatrix Theta_s, NumericVector omega_s, double lambda_s, double mu, double sigma, double mu_omega, double sigma_omega, double a, double b, int iter, int burn) {
  int N = id_start.size();
  int n = z.n_rows;
  int Q = Theta_s.nrow();
  int M = 1;
  int i, ii, q, qq, qqq, qqqq, count, k;
  int count_2 = 10;
  double tau = 0.1;
  double phi = 1;
  double hastings = 0;
  double accept_theta = 0;
  double accept_lambda = 0;
  double accept_omega = 0;
  double lambda = lambda_s;
  double lambda_temp;
  // double sum_temp;
  arma::vec z_temp(n);
  NumericMatrix theta_store(iter, Q*(Q + 1)/2);
  NumericVector lambda_store(iter);
  NumericMatrix Theta(Q, Q);
  NumericMatrix Theta_temp(Q, Q);
  NumericVector omega(Q);
  NumericVector omega_temp(Q);
  NumericMatrix omega_store(iter, Q);
  NumericVector zi(Q);
  NumericVector zi_temp(Q);
  
  // Initialization
  for(q = 0; q < Q; q++)
  {
    zi(q) = 0;
    zi_temp(q) = 0;
    omega(q) = omega_s(q);
    for(qq = 0; qq < Q; qq++)
    {
      Theta(q, qq) = Theta_s(q, qq);
    }
  }
  for(ii = 0; ii < n; ii++)
  {
    for(q = 0; q < Q; q++) {
      if(z(ii) == q + 1)
      {
        zi(q) = zi(q) + 1;
        break;
      }
    }
  }

  // MCMC
  for(i = 0; i < iter; i++)
  {
    // Update omega
    for(q = 0; q < Q - 1; q++)
    {
      for(qq = 0; qq < Q; qq++)
      {
        omega_temp(qq) = omega(qq);
      }
      omega_temp(q) = rnorm(1, omega(q), tau)(0);
      for(ii = 0; ii < n; ii++)
      {
        z_temp(ii) = z(ii);
      }
      for(ii = 0; ii < M; ii++)
      {
        for(k = 0; k < N; k++)
        {
          z_temp = model(z_temp, E, d, id_start(k), id_end(k), flag_start, flag_end, Theta, omega_temp, lambda);
        }
      }
      for(qq = 0; qq < Q; qq++)
      {
        zi_temp(qq) = 0;
      }
      for(ii = 0; ii < n; ii++)
      {
        for(qq = 0; qq < Q; qq++) 
        {
          if(z_temp(ii) == qq + 1)
          {
            zi_temp(qq)++;
            break;
          }
        }
      }
      hastings = -logenergy(z_temp, zi_temp, E, d, dd, Theta, omega, lambda) + logenergy(z, zi, E, d, dd, Theta, omega, lambda) - logenergy(z, zi, E, d, dd, Theta, omega_temp, lambda) + logenergy(z_temp, zi_temp, E, d, dd, Theta, omega_temp, lambda);
      hastings = hastings - (omega_temp(q) - mu_omega)*(omega_temp(q) - mu_omega)/2/sigma_omega/sigma_omega + (omega(q) - mu_omega)*(omega(q) - mu_omega)/2/sigma_omega/sigma_omega;
      if (hastings >= log(double(rand()%10001)/10000))
      {
        omega(q) = omega_temp(q);
        if (i > burn) {
          accept_omega++;
        }
      }
    }

    // Update Theta
    // sum_temp = 0;
    for(q = 0; q < Q - 1; q++)
    {
      for(qq = q; qq < Q; qq++)
      {
        //if(q != Q - 1 || qq != Q - 1)
        {
          for(qqq = 0; qqq < Q; qqq++)
          {
            for (qqqq = 0; qqqq < Q; qqqq++)
            {
              Theta_temp(qqq, qqqq) = Theta(qqq, qqqq);
            }
          }
          Theta_temp(q, qq) = rnorm(1, Theta(q, qq), tau)(0);
          Theta_temp(qq, q) = Theta_temp(q, qq);
          // Theta_temp(Q - 1, Q - 1) = Theta(Q - 1, Q - 1) + (Theta(q, qq) - Theta_temp(q, qq));
          for(ii = 0; ii < n; ii++)
          {
            z_temp(ii) = z(ii);
          }
          for(ii = 0; ii < M; ii++)
          {
            for(k = 0; k < N; k++)
            {
              z_temp = model(z_temp, E, d, id_start(k), id_end(k), flag_start, flag_end, Theta_temp, omega, lambda);
            }
          }
          for(qqq = 0; qqq < Q; qqq++)
          {
            zi_temp(qqq) = 0;
          }
          for(ii = 0; ii < n; ii++)
          {
            for(qqq = 0; qqq < Q; qqq++) 
            {
              if(z_temp(ii) == qqq + 1)
              {
                zi_temp(qqq)++;
                break;
              }
            }
          }
          hastings = -logenergy(z_temp, zi_temp, E, d, dd, Theta, omega, lambda) + logenergy(z, zi, E, d, dd, Theta, omega, lambda) - logenergy(z, zi, E, d, dd, Theta_temp, omega, lambda) + logenergy(z_temp, zi_temp, E, d, dd, Theta_temp, omega, lambda);
          hastings = hastings - (Theta_temp(q, qq) - mu)*(Theta_temp(q, qq) - mu)/2/sigma/sigma + (Theta(q, qq) - mu)*(Theta(q, qq) - mu)/2/sigma/sigma;
          // hastings = hastings - (Theta_temp(Q - 1, Q - 1) - mu)*(Theta_temp(Q - 1,  Q - 1) - mu)/2/sigma/sigma + (Theta(Q - 1, Q - 1) - mu)*(Theta(Q - 1, Q - 1) - mu)/2/sigma/sigma;
          if (hastings >= log(double(rand()%10001)/10000))
          {
            Theta(q, qq) = Theta_temp(q, qq);
            Theta(qq, q) = Theta(q, qq);
            // Theta(Q - 1, Q - 1) = Theta_temp(Q - 1, Q - 1);
            if (i > burn) {
              accept_theta++;
            }
          }
        }
        //sum_temp = sum_temp + Theta(q, qq);
      }
    }
    //Theta(Q - 1, Q - 1) = -sum_temp;
    
    // Updata lambda
    do {
      lambda_temp = rgamma(1, lambda*lambda/phi, phi/lambda)(0);
    } while (lambda_temp < 0.1);
    hastings = (a - 1)*(log(lambda_temp) - log(lambda)) - b*(lambda_temp - lambda);
    for(ii = 0; ii < n; ii++)
    {
      z_temp(ii) = z(ii);
    }
    for (ii = 0; ii < M; ii++)
    {
      for(k = 0; k < N; k++)
      {
        z_temp = model(z_temp, E, d, id_start(k), id_end(k), flag_start, flag_end, Theta, omega, lambda_temp);
      }
    }
    for(q = 0; q < Q; q++)
    {
      zi_temp(q) = 0;
    }
    for(ii = 0; ii < n; ii++)
    {
      for(q = 0; q < Q; q++) 
      {
        if(z_temp(ii) == q + 1)
        {
          zi_temp(q)++;
          break;
        }
      }
    }
    hastings = hastings - logenergy(z_temp, zi_temp, E, d, dd, Theta, omega, lambda) + logenergy(z, zi, E, d, dd, Theta, omega, lambda) - logenergy(z, zi, E, d, dd, Theta, omega, lambda_temp) + logenergy(z_temp, zi_temp, E, d, dd, Theta, omega, lambda_temp);
    if (hastings >= log(double(rand()%10001)/10000))
    {
      //Rcout<<lambda_temp<<"\n";
      lambda = lambda_temp;
      if (i > burn) {
        accept_lambda++;
      }
    }

    // Monitor the process
    if (i*100/iter == count_2)
    {
      Rcout<<count_2<< "% has been done\n";
      count_2 = count_2 + 10;
    }
    count = 0;
    for(q = 0; q < Q; q++)
    {
      omega_store(i, q) = omega(q);
      for(qq = q; qq < Q; qq++)
      {
        theta_store(i, count) = Theta(q, qq);
        count++;
      }
    }
    lambda_store(i) = lambda;
  }
  accept_omega = accept_omega/(iter - burn)/(Q - 1);
  accept_theta = accept_theta/(iter - burn)/(Q*(Q + 1)/2 - 1);
  accept_lambda = accept_lambda/(iter - burn);
  return Rcpp::List::create(Rcpp::Named("omega") = omega_store, Rcpp::Named("Theta") = Theta, Rcpp::Named("theta") = theta_store, Rcpp::Named("lambda") = lambda_store, Rcpp::Named("accept_theta") = accept_theta, Rcpp::Named("accept_omega") = accept_omega, Rcpp::Named("accept_lambda") = accept_lambda);
}

// [[Rcpp::export]]
Rcpp::List dist_list(NumericVector x, NumericVector y, double c) {
  int n = x.size();
  int i,j;
  int count = 0;
  double temp;
  IntegerMatrix net_temp(n*n, 2);
  NumericVector dist_temp(n*n);
  LogicalVector dup_temp(n*n);
  IntegerVector flag_start(n);
  IntegerVector flag_end(n);
  for(i = 0; i < n; i++)
  {
    flag_start(i) = count + 1;
    for(j = 0; j < n; j ++)
    {
      if(i != j && std::abs(x(i) - x(j)) <= c && std::abs(y(i) - y(j)) <= c)
      {
        temp = (x(i) - x(j))*(x(i) - x(j)) + (y(i) - y(j))*(y(i) - y(j));
        if(temp <= c*c)
        {
          net_temp(count, 0) = i + 1;
          net_temp(count, 1) = j + 1;
          dist_temp(count) = sqrt(temp);
          if(i > j)
          {
            dup_temp(count) = TRUE;
          }
          else
          {
            dup_temp(count) = FALSE;
          }
          count++;
        }
      }
    }
    flag_end(i) = count;
  }
  IntegerMatrix net(count, 2);
  NumericVector dist(count);
  LogicalVector dup(count);
  for (i = 0; i < count; i++)
  {
    for (j = 0; j < 2; j++)
    {
      net(i, j) = net_temp(i, j);
    }
    dist(i) = dist_temp(i);
    dup(i) = dup_temp(i);
  }
  return Rcpp::List::create(Rcpp::Named("edge") = net, Rcpp::Named("distance") = dist, Rcpp::Named("duplicate") = dup, Rcpp::Named("flag_start") = flag_start, Rcpp::Named("flag_end") = flag_end);
}

// [[Rcpp::export]]
arma::vec model(arma::vec z, IntegerMatrix E, NumericVector d, int id_start, int id_end, IntegerVector flag_start, IntegerVector flag_end, NumericMatrix Theta, NumericVector omega, double lambda) {
  int Q = Theta.nrow();
  NumericVector prob_temp(Q);
  IntegerVector state(Q);
  int i, j, q;
  double temp = 0;
  for(q = 0; q < Q; q++)
  {
    state(q) = q;
  }
  for(i = id_start - 1; i < id_end; i++)
  {
    if(flag_end(i) >= flag_start(i))
    {
      for(q = 0; q < Q; q++)
      {
        prob_temp(q) = 0;
      }
      for(j = flag_start(i) - 1; j < flag_end(i); j++)
      {
        for(q = 0; q < Q; q++)
        {
          prob_temp(q) = prob_temp(q) - Theta(q, z(E(j, 1) - 1) - 1)*exp(-lambda*d(j));
        }
      }
      for(q = 0; q < Q; q++)
      {
        prob_temp(q) = exp(prob_temp(q) - omega(q));
      }
      temp = 0;
      for (q = 0; q < Q; q++)
      {
        temp = temp + prob_temp(q);
      }
      for (q = 0; q < Q; q++)
      {
        prob_temp(q) = prob_temp(q)/temp;
      }
      z(i) = RcppArmadillo::sample(state, 1, true, prob_temp)(0) + 1;
    }
    else
    {
      for(q = 0; q < Q; q++)
      {
        prob_temp(q) = exp(-omega(q));
      }
      temp = 0;
      for (q = 0; q < Q; q++)
      {
        temp = temp + prob_temp(q);
      }
      for (q = 0; q < Q; q++)
      {
        prob_temp(q) = prob_temp(q)/temp;
      }
      z(i) = RcppArmadillo::sample(state, 1, true, prob_temp)(0) + 1;
    }
  }
  return z;
}

// [[Rcpp::export]]
double logenergy(arma::vec z, NumericVector zi, IntegerMatrix E, NumericVector d, LogicalVector dd, NumericMatrix Theta, NumericVector omega, double lambda) {
  double energy = 0;
  int m = E.nrow();
  int Q = zi.size();
  int j, q;
  for(j = 0; j < m; j++)
  {
    if(!dd(j))
    {
      energy = energy + Theta(z(E(j, 0) - 1) - 1, z(E(j, 1) - 1) - 1)*exp(-lambda*d(j));
    }
  }
  for(q = 0; q < Q; q++)
  {
    energy = energy + zi(q)*omega(q);
  }
  return energy;
}