# BayesMarkInteractionModel
This repository is used for access the performance of the Bayesian mark interaction model, which was proposed in the submitted manuscript titled "A Bayesian mark interaction model for analysis of tumor pathology images." Before running the code, please install two R packages: spatstat and Rcpp.

Firstly, please run "data_loader.R" to generate simulated data or load real data. This is a required step.

Secondly (optional), run "mcf_analyzer.R" to plot the data and the corresponding mark connection functions powered by R package spatstat.

Lastly, run "model_fitting.R" to fit the proposed model and obtain the preliminary results, such as runtime, estimated model parameters, and mark interaction function plots. 

*Note 1, the munuscript can be downloaded from https://arxiv.org/abs/1802.08308.

*Note 2, the notations in the code and data are followed the notations in the manuscript.

*Note 3, the results obtained by running the code in this repository may not be exact the same as the results reported in the manuscript, because we reported the results pooled from multiple MCMC chains in the manuscript.

*Note 4, for large datasets and large numbers of c, it takes longer time. Please first try small datasets such as "amacrine" and "betacells", or run the MCMC algorithm with small numbers of iterations (e.g. 5,000) and small value of c (e.g. 0.05).
