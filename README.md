# BayesMarkInteractionModel
This repository is used for accessing the performance of the Bayesian mark interaction model, which was proposed in the paper titled "A Bayesian mark interaction model for analysis of tumor pathology images" published on the Annals of Applied Statistics. Before running the code, please install two R packages: spatstat and Rcpp.

Firstly, run "data_loader.R" to generate simulated data or load real data. This is a required step.

Secondly (optional), run "mcf_analyzer.R" to plot the data and the corresponding mark connection functions powered by R package spatstat.

Lastly, run "model_fitting.R" to fit the proposed model and obtain the preliminary results, such as runtime, estimated model parameters, and mark interaction function plots. 

*Note 1, the munuscript can be downloaded from https://projecteuclid.org/journals/annals-of-applied-statistics/volume-13/issue-3/A-Bayesian-mark-interaction-model-for-analysis-of-tumor-pathology/10.1214/19-AOAS1254.full.

*Note 2, the notations in the code and data follow the notations in the manuscript.

*Note 3, the results obtained by running the code in this repository may not be exactly the same as the results reported in the manuscript, because we reported the results by pooling multiple MCMC chains in the manuscript.

*Note 4, for large dataset and large number of c, it takes longer time. Please first try small dataset such as "amacrine" and "betacells", or run the MCMC algorithm with small number of iterations (e.g. 5,000) and small value of c (e.g. 0.05).
