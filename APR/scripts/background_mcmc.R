library(foreach)
library(doSNOW)
source('igp_functions.R')
res = dir_mcmc('../data/data_for_mcmc', iter=3e4, lower = 3)