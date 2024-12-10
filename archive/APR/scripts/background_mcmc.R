library(foreach)
library(doSNOW)
source('igp_functions.R')
res = dir_mcmc('../data/data_for_mcmc', iter=3e4, lower = c(2,5,2,3,1,1,2),
               out.name = 'dir_trunc.out')