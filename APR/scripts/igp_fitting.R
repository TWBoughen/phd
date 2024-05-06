library(foreach)
library(doParallel)
library(doSNOW)
library(futile.logger)

# loading data ------------------------------------------------------------

jazz_raw = read.csv('../data/jazz/out.arenas-jazz',sep='\t',header = F)[-1,]
jazz = data.frame(table(table(unlist(jazz_raw))))
names(jazz) = c('x','Freq')
jazz[,1] = as.numeric(jazz[,1])

data.list = list(jazz=jazz,jazz=jazz)
# fitting model -----------------------------------------------------------


mcmc.out.list = list()

init = list(alpha=1, shape=1, scale=5, threshold=10)
prior.params = list(phi = c(1,1), alpha=5,shape=10, scale=0.001)
init.cov = list(alpha=0.1, shapescale=matrix(c(0.003,-0.0003,
                                               -0.0003,1),ncol=2), threshold=0.1)

#multitail -N 1 -m 1 ./* -w --no-mark-change -d

cl = makePSOCKcluster(length(data.list))
registerDoSNOW(cl)
N=5e4
results <- foreach(i=1:length(data.list)) %dopar% {
  source('functions.R')
  source('new_functions.R')
  library(networkdata)
  library(mvtnorm)
  out = pli.mcmc(N,data.list[[i]],init,prior.params,init.cov, cov.period = 1e7, type=1, show=F, logging = F, log_id = names(data.list)[[i]])
  saveRDS(out, file=paste0('../data/mcmc.outputs/pli-',names(data.list)[[i]], '.rds'))
  return(out)
}
stopCluster(cl)


# -------------------------------------------------------------------------
source('functions.R')
source('new_functions.R')
library(networkdata)
library(mvtnorm)
N=5e4
out = pli.mcmc(N,jazz,init,prior.params,init.cov, cov.period = 1e7, type=1, show=F, logging = F, log_id = 'jazz')
saveRDS(out, file=paste0('../data/mcmc.outputs/pli-','jazz', '.rds'))


marrangeGrob(list(plot.pli.mcmc(jazz, out$states,burn.in=1e4, thin.by=20)),nrow=1,ncol=1)
