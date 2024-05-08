library(foreach)
library(doParallel)
library(doSNOW)
library(futile.logger)

# loading data ------------------------------------------------------------

jazz_raw = read.csv('../data/jazz/out.arenas-jazz',sep='\t',header = F)[-1,]
jazz = data.frame(table(table(unlist(jazz_raw))))
names(jazz) = c('x','Freq')
jazz[,1] = as.numeric(jazz[,1])

ip_raw = read.csv('../data/as-caida20071105/out.as-caida20071105', sep=' ', header=F)[-(1:2),]
ip = data.frame(table(table(unlist(ip_raw[,1:2]))))
names(ip) = c('x','Freq')
ip[,1] = as.numeric(ip[,1])


erd_raw = read.csv('../data/pajek-erdos/out.pajek-erdos', sep='\t', header=F)[-1,]
erd = data.frame(table(table(unlist(erd_raw[,1:2]))))
names(erd) = c('x','Freq')
erd[,1] = as.numeric(erd[,1])

pro_raw = read.csv('../data/reactome/out.reactome', sep=' ', header=F)[-1,]
pro = data.frame(table(table(unlist(pro_raw[,1:2]))))
names(pro) = c('x','Freq')
pro[,1] = as.numeric(pro[,1])






library(igraph)

sim_ba_g = sample_pa(1e5, power=1, m=1, directed=F)
ba = data.frame(table(degree(sim_ba_g)))
names(ba) = c('x','Freq')
ba[,1] = as.numeric(ba[,1])

sim_ua_g = sample_pa(1e4, power=0, m=1, directed=F)
ua = data.frame(table(degree(sim_ua_g)))
names(ua) = c('x','Freq')
ua[,1] = as.numeric(ua[,1])


data.list = list(jazz=jazz,ip=ip,ba=ba)
# fitting model -----------------------------------------------------------


mcmc.out.list = list()

init = list(alpha=1, shape=0.5, scale=5, threshold=1.2)
prior.params = list(phi = c(1,1), alpha=5,shape=10, scale=0.001)
init.cov = list(alpha=0.1, shapescale=matrix(c(0.0003,-0.0003,
                                               -0.0003,1),ncol=2), threshold=0.1)

#multitail -N 1 -m 1 ./* -w --no-mark-change -d
{
# cl = makePSOCKcluster(length(data.list))
# registerDoSNOW(cl)
# N=5e4
# results <- foreach(i=1:length(data.list)) %dopar% {
#   source('functions.R')
#   source('new_functions.R')
#   library(networkdata)
#   library(mvtnorm)
#   out = pli.mcmc(N,data.list[[i]],init,prior.params,init.cov, cov.period = 1e7, type=1, show=F, logging = F, log_id = names(data.list)[[i]])
#   saveRDS(out, file=paste0('../data/mcmc.outputs/pli-',names(data.list)[[i]], '.rds'))
#   return(out)
# }
# stopCluster(cl)
# 
}
# -------------------------------------------------------------------------
source('functions.R')
source('new_functions.R')
library(networkdata)
library(mvtnorm)
N=1e4
out = pli.mcmc(N,ua,init,prior.params,init.cov, cov.period = 1e7, type=1, show=T, logging = F, log_id = 'ua')
saveRDS(out, file=paste0('../data/mcmc.outputs/pli-','ua', '.rds'))
marrangeGrob(list(plot.pli.mcmc(ua, out$states,burn.in=1e3, thin.by=5, type=1)),nrow=1,ncol=1)

out = readRDS('../data/mcmc.outputs/pli-ua.rds')




