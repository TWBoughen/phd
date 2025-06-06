library(foreach)
library(doParallel)
library(doSNOW)
# library(futile.logger)

# loading data ------------------------------------------------------------

jazz_raw = read.csv('../data/data_for_mcmc/out.jazz',sep='\t',header = F)[-1,]
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
library(wdnet)
set.seed(Sys.time())
#uniform attachment
ctrl_ua = rpa_control_preference(ftype='customized', pref = '1')
init_net = list(
  edgelist = matrix(c(1, 2), nrow = 1),
  edgeweight = 1,
  directed =FALSE
)
net = wdnet_to_igraph(rpanet(1e5, initial.network = init_net, ctrl_ua))

ua = data.frame(table(degree(net)))
names(ua) = c('x','Freq')
ua[,1] = as.numeric(ua[,1])
#Ba model
ctrl_ba = rpa_control_preference(ftype='customized', pref = 's', spref='outs', tpref='ins')
init_net = list(
  edgelist = matrix(c(1, 1,
                      1,1), nrow = 2),
  edgeweight = c(1,1),
  directed =FALSE
)


# -------------------------------------------------------------------------
library(igraph)
net = wdnet_to_igraph(rpanet(1e5, initial.network = init_net, ctrl_ba))

ba = data.frame(table(degree(net)))
names(ba) = c('x','Freq')
ba[,1] = as.numeric(ba[,1])
x = unique(ba[,1])
plot(unique(ba[,1]), 1-cumsum(ba[,2])/sum(ba[,2]), log='xy', pch=20)
lines(x, 2*x^-3)

plot(ua[,1], ua[,2]/sum(ua[,2]), log='xy', pch=20)

for(i in 1:20){
  net = wdnet_to_igraph(rpanet(1e5, initial.network = init_net, ctrl_ua))
  # net = sample_pa(1e5)
  ua = data.frame(table(degree(net)))
  names(ua) = c('x','Freq')
  ua[,1] = as.numeric(ua[,1])
  points(ua[,1], ua[,2]/sum(ua[,2]), pch=20)
}

e=exp(1)
x = 1:max(ua[,1])

lines(x,x^-3, col='red', lwd=2)




dat=ua
x =  1:max(dat[,1])
y = 1-cumsum(dat[,2])/sum(dat[,2])
plot(dat[,1],c(1,y[-length(y)]) , log='xy', pch=20, type='b')


# -------------------------------------------------------------------------




# -----------------------------------------------------------------------




res = dir_mcmc('../data/data_for_mcmc', iter=3e4)

res=readRDS('../data/mcmc.outputs/dir_trunc.out.rds')
n =ncol(res)
mx_dat = 0
min_prob = 1
dev.new()
for(i in 1:n){
  d = res[,i]$dat
  mx_dat = max(mx_dat, max(d[,1]))
  min_prob = min(min_prob, d[nrow(d),2]/sum(d[,2]))
}
dev.new()
par(mfrow=c(1,1))
for(i in 1:n){
  message(i)
  mcmc_plot(res[,i]$dat, as.data.frame(res[,i]$res_thinned,
                                       ylim=c(min_prob, 1), xlim=c(1,mx_dat)))
}

dev.new(wdith=10, height=10, unit='in', res=72)

# -------------------------------------------------------------------------



dev.new()
dat_list=konect_to_df_dir('../data/data_for_mcmc')

G = barabasi.game(1e5, m=1, directed=F, zero.appeal=0.000000000000000001)
x = degree(G, mode='in')
x=x[x>0]
sim = as.data.frame(table(x))
sim[,1] = as.numeric(as.character(sim[,1]))
write.csv(as_data_frame(G,what='edges'), file='../data/data_for_mcmc/out.BAsim', row.names = F)


G = barabasi.game(1e5,power=0, m=1, directed=F, zero.appeal=0.000000000000000001)
write.csv(as_data_frame(G,what='edges'), file='../data/data_for_mcmc/out.UAsim', row.names = F)



res_ua = mix_mcmc_wrap(dat_list[[3]][-1,], iter=3e4, burn.in=3e6, update_period = 1e2, plotting=T)

dev.new()
mcmc_plot(res_meta$dat, res_meta$res[-(1:1e4),])

plot(as.mcmc(res6$res))





# fitting model -----------------------------------------------------------

# 
# mcmc.out.list = list()
# 
# init = list(alpha=1, shape=0.5, scale=5, threshold=1.2)
# prior.params = list(phi = c(1,1), alpha=5,shape=10, scale=0.001)
# init.cov = list(alpha=0.1, shapescale=matrix(c(0.0003,-0.0003,
#                                                -0.0003,1),ncol=2), threshold=0.1)
# 
# #multitail -N 1 -m 1 ./* -w --no-mark-change -d
# {
# # cl = makePSOCKcluster(length(data.list))
# # registerDoSNOW(cl)
# # N=5e4
# # results <- foreach(i=1:length(data.list)) %dopar% {
# #   source('functions.R')
# #   source('new_functions.R')
# #   library(networkdata)
# #   library(mvtnorm)
# #   out = pli.mcmc(N,data.list[[i]],init,prior.params,init.cov, cov.period = 1e7, type=1, show=F, logging = F, log_id = names(data.list)[[i]])
# #   saveRDS(out, file=paste0('../data/mcmc.outputs/pli-',names(data.list)[[i]], '.rds'))
# #   return(out)
# # }
# # stopCluster(cl)
# # 
# }
# # -------------------------------------------------------------------------
# source('functions.R')
# source('new_functions.R')
# library(networkdata)
# library(mvtnorm)
# N=1e4
# out = pli.mcmc(N,ua,init,prior.params,init.cov, cov.period = 1e7, type=1, show=T, logging = F, log_id = 'ua')
# saveRDS(out, file=paste0('../data/mcmc.outputs/pli-','ua', '.rds'))
# marrangeGrob(list(plot.pli.mcmc(ua, out$states,burn.in=1e3, thin.by=5, type=1)),nrow=1,ncol=1)
# 
# out = readRDS('../data/mcmc.outputs/pli-ua.rds')
# 
# # -------------------------------------------------------------------------
# 
# 
# gp_ll = function(y, u, shape, scale){
#   x = y[y>0]
#   return(-sum(dgpd(x,loc=0, scale=scale, shape=shape, log=T)))
#   
# }
# 
# 
# gp_ll_optim = function(shapescale, y, u){
#   return(gp_ll(y,u,shapescale[1], shapescale[2]+shapescale[1]*u))
# }
# 
# dat=as.numeric(table(unlist(ip_raw[,1:2])))
# vals = data.frame(shape=numeric(1), scale=numeric(1))
# vs = seq(1.5e4, 1.9e4, by=1e2)
# for(v in vs){
#   print(v)
#   vals=rbind(vals, optim(c(1,1), gp_ll_optim,y=dat, u=v, method='L-BFGS-B', lower=c(0, 0.01), upper=c(Inf,Inf ))$par)
# }
# vals=vals[-1,]
# 
# plot(vs,vals$shape, type='b')
# 
# 
# 
# library(evmix)
# 
# par(mfrow=c(1,1))
# mrl.plot()
# 


