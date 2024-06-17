library(scales)



# -------------------------------------------------------------------------
#Misc functions

dhms <- function(t){
  paste(t %/% (60*60*24) 
        ,paste(formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0")
               ,formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0")
               ,formatC(t %% 60, width = 2, format = "d", flag = "0")
               ,sep = ":"
        )
  )
}

auto.mfrow <- function(nplots, setup=TRUE) {
  
  if(setup) {
    if(nplots <= 3) par(mfrow=c(1, nplots))
    else if(nplots <= 4)  par(mfrow=c(2,2))
    else if(nplots <= 6)  par(mfrow=c(2,3))
    else if(nplots <= 9)  par(mfrow=c(3,3))
    else if(nplots <= 12) par(mfrow=c(3,4))
    else if(nplots <= 16) par(mfrow=c(4,4))
    else if(nplots <= 20) par(mfrow=c(4,5))
    else if(nplots <= 25) par(mfrow=c(5,5))
    else if(nplots <= 30) par(mfrow=c(5,6))
    else if(nplots <= 36) par(mfrow=c(6,6))
    else if(nplots <= 42) par(mfrow=c(6,7))
    else if(nplots <= 49) par(mfrow=c(7,7))
    else if(nplots <= 56) par(mfrow=c(7,8))
    else if(nplots <= 64) par(mfrow=c(8,8))
    else {
      stop("Too many plots")
    }
  }
  else {
    nblankplots <- par("mfrow")[1] * par("mfrow")[2] - nplots
    if(nblankplots > 0)
      for(i in 1:nblankplots)
        plot_blank()
  }
}
# -------------------------------------------------------------------------

#These functions convert the data outputted by knoect.cc to the right format

konect_to_df = function(path){
  dat_raw = data.table::fread(path)
  dat = data.frame(table(table(unlist(dat_raw[,1:2]))))
  names(dat) = c('x','Freq')
  dat[,1] = as.numeric(dat[,1])
  return(dat)
}

konect_to_df_dir = function(dir,lower=2, return_names=F){
  dat_list = list()
  files = list.files(dir)
  if(length(lower)==1){
    lower=rep(lower, length(files))
  }
  for(i in 1:length(files) ){
    dat_list[[i]] = konect_to_df(paste0(dir,'/', files[i]))
    dat_list[[i]]  = dat_list[[i]][dat_list[[i]][,1]>=lower[i],]
  }
  if(return_names){
    return(list(dat = dat_list, names=files))
  }
  return(dat_list)
}


# -------------------------------------------------------------------------
#This function runs the mcmc on each data set in a directory in parallel.
dir_mcmc = function(dir, iter = 3e4, out.dir = '../data/mcmc.outputs', out.name='dir_out', lower=2){
  message('Start Time: ', Sys.time())
  dat_list = konect_to_df_dir(dir,lower=lower)
  cl = makeCluster(6)
  res = parSapply(cl, dat_list, FUN = mix_mcmc_wrap,iter=iter, burn.in=iter*0.1)
  stopCluster(cl)
  saveRDS(res, paste0(out.dir,'/',out.name,'.rds'))
  message('End Time: ', Sys.time())
  return(res)
}


#Density of the mixture
d_mix = function(x,v,phi,a,xi,sig){
  out = x*0
  if(sum(x<=v)>=1){
    out[x<=v] = (1-phi) * x[x<=v]^-(a+1) / sum((min(x):v)^-(a+1))
  }
  if(sum(x>v)>=1){
    p1 = pmax(x[x>v]*0, 1 + xi * (x[x>v]+1-v)/(sig+xi*v))
    p2 = pmax(x[x>v]*0, 1 + xi * (x[x>v]-v)/(sig+xi*v))
    pows2 = pows1 = x[x>v]*0 + 1
    pows1[p1>0] = -1/xi
    pows2[p2>0] = -1/xi
    out[x>v] =  phi*( p2^pows2 - p1^pows1 )
  }
  return(out)
}
#Cumulative density
#could improve this by putting in the actual formula instead of sim. 
p_mix = Vectorize(function(x,v,phi,a,xi,sig, lower=1){
  # return(sum(d_mix(lower:x,v,phi,a,xi,sig)))
  if(x<=v){
    return((1-phi)*sum((lower:x)^(-a-1))/sum((lower:v)^(-a-1)))
  }
  else if(x>v){
    p1 =  max(0,1+xi * (x+1-v)/(sig+xi*v))
    p2 = max(0,1+xi * (v+1-v)/(sig+xi*v))
    pows = c(1,1)
    pows[c(p1,p2)>0] = -1/xi
    return(1-phi*p1^pows[1])
  }
  else{
    return(0)
  }
}, vectorize.args = c('x', 'v', 'phi', 'a', 'xi', 'sig'))

p_mix_apply = function(pars, x){
  return(p_mix(x,pars[5], pars[1], pars[2], pars[3], pars[4], lower=min(x)))
}
#log-likelihood function
ll_mix = function(dat,v,phi,a,xi,sig){
  return(sum(dat[,2]*log(d_mix(dat[,1],v,phi,a,xi,sig))))
}
#log prior function
lpri = function(phi,a,xi,sig,params){
  return(dbeta(phi,params$phi[1], params$phi[2], log = T) +
           max(dgamma(a,1,rate=params$alpha,log=T), -1e-100) +
           dnorm(xi,0.5,sd=params$shape, log=T)+
           dgamma(sig, 1, rate=params$scale, log=T))
}
#log posterior
lpos = function(dat,v,phi,a,xi,sig,params){
  return(ll_mix(dat,v,phi,a,xi,sig) + lpri(phi,a,xi,sig,params))
}
#log posterior ratio
lpos_rat = function(dat, prop, cur, params){
  return(lpos(dat,prop$v,prop$phi,prop$a,prop$xi,prop$sig,params)-
           lpos(dat,cur$v,cur$phi,cur$a,cur$xi,cur$sig,params))
}
#functions for proposing phi (obsolote since using gibbs)
phi_prop = function(cur, var){
  n = cur*(1-cur)/var
  be_a = cur*n
  be_b = (1-cur)*n + 0.01
  return(rbeta(1,be_a,be_b))
}
phi_prop_p = function(cur,prop,var){
  n = cur*(1-cur)/var
  be_a = cur*n
  be_b = (1-cur)*n+ 0.01
  return(dbeta(prop,be_a,be_b,log=T))
}
#function to propose alpha
alpha_prop = function(cur,var){
  return(rnorm(1,cur,sqrt(var)))
}
#function to propose shape & scale
shapescale_prop = function(cur, var_mat){
  return(mvtnorm::rmvnorm(1,cur,sigma=var_mat))
}
#function to convert phi to a threshold
phi_to_u = function(dat, phi){
  lim = cumsum(dat[,2])/sum(dat[,2])
  if((1-phi)<min(lim[-1])){
    return(min(dat[,1]))
  }
  if((1-phi)>max(rev(lim)[-1])){
    return(max(dat[,1]))
  }
  return(head(dat[lim>= (1-phi),1],1))
}
#function for one step of mcmc
mix_mcmc_onestep = function(dat, init, params, covs, debug=F){
  
  acc.states = data.frame(phi=numeric(1),a=numeric(1),xi=numeric(1),sig=numeric(1),v=numeric(1))
  acc.states[1,] = c(init,phi_to_u(dat,init$phi))
  
  ####################################### phi steps
  a = tail(acc.states$a, 1)
  xi = tail(acc.states$xi, 1)
  sig = tail(acc.states$sig, 1)
  phis = 1 - cumsum(dat[,2])/sum(dat[,2])
  vs = unique(dat[,1])
  vs = vs[phis>0]
  phis = phis[phis>0]
  poss = c()
  for(i in 1:length(phis)){
    poss = c(poss, lpos(dat, vs[i], phis[i], a, xi, sig, params))
  }
  ws = poss-quantile(poss, probs=0.95)
  # ws[ws<0] = 0
  
  
  # print(sum(is.na(exp(ws^0.6))))
  if(runif(1)<0.95){
    ps = exp(ws)
    ps[ps==Inf] = 10^7
    phi_star = sample(phis[!is.nan(ps)], 1, prob=ps[!is.nan(ps)])
  }else{
    ps = exp(ws^0.001)
    ps[ps==Inf] = 10^7
    phi_star = sample(phis[!is.nan(ps)], 1, prob=ps[!is.nan(ps)]^0.001)
  }
  
  prop_state = tail(acc.states,1)
  prop_state$phi=phi_star
  prop_state$v = phi_to_u(dat, phi_star)
  acc.states = rbind(acc.states, prop_state)
  # mu = acc.states$phi[nrow(acc.states)]
  # be_a = ((1-mu)/covs$phi - 1/mu)*mu^2
  # be_b = be_a*(1/mu - 1)
  # phi_star = phi_prop(tail(acc.states$phi,1),covs$phi)
  # prop_state = tail(acc.states,1)
  # prop_state$phi = phi_star
  # prop_state$v = phi_to_u(dat, prop_state$phi)
  # if(debug) message('phi: ',prop_state)
  # 
  # #accepting/rejecting
  # logA = min(0,lpos_rat(dat, prop_state, tail(acc.states,1),params) + phi_prop_p(prop_state$phi,tail(acc.states$phi,1),covs$phi) - 
  #   phi_prop_p(tail(acc.states$phi,1),prop_state$phi,covs$phi))
  # if(log(runif(1))<logA){
  #   acc.states = rbind(acc.states, prop_state)
  # }
 
  ##########################(############## alpha steps
  alpha_star = rnorm(1, tail(acc.states$a, 1), sqrt(covs$alpha))
  prop_state = tail(acc.states, 1)
  prop_state$a = min(10^2, alpha_star)
  if(debug) message('alpha: ', prop_state)
  #accepting/rejecting
  logA = min(0,lpos_rat(dat, prop_state, tail(acc.states,1),params))
  
  if((!is.nan(logA)& log(runif(1))<logA | runif(1)<0.00) & prop_state$a>0){
    acc.states = rbind(acc.states, prop_state)
  }
  
  ######################################## shape & scale steps
  shapescale_star = mvtnorm::rmvnorm(1, unlist(tail(acc.states[,c('xi', 'sig')], 1)), sigma=covs$shapescale)
  prop_state = tail(acc.states, 1)
  prop_state$xi = shapescale_star[1]
  prop_state$sig = shapescale_star[2]
  if(debug) message('shapescale: ',prop_state)
  #accepting/rejecting
  if((prop_state$sig>0| runif(1)<0.0) & prop_state$xi>(-prop_state$sig/max(dat[,1])) ){
    logA = min(0,lpos_rat(dat, prop_state, tail(acc.states,1),params))
    if(!is.nan(logA) & log(runif(1))<logA ){
      acc.states = rbind(acc.states, prop_state)
    }
  }
  return(tail(acc.states,1))
}
#mcmc function
mix_mcmc = function(iter, dat, init, params, covs,update_period = 100, debug=F,plotting = F){
  start_time = Sys.time()
  last_time = Sys.time()
  out = mix_mcmc_onestep(dat, init, params, covs)
  for(i in 2:iter){
    out = rbind(out, mix_mcmc_onestep(dat, tail(out,1)[,-ncol(out)], params, covs, debug=debug))
    if(i>10*update_period & i%%update_period==0){
      temp = covs$shapescale
      covs$shapescale =0.1 + 2.38^2 * cov(tail(out[,c('xi', 'sig')], update_period*4)[!duplicated(tail(out[,c('xi', 'sig')], update_period*4)),])/2
      if(sum(is.nan(covs$shapescale))>0){
        covs$shapescale=temp
      }
      covs$alpha = 0.1+ 2.38^2 * var(tail(out$a,update_period*4))/2
    }
    if(i %% update_period == 0){
     
      
      recent = tail(out,1)
      message('---------------------------------')
      message('Iteration: ',i)
      message('Threshold: ',round(recent$v, 3))
      message('Phi: ',round(recent$phi, 7))
      message('Alpha:', round(recent$a, 3))
      message('Shape:', round(recent$xi, 3))
      message('Scale:', round(recent$sig, 3))
      message('Time since last update: ',dhms(as.numeric(Sys.time())- as.numeric(last_time)))
      message('Time since start: ', dhms(as.numeric(Sys.time())- as.numeric(start_time)))
      last_time = Sys.time()
      
      
      if(plotting){
        par(mfrow=c(2,3))
        plot(out$phi,type='l', log='y')
        
        plot(out$a, type='l')
        lines(cummean(out$a), col='red')
        
        plot(out$sig, type='l')
        lines(cummean(out$sig), col='red')
        # if(cov_update) abline(v=i, lty=2, col='green')
        plot(out$xi, type='l',main=sum(is.na(out$xi)))
        lines(cummean(out$xi), col='red')
        # if(cov_update) abline(v=i, lty=2, col='green')
        plot(out$sig, out$xi)
        
        
        legend('topleft', legend = diag(covs$shapescale))
        
        # plot(dat[,1],dat[,2]/sum(dat[,2]), pch=20, log='xy', col='grey')
        # lines(dat[,1], d_mix(dat[,1],phi_to_u(dat, mean(tail(out$phi,5*update_period))),
        #                        mean(tail(out$phi,5*update_period)),
        #                        mean(tail(out$a,5*update_period)),
        #                        mean(tail(out$xi,5*update_period)),
        #                        mean(tail(out$sig,5*update_period))), col='blue')
        
        plot(dat[,1],1-cumsum(dat[,2])/sum(dat[,2]), log='xy', col='blue',type='l')
        # abline(v = phi_to_u(dat, mean(tail(out$phi,5*update_period))), col='red')
        # lines(dat[,1], 1-p_mix(dat[,1],phi_to_u(dat, mean(tail(out$phi,5*update_period))),
        #                        mean(tail(out$phi,5*update_period)),
        #                        mean(tail(out$a,5*update_period)),
        #                        mean(tail(out$xi,5*update_period)),
        #                        mean(tail(out$sig,5*update_period)), lower=min(dat[,1])), col='blue')
        for(i in 1:update_period){
          iter_dat = tail(out,i)[1,]
          abline(v= phi_to_u(dat, iter_dat$phi), col=alpha('darkgreen',0.05))
          lines(dat[,1], 1-p_mix(dat[,1], phi_to_u(dat, iter_dat$phi),
                                 iter_dat$phi,
                                 iter_dat$a,
                                 iter_dat$xi,
                                 iter_dat$sig, lower=min(dat[,1])), col=alpha('red', 0.05))
        }
        lines(dat[,1],1-cumsum(dat[,2])/sum(dat[,2]), col='blue')
        # lines(dat[,1], 1-p_mix(dat[,1],phi_to_u(dat, quantile(mean(tail(out$phi,5*update_period)),0.95)),
        #                        quantile(tail(out$phi,5*update_period),0.95),
        #       quantile(mean(tail(out$a,5*update_period)),0.95),
        #       quantile(mean(tail(out$xi,5*update_period)),0.95),
        #       quantile(mean(tail(out$sig,5*update_period)),0.95), lower=min(dat[,1])), col='blue',lty=2)
        # lines(dat[,1], 1-p_mix(dat[,1],phi_to_u(dat, quantile(mean(tail(out$phi,5*update_period)),0.05)),
        #                        quantile(tail(out$phi,5*update_period),0.05),
        #                        quantile(mean(tail(out$a,5*update_period)),0.05),
        #                        quantile(mean(tail(out$xi,5*update_period)),0.05),
        #                        quantile(mean(tail(out$sig,5*update_period)),0.05), lower=min(dat[,1])), col='blue',lty=2)
        
      
      }
    }
  }
  return(out)
}
library(dplyr)
library(coda)
# -------------------------------------------------------------------------

#mcmc wrapper
mix_mcmc_wrap <- function(dat,iter=10000,  burn.in=1000, thin.by=10, update_period = 1000, debug=F, plotting=F) {
  source('igp_functions.R')
  init = data.frame(phi=0.05, a=2, xi=0.05, sig=10)
  params = list(
    phi=c(1,1),
    alpha=0.0001,
    shape=10,
    scale=0.01
  )
  covs = list(
    phi=0.001,
    alpha=0.0001,
    shapescale = t(matrix(ncol=2,
          c(0.02,-sqrt(0.02*10)/max(dat[,1]),
            -sqrt(0.02*10)/max(dat[,1]),10)
      ))
  )
  res = mix_mcmc(iter, dat, init, params, covs, update_period = update_period, debug=debug, plotting=plotting)
  #thinning and burning
  if(burn.in>nrow(res)){
    res_thinburn = res
  }else{
    res_burn = res[-(1:burn.in), ]
    res_thinburn = res_burn[seq(1,nrow(res_burn), by=thin.by),]
  }
  return(list(dat = dat,res_thinned = res_thinburn, res=res))
}
#
mcmc_plot = function(dat, res_thinburn,ylim=c(1e-5, 1), xlim=c(1, max(dat[,1])), mar=c(1,1,1,1),xaxt=NULL, yaxt=NULL){
  x = unique(dat[,1])
  cmfs = apply(res_thinburn, 1, p_mix_apply, x=x)
  cmfs_mean = apply(cmfs, 1, mean)
  cmfs_95 = apply(cmfs, 1, quantile, prob=0.975)
  cmfs_05 = apply(cmfs, 1, quantile, prob=0.025)
  #plotting
  par(mar=mar)
  Fk = cumsum(dat[,2])/sum(dat[,2])
  plot(dat[,1], 1- c(0,Fk[-length(Fk)]), log='xy', pch=20,
       xlab='Degree', ylab='Survival', xlim=xlim, ylim=ylim, xaxt=xaxt, yaxt=yaxt,las=1)
  abline(v=res_thinburn$v, col = alpha('blue', alpha=0.004))
  polygon(c(x, rev(x)), c(1- c(0,cmfs_95[-length(cmfs_95)]), rev(1- c(0,cmfs_05[-length(cmfs_05)]))), col=alpha('gray', 0.5), lty=0)
  lines(x, 1-c(0,cmfs_mean[-length(cmfs_mean)]), col='red')
}


# -------------------------------------------------------------------------
# set.seed()
# burn.in = 1e3
# thin.by = 5
# ba_res = mix_mcmc_wrap(ba, iter=5e4, update_period = 100, burn.in=1e4, thin.by = 5)
# ba_res_n1 = mix_mcmc_wrap(ba[-1,], iter=2e4, update_period = 100, burn.in=1e4, thin.by = 5)
# 
# ua_res = mix_mcmc_wrap(ua, iter=3e4, update_period = 100)
# 
# 
# ip_res = mix_mcmc_wrap(ip, update_period = 1e2)
# 
# plot(ip[,1], 1-cumsum(ip[,2])/sum(ip[,2]), log='xy', pch=20)
# 
# erd_res = mix_mcmc_wrap(erd, iter =1e4,update_period = 100)
# 
# jazz_res = mix_mcmc_wrap(jazz, 1e4, burn.in=2e3,update_period = 1e2, debug = F)
# pro_res = mix_mcmc_wrap(pro,iter = 3e4, update_period = 1e2)
# 
# # -------------------------------------------------------------------------
# dat=ba
# x = unique(ba[,1])
# res_thinburn = ba_res
# 
# cmfs = apply(res_thinburn, 1, p_mix_apply, x=x)
# cmfs_mean = apply(cmfs, 1, mean)
# cmfs_95 = apply(cmfs, 1, quantile, prob=0.95)
# cmfs_05 = apply(cmfs, 1, quantile, prob=0.05)
# #plotting
# par(mfrow=c(1,1))
# plot(dat[,1], 1- cumsum(dat[,2])/sum(dat[,2]), log='xy', pch=20)
# abline(v=res_thinburn$v, col = alpha('blue', alpha=0.004))
# polygon(c(x, rev(x)), c(1- cmfs_95, rev(1- cmfs_05)), col=alpha('gray', 0.5), lty=0)
# lines(x, 1-cmfs_mean, col='red')
# 
# 
# 
# 
# # -------------------------------------------------------------------------
# 
# set.seed(Sys.time())
# start_adj = diag(1)
# start_g = graph_from_adjacency_matrix(start_adj, mode='undirected')
# net = sample_pa(1e5, power=1.3, m=1, directed = F, zero.appeal = 0, start.graph = start_g)
# ba = data.frame(table(degree(net)))
# names(ba) = c('x','Freq')
# ba[,1] = as.numeric(ba[,1])
# plot(ba, log='xy')
# 
# plot(unique(ba[,1]), 1-cumsum(ba[,2])/sum(ba[,2]), log='xy', pch=20)
# 
# 
# # -------------------------------------------------------------------------
# 
# harvard_res=mix_mcmc_wrap(harvard, iter=3e4, burn.in=1e4, thin.by=5)
# ip_res = mix_mcmc_wrap(ip)
# 
# 
