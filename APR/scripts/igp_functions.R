library(scales)

dhms <- function(t){
  paste(t %/% (60*60*24) 
        ,paste(formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0")
               ,formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0")
               ,formatC(t %% 60, width = 2, format = "d", flag = "0")
               ,sep = ":"
        )
  )
}



d_mix = function(x,v,phi,a,xi,sig){
  out = x*0
  if(sum(x<=v)>=1){
    out[x<=v] = (1-phi) * x[x<=v]^-(a+1) / sum((1:v)^-(a+1))
  }
  if(sum(x>v)>=1){
    p1 = pmax(x[x>v]*0, 1 + xi * (x[x>v]+1-v)/(sig+xi*v))
    p2 = pmax(x[x>v]*0, 1 + xi * (x[x>v]-v)/(sig+xi*v))
    pows2 = pows1 = x[x>v]*0 + 1
    pows1[p1>0] = -1/xi
    pows2[p2>0] = -1/xi
    out[x>v] =  phi*( - p1^pows1 + p2^pows2)
  }
  return(out)
}
p_mix = Vectorize(function(x,v,phi,a,xi,sig){
  return(sum(d_mix(1:x,v,phi,a,xi,sig)))
}, vectorize.args = c('x', 'v', 'phi', 'a', 'xi', 'sig'))
p_mix_apply = function(pars, x){
  return(p_mix(x,pars[5], pars[1], pars[2], pars[3], pars[4]))
}
ll_mix = function(dat,v,phi,a,xi,sig){
  return(sum(dat[,2]*log(d_mix(dat[,1],v,phi,a,xi,sig))))
}
lpri = function(phi,a,xi,sig,params){
  return(dbeta(phi,params$phi[1], params$phi[2], log = T) +
           max(dgamma(a,1,rate=params$alpha,log=T), -1e-100) +
           dnorm(xi,0.5,sd=params$shape, log=T)+
           dgamma(sig, 1, rate=params$scale, log=T))
}
lpos = function(dat,v,phi,a,xi,sig,params){
  return(ll_mix(dat,v,phi,a,xi,sig) + lpri(phi,a,xi,sig,params))
}
lpos_rat = function(dat, prop, cur, params){
  return(lpos(dat,prop$v,prop$phi,prop$a,prop$xi,prop$sig,params)-
           lpos(dat,cur$v,cur$phi,cur$a,cur$xi,cur$sig,params))
}
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
alpha_prop = function(cur,var){
  return(rnorm(1,cur,sqrt(var)))
}
shapescale_prop = function(cur, var_mat){
  return(mvtnorm::rmvnorm(1,cur,sigma=var_mat))
}
phi_to_u = function(dat, phi){
  lim = cumsum(dat[,2])/sum(dat[,2])
  if((1-phi)<min(lim[-1])){
    return(min(dat[,1]))
  }
  if((1-phi)>max(rev(lim)[-1])){
    return(max(dat[,1]))
  }
  return(tail(dat[lim< (1-phi),1],1))
}
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
  ws = poss-quantile(poss, probs=0.05)
  ws[ws<0] = 0
  
  
  
  phi_star = sample(phis, 1, prob=exp(ws^0.2))
  v_star = sample(vs, 1, prob=ws)
  prop_state = tail(acc.states,1)
  prop_state$phi=phi_star
  prop_state$v = v_star
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
  prop_state$a = alpha_star
  if(debug) message('alpha: ', prop_state)
  #accepting/rejecting
  logA = min(0,lpos_rat(dat, prop_state, tail(acc.states,1),params))
  if((log(runif(1))<logA | runif(1)<0.05) & prop_state$a>0){
    acc.states = rbind(acc.states, prop_state)
  }
  
  ######################################## shape & scale steps
  shapescale_star = mvtnorm::rmvnorm(1, unlist(tail(acc.states[,c('xi', 'sig')], 1)), sigma=covs$shapescale)
  prop_state = tail(acc.states, 1)
  prop_state$xi = shapescale_star[1]
  prop_state$sig = shapescale_star[2]
  if(debug) message('shapescale: ',prop_state)
  #accepting/rejecting
  if((prop_state$sig>0| runif(1)<0.05) & prop_state$xi>(-prop_state$sig/prop_state$v) ){
    logA = min(0,lpos_rat(dat, prop_state, tail(acc.states,1),params))
    if(log(runif(1))<logA ){
      acc.states = rbind(acc.states, prop_state)
    }
  }
  return(acc.states)
}

mix_mcmc = function(iter, dat, init, params, covs,update_period = 100, debug=F){
  start_time = Sys.time()
  last_time = Sys.time()
  out = mix_mcmc_onestep(dat, init, params, covs)
  for(i in 2:iter){
    out = rbind(out, mix_mcmc_onestep(dat, tail(out,1)[,-ncol(out)], params, covs, debug=debug))
    if(i %% update_period == 0){
      covs$shapescale = 2.38^2 * cov(out[,c('xi', 'sig')][!duplicated(out[,c('xi', 'sig')]),])/2
      covs$alpha = 2.38^2 * var(out$a)
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
      
      plotting = F
      if(plotting){
        par(mfrow=c(2,3))
        hist(tail(out$v,5*update_period), type='l',breaks = max(dat[,1]), border=F)
        plot(out$a, type='l', log='y')
        lines(cummean(out$a), col='red')
        
        plot(out$sig, type='l')
        lines(cummean(out$sig), col='red')
        plot(out$xi, type='l',main=sum(is.na(out$xi)))
        lines(cummean(out$xi), col='red')
        
        plot(dat[,1],dat[,2]/sum(dat[,2]), pch=20, log='xy')
        lines(dat[,1], d_mix(dat[,1],phi_to_u(dat, mean(tail(out$phi,5*update_period))),
                               mean(tail(out$phi,5*update_period)),
                               mean(tail(out$a,5*update_period)),
                               mean(tail(out$xi,5*update_period)),
                               mean(tail(out$sig,5*update_period))), col='blue')
        
        plot(dat[,1],1-cumsum(dat[,2])/sum(dat[,2]), pch=20, log='xy')
        abline(v = phi_to_u(dat, mean(tail(out$phi,5*update_period))), col='red')
        lines(dat[,1], 1-p_mix(dat[,1],phi_to_u(dat, mean(tail(out$phi,5*update_period))),
                               mean(tail(out$phi,5*update_period)),
                               mean(tail(out$a,5*update_period)),
                               mean(tail(out$xi,5*update_period)),
                               mean(tail(out$sig,5*update_period))), col='blue')
      
      }
    }
  }
  return(out)
}
library(dplyr)
library(coda)
# -------------------------------------------------------------------------


mix_mcmc_wrap <- function(dat,iter=10000,  burn.in=1000, thin.by=10, update_period = 1000, debug=F) {
  source('igp_functions.R')
  init = data.frame(phi=0.1, a=3, xi=0.1, sig=1)
  params = list(
    phi=c(1,1),
    alpha=0.0001,
    shape=10,
    scale=0.01
  )
  covs = list(
    phi=0.001,
    alpha=0.01,
    shapescale = t(matrix(ncol=2,
          c(1,-0.4,
            -0.4,1)
      ))
  )
  res = mix_mcmc(iter, dat, init, params, covs, update_period = update_period, debug=debug)
  #thinning and burning
  res_burn = res[-(1:burn.in), ]
  res_thinburn = res_burn[seq(1,nrow(res_burn), by=thin.by),]
  return(list(dat = dat,res = res_thinburn))
}

mcmc_plot = function(dat, res_thinburn){
  x = unique(dat[,1])
  cmfs = apply(res_thinburn, 1, p_mix_apply, x=x)
  cmfs_mean = apply(cmfs, 1, mean)
  cmfs_95 = apply(cmfs, 1, quantile, prob=0.95)
  cmfs_05 = apply(cmfs, 1, quantile, prob=0.05)
  #plotting
  par(mfrow=c(1,1))
  plot(dat[,1], 1- cumsum(dat[,2])/sum(dat[,2]), log='xy', pch=20)
  abline(v=res_thinburn$v, col = alpha('blue', alpha=0.004))
  polygon(c(x, rev(x)), c(1- cmfs_95, rev(1- cmfs_05)), col=alpha('gray', 0.5), lty=0)
  lines(x, 1-cmfs_mean, col='red')
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
