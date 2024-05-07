d_mix = function(x,v,phi,a,xi,sig){
  out = x*0
  out[x<=v] = (1-phi) * x[x<=v]^-(a+1) / sum((1:v)^-(a+1))
  p1 = pmax(x[x>v]*0, 1 + xi * (x[x>v]+1-v)/(sig+xi*v))
  p2 = pmax(x[x>v]*0, 1 + xi * (x[x>v]-v)/(sig+xi*v))
  pows2 = pows1 = x[x>v]*0 + 1
  pows1[p1>0] = -1/xi
  pows2[p2>0] = -1/xi
  out[x>v] =  phi*( - p1^pows1 + p2^pows2)
  return(out)
}
ll_mix = function(dat,v,phi,a,xi,sig){
  return(sum(dat[,2]*log(d_mix(dat[,1],v,phi,a,xi,sig))))
}
lpri = function(phi,a,xi,sig,params){
  return(dbeta(phi,params$phi[1], params$phi[2], log = T) +
           dnorm(a,0,sd=params$alpha,log=T) +
           dnorm(xi,0,sd=params$shape, log=T)+
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
  be_a = cur^2 * ((1-cur)/var - 1/cur)
  be_b = be_a*(1/cur - 1)
  return(rbeta(1,be_a,be_b))
}
phi_prop_p = function(cur,prop,var){
  be_a = cur^2 * ((1-cur)/var - 1/cur)
  be_b = be_a*(1/cur - 1)
  return(dbeta(prop,be_a,be_b,log=T))
}
alpha_prop = function(cur,var){
  return(rnorm(1,cur,sqrt(var)))
}
shapescale_prop = function(cur, var_mat){
  return(mvtnorm::rmvnorm(1,cur,sigma=var_mat))
}
phi_to_u = function(dat, phi){
  return(tail(dat[cumsum(dat[,2])/sum(dat[,2]) < (1-phi),1],1))
}
mix_mcmc = function(dat,iter, init, params, covs){
  acc.states = data.frame(phi=numeric(1),a=numeric(1),xi=numeric(1),sig=numeric(1),v=numeric(1))
  acc.states[1,] = c(init,phi_to_u(dat,init$phi))
  #######################################phi steps
  mu = acc.states$phi[nrow(acc.states)]
  be_a = ((1-mu)/covs$phi^2 - 1/mu)*mu^2
  be_b = be_a*(1/mu - 1)
  phi_star = phi_prop(tail(acc.states$phi,1),covs$phi)
  prop_state = tail(acc.states,1)
  prop_state$phi = phi_star
  prop_state$v = phi_to_u(dat, prop_state$phi)
  #accepting/rejecting
  logA = min(0,lpos_rat(dat, prop_state, tail(acc.states,1),params) + phi_prop_p(prop_state$phi,tail(acc.states$phi,1),covs$phi) - 
    phi_prop_p(tail(acc.states$phi,1),prop_state$phi,covs$phi))
  if(log(runif(1))<logA){
    acc.states = rbind(acc.states, prop_state)
  }
  
}







