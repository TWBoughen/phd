# Packages ----------------------------------------------------------------
require(VGAM)
require(gsl)
require(mvtnorm)

# Distribution Functions --------------------------------------------------
g = function(x,a){
  return(x^-(1+a) / gsl::zeta(1+a))
}
G = function(x,a){
  return(1-gsl::hzeta(1+a, x+1)/gsl::zeta(1+a))
}
h = function(x,xi,sig){
  return(
    VGAM::pgpd(x,scale=sig, shape=xi) - VGAM::pgpd(x-1,scale=sig, shape=xi)
  )
}
H = function(x,xi,sig){
  return(VGAM::pgpd(x,scale=sig,shape=xi))
}
dmix = Vectorize(function(x,v,a,xi,sig,log=T,type=1,phi=NULL){
  # type in {1,2,3}
  # 1: constrained power law
  # 2: free parameter phi (use this for empirical phi)
  if(type==2 & is.null(phi)){
    message('phi must be provided')
    return(1)
  }
  if(type==2){
    A = log(1-phi)-log(G(v,a))
    B = log(phi)
  }else if(type==1){
    gam = log(g(v+1,a))-log(h(1,xi,sig+xi*v))-log(G(v+1,a))
    phi = exp(gam)/(1+exp(gam))
    A = log(1-phi)-log(G(v,a))
    B = log(phi)
  }else return(1)
  if(x<=v){
    out = log(g(x,a)) + A
  }
  else{
    out = log(h(x-v,xi,sig+xi*v)) + B
  }
  
  if(!log){
    print(phi)
    return(exp(out))
  }
  return(out)
}, 'x')
# Priors ------------------------------------------------------------------
v_prior = function(v, upper=NULL, dampen = 1){
  if(is.null(upper)){
    return(1)
  }else{
    return(dampen*dbinom(v, upper, 0.5, log=T))
  }
}
a_prior = function(a){
  return(dgamma(a, 1, 1e-1,log=T))
}
xi_prior = function(xi){
  return(dnorm(xi, 0, 50,log=T))
}
sig_prior = function(sig){
  return(dgamma(sig, 1, 1e-1,log=T))
}
phi_prior = function(phi){
  return(dbeta(phi, 1, 10, log=T))
}
# Posterior ---------------------------------------------------------------
# data in degree,frequency format
mix_posterior = function(vals,dat, mode='constrained', nms = NULL){
  # mode: constrained,free, or empirical
  if(!is.null(nms)){
    vals = as.list(vals)
    names(vals) = nms
    vals = as.data.frame(vals)
  }
  prior_sum = v_prior(vals$v, upper=max(dat[,1])) +
              a_prior(vals$a) +
              xi_prior(vals$xi) +
              sig_prior(vals$sig)
  if(mode == 'free'){
    prior_sum = prior_sum + phi_prior(vals$phi)
    likelihood = sum(dat[,2]*dmix(dat[,1], vals$v, vals$a, vals$xi, vals$sig, type=2, phi=vals$phi))
  }
  else if(mode=='empirical'){
    likelihood = sum(dat[,2]*dmix(dat[,1], vals$v, vals$a, vals$xi, vals$sig, type=2, phi=sum(dat[dat[,1]>v,2]) / sum(dat[,2])))
  }
  else{
    likelihood = sum(dat[,2]*dmix(dat[,1], vals$v, vals$a, vals$xi, vals$sig))
  }
  if(is.nan(likelihood+prior_sum)){
    return(-Inf)
  }
  return(likelihood+prior_sum)
}
# Steps -------------------------------------------------------------------
threshold_gibbs_step = function(vals, dat, mode='constrained', jumpsize=1e4){
  probs = c()
  df = vals
  vs_to_sample = unique(dat[,1])[abs(unique(dat[,1])-vals$v)<=jumpsize]
  # vs_to_sample = unique(dat[,1])
  for(v in vs_to_sample){
    temp = vals
    temp$v = v
    df = rbind(df, temp)
  }
  df = df[-1,]
  probs = 1e-1*apply(df, 1, FUN=mix_posterior,dat=dat,nms = names(vals),mode=mode)
  # probs=  unique(dat[,1])[abs(unique(dat[,1])-vals$v)<=jumpsize]
  # probs=  unique(dat[,1])
  # probs= probs/probs
  vals$v = sample(vs_to_sample, 1, prob=exp(probs - max(probs))) 
  # plot(vs_to_sample,exp(probs - max(probs)))
  return(vals)
}
threshold_metropolis_step = function(vals, dat, mode='constrained', jumpsize=5){
  temp = vals
  temp$v = rbinom(1, max(dat[,1]), vals$v/max(dat[,1]))
  log_prob = min(0, mix_posterior(temp, dat, mode) - mix_posterior(vals, dat, mode)
                 +dbinom(vals$v, max(dat[,1]),temp$v/max(dat[,1]),log=T)-
                   dbinom(temp$v, max(dat[,1]),vals$v/max(dat[,1]),log=T))
  if(log(runif(1))<log_prob){
    return(temp)
  }else{
    return(vals)
  }
}
alpha_metropolis_step = function(vals, dat, mode='constrained',sd=0.05){
  temp = vals
  a_star = rnorm(1, vals$a, sd)
  temp$a = a_star
  if(a_star<=0){
    return(vals)
  }
  log_prob = min(0, mix_posterior(temp, dat, mode) - mix_posterior(vals, dat, mode))
  # print(log_prob)
  if(log(runif(1))<log_prob){
    return(temp)
  }else{
    return(vals)
  }
  
}
xi_sig_metropolis_step = function(vals, dat, mode='constrained', corr_mat=diag(0.1,2)){
  temp=vals
  xi_star = rnorm(1,vals$xi, sd=0.0)
  sig_star = rnorm(1, vals$sig, sd=0.1)
  comb = rmvnorm(1, mean=c(vals$xi, vals$sig), sigma=corr_mat)
  temp$xi = comb[1]
  temp$sig = comb[2]
  if(temp$sig<=0){
    return(vals)
  }
  log_prob = min(0, mix_posterior(temp, dat, mode) - mix_posterior(vals, dat, mode))
  if(is.nan(log_prob)){
    return(vals)
  }
  if(log(runif(1))<log_prob){
    return(temp)
  }else{
    return(vals)
  }
}
acc_rej = function(temp,vals, dat, mode='constrained'){
  log_prob = min(0, mix_posterior(temp, dat, mode) - mix_posterior(vals, dat, mode))
  if(is.nan(log_prob)){
    return(vals)
  }
  if(log(runif(1))<log_prob){
    return(temp)
  }else{
    return(vals)
  }
}






df = rbind(vals, vals+1)
mix_posterior(c(1,2,3,4), dat, nms = names(vals))
apply(df, 1, FUN = mix_posterior, dat=dat, nms = names(vals))
threshold_gibbs_step(vals, dat)
# -------------------------------------------------------------------------
a_ls = c()
v_ls = c()
vals=data.frame(v = 40, a=1, xi=1, sig=40)
vals$v = 40
vals$a = 2
vals$xi = 0.5
vals$sig=10

init = vals



# dat=dat_list[[2]]
mcmc_wrap2 <- function(init, dat, n, mode='constrained', plot=F) {
  vals= init
  df = as.data.frame(vals)
  # ps = c()
  if(plot){
    dev.new()
    par(mfrow=c(1,1))
  }
 
  for(i in 1:n){
    vals = threshold_metropolis_step(vals, dat,jumpsize=15,mode=mode)
    if(i>1e6){
      vals = alpha_metropolis_step(vals, dat,sd = 2.38*sd(tail(df$a, 100)))
      vals = xi_sig_metropolis_step(vals, dat, corr_mat = 2.38^2* cov(tail(df[,c('xi','sig')],100)))
    }else{
      vals = alpha_metropolis_step(vals, dat,mode=mode)
      vals = xi_sig_metropolis_step(vals, dat,mode=mode)
    }
    
    df = rbind(df, vals)
    message(i,'/v=',vals$v,
            '/alpha=',round(vals$a,2),
            '/xi=',round(vals$xi, 2), 
            '/sig=', round(vals$sig,2))
   # ps = c(ps, mix_posterior(vals, dat))
    if(plot){
      plot(tail(1/df$a,30),tail(df$xi,30), type='l', xlim=c(0,max(1/df$a)), ylim=c(-1,max(df$xi)), col='red')
      points(1/init$a,init$xi, col='blue', pch=5)
      points(1/df$a, df$xi, col=alpha('black', 0.05),pch=20)
      abline(a=0,b=1, lty=2)
    }
  }
  return(list(dat = dat, init=init, res = df))
}




par(mfrow=c(1,1))

hist(df$v)
plot(df$a, type='l')
plot(df$xi, type='l')
plot(df$sig, type='l')
plot(df$xi, 1/df$a, type='l')

plot(1/df$a, df$xi, type='p')
abline(a=0, b=1)

x = 1:200
df2 = df[-(1:200),]
v = df2$v[length(df2$a)]
a= df2$a[length(df2$a)]
xi=df2$xi[length(df2$a)]
sig=df2$sig[length(df2$a)]
y = dmix(x,v,a,xi ,sig , log=F)

plot(dat$x,1-(cumsum(dat$Freq)/sum(dat$Freq)), log='xy')
lines(x,1-cumsum(y))

# -------------------------------------------------------------------------
xi = seq(0,10,length.out=1e3)
alpha = seq(0,5,length.out=1e3)


join_prior = function(xi,a){
  # xi=xi_a[1]
  # a=x_a[2]
  return(xi_prior(xi)*a_prior(a))
}
z = outer(xi,alpha,FUN=join_prior)




# -------------------------------------------------------------------------
library(scales)
plot(dat$x, dat$Freq/sum(dat$Freq), log='xy')
abline(v=df$v, col=scales::alpha('blue', 0.005) )
lines(dat$x, dmix(dat$x,9,1.313, 0.29,11.747,log=F))
vals$v = 10
vals$a = 1
vals$xi = 1
vals$sig=1
df2 = df[-(1:2e3),]


plot(1/df2$a, df2$xi, xlim=c(0,2),ylim=c(0,2),pch=20, col=scales::alpha('grey',0.05),
     xlab='Implied tail heaviness', ylab='Tail heaviness')
abline(a=0,b=1,lty=2)
# contour(MASS::kde2d(1/df2$a,df2$xi,n=1000,lims = c(range(1/df2$a)+c(-1,1), range(df2$xi)+c(-1,1))), add=T, col='blue',
        # nlevels=7)


# -------------------------------------------------------------------------

df1 <-
  crandep::cran_dependencies |>
  filter(type == "imports", reverse == TRUE) |>
  count(from, name = "x") |>
  count(x, name = "count")
names(df1) = c('x','Freq')
dat=as.data.frame(df1)


# -------------------------------------------------------------------------
sim_pa_wrapper <- function(seed, n, m, zero.appeal) {
  set.seed(seed)
  igraph::sample_pa(n = n, m = m, zero.appeal = zero.appeal) |>
    igraph::degree(mode = "all") |>
    tibble::tibble(x = _) |>
    count(x, name = "count")
}
data.seed <- 4839L ####
n <- 3589 ####
m <- 1L ####
zero.appeal <- 0.1 ####
## simulate
df0 <- sim_pa_wrapper(data.seed, n, m, zero.appeal) |> filter(x > 0L)
dat = as.data.frame(df0)

# -------------------------------------------------------------------------


vals=data.frame(v = 40, a=1, xi=1, sig=40)
vals$v = 40
vals$a = 2
vals$xi = 0.5
vals$sig=10

set.seed(123L)
mcmc_wrap2(vals, as.data.frame(df1), 1e4, plot=T)














