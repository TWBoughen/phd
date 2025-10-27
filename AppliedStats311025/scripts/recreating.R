sample_rnd = function(n, pref, gamma_g,quiet = T){
  degs = c(0)
  gammas = c(gamma_g(1))
  prefs = c(pref(degs[1], gammas[1]))
  for(i in 2:n){
    if(!quiet) message(i)
    selected = sample(1:(i-1), 1, prob = prefs)
    degs[selected] = degs[selected]+1
    prefs[selected] = pref(degs[selected], gammas[selected])
    degs = c(degs, 0)
    gammas = c(gammas, gamma_g(1))
    prefs  = c(prefs, pref(degs[i], gammas[i]))
  }
  return(data.frame(degree=degs, gamma = gammas, pref = prefs))
}


# -------------------------------------------------------------------------
library(ggplot2)

k0 = 20
a = 1
eps = 1


pref = function(x, gamma){
  return(ifelse(x<=k0, x^a+gamma[,2], k0^a + eps + gamma[,1]*(x-k0)))
}
gamma_g = function(n){
  return(rgamma(n*2, 1, 10))
}

n=1e4
G = sample_rnd(n, pref, gamma_g, quiet = F)
degs =  G$degree
par(mfrow=c(2,2))
plot(twbfn::deg_dist(degs[degs>0]), log='xy',pch=20,ylim = c(1/n,1), xlim=c(1,n+1))
plot(twbfn::deg_surv(degs[degs>0]), log='xy',pch=20,ylim = c(1/n,1), xlim=c(1,n+1))
# plot(density(G$gamma))
plot(G$degree+1, G$pref,pch=20, col = scales::alpha('blue', 0.3), log='xy', ylim=c(0.01,1e3))

# mcmc network ------------------------------------------------------------

source('funcs.R')

fits = readRDS('modelfits.rds')


# -------------------------------------------------------------------------

j = 1
mcmcdat = fits[[j]]$smps

N = 1e4
par_rows = sample(1:nrow(mcmcdat),N,replace=T)
par_mat = mcmcdat[par_rows, ]
par_vec = apply(par_mat,2,median)
par_mat_mn = t(matrix(rep(par_vec, N), ncol=N))

mcmc_network = function(N, par_mat){
  pref = function(x,pars){
    pars = as.numeric(pars)
    return(ifelse(x<=pars[3], x^pars[1] + pars[2], pars[3]^pars[1] + pars[2] + pars[4]*(x-pars[3])))
  }
  degs = c(0)
  prefs = c(pref(degs[1], par_mat[1,]))
  for(i in 2:N){
    if(!quiet) message(i)
    selected = sample(1:(i-1), 1, prob = prefs)
    degs[selected] = degs[selected]+1
    prefs[selected] = pref(degs[selected], par_mat[selected,])
    degs = c(degs, 0)
    prefs  = c(prefs, pref(degs[i], par_mat[i,]))
  }
  return(list(degs = degs, prefs=prefs))
}

degs = mcmc_network(N, par_mat)$degs
degs_mn = mcmc_network(N,par_mat_mn)$degs
g = function(x){
  return(x^par_vec[1] + par_vec[2])
}
dat_list = readRDS('dat_list.rds')
par(mfrow = c(2,2))



x = 1:1e3
lambda = find_lambda2(g, par_vec[4], par_vec[3])
y = S(x, g, lambda, par_vec[3], par_vec[4])


plot(twbfn::deg_surv(counts_to_degs(dat_list[[j]])+1), log='xy',pch=20,ylim=c(1/sum(dat_list[[j]][,2]),1))
points(twbfn::deg_surv(degs),pch=20, col='orange')
lines(x, y)


plot(twbfn::deg_surv(counts_to_degs(fits[[j]]$dat)+1), log='xy',pch=20,ylim=c(1/sum(dat_list[[j]][,2]),1))
lines(x, y/S(min(fits[[j]]$dat[,1])-1,g, lambda, par_vec[3], par_vec[4]))
points(twbfn::deg_surv(degs[degs>=min(fits[[j]]$dat[,1])]),pch=20, col='orange')
lines(fits[[j]]$dat[,1], fits[[j]]$surv$est, lty=2, col='red')

plot(twbfn::deg_surv(counts_to_degs(dat_list[[j]])+1), log='xy',pch=20,ylim=c(1/sum(dat_list[[j]][,2]),1))
lines(x, y)
points(twbfn::deg_surv(degs_mn),pch=20, col='orange')

plot(twbfn::deg_surv(counts_to_degs(fits[[j]]$dat)+1), log='xy',pch=20,ylim=c(1/sum(dat_list[[j]][,2]),1))
lines(x, y/S(min(fits[[j]]$dat[,1])-1,g, lambda, par_vec[3], par_vec[4]))
points(twbfn::deg_surv(degs_mn[degs_mn>=min(fits[[j]]$dat[,1])]),pch=20, col='orange')





















































































































