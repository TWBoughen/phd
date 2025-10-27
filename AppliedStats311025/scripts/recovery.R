f = Vectorize(function(x, g, lambda, k0, b){
  if(x>=k0){
    k0s = 0:(k0-1)
    phi_g = lambda*gamma((lambda+g(k0))/b)*prod(g(k0s)/(lambda+g(k0s)))/
      (b*gamma(lambda/b + 1)*gamma(g(k0)/b))
    return(phi_g*beta(g(k0)/b + x-k0, lambda/b + 1))
  }
  if(x>0){
    ks = 0:(x-1)
    return(lambda/(lambda+g(x))*prod(g(ks)/(lambda+g(ks))))
  }
  if(x==0){
    return(lambda/(lambda+g(x)))
  }
}, vectorize.args = 'x')
S = Vectorize(function(x, g, lambda, k0, b){
  if(x>=k0){
    k0s = 0:(k0-1)
    
    phi_g = lgamma((lambda+g(k0))/b)  + sum(log(g(k0s)/(lambda+g(k0s))))-
      (lgamma(lambda/b)+lgamma(g(k0)/b))
    
    return(exp(phi_g)*beta(g(k0)/b + x-k0+1, lambda/b))
  }
  if(x>0){
    ks = 0:(x)
    return(prod(g(ks)/(lambda+g(ks))))
  }
  if(x==0){
    return(g(x)/(lambda+g(x)))
  }
  
},vectorize.args = 'x')
sample_gpa_tree = function(n, g,b,k0,m=1){
  degs = c(0)
  prefs = c(g(0))
  for(i in 1:n){
    # message(i)
    selected = sample(1:length(degs), min(m, length(degs)), prob=prefs)
    degs = c(degs, 0)
    prefs = c(prefs, g(0))
    degs[selected] = degs[selected]+1
    prefs[selected] = ifelse(degs[selected]<k0, g(degs[selected]), g(k0) + b*(degs[selected]-k0))
  }
  return(degs)
}
rho = Vectorize(function(lambda, g, b, k0){
  n1 = 1:k0
  n2 = 0:(k0-1)
  out1 = sum(exp(cumsum(log(1-lambda/(lambda+g(n1-1))))))
  out2 = exp(sum(log(1-lambda/(lambda+g(n2)))))
  return(out1 + g(k0)*out2/(lambda-b))
}, vectorize.args = 'lambda')
rho_optim = function(lambda, g, b, k0){
  return(rho(lambda, g, b, k0)-1)
}
find_lambda2 = function(g, b, k0){
  return(uniroot(rho_optim, c(b,100), g=g,b=b,k0=k0, extendInt = 'yes')$root)
}
polylincurve = function(x, a,eps, b, k0){
  g01 = function(x){
    return(x^a+eps)
  }
  l01 = find_lambda2(g01,b, k0)
  y01 = S(x, g01,l01, k0,b)
  return(geom_line(aes(x=x+1, y=y01)))
}

# -------------------------------------------------------------------------


p = 3
plot(which(accepted==0),par_mat[accepted==0,p], col='red', cex=0.4,pch=20)
lines(which(accepted==1), par_mat[accepted==1,p])

plot(twbfn::deg_surv(degs), log='xy', pch=20, xlim=c(1,1e4))
abline(v=post_pars[,3], col = scales::alpha('blue', 0.1), lwd=2)
points(twbfn::deg_surv(degs),pch=20)
# -------------------------------------------------------------------------
plot(density(post_pars$a))
abline(v=a)
# -------------------------------------------------------------------------
a = 1.3
eps = 1
b = 0.5
k0 = 30
lambda  = find_polylambda(c(a,eps,k0,b))


g = function(x){
  return(x^a + eps)
}
#######samplepolylin(N, a, eps, b, k0)
degs = samplepolylin(1e4, a, eps, b, k0)
# plot(twbfn::deg_surv(degs), log='xy', pch=20)
dat = twbfn::deg_count(degs)

library(mev)
u = k0
fit2 = fit.gpd(degs, threshold = u, MCMC=T,method='zhang')
# x = k0:max(degs)
# y = pgp(x, k0, fit2$estimate[1], fit2$estimate[2], lower.tail = F)
# lines(x, y*sum(degs>=k0)/length(degs))


x = 0:max(degs)
y_est  = pgp(u:max(degs), u, fit2$estimate[1], fit2$estimate[2], lower.tail = F)
ggplot() + polylincurve(x, a, eps, b, k0) +
  geom_line(aes(x=u:max(degs), y=y_est*sum(degs>=u)/length(degs)), color='red', linetype='dashed',lwd=1)+
  geom_point(data=twbfn::deg_surv(degs+1), aes(x=degree, y=surv))+
  scale_x_log10() + scale_y_log10()

c(k0^a + eps+b,b)/lambda
fit2$post.mean

# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
a = 1.3
eps = 1
b = 0.5
k0 = 30
find_lambda2 = function(g, b, k0){
  return(uniroot(rho_optim, c(b,100), g=g,b=b,k0=k0, extendInt = 'yes')$root)
}
polylin = function(a,eps){
  return(function(x){
    return(x^a+eps)
  })
}
library(mvtnorm)
degs = samplepolylin(1e3, a, eps, b, k0)
dat = twbfn::deg_count(degs)
#pars: a, eps, k0, b
llh_x = Vectorize(function(x, pars,lambda){
  a = pars[1]
  eps = pars[2]
  k0 = pars[3]
  b = pars[4]
  if(x==0){
    return(log(lambda) - log(lambda + eps))
  }else if(x<k0){
    return(log(lambda) - log(lambda + eps + x^a) +
             sum(log(1 - lambda/(lambda + eps + (0:(x-1))^a))))
  }else{
    return(sum(log(1 - lambda/(lambda + eps + (0:(k0-1))^a))) +
             lbeta(x-k0 + (k0^a + eps)/b,1 + lambda/b) - lbeta((k0^a + eps)/b,lambda/b))
  }
},vectorize.args='x')
llh = function(dat, pars, lambda){
  if(any(pars<0)){
    return(-Inf)
  }
  return(sum(dat[,2]*llh_x(dat[,1], pars, lambda)))
}
prior = function(pars,xmax = 1e3){
  return(
    dgamma(pars[1],1,0.01,log=T)+
      dgamma(pars[2],1,0.01,log=T)+
      dunif(pars[3],1,xmax-1,log=T)+
      dgamma(pars[4],1,0.01,log=T)
  )
}
posterior = function(dat, pars, lambda){
  return(llh(dat, pars, lambda) + prior(pars,xmax = max(dat[,1])))
}
posterior_k0 = Vectorize(function(dat, k0, pars, lambda){
  if(k0<=min(dat[,1])){
    return(-Inf)
  }
  if(k0>=max(dat[,1])){
    return(-Inf)
  }
  return(posterior(dat, c(pars[1],pars[2],k0,pars[4]), lambda))
},vectorize.args = 'k0')



# -------------------------------------------------------------------------
library(coda)
a = 1.5
b = 1.5
eps=0.1
k0=25
out = recover_params(a, eps, k0, b, quiet=F,stop_at=1e3,network.size = 1e4,N=1e3)
plot(out$output$pars$a,type='l')
thin.by = 1
mcmcout = coda::as.mcmc(out$output$pars[seq(1,nrow(out$output$pars),by=thin.by),c(1,2,4)])
plot(mcmcout)
effectiveSize(mcmcout)

plot(out$output$pars$a,out$output$pars$eps)

# -------------------------------------------------------------------------

library(ggplot2)

source('funcs.R')


as = c(.5,1,1.5)
bs = c(.1,.5,1,1.5)
eps = c(0.1,0.5,1)
k0s = 20

pars = expand.grid(as,eps,k0s,bs)
recovery_list = list()
names(pars) = c('a','eps','k0','b')
for(i in 28:nrow(pars)){
  message(i,'/',nrow(pars))
  recovery_list[[i]] = recover_params(pars[i,1], pars[i,2], pars[i,4],pars[i,3],quiet=F,n.iter=2e4,network.size = 1e5)
  recovery_list[[i]]$output$pars$true_a = pars[i,1]
  recovery_list[[i]]$output$pars$true_eps = pars[i,2]
  recovery_list[[i]]$output$pars$true_k0 = pars[i,3]
  recovery_list[[i]]$output$pars$true_b = pars[i,4]
}


plot(twbfn::deg_surv(recovery_list[[20]]$degs), log='xy')

pars[28,]
plot(recovery_list[[28]]$mcmc$mcmc$pars[1:10,4],type='l')

recovery_list[[28]]$mcmc$mcmc$sig

saveRDS(recovery_list, 'recovery_dat.rds')

saveRDS(pars, 'recovery_pars.rds')

recovery_list  =readRDS('recovery_dat.rds')
recover_pars = readRDS('recovery_pars.rds')
thin.by = 5
full_pars = recovery_list[[1]]$mcmc$smps
for(i in 2:length(recovery_list)){
  full_pars = rbind(full_pars, recovery_list[[i]]$mcmc$smps)
}

library(latex2exp)
library(ggridges)
labeller = label_bquote(cols = `epsilon`==.(true_eps),rows = `alpha`==.(true_a))
pa = ggplot(data=full_pars) + geom_density_ridges(aes(x=a,y=as.character(true_b), fill = as.character(true_b))) + geom_point(aes(x=true_a, y=as.character(true_b)))+xlim(0,2)+
  xlab(TeX('\\alpha'))+ylab(TeX(''))+labs(fill=TeX('\\beta'))+
  facet_grid(true_a~true_eps ,labeller=labeller, scales='fixed')+ theme(aspect.ratio = .5,axis.title.y=element_blank(),
                                                                       axis.text.y=element_blank(),
                                                                       axis.ticks.y=element_blank()) + theme_bw()

labeller = label_bquote(cols = `epsilon`==.(true_eps),rows = `alpha`==.(true_a))
peps = ggplot(data=full_pars) + geom_density_ridges(aes(x=eps,y=as.character(true_b), fill = as.character(true_b))) +geom_point(aes(x=true_eps, y=as.character(true_b)))+
  xlab(TeX('\\epsilon'))+ylab(TeX(''))+labs(fill=TeX('\\beta'))+
  facet_grid(true_a~true_eps ,labeller=labeller, scales='fixed')+ theme(aspect.ratio = .5,axis.title.y=element_blank(),
                                                                       axis.text.y=element_blank(),
                                                                       axis.ticks.y=element_blank()) + theme_bw()

labeller = label_bquote(cols = `epsilon`==.(true_eps),rows = `alpha`==.(true_a))
pk0 = ggplot(data=full_pars) + stat_binline(aes(x=k0,y=as.character(true_b), fill = as.character(true_b)),binwidth = 1) +geom_point(aes(x=true_k0, y=as.character(true_b)))+
  xlab(TeX('$k_0$'))+ylab(TeX(''))+labs(fill=TeX('\\beta'))+xlim(0,200)+
  facet_grid(true_a~true_eps ,labeller=labeller, scales='fixed')+ theme(aspect.ratio = .5,axis.title.y=element_blank(),
                                                                       axis.text.y=element_blank(),
                                                                       axis.ticks.y=element_blank()) + theme_bw()

labeller = label_bquote(cols = `epsilon`==.(true_eps),rows = `alpha`==.(true_a))
pb=ggplot(data=full_pars) + geom_density_ridges(aes(x=b,y=as.character(true_b), fill = as.character(true_b))) +geom_point(aes(x=true_b, y=as.character(true_b)))+
  xlab(TeX('\\beta'))+ylab('')+labs(fill=TeX('\\beta'))+
  facet_grid(true_a~true_eps ,labeller=labeller, scales='fixed') + theme(aspect.ratio = .5,axis.title.y=element_blank(),
                                                                        axis.text.y=element_blank(),
                                                                        axis.ticks.y=element_blank()) + theme_bw()
full_plot = ggpubr::ggarrange(pa,peps,pk0,pb,common.legend = T, nrow=2,ncol=2,legend='right')

plot(full_plot)


saveRDS(full_plot, 'mcmc_plot.rds')


ggplot(data=full_pars) + geom_point(aes(x=a, y=eps, colour = as.character(true_b)),alpha=0.1) + geom_point(aes(x=true_a,y=true_eps)) + facet_grid(true_a~true_eps)

plot(recovery_list[[1]]$output$pars$a, type='l')

recovery_list[[1]]$output$sig
# -------------------------------------------------------------------------

pref = function(x,a,eps, k0, b){
  return(ifelse(x<k0, x^a + eps, k0^a  +eps + b*(x-k0)))
}
j=36
mcdat = recovery_list[[j]]$output$pars
degs = sort(unique(recovery_list[[j]]$degs))
pref_mat = matrix(ncol=length(degs), nrow=nrow(mcdat))
for(i in 1:nrow(mcdat)){
  print(i)
  pref_mat[i,] = pref(degs, mcdat$a[i], mcdat$eps[i],mcdat$k0[i], mcdat$b[i])/sum(pref(degs, mcdat$a[i], mcdat$eps[i],mcdat$k0[i], mcdat$b[i]))
}

pref975 = apply(pref_mat, 2, quantile, prob=0.975)
pref025 = apply(pref_mat, 2, quantile, prob=0.025)
pref50 = apply(pref_mat, 2, quantile, prob=0.5)
plot(degs, pref50,pch=20,type='l', log='xy')
lines(degs,pref975,lty=2)
lines(degs, pref025,lty=2)
lines(degs, pref(degs, recover_pars$a[j],recover_pars$eps[j],recover_pars$k0[j],recover_pars$b[j])/
        sum(pref(degs, recover_pars$a[j],recover_pars$eps[j],recover_pars$k0[j],recover_pars$b[j])), col='red')



# -------------------------------------------------------------------------



out2 = out

pars1 = out$output$pars
pars2 = out2$output$pars

pars1$true_b = 0.5
pars2$true_b = 1
pars_full = rbind(pars1, pars2)

library(ggridges)
library(gridExtra)
ggplot(data=pars_full,aes(x=a, y=as.character(true_b),fill=as.character(true_b))) + geom_density_ridges()
ggplot(data=pars_full,aes(x=eps, y=as.character(true_b),fill=as.character(true_b))) + geom_density_ridges()
ggplot(data=pars_full,aes(x=k0, y=as.character(true_b),fill=as.character(true_b))) + geom_density_ridges()
ggplot(data=pars_full,aes(x=b, y=as.character(true_b),fill=as.character(true_b))) + geom_density_ridges()

grid.arrange(pa,peps,pk0,pb)


# -------------------------------------------------------------------------


library(grid)
library(gridExtra)

rplots = list()
for(i in 1:length(recovery_list)){
  df = recovery_list[[i]]$output$pars
  rplots[[i]] = ggplot() + geom_boxplot(data=stack(df),aes(y=values)) +
    geom_hline(data=stack(pars[i,]),aes(yintercept=values), colour='red',linetype='dashed',lwd=1.5)  +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          aspect.ratio = 1)+
    facet_wrap(ind~., scales='free') + ggtitle(paste0('a=',pars[i,1],', eps=',pars[i,2],',k0=',pars[i,3],',b=',pars[i,4]))
    
}
grid.arrange(grobs = rplots,nrow=2)


surfit_list = list()

for(j in 1:length(recovery_list)){
  print(j)
  x = sort(unique(recovery_list[[j]]$degs))
  y_mat = numeric(length(x))
  smps = recovery_list[[j]]$output$pars[-(1:1e3),]
  ls = c()
  for(i in 1:nrow(smps)){
    est = as.numeric(smps[i,])
    lambda = find_lambda2(polylin(est[1], est[2]),est[4],est[3])
    ls = c(ls,lambda)
    y_mat = rbind(y_mat,S(x, polylin(est[1], est[2]),lambda, est[3],est[4]))
  }
  y_mat =y_mat[-1,]
  y_975 = apply(y_mat, 2, quantile, prob=0.975,na.rm=T) 
  y_025 = apply(y_mat, 2, quantile, prob=0.025,na.rm=T)
  y_50 = apply(y_mat, 2, quantile, prob=0.5,na.rm=T) 
  
  surfit_list[[j]] = ggplot() + geom_point(data=twbfn::deg_surv(recovery_list[[j]]$degs), aes(x=degree+1, y=surv)) +
    geom_line(data=NULL, aes(x=!!x+1, y=!!y_975), colour = 'red', linetype='dashed')+
    geom_line(data=NULL, aes(x=!!x+1, y=!!y_50),colour='red')+
    geom_line(data=NULL, aes(x=!!x+1, y=!!y_025), colour = 'red', linetype='dashed')+
    scale_x_log10()  +scale_y_log10(limits = c(1/length(recovery_list[[j]]$degs),1))+theme_bw() + theme(aspect.ratio = 1)
}
grid.arrange(grobs = surfit_list,nrow=2)



y_mat = numeric(length(x))
  smps = recovery_list[[j]]$output$pars[-(1:1e3),]
  ls = c()
  for(i in 1:nrow(smps)){
    est = as.numeric(smps[i,])
    lambda = find_lambda2(polylin(est[1], est[2]),est[4],est[3])
    ls = c(ls,lambda)
    y_mat = rbind(y_mat,S(x, polylin(est[1], est[2]),lambda, est[3],est[4]))
  }
  y_mat =y_mat[-1,]
  y_975 = apply(y_mat, 2, quantile, prob=0.975,na.rm=T) 
  y_025 = apply(y_mat, 2, quantile, prob=0.025,na.rm=T)
  y_50 = apply(y_mat, 2, quantile, prob=0.5,na.rm=T) 
  
  surfit_list[[j]] = ggplot() + geom_point(data=twbfn::deg_surv(recovery_list[[j]]$degs), aes(x=degree+1, y=surv)) +
    geom_line(data=NULL, aes(x=!!x+1, y=!!y_975), colour = 'red', linetype='dashed')+
    geom_line(data=NULL, aes(x=!!x+1, y=!!y_50),colour='red')+
    geom_line(data=NULL, aes(x=!!x+1, y=!!y_025), colour = 'red', linetype='dashed')+
    scale_x_log10()  +scale_y_log10(limits = c(1/length(recovery_list[[j]]$degs),1))+theme_bw() + theme(aspect.ratio = 1)





















































































































# -------------------------------------------------------------------------
counts_to_degs = function(dat){
  return(rep(dat[,1],dat[,2]))
}


full = read.csv('degs.csv')
unique(full$name)
dat = full[full$name=="as-caida20071105",-1]
plot(dat, log='xy')

dat = dat[dat[,1]>=5,]
dat[,1] = dat[,1]-5
degs = counts_to_degs(dat)
plot(twbfn::deg_surv(degs), log='xy',ylim=c(1/length(degs),1))

# -------------------------------------------------------------------------
out = polylin_mcmc(1e5, dat, stop_at = 1e4,quiet = F)


smps = out$pars
ls = c()
x = sort(unique(degs))
y_mat = numeric(length(x))
for(i in 1:nrow(smps)){
  print(i)
  est = as.numeric(smps[i,])
  lambda = find_lambda2(polylin(est[1], est[2]),est[4],est[3])
  ls = c(ls,lambda)
  y_mat = rbind(y_mat,S(x, polylin(est[1], est[2]),lambda, est[3],est[4]))
}
y_mat =y_mat[-1,]
y_975 = apply(y_mat, 2, quantile, prob=0.975,na.rm=T) 
y_025 = apply(y_mat, 2, quantile, prob=0.025,na.rm=T)
y_50 = apply(y_mat, 2, quantile, prob=0.5,na.rm=T) 

ggplot() + geom_point(data=twbfn::deg_surv(degs), aes(x=degree, y=surv)) +
  geom_line(data=NULL, aes(x=x+1, y=y_975), colour = 'red', linetype='dashed')+
  geom_line(data=NULL, aes(x=x+1, y=y_50),colour='red')+
  geom_line(data=NULL, aes(x=x+1, y=y_025), colour = 'red', linetype='dashed')+
  scale_x_log10()  +scale_y_log10(limits = c(1/length(degs),1))+theme_bw() + theme(aspect.ratio = 1)













# -------------------------------------------------------------------------

as = c(.5,1,1.5)
bs = c(.1,.5,1,1.5)
eps = c(.001,.05,.5)
k0s = 50
x = k0s:500
pars = expand.grid(as,bs,eps,k0s,x)
names(pars) = c('a','b','eps','k0','x')
lambdas = numeric(nrow(pars))
for(i in 1:nrow(pars)){
  lambdas[i] = find_lambda2(polylin(pars$a[i], pars$eps[i]), pars$b[i], pars$k0[i])
}
pars = cbind(pars,lambdas)
names(pars) = c('a', 'b', 'eps', 'k0','x','lambda')
surv =numeric(nrow(pars))
gp_surv =numeric(nrow(pars))
for(i in 1:nrow(pars)){
  surv[i] = S(pars$x[i], polylin(pars$a[i], pars$eps[i]),pars$lambda[i], pars$k0[i],pars$b[i])/
    S(pars$k0[i], polylin(pars$a[i], pars$eps[i]),pars$lambda[i], pars$k0[i],pars$b[i])
  gp_surv[i] = mev::pgp(pars$x[i],pars$k0[i]-1,(pars$k0[i]^pars$a[i] +pars$eps[i] )/pars$lambda[i],pars$b[i]/pars$lambda[i],lower.tail = F)
}
pars = cbind(pars,surv,gp_surv)
names(pars) = c('a', 'b', 'eps', 'k0','x','lambda','surv','gp_surv')

library(latex2exp)
labeller = label_bquote(cols = `epsilon`==.(eps),rows = `alpha`==.(a))

ggplot(data = pars) + geom_line(aes(x=(x),y=surv, linetype=as.character(b),colour =as.character(b)), lwd=1) +
  geom_line(data = pars[pars$x>=pars$k0,],aes(x=(x+1),y=gp_surv, linetype=as.character(b)),colour ='gray')+
  scale_x_log10() + scale_y_log10()+
  theme(aspect.ratio = 1/2) + theme_bw() + xlab('Total Degree')  +ylab('Survival') + labs(linetype=TeX('\\beta'), colour =TeX('\\beta')) +
  facet_grid(a~eps,labeller = labeller,scales='free')


# ggplot(data = pars) + geom_line(aes(x=(x+1),y=(gp_surv-surv)/surv, linetype=as.character(b),colour =as.character(b)), lwd=1)+
#   theme(aspect.ratio = 1/2) + theme_bw() + xlab('Total Degree')  +ylab('Survival') + labs(linetype=TeX('\\beta'), colour =TeX('\\beta')) + 
#   facet_grid(a~eps,labeller = labeller)

# -------------------------------------------------------------------------

rho_optim_ba= Vectorize(function(a,eps, b, k0){
  return(abs(rho(2*b, polylin(a, eps), b, k0)-1))
},vectorize.args = 'a')
find_a_ba = Vectorize(function(b, eps, k0){
  out = optimise(rho_optim_ba, c(0.00001,3), b=b, eps=eps, k0=k0)$minimum
  return(out)
}, vectorize.args = 'b')
N =50
as = seq(0,2,l=N+1)[-1]
bs = seq(0,2,l=N+1)[-1]
eps = c(0.001,0.1,1)
k0 = c(25,100)
pars = expand.grid(as,bs,eps,k0)
names(pars) = c('a', 'b', 'eps', 'k0')
lambdas = numeric(nrow(pars))
a_for_ba = numeric(nrow(pars))
for(i in 1:nrow(pars)){
  message(i,'/',nrow(pars))
  lambdas[i] = find_lambda2(polylin(pars$a[i], pars$eps[i]), pars$b[i], pars$k0[i])
  a_for_ba[i] = find_a_ba(pars$b[i], pars$eps[i], pars$k0[i])
}
pars = cbind(pars,lambdas,a_for_ba)
names(pars) = c('a', 'b', 'eps', 'k0','lambda','ba')
labeller = label_bquote(cols = `epsilon`==.(eps),rows = ~k[0]==.(k0))

ggplot(pars) + geom_raster(aes(x=a,y=b,fill=b/lambda)) +
  geom_line(aes(x=ba, y=b),linetype='dashed', colour='red',lwd=1)+geom_point(aes(x=1,y=1))+
  scale_fill_paletteer_c(palette='grDevices::Blues',limits=c(0,1),direction=-1)+theme(aspect.ratio = 1)+ylim(min(bs), max(bs))+xlim(min(as), max(as))+
  facet_grid(k0~eps,labeller = labeller) + labs(fill=TeX('\\xi')) + xlab(TeX('\\alpha')) + ylab(TeX('\\beta'))

find_a_ba(0.1, 1, 25)


# -------------------------------------------------------------------------
sample_gpa_tree = function(n, g,b,k0,m=1){
  degs = c(0)
  prefs = c(g(0))
  for(i in 1:n){
    # message(i)
    selected = sample(1:length(degs), min(m, length(degs)), prob=prefs,replace = T)
    degs = c(degs, 0)
    prefs = c(prefs, g(0))
    degs[selected] = degs[selected]+1
    prefs[selected] = ifelse(degs[selected]<k0, g(degs[selected]), g(k0) + b*(degs[selected]-k0))
  }
  return(degs)
}

samplepolylin2 = function(n,a,eps,b,k0,m=1){
  g01 = function(x){
    return((x+m)^a+eps-m^a)
  }
  degs = sample_gpa_tree(n,g01,b,k0,m)
  return(degs)
}

n=1e5
a=1
b=0.1
eps=0.1
k0=10
degs = samplepolylin2(n, a, eps, b,k0,m=1)
plot(twbfn::deg_surv(degs+1), log='xy',pch=20,ylim=c(1/n,1),xlim=c(1,n))
for(i in 2:5){
  degs = samplepolylin2(n, a, eps, b, k0,m=i)
  points(twbfn::deg_surv(degs+i), col=i,pch=20)
}
 










































