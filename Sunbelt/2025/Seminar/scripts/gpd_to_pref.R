logdiff <- function(l1, l2) { l1 + VGAM::log1mexp(l1-l2); }
digp_log = function(x, shape, scale, threshold){
  out=x
  v=threshold
  sig=scale
  xi=shape
  p1 = pmax(x[x>v]*0, 1 + xi * (x[x>v]+1-v)/sig)
  p2 = pmax(x[x>v]*0, 1 + xi * (x[x>v]-v)/sig)
  pows2 = pows1 = x[x>v]*0 + 1
  pows1[p1>0] = -1/xi
  pows2[p2>0] = -1/xi
  # out[x>v] =  ( p2^pows2 - p1^pows1)
  out[x>v] = logdiff(p2^pows2,p1^pows1)
  out[x<=v] = 0
  return(out)
}
igp_ll_neg = function(pars,dat,threshold=0){
  shape = pars[1]
  scale = pars[2]
  if(shape==0){
    return(1e9)
  }
  if(shape<0 & max(dat[,1]+1)<(threshold - scale/shape)){
    return(1e9)
  }
  out = -sum(dat[,2]*digp_log(dat[,1], shape, scale, threshold))
  if(is.infinite(out)){
    return(1e9)
  }else{
    return(out)
  }
}
fit.igp = function(dat, threshold=0){
  return(optim(c(0.5, 1), fn = igp_ll_neg,dat=dat[dat[,1]>threshold,], threshold=threshold,
               lower=c(-Inf,0.001), upper=c(Inf, Inf),method='L-BFGS-B'))
}

ft = fit.igp(fit$dat,10)


degs = counts_to_degs(fit$dat)


p = digp(degs[degs>10], -0.001, 1, 10)

max(p)
min(p)


which.min(p)

degs[degs>10][2645]


# -------------------------------------------------------------------------


igpd_rho = function(pars, shape, scale, threshold, exceedance){
  p2 = shape*exceedance / (1-shape)
  x = 1:(threshold+1)
  p1 = sum(exp(cumsum(
    log(x^pars[1] - (threshold+1)^pars[1] + pars[2]*scale)-
      log(x^pars[1] - (threshold+1)^pars[1] + pars[2]*(scale+1))
  )))
  return(p1 + p2)
}

igpd_rho2 = function(pars, shape, scale, threshold, exceedance){
  p2 = shape*exceedance / (1-shape)
  
  
  x = 0:(threshold)
  
  alpha = log(pars[2]*scale-shape)/log(threshold+1)
  
  p1 = sum(cumprod(
    (x^alpha +pars[1] )/
      (x^alpha + pars[1] + pars[2])
  ))
  return(p1 + p2)
}


library(ggplot2)
# -------------------------------------------------------------------------


shape = 0.2
scale = 10
threshold = 100
exceedance = 0.05


n=250
as  = seq(0,2,l=n+1)[-1]
ls = seq(0,3,l=n+1)[-1]
pars = data.frame(expand.grid(as,ls))
names(pars) = c('alpha', 'lambda')

pars$rho = apply(pars, 1, igpd_rho,
      shape=shape, scale=scale, threshold=threshold, exceedance=exceedance)
tol = 0.1
pars$close = as.numeric(abs(pars$rho-1)<tol)


ggplot(data=pars, aes(x=alpha, y=lambda)) + 
  geom_raster(aes(fill=close))


# -------------------------------------------------------------------------
source('funcs.R')
degs = counts_to_degs(fits[[3]]$dat)



a=1
eps=1
b=0.5
k0=25

degs = samplepolylin(1e5, a, eps, b, k0)

plot(twbfn::deg_surv(degs), log='xy')


threshold = k0-1
gp = mev::fit.gpd(degs, threshold)


igp = fit.igp(twbfn::deg_count(degs),threshold)
shape = igp$par[1]
scale = igp$par[2]
exceedance = sum(degs>threshold+1)/length(degs)


# -------------------------------------------------------------------------
# 
shape = gp$estimate[2]
scale = gp$estimate[1]


n=250
epss  = seq(0,2,l=n+1)[-1]
ls = seq(0,5,l=n+1)[-1]
pars = data.frame(expand.grid(epss,ls))
names(pars) = c('epsilon', 'lambda')

pars$rho = apply(pars, 1, igpd_rho2,
                 shape=shape, scale=scale, threshold=threshold, exceedance=exceedance)
tol = 0.1
pars$close = as.numeric(abs(pars$rho-1)<tol)

ggplot(data=pars, aes(x=epsilon, y=lambda)) + 
  geom_raster(aes(fill=close))


# full_pars = pars[pars$epsilon==eps[20],]
full_pars = pars
full_pars$beta = full_pars$lambda*shape
full_pars$k0 = threshold+1
full_pars$alpha = log(full_pars$lambda*scale-shape)/log(threshold+1)

full_pars$close[is.na(full_pars$close)]=0


ggplot(data=full_pars[full_pars$close==1,], aes(x=alpha, y=beta)) + 
  geom_point() +
  geom_point(x=a,y=b, colour='red')+
  xlim(0,3) + ylim(0,3)


ggplot(data=full_pars[full_pars$close==1,], aes(x=alpha, y=epsilon)) + 
  geom_point() +
  geom_point(x=a,y=eps, colour='red')+
  xlim(0,2) + ylim(0,2)

close_pars = full_pars[full_pars$close==1,]

pref = function(x,a,eps, k0, b){
  return(ifelse(x<k0, (x)^a + eps, (k0)^a  +eps + b*(x-k0)))
}


counts = twbfn::deg_count(degs)
x = counts[,1]

pref_mat = matrix(ncol=length(x), nrow = nrow(close_pars))

for(r in 1:nrow(close_pars)){
  print(r)
  pref_mat[r,] = pref(x, close_pars$alpha[r], close_pars$epsilon[r],
                      close_pars$k0[r], close_pars$beta[r])
}

for(r in 1:nrow(close_pars)){
  if(r==1){
    plot(x, pref_mat[r,]/sum(counts[,2]*pref_mat[r,]),pch=20,log='xy')
  }else{
    points(x, pref_mat[r,]/sum(counts[,2]*pref_mat[r,]),pch=20)
  }
}

real_pref = pref(x,a, eps, k0, b)
points(x,real_pref/sum(counts[,2]*real_pref), col='red',pch=20)









# -------------------------------------------------------------------------

























# -------------------------------------------------------------------------
# 
# shape = gp$estimate[2]
# scale = gp$estimate[1]
# exceedance = gp$pat
# 
# n=200
# epss  = seq(0,2,l=n+1)[-1]
# ls = seq(0,5,l=n+1)[-1]
# pars = data.frame(expand.grid(epss,ls))
# names(pars) = c('alpha', 'lambda')
# 
# pars$rho = apply(pars, 1, igpd_rho,
#                  shape=shape, scale=scale, threshold=threshold, exceedance=exceedance)
# tol = 0.1
# pars$close = as.numeric(abs(pars$rho-1)<tol)
# 
# ggplot(data=pars, aes(x=alpha, y=lambda)) + 
#   geom_raster(aes(fill=close))
# 
# 
# # full_pars = pars[pars$epsilon==eps[20],]
# full_pars = pars
# full_pars$beta = full_pars$lambda*shape
# full_pars$k0 = threshold+1
# full_pars$epsilon = full_pars$lambda * scale - (threshold+1)^full_pars$alpha
# 
# full_pars$close[is.na(full_pars$close)]=0
# 
# 
# ggplot(data=full_pars[full_pars$close==1,], aes(x=alpha, y=beta)) + 
#   geom_point() + xlim(0,2) + ylim(0,2)
# 
# 
# ggplot(data=full_pars[full_pars$close==1,], aes(x=alpha, y=epsilon)) + 
#   geom_point() + xlim(0,2) + ylim(0,2)
# 
# close_pars = full_pars[full_pars$close==1,]
# 
# pref = function(x,a,eps, k0, b){
#   return(ifelse(x<k0, (x)^a + eps, (k0)^a  +eps + b*(x-k0)))
# }
# 
# 
# counts = twbfn::deg_count(degs)
# x = counts[,1]
# 
# pref_mat = matrix(ncol=length(x), nrow = nrow(close_pars))
# 
# for(r in 1:nrow(close_pars)){
#   print(r)
#   pref_mat[r,] = pref(x, close_pars$alpha[r], close_pars$epsilon[r],
#                       close_pars$k0[r], close_pars$beta[r])
# }
# 
# for(r in 1:nrow(close_pars)){
#   if(r==1){
#     plot(x, pref_mat[r,]/sum(counts[,2]*pref_mat[r,]),pch=20, log='xy')
#   }else{
#     points(x, pref_mat[r,]/sum(counts[,2]*pref_mat[r,]),pch=20)
#   }
# }
# 
# real_pref = pref(x,a, eps, k0, b)
# points(x,real_pref/sum(counts[,2]*real_pref), col='red',pch=20)





































































































