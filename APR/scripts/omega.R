library(fishualize)
require(ggplot2)
require(reshape2)
g = function(s,a=1){
  return((s+1)^a)
}
f = Vectorize(function(k,lambda,g,a=1){
  if(k==0){
    return(lambda/(lambda+g(k,a)))
  }
  ks = 0:(k-1)
  fracs = 1 - (lambda/(lambda+g(ks,a)))
  return(lambda/(lambda+g(k,a)) * prod(fracs))
}, vectorize.args = 'k')


f_rat = Vectorize(function(k,t,lambda,g,a=1){
  return(log(f(t*k, lambda, g,a))-log(f(k,lambda, g,a)))
}, vectorize.args = 'k')

CDF = Vectorize(function(k,lambda,g,a=1){
  return(sum(f(0:(k-1), lambda, g,a)))
}, vectorize.args='k')
CCDF_rat = Vectorize(function(k,t,lambda, g,a=1){
  return((1-CDF(t*k, lambda, g,a))/(1-CDF(k, lambda, g,a))) 
}, vectorize.args = 'k')

CCDF_rat.vec = function(k,t,lambda,g,a=1){
  o =  uniroot(cond.b, interval=c(0,30), g=g,a=a)
  lambda =o$root
  # print(o)
  return(return(CCDF_rat(k,t,lambda,g,a)))
}
omega = Vectorize(function(k,lambda,g,a=1){
  p1 = log((1-CDF(k+1,lambda,g,a))/(1-CDF(k+2,lambda,g,a)))^-1
  p2 = log((1-CDF(k,lambda,g,a))/(1-CDF(k+1,lambda,g,a)))^-1
  return(p1-p2)
}, vectorize.args = c('k','a'))
omega.vec = function(k,g,a=1){
  lambda = uniroot(cond.b, interval=c(0,30), g=g,a=a)$root
  return(omega(k,lambda,g,a))
}

cond.b = Vectorize(function(lambda, g,a=1){
  vals = 0:1e4
  fracs = g(vals,a)/(lambda+g(vals,a))
  out = 0
  for(n in 1:length(vals)){
    eps = prod(fracs[1:n])
    # print(eps)
    out = out + eps
    if(eps<1e-4){
      return(out-1)
    }
  }
  return(out-1)
}, vectorize.args = 'lambda')












