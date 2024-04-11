
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
  print(o)
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


# omega.vec(1:100, g)



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



# Omega -------------------------------------------------------------------

g = function(s,a=1){
  return((s+1)^a)
}

x=1:200
n=6
a.list = rev(c(seq(0.7,1,length.out=n),seq(1,1.3,length.out=n)))
a.list = unique(a.list)
pal = fish(2*n-1,option = "Thalassoma_hardwicke")
df = data.frame(k=x)
df = cbind(df,omega.vec(x, g, a.list[1]))
plot(x,df[,2], type='l', ylim = c(5e-2,max(omega.vec(x, g, a.list[1]))),log='xy',
     ylab = 'omega', xlab = 'k', col=pal[1])
i=3
for(a in a.list[-1]){
  df = cbind(df,omega.vec(x, g, a))
  # lines(x,df[,i],col=pal[which(a.list==a)])
  i=i+1
}
lines(x,omega.vec(x,g, 1),col=1, lty=2)

names(df)[-1] = round(a.list,3)

require(ggplot2)
require(reshape2)
require(fishualize)
df = melt(df, id.vars = 'k', variable.name = 'a')



pal = fish(2*n-1,option = "Bodianus_rufus")
ggplot(df, aes(k,value))+geom_line(aes(colour = a),linewidth=1) + scale_x_log10() + scale_y_log10() + 
  scale_color_manual(values=pal)

# -------------------------------------------------------------------------
require(latex2exp)

omega.outer = function(x,a,g, lambda=2){
  return(omega.vec(x, g, a))
}
n = 100
x = 1:(2*n)
a.list = 10^seq(-1,0.5, length.out=2*n)
z = (outer(x,a.list, omega.outer, g=g))
data = expand.grid(x,a.list)
n.cut = 50
levelplot(log10(z),col.regions =  magma(n.cut+1), column.values = a.list, row.values = x,
          xlab='k', ylab='a',aspect = 1,ylim=range(a.list),xlim=range(x),
          contour=F,cuts=n.cut)
#heatmap of omega_k for different k values, showing that it diverges for a>1 and converges to 0 for a< 1,
#and converges a non-zero constant for a=1.

#should include plot of a=1 on its own.

# z[z==Inf] = max(z[z!=Inf]) 
# z[z==-Inf] = min(z[z!=-Inf]) 
# plot3d(x,a.list,z)


# -------------------------------------------------------------------------


gs = list()
as = c(0.5,1,1.1,1.2,1.3,1.5,1.8)

# fun1 <- function(n) {
#   force(n)
#   fun2 <- function(x,a) {
#     x^n
#   }
#   return(fun2)
# }
# for(i in 1:length(as)){
#   p = as[i]
#   gs[[i]] = fun1(p)
# }

for(i in 1:length(as)){
  print(CCDF_rat(1e3, 2, 2, g,as[i]))
}



CCDF_rat.vec(k,2,2,g,1.3)


# CCDF_ratio --------------------------------------------------------------


g = function(s,a=1){
  return((s+1)^a)
}
t=5
k = round(10^seq(1,3,length.out=20))
n=25
a.list = rev(seq(1,2,length.out=2*n-1))
a.list = unique(a.list)
pal = fish(2*n-1,option = "Thalassoma_hardwicke")
df = data.frame(k=k)
df = cbind(df,CCDF_rat.vec(k,t,2,g,a.list[1]))
# plot(df[,1],df[,2], type='l',log='xy',
     # ylab = 'omega', xlab = 'k', col=pal[1])
for(a in a.list[-1]){
  df = cbind(df,CCDF_rat.vec(k,t,2,g,a))
  # lines(x,df[,i],col=pal[which(a.list==a)])
}
# lines(x,CCDF_rat.vec(k,t,2,g,1),col=1, lty=2)

names(df)[-1] = round(a.list,3)

require(ggplot2)
require(reshape2)
require(fishualize)
df = melt(df, id.vars = 'k', variable.name = 'a')



pal = fish(2*n-1,option = "Bodianus_rufus")
ggplot(df, aes(k,value))+geom_line(aes(colour = a),linewidth=1) + scale_x_log10()+ 
  scale_color_manual(values=pal)

#seems to imply that CCDF is slowly varying for any a>1


# -------------------------------------------------------------------------


g = function(s,a=1){
  return((s+1)^a)
}
t= seq(2,10)
k = 1e3
n=6
a.list = rev(c(seq(0.7,1,length.out=n),seq(1,1.3,length.out=n)))
a.list = unique(a.list)
pal = fish(2*n-1,option = "Thalassoma_hardwicke")
df = data.frame(t=t)
df = cbind(df,CCDF_rat.vec(k,t,2,g,a.list[1]))
# plot(df[,1],df[,2], type='l',log='xy',
# ylab = 'omega', xlab = 'k', col=pal[1])
for(a in a.list[-1]){
  df = cbind(df,CCDF_rat.vec(k,t,2,g,a))
  # lines(x,df[,i],col=pal[which(a.list==a)])
}
# lines(x,CCDF_rat.vec(k,t,2,g,1),col=1, lty=2)

names(df)[-1] = round(a.list,3)

require(ggplot2)
require(reshape2)
require(fishualize)
df = melt(df, id.vars = 't', variable.name = 'a')



pal = fish(2*n-1,option = "Bodianus_rufus")
ggplot(df, aes(t,value))+geom_line(aes(colour = a),linewidth=1) + scale_x_log10() + scale_y_log10() + 
  scale_color_manual(values=pal)












