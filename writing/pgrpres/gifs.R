library(igraph)
library(ggraph)
library(animation)

m=2
N=20
G = sample_pa(N, power=1, m=m)
cds = as.data.frame(layout_with_kk(G))
names(cds) = c('x','y')

plot_list = list()


plot_list[[1]] = ggraph(subgraph(G,1), layout=cds[1,]) + geom_node_point(size=6, color='red') + theme_void()


ne = length(E(subgraph(G,1:1)))
for(i in 2:N){
  print(ne)
  ecol = ( (1:length(E(subgraph(G,1:i)))) > (ne)  ) + 1
  ecol[ecol==1] = '#000'
  ecol[ecol==2] = '#f00'
  plot_list[[i]] = ggraph(subgraph(G,1:i), layout=cds[1:i,]) +
    geom_node_point(size=6,colour = 1 + ((1:i)==i)) +
    geom_edge_link0(color = ecol, arrow = arrow(length=unit(4,'mm'),type='closed'))+
    lims(x=range(cds$x), y=range(cds$y)) + theme_void()
  ne = length(E(subgraph(G,1:i)))
}


saveHTML({
  for(i in 1:N){
    print(plot_list[[i]])
  }
}, img.name='ba.html', verbose=F)
# -------------------------------------------------------------------------

N=20
G = twbfn::sample_ppa(N,a=1.5, b=1,n0=5)
G = igraph::delete.edges(G,1)


cds = as.data.frame(layout_with_kk(G))
names(cds) = c('x','y')

plot_list = list()


plot_list[[1]] = ggraph(subgraph(G,1), layout=cds[1,]) + geom_node_point(size=6, color='red') + theme_void()


ne = length(E(subgraph(G,1:1)))
for(i in 2:N){
  print(ne)
  ecol = ( (1:length(E(subgraph(G,1:i)))) > (ne)  ) + 1
  ecol[ecol==1] = '#000'
  ecol[ecol==2] = '#f00'
  plot_list[[i]] = ggraph(subgraph(G,1:i), layout=cds[1:i,]) +
    geom_node_point(size=6,colour = 1 + ((1:i)==i)) +
    geom_edge_link0(color = ecol, arrow = arrow(length=unit(4,'mm'),type='closed'))+
    lims(x=range(cds$x), y=range(cds$y)) + theme_void()
  ne = length(E(subgraph(G,1:i)))
}


saveHTML({
  for(i in 1:N){
    print(plot_list[[i]])
  }
}, img.name='gpa',htmlfile='gpa.html', verbose=F)

# -------------------------------------------------------------------------

x = 1:100
lam = twbfn::find_lambda_cpp(20,1.2,1,0)
y = exp(twbfn::pppa(x,20,1.2,1,0,lam))

n0s = seq(10,50,by=5)
as = seq(0.5,1.5,length.out=7)
bs = seq(0.5,1.5,length.out=7)

pars = expand.grid(n0s,as,bs)
names(pars) = c('n0','a','b')
pars$lambda = 0
for(r in 1:nrow(pars)){
  print(r)
  pars$lambda[r] = twbfn::find_lambda_cpp(pars$n0[r],pars$a[r], pars$b[r], 0)
}

ys = matrix(nrow=nrow(pars), ncol=length(x))
for(i in 1:nrow(ys)){
  print(i)
  ys[i,] = exp(twbfn::pppa(x,pars$n0[i], pars$a[i], pars$b[i],0, pars$lambda[i]))
}

# -------------------------------------------------------------------------
library(dplyr)
library(data.table)
dat = cbind(pars, ys)
data2 = melt(ys, id.vars='x')
names(data2) = c('parid', 'x','value')
data2$n0 = 0
data2$a=0
data2$b=0
data2$lambda=0
for(r in 1:nrow(data2)){
  print(r)
  data2$n0[r] = pars$n0[data2$parid[r]]
  data2$a[r] = pars$a[data2$parid[r]]
  data2$b[r] = pars$b[data2$parid[r]]
  data2$lambda[r] = pars$lambda[data2$parid[r]]
}

# -------------------------------------------------------------------------
head(data2)


# devtools::install_github('thomasp85/gganimate') 

library(gganimate)


bplot = ggplot(data=data2[data2$n0==10 & data2$a==4/3,]) +
  geom_point(aes(x=x, y=value), colour='darkred') +geom_vline(xintercept = 10, linetype=2)+
  scale_x_log10()+scale_y_log10(limits=c(1e-4,1))+ggtitle('a=1, b varying, T=10')+
  transition_time(b)

aplot = ggplot(data=data2[data2$n0==10 & data2$b==0.5,]) +
  geom_point(aes(x=x, y=value), colour='darkred') +geom_vline(aes(xintercept = 10), linetype=2)+
  scale_x_log10()+scale_y_log10(limits=c(1e-4,1))+ggtitle('a varying, b=1, T=10')+
  transition_time(a)

n0plot = ggplot(data=data2[data2$a==4/3 & data2$b==0.5,]) +
  geom_point(aes(x=x, y=value), colour='darkred') + geom_vline(aes(xintercept=n0),linetype=2)+
  scale_x_log10()+scale_y_log10(limits=c(1e-4,1))+ggtitle('a=4/3, b=0.5, T varying')+
  transition_time(n0)

# 
# anim_save('aplot.gif',animation = aplot, rewind=T)
# anim_save('bplot.gif',animation = bplot, rewind=T)
anim_save('n0plot.gif',animation = n0plot, rewind=T)

# -------------------------------------------------------------------------

degs = read.csv('degs.csv')


par(mfrow=c(6,6))
plist = list()
for(i in 1:length(unique(degs$name))){
  nm = unique(degs$name)[i]
  dat = degs[degs$name==nm,][2:3]
  plist[[i]] = ggplot(data=dat) + geom_point(aes(x=x, y=1-cumsum(count)/sum(count))) +ggtitle(nm)+
    scale_x_log10()+scale_y_log10()
}

gridExtra::grid.arrange(grobs=plist)

nm = 'as20000102'
dat = degs[degs$name==nm,][2:3]
fit = twbfn::ppa_mle(as.matrix(dat))


x=unique(dat$x)
lambda = twbfn::find_lambda_cpp(fit$par[1], fit$par[2], fit$par[3], 0)
y= exp(twbfn::pppa(x,fit$par[1], fit$par[2], fit$par[3],0,lambda))
fitdat = data.frame(x=x,y=y)


G = twbfn::sample_ppa(1e4, 1.4, 0.1, 50)





mcmc = twbfn::ppa_mcmc(as.matrix(dat), 10, 1,1,2)
library(parallel)
twbfn::mcmc_surv_plot(mcmc, burn_in = 1e4)

ggplot(data=dat) + geom_point(aes(x=x, y=1-cumsum(count)/sum(count))) +
  geom_line(data=fitdat, aes(x=x,y=y))+
  ggtitle(nm)+
  scale_x_log10()+scale_y_log10()

# -------------------------------------------------------------------------
a=1.2
b=0.5
n0=20
smplr = function(n){
  return(rep(3, n))
}
library(wdnet)
G = twbfn::sample_ppa(5e3, a, b, n0,smplr=smplr, quiet=F)

df = twbfn::deg_surv(degree(G))
lambda = twbfn::find_lambda_cpp(n0,a,b,0)
x = 1:100
y = exp(twbfn::pppa(x, n0,a,b,0,lambda))
ggplot(data=df)  +geom_point(aes(x=degree,y=surv)) +
  geom_line(data=data.frame(x=x,y=y), aes(x=x, y=y))+
  scale_x_log10() + scale_y_log10(limits=c(1e-4,1))

# -------------------------------------------------------------------------
































x = 1:100
lam = twbfn::find_lambda_cpp(20,1.2,1,0)
y = exp(twbfn::pppa(x,20,1.2,1,0,lam))

n0s = 100
as = seq(0.5,1,length.out=7)
bs = 1

pars = expand.grid(n0s,as,bs)
names(pars) = c('n0','a','b')
pars$lambda = 0
for(r in 1:nrow(pars)){
  print(r)
  pars$lambda[r] = twbfn::find_lambda_cpp(pars$n0[r],pars$a[r], pars$b[r], 0)
}

ys = matrix(nrow=nrow(pars), ncol=length(x))
for(i in 1:nrow(ys)){
  print(i)
  ys[i,] = exp(twbfn::pppa(x,pars$n0[i], pars$a[i], pars$b[i],0, pars$lambda[i]))
}

# -------------------------------------------------------------------------
library(dplyr)
library(data.table)
dat = cbind(pars, ys)
data2 = melt(ys, id.vars='x')
names(data2) = c('parid', 'x','value')
data2$n0 = 0
data2$a=0
data2$b=0
data2$lambda=0
for(r in 1:nrow(data2)){
  print(r)
  data2$n0[r] = pars$n0[data2$parid[r]]
  data2$a[r] = pars$a[data2$parid[r]]
  data2$b[r] = pars$b[data2$parid[r]]
  data2$lambda[r] = pars$lambda[data2$parid[r]]
}

# -------------------------------------------------------------------------
head(data2)


# devtools::install_github('thomasp85/gganimate') 

library(gganimate)



aplot_2 = ggplot(data=data2) +
  geom_point(aes(x=x, y=value), colour='darkred') +
  scale_x_log10()+scale_y_log10(limits=c(1e-4,1))+ggtitle(label = '{round(frame_time,2)}')+
  transition_time(a)



anim_save('aplot_2.gif',animation = aplot_2, rewind=T)






# -------------------------------------------------------------------------






degs = read.csv('degs.csv')






























