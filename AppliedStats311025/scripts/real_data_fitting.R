source('funcs.R')
library(ggplot2)
data_raw = read.csv('degs.csv')
nms = unique(data_raw$name)
selected = c(2,3,5,6,7,8,9,11,12,13,17,22)
fits = list()
plots = list()
for(i in 1:length(selected)){
  print(i)
  fits[[i]] = fit_model(data_raw[data_raw$name==nms[selected[i]],2:3],N=5e4,quiet=F,
                        trunc.at = 5,k0_jump = 1)
}


saveRDS(fits, 'real_mcmc.rds')


for(i in 1:length(fits)){
  plots[[i]] = ggplot() + scale_x_log10() + scale_y_log10() +
    geom_point(data = twbfn::deg_surv(counts_to_degs(fits[[i]]$dat)),aes(x=degree, y=surv))+
    geom_line(aes(x=!!fits[[i]]$dat[,1],y=!!fits[[i]]$surv$est),colour='red') + 
    geom_line(aes(x=!!fits[[i]]$dat[,1],y= !!fits[[i]]$surv$CI[1,]),linetype=2, colour='red') + 
    geom_line(aes(x=!!fits[[i]]$dat[,1],y= !!fits[[i]]$surv$CI[2,]),,linetype=2,colour='red')+
    ggtitle(nms[selected[i]])
}


real_survs = ggpubr::ggarrange(plotlist = plots)

saveRDS(real_survs,'real_survs.rds')




real_survs
# igpd --------------------------------------------------------------------


fits = readRDS('real_mcmc.rds')

gpfit = mev::fit.gpd(counts_to_degs(fits[[1]]$dat),threshold=6)

scale = gpfit$estimate[1]
shape = gpfit$estimate[2]

x = fits[[1]]$dat[,1]
y = mev::pgp(x,loc=6,scale=scale,shape=shape,lower.tail = F)
degs = counts_to_degs(fits[[1]]$dat)

plot(x,y,log='xy',type='l')
points(twbfn::deg_surv(degs[degs>6]))

# -------------------------------------------------------------------------

data_raw = read.csv('degs.csv')
nms = unique(data_raw$name)
selected = c(2,3,5,6,7,8,9,11,12,13,17,22)

nm <- nms[selected[1]]

dat = data_raw[data_raw$name==nm,][,2:3]


degs = counts_to_degs(dat)


plot(twbfn::deg_surv(degs), log='xy')


fit = fit.gpd(degs, threshold=15)


# -------------------------------------------------------------------------


fits = readRDS('modelfits.rds')



plt = readRDS('real_survs.rds')




library(coda)
plot(as.mcmc(fits[[10]]$smps), type='l')


# gpd  vs model -----------------------------------------------------------


plts = list()

for(j in 1:length(fits)){
  print(j)
  fit = fits[[j]]
  fit$mcmc$pars  = as.data.frame(fit$mcmc$pars)
  names(fit$mcmc$pars) = c('a','eps','k0','b')
  N = nrow(fit$mcmc$pars)
  shapes<- numeric()
  xs <- numeric()
  
  for(i in 1:(nrow(fit$dat)-3)){
    if(fit$dat[i,1]< min(fit$mcmc$pars$k0[seq(5e3, N, by=10)])){
      next;
    }
    if(fit$dat[i,1]>max(fit$mcmc$pars$k0[seq(5e3, N, by=10)] )){
      break;
    }
    igp = mev::fit.gpd(counts_to_degs(fit$dat),threshold=fit$dat[i,1])
    shapes = c(shapes,igp$estimate[2])
    xs = c(xs,fit$dat[i,1])
  }
  plts[[j]] = ggplot() +
    geom_jitter(aes(x=!!fit$mcmc$pars$k0[seq(5e3, N, by=10)] ,y=!!(fit$mcmc$pars$b[seq(5e3, N, by=10)] /fit$mcmc$lambdas[seq(5e3, N, by=10)] )),
                colour='blue', alpha=0.025) + 
    geom_point(aes(x=!!xs[xs<=max(fit$mcmc$pars$k0[seq(5e3, N, by=10)] )],y=!!shapes[xs<=max(fit$mcmc$pars$k0[seq(5e3, N, by=10)] )]), colour='red') +
    ylab('Shape') + xlab('k0') + ggtitle(ggtitle(nms[selected[j]])) + ylim(-0.5,1) + xlim(min(fit$dat$x), max(xs))
}




igp_plt = ggpubr::ggarrange(plotlist = plts,nrow=3,ncol=4)
plot(igp_plt)


saveRDS(igp_plt, 'gp_plt.rds')




df = data.frame(k0 = fit$mcmc$pars$k0,shape = fit$mcmc$pars$b/fit$mcmc$lambdas)

df$real = shapes[match(df$k0,xs)]


plt = readRDS('real_survs.rds')

# -------------------------------------------------------------------------

mixfits = readRDS('../results/mixfits.rds')
fits = readRDS('../results/real_mcmc.rds')
nms = readRDS('../results/nms.rds')
plts = list()
for(j in 1:length(fits)){
  fit = fits[[j]]
  fit$mcmc$pars  = as.data.frame(fit$mcmc$pars)
  names(fit$mcmc$pars) = c('a','eps','k0','b')
  N = nrow(fit$mcmc$pars)
  plts[[j]] = ggplot() +
    geom_boxplot(aes(y=!!mixfits[[j]]$pars$shape,x='Zipf-IGP', fill='Zipf-IGP'))+
    geom_boxplot(aes(x='GPA',y=!!(fit$mcmc$pars$b[seq(1e4, N, by=10)] /fit$mcmc$lambdas[seq(1e4, N, by=10)]),fill='GPA')) +
    labs(fill='Model', x='',y='Shape')+ylim(0,1) +ggtitle(nms[j])
}

shape_plt = ggpubr::ggarrange(plotlist = plts, common.legend = T, legend = 'right')
saveRDS(shape_plt, 'shape_plt.rds')


comp_plot

# -------------------------------------------------------------------------
library(ggplot2)



fits = readRDS('fits.rds')
nms = readRDS('nms.rds')
plts = list()


for(i in 1:length(fits)){
  plts[[i]] = ggplot()  +
    geom_point(aes(x=!!fits[[i]]$dat$x, y=!!fits[[i]]$PA$est),size=0.9)+
    geom_line(aes(x=!!fits[[i]]$dat$x, y=!!fits[[i]]$PA$CI[1,]), lty=2)+
    geom_line(aes(x=!!fits[[i]]$dat$x, y=!!fits[[i]]$PA$CI[2,]), lty=2)+
    scale_x_log10() + scale_y_log10() + ggtitle(nms[i])+
    xlab('Degree') + ylab('P') + theme_bw()
}

paplot = ggpubr::ggarrange(plotlist = plts)

compplot =  readRDS('comp_plot.rds')

dev.new()
paplot


# -------------------------------------------------------------------------


fits = readRDS('fits.rds')
nms = readRDS('nms.rds')

plts = list()



pref = function(pars, x){
  pars = unlist(pars)
  out = x
  out[x<pars[3]] = x[x<pars[3]]^pars[1] + pars[2]
  out[x>=pars[3]] = pars[3]^pars[1] + pars[2] + pars[4]*(x[x>=pars[3]]-pars[3])
  return(out)
}

for(i in 1:length(fits)){
  message(i)
  x = 1:max(fits[[i]]$dat$x)
  prefmat = apply(as.matrix(fits[[i]]$smps), 1, pref, x=x)
  prefCI = apply(prefmat, 1, quantile, prob = c(0.025, 0.5, 0.975))
  
  plts[[i]] = ggplot() + geom_line(aes(x = !!x, y= !!prefCI[1,]),lty=2) +
    geom_line(aes(x = !!x, y= !!prefCI[2,])) +
    geom_line(aes(x = !!x, y= !!prefCI[3,]),lty=2) +
     ggtitle(nms[i])
  
}

ggpubr::ggarrange(plotlist = plts)
# -------------------------------------------------------------------------

sample_tree = function(n, g, m=1,quiet=T){
  df = data.frame(deg = c(0),
                  count = c(1),
                  pref = c(g(0)))
  for(i in 1:n){
    if(!quiet){
      message(i)
    }
    selected = sample(1:nrow(df), min(m, sum(df$count)),prob = df$pref*df$count,replace=T)
    for(j in 1:length(selected)){
      df$count[selected[j]] = df$count[selected[j]]-1
      if((df$deg[selected[j]]+1) %in% df$deg){
        df$count[df$deg==df$deg[selected[j]]+1] = df$count[df$deg==df$deg[selected[j]]+1]+1
      }else{
        df = rbind(df, c(df$deg[selected[j]]+1, 1,
                         g(df$deg[selected[j]]+1)))
      }
      df$count[1] = df$count[1] + 1
      df = df[df$count!=0,]
    }
  }
  return(df)
}






i = 6
x = 0:1e3
prefmat = apply(as.matrix(fits[[i]]$smps), 1, pref, x=x)
prefCI = apply(prefmat, 1, quantile, prob = c(0.025, 0.5, 0.975))
pas  = prefCI[2,]
g = function(x){
  return(pas[x+1])
}
datdegs = counts_to_degs(dat_list[[i]])

# -------------------------------------------------------------------------


n = length(datdegs)
G = sample_tree(n, g, m=1, quiet=F)
degs = counts_to_degs(G[,1:2])
trunc.at=4
plot(twbfn::deg_surv(datdegs[datdegs>trunc.at]), log='xy',pch=20,ylim=c(1/n, 1))
points(twbfn::deg_surv(degs[degs>trunc.at]))

# -------------------------------------------------------------------------

g = function(x){
  return(x+1)
}
newsample = function(n, g, A, p){
  degs = c(0)
  prefs = c(g(0))
  types = c(0)
  for(i in 2:n){
    selected = sample(1:(i-1), 1, prob = prefs)
    degs[selected] = degs[selected]+1
    prefs[selected] = (1-types[selected])*g(degs[selected]) + types[selected]*A*g(degs[selected])
    degs = c(degs, 0)
    newtype = rbinom(1, 1, p)
    types = c(types, newtype)
    prefs = c(prefs, (1-newtype)*g(0) + newtype*A*g(0))
  }
  return(degs)
}

n = 1e4
A = 10
p=0.05
ps = seq(0,1,l=10)

pal = RColorBrewer::brewer.pal()


for(p in ps){
  message(p)
  degs = newsample(n, g, A, p)
  if(p==ps[1]){
    plot(twbfn::deg_surv(degs), log='xy', pch=20, col=which(ps==p))
  }else{
    points(twbfn::deg_surv(degs), pch=20, col=which(ps==p))
  }
}

legend('topright', legend = round(ps,2), fill=1:length(ps))


# -------------------------------------------------------------------------


source('funcs.R')


Smix  = function(x,p,eps,k0,b,lambda){
  g = function(x){
    return(x+eps)
  }
  out = x
  out[x<=k0] = S(x[x<=k0], g, lambda, k0, 1)
  if(b>0){
    out[x>k0] = (1-p)*S(x[x>k0], g, lambda, k0, 1) +
      p*S(x[x>k0], g, lambda, k0, b)
  }else{
    out[x>k0] = (1-p)*S(x[x>k0], g, lambda, k0, 1) +
      p*S(k0, g, lambda, k0, 1)*((k0+eps)/(lambda+k0+eps))^(x[x>k0]-k0)
  }
  
  return(out)
}

rhomix = function(lambda, p, eps, k0,b){
  g = function(x){
    return(x+eps)
  }
  return(
    rho(lambda, g, 1, k0)*(1-p) + rho(lambda, g, b, k0)*p
  )
}

rhomix_optim = function(lambda, p,eps,k0,b){
  return(rhomix(lambda, p, eps, k0,b)-1)
}

find_lambda_mix = function(p, eps, k0,b){
  out = try(uniroot(rhomix_optim, c(b,100), p=p,eps=eps,k0=k0,b=b, extendInt = 'yes')$root,silent=F)
  if(is.character(out)){
    return(NULL)
  }
  return(out)
}


rhomix()

# -------------------------------------------------------------------------



x = 0:1000
eps = 1
p = 0.9
k0 = 20
b=0

lambda = find_lambda_mix(p, eps, k0,b)
y = Smix(x, p, eps, k0,b, lambda)
plot(x+1, y, log='xy',type = 'l',ylim=c(min(y),1))
abline(v = k0, lty=2)



























  
  
  
  
  
  
  
  
  
  
  
  
  
  














