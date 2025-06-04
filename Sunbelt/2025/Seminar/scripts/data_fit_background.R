source('funcs.R')
library(ggplot2)
library(crandep)
library(data.table)
data_raw = read.csv('degs.csv')
nms = unique(data_raw$name)
selected = c(2,3,5,6,7,8,9,11,12,13,17,22)
fits = list()
mixfits = list()
plots = list()
dat_list = list()
trunc.at = 3



data.files = list.files('./data') 
nms = tstrsplit(data.files,'out.',keep=2)[[1]]

saveRDS(nms, 'nms.rds')


i=4
message(nms[i])
df  =as.data.frame(table(table(fread(paste0('./data/', data.files[i]))[,2])))
df[,1] = as.numeric(as.character(df[,1]))
dat_list[[i]] = as.matrix(df)
fits[[i]] = fit_model(dat_list[[i]],N=1e5,quiet=F,
                      burn_in=1e4,thin.by=5,trunc.at = 4,k0_jump = 2,lookback = 1e3, cov_scale=0.05)



for(i in 2:2){
  message(nms[i])
  df  =as.data.frame(table(table(fread(paste0('./data/', data.files[i]))[,2])))
  df[,1] = as.numeric(as.character(df[,1]))
  dat_list[[i]] = as.matrix(df)
  fits[[i]] = fit_model(dat_list[[i]],N=1e5,quiet=F,
                        burn_in=1e4,thin.by=5,trunc.at = trunc.at,k0_jump = 2,lookback = 1e3, cov_scale=0.05)
}


saveRDS(dat_list, 'dat_list.rds')

dat_list = readRDS('dat_list.rds')

plots = list()

for(i in 1:length(fits)){
  udf = as.data.frame(table(fits[[i]]$smps$k0))
  udf$Var1 = as.numeric(as.character(udf$Var1))
  sdf = twbfn::deg_surv(counts_to_degs(dat_list[[i]]))
  ptrunc = sdf$surv[sdf$degree==min(fits[[i]]$dat[,1])]
  plots[[i]] = ggplot() + scale_x_log10(limits = c(1, max(dat_list[[i]][,1]))) + scale_y_log10(limits=c(1/sum(dat_list[[i]][,2]),1)) +
    geom_vline(aes(xintercept = jitter(!!fits[[i]]$smps$k0[seq(1,length(fits[[i]]$smps$k0),by=20)]),2), colour='blue', alpha=0.1)+
    geom_point(data = twbfn::deg_surv(counts_to_degs(dat_list[[i]])),aes(x=degree, y=surv))+
    geom_line(aes(x=!!fits[[i]]$dat[,1],y=!!fits[[i]]$surv$est*!!ptrunc/fits[[i]]$surv$est[1]), colour='red') + 
    geom_line(aes(x=!!fits[[i]]$dat[,1],y= !!fits[[i]]$surv$CI[1,]*!!ptrunc/fits[[i]]$surv$CI[1,1]),linetype=2, colour='red') + 
    geom_line(aes(x=!!fits[[i]]$dat[,1],y= !!fits[[i]]$surv$CI[2,]*!!ptrunc/fits[[i]]$surv$CI[2,1]),linetype=2, colour='red')+
    ggtitle(nms[i]) + theme(legend.position = 'NULL',aspect.ratio = 1)
}

tmp = fits[[4]]


fits = readRDS('modelfits.rds')

fits[[4]] = tmp

saveRDS(fits, 'modelfits.rds')


# -------------------------------------------------------------------------


# saveRDS(fits,'modelfits.rds')
mixfits = readRDS('mixfits.rds')
dat_list = readRDS('dat_list.rds')
real_survs = ggpubr::ggarrange(plotlist = plots)
real_survs
saveRDS(real_survs,'real_survs.rds')
igpplots = list()

for(i in 1:length(selected)){
  print(i)
  df = as.data.frame(dat_list[[i]])
  names(df) = c('x', 'count')
  df = df[df[1,]>trunc.at,]
  mixfits[[i]] = crandep::mcmc_mix2_wrapper(df,seed=1234L,burn=1e4L)
}


print(i)
i = 2
df = as.data.frame(dat_list[[i]])
names(df) = c('x', 'count')
df = df[df[1,]>trunc.at,]
mixfits[[i]] = crandep::mcmc_mix2_wrapper(df,seed=123L,burn=1e4L)




# -------------------------------------------------------------------------


mixfits  =readRDS('mixfits.rds')

for(i in 1:length(mixfits)){
  udf = as.data.frame(table(mixfits[[i]]$pars$u))
  udf$Var1 = as.numeric(as.character(udf$Var1))
  sdf = twbfn::deg_surv(counts_to_degs(dat_list[[i]]))
  ptrunc = sdf$surv[sdf$degree==min(mixfits[[i]]$data[,1])]
  igpplots[[i]] = ggplot() +
    geom_vline(data=NULL,aes(xintercept = jitter(mixfits[[i]]$pars$u[seq(1,length(mixfits[[i]]$pars$u),by=20)])), colour='blue', alpha=0.5)+
    geom_point(data = twbfn::deg_surv(counts_to_degs(dat_list[[i]])),aes(x=degree, y=surv))+
    geom_line(data = mixfits[[i]]$fitted,aes(x=x, y=S_975*!!ptrunc/S_975[1]), lty=2, colour='red')+
    geom_line(data = mixfits[[i]]$fitted,aes(x=x, y=S_025*!!ptrunc/S_025[1]),lty=2, colour='red')+
    geom_line(data = mixfits[[i]]$fitted,aes(x=x, y=S_med*!!ptrunc/S_med[1]), colour='red')+
    scale_x_log10() + scale_y_log10() + ggtitle(nms[i]) + theme(aspect.ratio = 1)
}

mix_survs = ggpubr::ggarrange(plotlist = igpplots)
mix_survs

saveRDS(mix_survs,'mix_survs.rds')

saveRDS(mixfits,'mixfits.rds')


mixplt = readRDS('mix_survs.rds')

comp_plot = ggpubr::ggarrange(mixplt, real_survs,labels=c('Zipf-IGP', 'GPA'), ncol=1)

saveRDS(comp_plot, 'comp_plot.rds')
fits = readRDS('modelfits.rds')
plots = list()
for(i in 1:length(fits)){
  udf = as.data.frame(table(fits[[i]]$smps$k0))
  udf$Var1 = as.numeric(as.character(udf$Var1))
  sdf = twbfn::deg_surv(counts_to_degs(dat_list[[i]]))
  ptrunc = sdf$surv[sdf$degree==min(fits[[i]]$dat[,1])]
  ptrunc_mix = sdf$surv[sdf$degree==min(mixfits[[i]]$data[,1])]
  plots[[i]] = ggplot() + scale_x_log10(limits = c(1, max(dat_list[[i]][,1]))) + scale_y_log10(limits=c(1/sum(dat_list[[i]][,2]),1)) +
    geom_segment(aes(x = jitter(!!fits[[i]]$smps$k0[seq(1,length(fits[[i]]$smps$k0),by=20)]),yend=0.5,y=1, colour='GPA'), alpha=0.1)+
    geom_segment(data=NULL,aes(x = jitter(mixfits[[i]]$pars$u[seq(1,length(mixfits[[i]]$pars$u),by=20)]),yend=0.5,y=1, colour='Zipf-IGP'), alpha=0.5)+
    geom_point(data = twbfn::deg_surv(counts_to_degs(dat_list[[i]])),aes(x=degree, y=surv))+
    geom_line(aes(x=!!fits[[i]]$dat[,1],y=!!fits[[i]]$surv$est*!!ptrunc/fits[[i]]$surv$est[1], colour='GPA')) + 
    geom_line(aes(x=!!fits[[i]]$dat[,1],y= !!fits[[i]]$surv$CI[1,]*!!ptrunc/fits[[i]]$surv$CI[1,1], colour='GPA'),linetype=2) + 
    geom_line(aes(x=!!fits[[i]]$dat[,1],y= !!fits[[i]]$surv$CI[2,]*!!ptrunc/fits[[i]]$surv$CI[2,1], colour='GPA'),linetype=2)+
    geom_line(data = mixfits[[i]]$fitted,aes(x=x, y=S_975*!!ptrunc_mix/S_975[1], colour='Zipf-IGP'), lty=2)+
    geom_line(data = mixfits[[i]]$fitted,aes(x=x, y=S_025*!!ptrunc_mix/S_025[1], colour='Zipf-IGP'),lty=2)+
    geom_line(data = mixfits[[i]]$fitted,aes(x=x, y=S_med*!!ptrunc_mix/S_med[1], colour='Zipf-IGP'))+
    ggtitle(nms[i]) + theme_bw()+ theme(aspect.ratio = 1, axis.title.x = element_blank(), axis.title.y =element_blank())+ labs(colour='Model')
    
}



comp_plot = ggpubr::ggarrange(plotlist = plots,common.legend = T,label.x = 'Degree', label.y = 'Survival',legend = 'right')
comp_plot

# -------------------------------------------------------------------------

mixfits = readRDS('mixfits.rds')
plots = list()
for(i in 1:length(fits)){
  plots[[i]] = ggplot() + geom_boxplot(aes(y=!!fits[[i]]$mcmc$pars[-(1:2e4),4]/!!fits[[i]]$mcmc$lambdas[-(1:2e4)], x='GPA', fill='GPA')) +
    geom_boxplot(aes(y = !!mixfits[[i]]$pars$shape, x='Zipf-IGP', fill = 'Zipf-IGP')) + ylim(-1,1) + ggtitle(nms[i]) + xlab('')+ylab('')
}

shape_plot = ggpubr::ggarrange(plotlist = plots, common.legend = T,label.x='', label.y='\\xi', legend='right')
shape_plot
saveRDS(shape_plot, 'shape_plot.rds')

# PA plot -----------------------------------------------------------------

i = 1

pref = function(pars,x){
  return(ifelse(x<=pars[3], x^pars[1] + pars[2], pars[3]^pars[1] + pars[2] + pars[4]*(x-pars[3])))
}

plots  = list()
PA_overlay = ggplot()
for(i in 1:length(fits)){
  print(i)
  x = 1:1e3
  pref_mat = apply(fits[[i]]$smps, 1, pref, x = x)
  pref_CI = apply(pref_mat, 1, quantile, prob =c(0.025, 0.5, 0.975))
  
  plots[[i]] = ggplot() + geom_line(aes(x=!!x, y = !!pref_CI[1,]), lty=2)+
    geom_line(aes(x=!!x, y = !!pref_CI[2,]), lty=1)+
    geom_line(aes(x=!!x, y = !!pref_CI[3,]), lty=2) + xlab('Degree') + ylab('Preference')+
    theme(aspect.ratio = 1, axis.title.x = element_blank(), axis.title.y =element_blank()) + ggtitle(nms[i]) + scale_x_log10() + scale_y_log10()
  PA_overlay = PA_overlay + geom_line(aes(!!x, y =!! pref_CI[2,]/sum(!!pref_CI[2,]), colour=!!as.character(nms[i])))
}

PA_plot = ggpubr::ggarrange(plotlist = plots,common.legend = T,label.x = 'Degree', label.y='Preference')
PA_plot

saveRDS(PA_plot, 'PA_plot.rds')


PA_overlay = PA_overlay + scale_x_log10() + scale_y_log10()

saveRDS(PA_overlay, 'PA_overlay.rds')


# -------------------------------------------------------------------------

library(ggridges)

full_par_mat = data.frame(a=NA, eps=NA, k0=NA, b=NA, name=NA)

for(i in 1:length(fits)){
  df = cbind(fits[[i]]$smps, rep(nms[i], nrow(fits[[i]]$smps)))
  names(df) = c('a','eps', 'k0', 'b', 'name')
  full_par_mat = rbind(full_par_mat,df)
}
full_par_mat = full_par_mat[-1,]

meanmat = full_par_mat |> aggregate(.~name, FUN=mean)
aplot = ggplot() + geom_density_ridges2(data = full_par_mat, aes(x=a, y=name, fill=name))+theme(aspect.ratio = 1, axis.title.x = element_blank(), axis.title.y =element_blank()) + ggtitle('a')
epsplot =ggplot() + geom_density_ridges2(data = full_par_mat, aes(x=eps, y=name, fill=name))+theme(aspect.ratio = 1, axis.title.x = element_blank(), axis.title.y =element_blank())+ ggtitle('eps')
k0plot=ggplot() + geom_density_ridges2(data = full_par_mat, aes(x=k0, y=name, fill=name))+theme(aspect.ratio = 1, axis.title.x = element_blank(), axis.title.y =element_blank())+ ggtitle('k0')
bplot = ggplot() + geom_density_ridges2(data = full_par_mat, aes(x=b, y=name, fill=name))+theme(aspect.ratio = 1, axis.title.x = element_blank(), axis.title.y =element_blank())+ ggtitle('b')


parsplot = ggpubr::ggarrange(aplot, epsplot, k0plot, bplot, common.legend = T, legend='none')

parsplot
saveRDS(parsplot, 'pars_plot.rds')




ggsave('pars_plot.png', parsplot)



























































