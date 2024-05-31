
# -------------------------------------------------------------------------
library(plot3D)
library(ggplot2)
library(ggallin)
res=readRDS('../data/mcmc.outputs/dir_out.rds')


res2 = res[,2]$res
res2 = res2[-(1:1e4),]
library(MASS)
library(fishualize)
n=635
pal <- viridis_pal()(n)
dens = kde2d(res2$xi, res2$sig)
filled.contour(dens,col=pal,nlevels = 400)



ggplot(res2[,c('xi','sig')], aes(x=xi, y=sig,fill = ..level..)) +stat_density_2d(geom = "polygon")


# plot(density(df[df$from=='BAsim',]$xi))

