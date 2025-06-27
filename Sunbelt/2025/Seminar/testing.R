library(tidyverse)
library(igraph)
library(ggraph)
library(gifski)
library(gganimate)
library(ggforce)
library(qgraph)



g = sample_pa(1e1)


pal = paletteer::paletteer_c("grDevices::Heat 2", 1e1)

gplots = qgraph.animate(as_adjacency_matrix(g), color=pal, labels=F, plotGraphs = F, smooth = T)




# -------------------------------------------------------------------------

library(igraph)
library(network)
library(networkDynamic)
library(ndtv)
library(intergraph)
g_igraph <- sample_pa(20, directed = FALSE)
el = as_edgelist(g_igraph)
g_net <- network.initialize(20)

add.edges.active(g_net, tail = el[,2], head=el[,1], onset = 2:(nrow(el)+1), terminus = rep(Inf, nrow(el)))
activate.vertices(g_net, onset=1:20, terminus = rep(Inf,20))
slice.par<-list(start=1,end=20,interval=1,
                  aggregate.dur=1,rule='earliest')
compute.animation(g_net,slice.par=slice.par )
reconcile.vertex.activity(g_net)
reconcile.edge.activity(g_net)
render.d3movie(g_net, output.mode = 'inline', script.type = 'remoteSrc')

# -------------------------------------------------------------------------

library(igraph)
library(network)
library(networkDynamic)
library(ndtv)
library(intergraph)
g_igraph <- sample_pa(20, directed = FALSE)
el = as_edgelist(g_igraph)
g_net <- network.initialize(20)

add.edges.active(g_net, tail = el[,2], head=el[,1], onset = 2:(nrow(el)+1), terminus = rep(Inf, nrow(el)))
activate.vertices(g_net, onset=1:20, terminus = rep(Inf,20))
slice.par<-list(start=1,end=20,interval=1,
                aggregate.dur=1,rule='earliest')
compute.animation(g_net,slice.par=slice.par )
reconcile.vertex.activity(g_net)
reconcile.edge.activity(g_net)
render.d3movie(g_net, file.name='new.html',
               verbose = TRUE, output.mode = 'inline',launchBrowser = F, script.type = 'embedded')


# -------------------------------------------------------------------------
source('scripts/funcs.R')
fits = readRDS('results/modelfits.rds')

pref = function(pars,x){
  return(ifelse(x<=pars[3], x^pars[1] + pars[2], pars[3]^pars[1] + pars[2] + pars[4]*(x-pars[3])))
}

selected = c(1,8,10)
plts = list()
j=3
for(j in 1:length(selected)){
  i = selected[j]
  pars = fits[[i]]$smps
  x = fits[[i]]$dat$x
  
  plots  = list()
  PA_overlay = ggplot()
  
  pref_mat = apply(pars, 1, pref, x = x)
  pref_CI = apply(pref_mat, 1, quantile, prob =c(0.025, 0.5, 0.975))
  
  plts[[j]] = ggplot() +
    geom_line(aes(x=!!x, y=!!pref_CI[3,]), lty=2)+
    geom_line(aes(x=!!x, y=!!pref_CI[2,]))+
    geom_line(aes(x=!!x, y=!!pref_CI[1,]), lty=2)+
    theme(aspect.ratio = 1, axis.title.x = element_blank(), axis.title.y =element_blank(), legend.position = 'none',
          text = element_text(size = 20),plot.background = element_rect(fill="#f1f1f1", color = NA),
          panel.background = element_rect(fill="#f1f1f1", color = NA),
          plot.margin=grid::unit(c(0,0,0,0), "mm"))
}














