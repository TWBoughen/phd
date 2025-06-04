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

