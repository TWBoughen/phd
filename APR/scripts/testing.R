library(igraph)
library(visNetwork)

g = sample_pa(2, 1, 1, directed = F)
# -------------------------------------------------------------------------

old = as_edgelist(g)

g = sample_pa(3,1,m=2,directed = F, start.graph = g, zero.appeal = 0.0000001)


combined = 
