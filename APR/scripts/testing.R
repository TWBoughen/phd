library(igraph)
library(visNetwork)
library(wdnet)
g = sample_pa(1, 1, 1, directed = F)
i=1
# -------------------------------------------------------------------------
i=i+1
old = as_edgelist(g)
print((old))

g = sample_pa(i,1,m=1,directed = F, start.graph = g, zero.appeal = 0.0000001)
new_edges = as.data.frame(as_edgelist(g))
names(new_edges) = c('from', 'to')
new_edges$id = get.edge.ids(g,as.vector(t(as_edgelist(g))))*1000
new_edges$label = get.edge.ids(g,as.vector(t(as_edgelist(g))))*1000


# -------------------------------------------------------------------------

ctrl_ba = rpa_control_preference(ftype='customized', pref = 's', spref='outs', tpref='ins')+ rpa_control_scenario(alpha=1, beta=0, gamma=0, xi=0, rho=0, beta.loop=TRUE,source.first = TRUE)
init_net = list(
  edgelist = matrix(c(1, 2), nrow = 1),
  edgeweight = 1,
  directed =FALSE
) 

g1 = rpanet(1, init_net, control=ctrl_ba)
for(i in 1:1e5){
  print(i)
  g1 = rpanet(0, g1, control=ctrl_ba)
}
k = 1:200
m=1
pk_ba = 2*m*(m+1) / (k*(k+1)*(k+2))

g <- wdnet_to_igraph(g1)

x = degree(g)
# print(sort(unique(x)))
p = table(x)/ length(x)
print(names(p))
par(mar=c(2,2,0,0))
# plot(sort(unique(x)),as.numeric(p), log='xy', xlab='', ylab='')
q = c(1,(1-cumsum(as.numeric(p)/sum(as.numeric(p))))[-length(p)])
plot(names(p),q, log='xy', xlab='Degree (x)', ylab='1-F(x)')
lines(k, c(1,1-cumsum(pk_ba)[-length(pk_ba)]))
df = data.frame(p=p ,q=q)
degp = ggplot(data=df, aes(x=as.numeric(p.x), y=q)) + geom_point() + scale_x_log10() + scale_y_log10()
#-------------------------------------------------------------------------
library(ggnetwork)
library(GGally)
library(network)
library(sna)
library(intergraph)


par(mfrow=c(1,1))
netplot <-recordPlot() 



netplot <- ggnet2(g, size='degree', legend.position = '')

print(grid.arrange(degp, netplot, nrow=2))


source('igp_functions.R')
# -------------------------------------------------------------------------



v = 50
a=2
xi=2
sig=50
phi=1-p_mix(v,1e4,0,a,xi,sig)

q = p_mix(k,v,phi,a,xi,sig)



plot(k,1-q,log='xy', type='l')





