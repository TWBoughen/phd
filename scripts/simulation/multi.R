# -------------------------------------------------------------------------

library(Rcpp)
compute_mat_eig = function(theta, ntype,xmax=1e3){
  mat = matrix(nrow=ntype, ncol=ntype)
  for(i in 1:ntype){
    for(j in 1:ntype){
      mat[i,j] = compute_sum(theta, i, j, ntype, x_max)
    }
  }
  eig = eigen(mat)
  return(eig)
}
perron_optim = Vectorize(function(theta, ntype,xmax=1e3){
  mat = matrix(nrow=ntype, ncol=ntype)
  for(i in 1:ntype){
    for(j in 1:ntype){
      mat[i,j] = compute_sum(theta, i, j, ntype, x_max)
    }
  }
  eig = eigen(mat)
  return(abs(eig$values[1]-1))
},vectorize.args='theta')
set_pref = function(pref_string){
  w_ij_full_string = paste0("double w_ij(double x, int i, int j){ return ",pref_string,';}')
  contents = readLines('nopref.cpp')
  new_contents = c('#include <Rcpp.h>','using namespace Rcpp;','// [[Rcpp::export]]',w_ij_full_string,contents)
  writeLines(new_contents, 'mt_fun.cpp')
  sourceCpp('mt_fun.cpp')
}
{
# type_weight = function(i, degrees, types, ntype=2){
#   s  = 0
#   for(j in 1:length(degrees)){
#     s = s + w_ij(degrees[j],types[j],i )
#   }
#   return(s)
# }
# node_weight = function(degrees, types, ntype=2){
#   vec = c()
#   for(i in 1:length(degrees)){
#     vec=c(vec, w_i(degrees[i],types[i],ntype))
#   }
#   return(vec)
# }
}
get_weights = function(degrees, types, ntype=2){
  node_weights = c()
  type_weights = numeric(ntype)
  
  for(i in 1:length(degrees)){
    node_weights = c(node_weights,w_i(degrees[i],types[i],ntype))
    for(j in 1:ntype){
      type_weights[j] = type_weights[j] + w_ij(degrees[i],types[i],j )
    }
  }
  return(list(type = type_weights, node=node_weights))
}
sim_mt = function(n,ntype=2,quiet=T){
  types = c(1,2)
  degrees = c(1,1)
  for(i in 1:n){
    if(!quiet) message(i)
    #decide on type
    weights = get_weights(degrees, types, ntype)
    new_type = sample(1:ntype, 1, prob=weights$type)
    #decide what to attach to
    attach_to = sample(1:length(degrees),1,prob = weights$node)
    degrees[attach_to] = degrees[attach_to] + 1
    degrees = c(degrees,0)
    types = c(types, new_type)
  }
  return(list(d=degrees, t=types))
}

# NEED TO IMPLEMENT SIMULATION IN RCPP TO IMPROVE SPEED


# -------------------------------------------------------------------------

#w_ij_string must be written using valid Rcpp functions  

ntype <- 2
w_ij_string = "pow(i*x + 1.0,1.0)"
set_pref(w_ij_string)

opt = optim(3,perron_optim, ntype=ntype,method='Brent', lower = 0.1, upper=15)
alpha = opt$par
u = compute_mat_eig(alpha, ntype)$vectors[,1]

#conditional survivals
n=1000
S = matrix(ncol = n, nrow=ntype)
S_tot = matrix(ncol=n, nrow=1)
for(x in 0:(n-1)){
  for(i in 1:ntype){
    S[i,x+1] = S_i(x, alpha, i, ntype)
  }
  S_tot[x+1] = sum(u/sum(u) * S[,x+1])
}
plot(0:(n-1) + 1 , S_tot, log='xy', type='l', col=1, ylim=c(1e-6,1))
for(i in 1:ntype){
  lines(0:(n-1) + 1, S[i,], col=i+1)
}
legend('topright', legend = paste0(c(1:ntype),' : ',round(u /sum(u),2)), fill = 1:ntype + 1)
out = sim_mt(5e3, ntype, quiet=F)
points(twbfn::deg_surv(out$d + 1), pch=20, col=1)
points(twbfn::deg_surv(out$d[out$t==1] + 1),pch=20,col=2)
points(twbfn::deg_surv(out$d[out$t==2] + 1),pch=20,col=3)







































