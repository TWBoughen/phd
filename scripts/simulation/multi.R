# -------------------------------------------------------------------------

library(Rcpp)
sourceCpp('test.cpp')

theta <- 3
i <- 2
j <- 1
ntype <- 4
x_max <- 1000  # Truncation point; adjust based on convergence analysis

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

opt = optim(3,perron_optim, ntype=2,method='Brent', lower = 0.1, upper=15)
alpha = opt$par
u = compute_mat_eig(alpha, ntype)$vectors[,1]
y1 = y2 = y3 = c()
n = 1e3
for(x in 0:(n)){
  y1 = c(y1,I_i(x,alpha,1,2) * (u[1]/sum(u)))
  y2 = c(y2,I_i(x,alpha,2,2) * (u[2]/sum(u)))
  y3 = c(y2,I_i(x,alpha,3,2) * (u[3]/sum(u)))
}
y1 = alpha*y1
y2 = alpha*y3
# out = sim_mt(5e3, quiet=F)

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
# 
# points(twbfn::deg_surv(out$d + 1), pch=20, col=1)
# points(twbfn::deg_surv(out$d[out$t==1] + 1),pch=20,col=2)
# points(twbfn::deg_surv(out$d[out$t==2] + 1),pch=20,col=3)

#simulating --------------------------------------------------------------

type_weight = function(i, degrees, types, ntype=2){
  s  = 0
  for(j in 1:length(degrees)){
    s = s + w_ij(degrees[j],types[j],i )
  }
  return(s)
}
node_weight = function(degrees, types, ntype=2){
  vec = c()
  for(i in 1:length(degrees)){
    vec=c(vec, w_i(degrees[i],types[i],ntype))
  }
  return(vec)
}

get_weights = function(degrees, types, ntype=2){
  type_weights = numeric(ntype)
  node_weights = numeric(length(degrees))
  
  
  
  
   
}







sim_mt = function(n,ntype=2,quiet=T){
  types = c(1,2)
  degrees = c(1,1)
  for(i in 1:n){
    if(!quiet) message(i)
    #decide on type
    w_vec = c()
    for(type in 1:ntype){
      w_vec = c(w_vec, type_weight(type, degrees, types))
    }
    new_type = sample(1:ntype, 1, prob=w_vec)
    #decide what to attach to
    attach_to = sample(1:length(degrees),1,prob = node_weight(degrees, types))
    degrees[attach_to] = degrees[attach_to] + 1
    degrees = c(degrees,0)
    types = c(types, new_type)
  }
  return(list(d=degrees, t=types))
}




































