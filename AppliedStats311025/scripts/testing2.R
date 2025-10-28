source('phd/AppliedStats311025/scripts/funcs.R')


as = seq(0,2,l=50)
eps = seq(0.1,2,l=50)
k0 = 20
bs = seq(0.1,2, l=50)


df =expand.grid(as, eps, k0, bs)
df$lambda = 0

for(r in 1:nrow(df)){
  print(r)
  est = as.numeric(df[r,])
  df$lambda[r] = find_lambda2(polylin(est[1], est[2]),est[4],est[3])
}









