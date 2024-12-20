w_ij = function(x,i,j){
  return((x + j) ^ (1/i))
}
get_type_weight = function(weight_mat){
  return(colSums(weight_mat))
}
get_node_weight = function(weight_mat){
  return(rowSums(weight_mat))
}
sample_mt = function(n,ntype,init = NULL){
  if(is.null(init)){
    deg_vec = c(0)
    type_vec = c(1)
    weight_mat = matrix(nrow=1,ncol=ntype)
    weight_mat[1,] = w_ij(deg_vec[1],1:ntype,type_vec[1])
  }else{
    deg_vec = init$degree
    type_vec = init$type
    weight_mat = matrix(nrow=length(deg_vec),ncol=ntype)
    for(i in 1:length(deg_vec)){
      weight_mat[i,] = w_ij(deg_vec[i],1:ntype,type_vec[i])
    }
  }
  for(i in 1:n){
    new_type = sample(1:ntype, 1, prob=get_type_weight(weight_mat))
    select_node = sample(1:length(deg_vec), 1, prob=get_node_weight(weight_mat))
    deg_vec[select_node] = deg_vec[select_node] + 1
    type_vec = c(type_vec,new_type)
    deg_vec = c(deg_vec,0)
    weight_mat[select_node,] = w_ij(deg_vec[select_node],1:ntype, type_vec[select_node])
    weight_mat = rbind(weight_mat, w_ij(0,1:ntype, new_type))
  }
  return(list(degree = deg_vec, type=type_vec))
}























