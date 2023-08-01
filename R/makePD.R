#' Positive Definite (PD) approximation of a Matrix
#'
#' @description The function uses Eigen values and Eigen vectors to create a PD matrix
#'
#' @param mat Input a matrix
#'
#' @return A modified PD matrix
#' @export

makePD = function(mat){
  N = nrow(mat)
  HC = mat
  D = eigen(mat)
  E = D$values
  U = D$vectors
  v = as.numeric(E < 0)
  idx = which(v==1)
  m = sum(v) # number of negative values
  if(m > 0){
    S = sum(v*E)*2
    W = (S*S*100)+1
    P = min(abs(E[idx])) # smallest positive value
    for(i in idx){
      C = E[i]
      E[i] = P * (S-C)*(S-C)/W
    }
  }
  return(E)
}