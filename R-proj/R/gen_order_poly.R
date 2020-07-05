gen_order_poly <- function(path) {
  
  G = scan(path, character(), quote = "")
  G = as.numeric(G)
  G = matrix(G, nrow = sqrt(length(G)), byrow = TRUE)
  
  n = dim(G)[1]
  if (length(which(G %in% 1)) == 0) {
    P = gen_cube(n,'H')
    return(P)
  }
  
  E = matrix(0,0,2)
  A = matrix(0,0,n)
  for (i in 1:n) {
    
    r = G[i,]
    indices = which(r %in% 1)
    
    if (length(indices) == 0) {
      next
    }
    
    for (j in 1:length(indices)) {
      rowA = rep(0,n)
      rowA[i] = -1
      rowA[indices[j]] = 1
      A = rbind(A, c(rowA))
      E = rbind(E, c(i,indices[j]))
    }
    
  }
  A = rbind(rbind(diag(n), -diag(n)), A)
  P = OrderPolytope$new(E, A)
  return(P)
  
}