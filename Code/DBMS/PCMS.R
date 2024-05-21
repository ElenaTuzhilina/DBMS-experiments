
PCMS = function(Z, H){
  ED3 = rARPACK::eigs_sym(t(H) %*% Z %*% H, 3, which = "LA")
  U3 = ED3$vectors
  d3 = ED3$values
  rank = sum(d3 > 1e-12)
  d3 = pmax(d3, 0)
  Theta = U3 %*% diag(sqrt(d3))
  X = H %*% Theta
  return(list(Theta = Theta, loss = loss_PCMS(X, Z), rank = rank))
}

loss_PCMS = function(X, Z){
  obj = (Z - X %*% t(X))^2
  return(mean(obj))
}
