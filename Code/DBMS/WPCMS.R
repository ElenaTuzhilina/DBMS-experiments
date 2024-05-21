WPCMS = function(Z, H, W = matrix(1, nrow(Z), ncol(Z)), 
                 beta = -min(Z), Theta = matrix(rnorm(ncol(H) * 3), ncol(H), 3), 
                 update_beta = TRUE, eps = 1e-6, maxiter = 100, verbose = FALSE){
  #Initialize
  X = H %*% Theta
  X = scale(X, scale = FALSE, center = TRUE)
  obj = loss_WPCMS(X, Z, W, beta)
  
  delta = Inf
  iter = 0
  
  info = c(iter, 0, beta, obj, delta)
  info_names = c("iter", "rate", "beta", "objective", "delta")
  if(verbose) cat(paste(info_names, ":", info), "\n")
  
  #Iterate
  while(delta > eps && iter < maxiter){
    iter = iter + 1
    obj0 = obj
    
    #Line search
    find = line_search_WPCMS(X, Z, W, H, beta, obj0)
    
    #Update solution
    if(find$rate > 0) Theta = find$Theta
    X = H %*% Theta
    X = scale(X, scale = FALSE, center = TRUE)
    
    #Update beta
    if(update_beta) beta = update_param_WPCMS(X, Z, W)
    
    #Update loss
    obj = loss_WPCMS(X, Z, W, beta)
    delta = abs((obj0 - obj)/obj0)
    info = rbind(info, c(iter, find$rate, beta, obj, delta))
    
    if(verbose) cat(paste(info_names, ":", info[iter + 1,]), "\n")
  }
  info = data.frame(info)
  names(info) = info_names
  info$df = ncol(H)
  return(list(Theta = Theta, X = X, beta = beta, info = info, iter = iter, obj = obj))
}

loss_WPCMS = function(X, Z, W, beta){
  D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
  obj = W * (Z + beta - D^2)^2
  return(mean(obj))
}

deriv_WPCMS = function(X, Z, W, beta){
  D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
  G = W * (Z + beta - D^2)
  G_plus = diag(rowSums(G))
  return(G - G_plus)
}

line_search_WPCMS = function(X, Z, W, H, beta, obj0){
  S = X %*% t(X)
  G = deriv_WPCMS(X, Z, W, beta)
  
  rate = 2
  obj = Inf
  rank = 0
  
  while(obj0 < obj || rank < 3){
    rate = rate / 2
    pcms = PCMS(S - rate * G, H)
    obj = loss_WPCMS(H %*% pcms$Theta, Z, W, beta)
    rank = pcms$rank
    if(rate < 1e-20) return(list(Theta = NA, rate = 0))
  }
  return(list(Theta = pcms$Theta, rate = rate))
}

update_param_WPCMS = function(X, Z, W){
  D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
  beta = -sum(W * (Z - D^2))/sum(W)
  return(beta)
}