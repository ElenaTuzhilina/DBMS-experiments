DBMS = function(C, H, param, method, Theta = matrix(rnorm(ncol(H) * 3), ncol(H), 3), update_param = TRUE, eps_wpcms = 1e-6, maxiter = 100, verbose_wpcms = FALSE, eps_dbms = 1e-6, maxepoch = 100, verbose_dbms = FALSE){
  #Initialize
  X = H %*% Theta
  X = scale(X, scale = FALSE, center = TRUE)
  obj = method$loss(X, C, param)
  
  delta = Inf
  deltaX = Inf
  iter_total = 0
  epoch = 0
  
  info = c(epoch, 0, iter_total, unlist(param), obj, delta, deltaX)
  info_names = c("epoch", "rate", "iter_total", names(param), "objective", "delta", "deltaX")
  if(verbose_dbms) cat(paste(info_names, ":", info), "\n")
  
  #Iterate
  while(deltaX > eps_dbms && epoch < maxepoch){
    epoch = epoch + 1
    obj0 = obj
    X0 = X
    
    #SOA
    soa = SOA(X, method$deriv(X, C, param, "D2"))
    W = soa$W
    Z = soa$Z
    
    #WPCMS (beta = 0)
    wpcms = WPCMS(Z, H, W, 0, Theta, FALSE, eps_wpcms, maxiter, verbose_wpcms)
    iter_total = iter_total + wpcms$iter
    
    #Line search
    find = line_search(X, wpcms$X, C, H, param, obj, method$loss)
    
    #Update solution
    if(find$rate > 0) Theta = find$Theta
    X = H %*% Theta
    X = scale(X, scale = FALSE, center = TRUE)
    
    #Update parameters
    if(update_param) param = method$update_param(X, C, param, method, eps_dbms, verbose_dbms)
    
    #Update loss
    obj = method$loss(X, C, param)
    delta = abs((obj0 - obj)/obj0)
    deltaX = vegan::procrustes(X0, X, scale = TRUE, symmetric = TRUE)$ss
    info = rbind(info, c(epoch, find$rate, iter_total, unlist(param), obj, delta, deltaX))
    
    if(verbose_dbms) cat(paste(info_names, ":", info[epoch + 1,]), "\n")
  }
  info = data.frame(info)
  names(info) = info_names
  info$df = ncol(H) 
  return(list(Theta = Theta, X = X, param = param, info = info, epoch = epoch, iter_total = iter_total, obj = obj))
}

############################################

line_search = function(X0, X, C, H, param, obj0, loss){
  S0 = X0 %*% t(X0)
  S = X %*% t(X)
  
  rate = 2
  obj = Inf
  rank = 0
  
  while(obj0 < obj || rank < 3){
    rate = rate / 2
    pcms = PCMS((1 - rate) * S0 + rate * S, H)
    obj = loss(H %*% pcms$Theta, C, param)
    rank = pcms$rank
    if(rate < 1e-20) return(list(Theta = NA, rate = 0))
  }
  return(list(Theta = pcms$Theta, rate = rate))
}

############################################

SOA = function(X, deriv){
  D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
  G = deriv$G
  H = deriv$H
  W = H
  diag(W) = 0
  W = W/max(W)
  Z = D^2 - G/H
  Z[W == 0] = 0
  return(list(Z = Z, W = W))
}

