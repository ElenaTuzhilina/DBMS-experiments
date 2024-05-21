
newton = function(param, method, eps, verbose){
  obj = method$loss(param)
  delta = Inf
  while(delta > eps){
    der = method$deriv(param)
    rate = 1
    while(!method$cond(param - rate * der$G/der$H)) rate = rate/2
    param = param - rate * der$G/der$H
    obj0 = obj
    obj = method$loss(param)
    delta = abs((obj0 - obj)/obj0)
    if(verbose) cat("\n param:", param, "rate:", rate, "loss:", obj, "delta:", delta)
  }
  cat("\n")
  return(param)
}

sgd = function(param, method, eps, verbose){
  obj = method$loss(param)
  delta = Inf
  while(delta > eps){
    der = method$deriv(param)
    rate = 1
    while(!method$cond(param - rate * der$G)) rate = rate/2
    param = param - rate * der$G
    obj0 = obj
    obj = method$loss(param)
    delta = abs((obj0 - obj)/obj0)
    if(verbose) cat("\n param:", param, "rate:", rate, "loss:", obj, "delta:", delta)
  }
  cat("\n")
  return(param)
}

######################discrete#########################

pois = list(
  loss = function(X, C, param){
    beta = param$beta
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    logL = -D^2 + beta
    L = exp(logL)
    obj = L - C * logL
    return(mean(obj))
  },
  
  deriv = function(X, C, param, var){
    beta = param$beta
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    logL = -D^2 + beta
    L = exp(logL)
    H = L
    G = C - L
    if(var == "D2") return(list(G = G, H = H))
  },
  
  update_param = function(X, C, param, method, eps, verbose){
    cat("\n update beta")
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    param$beta = log(sum(C)/sum(exp(-D^2)))
    cat("\n\n")
    return(param)
  }
)

hpois = list(
  loss = function(X, C, param){
    beta = param$beta
    p = param$p
    
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    logL = -D^2 + beta
    L = exp(logL)
    
    obj = matrix(0, nrow(C), ncol(C))
    obj[C == 0] = -log(p)
    obj[C != 0] = -log(1 - p) + (L - C * logL + log(1 - exp(-L)))[C != 0]
    
    return(mean(obj))
  },
  
  deriv = function(X, C, param, var){
    beta = param$beta
    
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    logL = -D^2 + beta
    L = exp(logL)
    M = L/(1 - exp(-L))
    
    G = H = matrix(0, nrow(C), ncol(C))
    G[C != 0] = (C - M)[C != 0]
    H[C != 0] =  (M * (1 - M * exp(-L)))[C != 0]
    
    if(var == "D2") return(list(G = G, H = H))
    if(var == "beta") return(list(G = -sum(G), H = sum(H)))
  },
  
  update_param = function(X, C, param, method, eps, verbose){
    cat("\n update p")
    param$p = mean(C == 0)
    
    cat("\n\n update beta")
    method_newton = list(loss = function(beta) method$loss(X, C, list(p = param$p, beta = beta)),
                         deriv = function(beta) method$deriv(X, C, list(p = param$p, beta = beta), "beta"),
                         cond = function(beta) TRUE)
    param$beta = newton(param$beta, method_newton, eps, verbose)
    cat("\n\n")
    return(param)
  }
)


zipois = list(
  loss = function(X, C, param){
    p = param$p
    beta = param$beta
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    logL = -D^2 + beta
    L = exp(logL)
    a = p /(1 - p)
    obj = matrix(0, nrow(C), ncol(C))
    obj[C == 0] = -log(a + exp(-L))[C == 0] - log(1 - p)
    #obj[C == 0] = -log(a + 1/exp(L))[C == 0] - log(1 - p)
    obj[C != 0] = (L - C * logL)[C != 0] - log(1 - p)
    return(mean(obj))
  },
  deriv = function(X, C, param, var){
    p = param$p
    beta = param$beta
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    logL = -D^2 + beta
    L = exp(logL)
    expL = exp(L)
    a = p/(1 - p)  
    
    G = H = matrix(0, nrow(C), ncol(C))
    
    if(var %in% c("D2", "beta")){
      G[C == 0] = (-L / (a * expL + 1))[C == 0]
      G[C != 0] = (C - L)[C != 0]
      #H[C == 0] =  (L * (a * expL * (1 - L) + 1)/(a * expL + 1)^2)[C == 0]
      H[C == 0] =  (L * (a * (1 - L) + 1 / expL)/(a * expL + 1)/(a + 1/expL))[C == 0]
      H[C != 0] = L[C != 0]
      if(var == "D2") return(list(G = G, H = H))
      else return(list(G = -sum(G), H = sum(H)))
    }
    if(var == "p"){
      #G[C == 0] = (-(expL - 1) / (p * (expL - 1) + 1))[C == 0]
      G[C == 0] = (-1 / (p + 1/(expL - 1)))[C == 0]
      G[C != 0] = 1/(1 - p)
      #H[C == 0] =  ((expL - 1)^2 / (p * (expL - 1) + 1)^2)[C == 0]
      H[C == 0] =  (1 / (p + 1/(expL - 1))^2)[C == 0]
      H[C != 0] = 1/(1 - p)^2
      return(list(G = sum(G), H = sum(H)))
    }
  },
  update_param = function(X, C, param, method, eps, verbose){
    cat("\n update beta")
    method_newton = list(loss = function(beta) method$loss(X, C, list(p = param$p, beta = beta)),
                         deriv = function(beta) method$deriv(X, C, list(p = param$p, beta = beta), "beta"), 
                         cond = function(beta) TRUE)
    param$beta = newton(param$beta, method_newton, eps, verbose)
    
    cat("\n update p")
    method_newton = list(loss = function(p) method$loss(X, C, list(p = p, beta = param$beta)),
                         deriv = function(p) method$deriv(X, C, list(p = p, beta = param$beta), "p"), 
                         cond = function(p) p > 0 & p < 1)
    param$p = newton(param$p, method_newton, eps, verbose)
    cat("\n\n")
    return(param)
  }
)

nbin = list(
  loss = function(X, C, param){
    r = param$r
    beta = param$beta
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    logL = -D^2 + beta
    L = exp(logL)
    obj =  lgamma(r) - lgamma(C + r) - r * log(r)  + (C + r) * log(r + L) - C * logL
    return(mean(obj))
  },
  
  deriv = function(X, C, param, var){
    r = param$r
    beta = param$beta
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    logL = -D^2 + beta
    L = exp(logL)
    if(var %in% c("D2", "beta")){
      G = r * (C - L) / (L + r)
      H = r * L * (C + r)/(L + r)^2
      if(var == "D2") return(list(G = G, H = H))
      if(var == "beta") return(list(G = -sum(G), H = sum(H)))
    }
    if(var == "r"){
      G = digamma(r) - digamma(C + r) - log(r) + log(L + r) + (C - L)/(L + r)
      H = trigamma(r) - trigamma(C + r) - 1/r + 1/(L + r) - (C - L)/(L + r)^2
      return(list(G = sum(G), H = sum(H)))
    }
  },
  
  update_param = function(X, C, param, method, eps, verbose){
    method_newton = list(loss = function(beta) method$loss(X, C, list(r = param$r, beta = beta)),
                         deriv = function(beta) method$deriv(X, C, list(r = param$r, beta = beta), "beta"),
                         cond = function(beta) TRUE)
    cat("\n update beta")
    param$beta = newton(param$beta, method_newton, eps, verbose)
    
    method_newton = list(loss = function(r) method$loss(X, C, list(r = r, beta = param$beta)),
                         deriv = function(r) method$deriv(X, C, list(r = r, beta = param$beta), "r"), 
                         cond = function(r) r > 0)
    cat("\n update r")
    param$r = sgd(param$r, method_newton, eps, verbose)
    
    cat("\n\n")
    return(param)
  }
)

##################continuous########################

norm = list(
  loss = function(X, C, param){
    beta = param$beta
    sigma2 = param$sigma2
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    M = -D^2 + beta
    obj = log(2*pi*sigma2) + (C - M)^2/sigma2
    return(1/2*mean(obj))
  },
  
  deriv = function(X, C, param, var){
    beta = param$beta
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    M = -D^2 + beta
    H = matrix(1, nrow(X), nrow(X))
    G = -M + C
    if(var == "D2") return(list(G = G, H = H))
  },
  
  update_param = function(X, C, param, method, eps, verbose){
    cat("\n update beta")
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    param$beta = mean(D^2 + C)
    cat("\n update sigma")
    param$sigma2 = mean((C + D^2 - param$beta)^2)
    cat("\n\n")
    return(param)
  }
)

expon = list(
  loss = function(X, C, param){
    beta = param$beta
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    logL = D^2 - beta
    L = exp(logL)
    obj = C * L - logL
    return(mean(obj))
  },
  
  deriv = function(X, C, param, var){
    beta = param$beta
    D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
    logL = D^2 - beta
    L = exp(logL)
    H = C * L
    G = C * L - 1
    if(var == "D2") return(list(G = G, H = H))
    else return(list(G = -sum(G), H = sum(H)))
  },
  
  update_param = function(X, C, param, method, eps, verbose){
    method_newton = list(loss = function(beta) method$loss(X, C, list(beta = beta)),
                         deriv = function(beta) method$deriv(X, C, list(beta = beta), "beta"),
                         cond = function(beta) TRUE)
    cat("\n update beta")
    param$beta = newton(param$beta, method_newton, eps, verbose)

    cat("\n\n")
    return(param)
  }
)
