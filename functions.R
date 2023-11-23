library(mda)
library(ggplot2)
library(fields)
library(dplyr)
library(splines)
library(plotly)
library(vegan)
library(matlib)

#create orthogonal basis
load_H = function(df, index, orth = TRUE){
  n_knots = df - 2
  knots = unique(seq(from = min(index), to = max(index), length = n_knots))
  knots = knots[-c(1,n_knots)]
  H = bs(index, knots = knots, intercept = TRUE)
  if(orth) H = qr.Q(qr(H))
  return(H)
}

#compute Demmler-Reinch basis and K
DR = function(index){
  n = length(index)
  h = index[-1] - index[-length(index)]
  Q = matrix(0, n, n-2)
  diag(Q) = 1/h[-length(h)]
  diag(Q[-1,]) = - 1/h[-1] - 1/h[-length(h)]
  diag(Q[-(1:2),]) = 1/h[-1]
  R = matrix(0, n-2, n-2)
  diag(R) = (h[-length(h)] + h[-1])/3
  diag(R[,-1]) = diag(R[-1,]) = h[-c(1,length(h))]/6
  K = Q %*% solve(R) %*% t(Q)
  ED = eigen(K)
  s = ED$values
  H = ED$vectors[,order(s)]
  return(list(K = K, H = H, s = s))
}

#inner products to distances
StoD2 = function(S){
  D2 = 2*S
  D2 = scale(D2, center = diag(S), scale = FALSE)
  D2 = t(scale(t(D2) , center = diag(S), scale = FALSE))
  return(-D2)
}

#compute Lambda matrix
Lambda = function(X, beta){
  D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
  logL = -D^2 + beta
  L = exp(logL)
  return(L)
}

#compute expected counts for DBMS
expected = function(X, param, method_name){
  L = Lambda(X, param$beta)
  if(method_name %in% c("PoisMS", "NBMS")) E = L
  if(method_name == "HPoisMS") E = (1 - param$p) * L/(1 - exp(-L))
  if(method_name == "ZIPoisMS") E = (1 - param$p) * L
  return(E)
}

#compute different metrics for DBMS solution: lokglik, mse, correlation
dbms_evaluate = function(X, param, C, method, method_name){
  L = Lambda(X, param$beta)
  E = expected(X, param, method_name)
  data.frame(loglik = method$loss(X, C, param),
             corE = cor(c(C), c(E)),
             corL = cor(c(C), c(L)),
             spcorE = cor(c(C), c(E), method = "spearman"),
             spcorL = cor(c(C), c(L), method = "spearman"),
             strcorE = get.scc(C, E, resol = 50000, h = 5)$scc,
             strcorL = get.scc(C, L, resol = 50000, h = 5)$scc,
             mseE = mean((C-E)^2),
             mseL = mean((C-L)^2))
}

#find reconstruction curvature
curvature = function(X){
  angles = numeric(nrow(X) - 2)
  for(i in 3:nrow(X)){
    e1 = X[i - 1,] - X[i - 2,]
    e2 = X[i,] - X[i - 1,]
    angles[i] = angle(e1, e2)
  }
  angles
}

#align two reconstructions using procrustes
align = function(X, Y){
  X = scale(X, center = TRUE, scale = FALSE)
  Y = scale(Y, center = TRUE, scale = FALSE)
  flip = TRUE
  while(flip){
    Y = procrustes(X, Y, scale = FALSE)$Yrot
    refl = diag(sign(cor(X, Y)))
    Y = scale(Y, center = FALSE, scale = refl)
    flip = (sum(refl < 0) > 0)
  }
  return(list(X = X, Y = Y))
}

#align two reconstructions using procrustes and pre-smoothing

smooth_X = function(X, width){
  Xsmooth = data.frame(X) %>% 
    mutate_all(~rollapply(., width = width, FUN = function(x) mean(x, na.rm=TRUE), by = 1, by.column = TRUE, partial = TRUE, fill = NA, align = "center")) %>%
    as.matrix()
}

smooth_align = function(X, Y, width, smoothX = TRUE, smoothY = TRUE){
  if(smoothX) Xsmooth = smooth_X(X, width)
  else Xsmooth = X
  indexX = which(rowSums(is.na(Xsmooth)) == 0)
  if(smoothY) Ysmooth = smooth_X(Y, width)
  else Ysmooth = Y
  indexY = which(rowSums(is.na(Ysmooth)) == 0)
  index = intersect(indexX, indexY)
  Xsmooth = Xsmooth[index,]
  Ysmooth = Ysmooth[index,]
  mX = colMeans(Xsmooth)
  mY = colMeans(Ysmooth)
  Xsmooth = scale(Xsmooth, center = mX, scale = FALSE)
  Ysmooth = scale(Ysmooth, center = mY, scale = FALSE)
  X = scale(X, center = mX, scale = FALSE)
  Y = scale(Y, center = mY, scale = FALSE)
  flip = TRUE
  while(flip){
    proc = procrustes(Xsmooth, Ysmooth, scale = TRUE)
    R = proc$rotation
    Ysmooth = Ysmooth %*% R
    Y = Y %*% R
    refl = sign(diag(cor(Xsmooth, Ysmooth)))
    Ysmooth = scale(Ysmooth, center = FALSE, scale = refl)
    Y = scale(Y, center = FALSE, scale = refl)
    flip = (sum(refl < 0) > 0)
  }
  return(list(X = X, Y = Y, ss = mean((X - Y)^2, na.rm = T), ssmooth = mean((Xsmooth - Ysmooth)^2)))
}

#reduce resolution of HiC data
low_resolution = function(C, n_avg){    
  n = ncol(C)
  size = ceiling(n/n_avg)
  split = (1:n) %/% size + 1
  C_avg = matrix(0, n_avg, n_avg)
  for(i in 1:n_avg){
    for(j in 1:n_avg){
      C_avg[i,j] = sum(C[split == i, split == j])
    }    
  }
  return(C_avg)
}

#plot reconstruction
visualize = function(X, type = 'projection', index = 1:nrow(X), title = NULL){
  n = nrow(X)
  before_centromere = which(index < (n * 0.45))
  after_centromere = which(index >= (n * 0.45))
  col = c(rep('orange', length(before_centromere)), rep('darkturquoise', length(after_centromere)))
  colnames(X) = c('x', 'y', 'z')
  par(mfrow = c(1,1), oma = c(0, 0, 2, 0))
  panelf = function(x, y){
    col = c(rep('orange', length(before_centromere)), rep('darkturquoise', length(after_centromere)))
    points(x, y, pch = 19, cex = 1, col = col)
    lines(x, y, col = 'orange', lwd = 2)
    lines(x[after_centromere], y[after_centromere], col = 'darkturquoise', lwd = 2)
  }
  if(type == 'projection') return(pairs(X, panel = panelf, cex.labels = 5, main = title))
  
  if(type == '3D'){
    ax = list(title = " ", nticks = 8, range = c(min(X[,1]), max(X[,1])), showticklabels = FALSE)
    ay = list(title = " ", nticks = 8, range = c(min(X[,2]), max(X[,2])), showticklabels = FALSE)
    az = list(title = " ", nticks = 8, range = c(min(X[,3]), max(X[,3])), showticklabels = FALSE)
    return(plot_ly(x = X[,1], y = X[,2], z = X[,3], type = 'scatter3d', mode = 'lines+markers',
                   line = list(width = 6, color = col), marker = list(size = 3.5, color = col)) %>% 
             layout(scene = list(xaxis = ax, yaxis = ay, zaxis = az)))
    
  }
}

#plot multiple reconstructions together
plot_3D = function(Xs, color, name){
  if(!is.list(Xs)){
    nconf = 1
    X = Xs
  } 
  else{
    nconf = length(Xs)
    X = Xs[[1]]
  } 
  ax = list(title = " ", nticks = 8, showticklabels = FALSE)
  ay = list(title = " ", nticks = 8, showticklabels = FALSE)
  az = list(title = " ", nticks = 8, showticklabels = FALSE)
  plt = plot_ly(x = X[,1], y = X[,2], z = X[,3], type = 'scatter3d', mode = 'lines+markers', name = name[1],
                line = list(width = 6, color = color[1]), marker = list(size = 3.5, color = color[1])) %>% 
    layout(title = title, scene = list(xaxis = ax, yaxis = ay, zaxis = az))
  if(nconf > 1) for(i in 2:nconf) plt = plt %>% add_trace(x = Xs[[i]][,1], y = Xs[[i]][,2], z = Xs[[i]][,3], mode = 'lines+markers', name = name[i], line = list(width = 6, color = color[i]), 
                                                          marker = list(size = 3.5, color = color[i]))
  return(plt)
}

#plot X,Y,Z separately
plot_flat = function(X){
  Xflat = melt(X)
  colnames(Xflat) = c("index", "coordinate", "value")
  plt = Xflat %>%
    mutate(index = as.numeric(index)) %>%
    ggplot(aes(index, value, color = coordinate))+
    geom_line()+
    facet_grid(rows = vars(coordinate))+
    theme(legend.position = "none")
  return(plt)
}

#emulate ggplot colors
gg_color_hue = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[seq(n, 1, -1)]
}

cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


