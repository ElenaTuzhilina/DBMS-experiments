source("PoisMS/PCMS.R")
source("PoisMS/WPCMS.R")
source("PoisMS/PoisMS.R")
source("functions.R")
library(matlib)
library(htmlwidgets)
library(tidyr)
library(dplyr)
library(plotly)
library(reshape2)
library(GGally)
cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

############# Load data ####################

chr = 21
start = 15200000
end = 48100000
ext = ".RAWobserved"

load_C = function(res){
  filename = paste0("Data/IMR90/", res, "kb_resolution_intrachromosomal/chr", chr, "/MAPQGE30/chr", chr, "_", res, "kb", ext)
  data = read.csv(filename, header = FALSE, sep = "\t")
  colnames(data) = c("i", "j", "count")
  data = data %>% filter(i >= start & j >= start) %>%
    filter(i <= end & j <= end)
  data = data %>% mutate(i = i/res/1000 + 1, j = j/res/1000 + 1)
  n = max(c(data$i, data$j))
  C = matrix(0, n, n)
  C[cbind(data$i, data$j)] = data$count
  C[cbind(data$j, data$i)] = data$count
  #filter
  index = which(colSums(C) != 0)
  C = C[index, index]
  n = ncol(C)
  saveRDS(list(n = n, C = C, index = index, index_kb = (index-1)*res), paste0("Results/info/conformations_chr", chr,"_res", res, ".rds"))
  #image.plot(log(C+1))
  return(list(C = C, index = index))
}


############# Fit SPoisMS and calculate df for different resolutions ####################

loglambdas = seq(4, -1, -1)
ress = c(100, 50, 25, 10)
ks = c(150, 300, 600, 1500)


degrees = function(X, s, lambda){
  rank = ncol(X)
  dfs = rep(0, rank)
  for(i in 1:rank) dfs[i] = sum(1/(1 + lambda * s / sum(X[,i]^2)))
  mean(dfs)
}

for(j in 1:length(ress)){
  res = ress[j]
  cat("\nres:", res)
  load = load_C(res)
  C = load$C
  n = nrow(C)
  index = load$index
  index_kb = (index-1) * res
  
  dr = DR(index)
  K = dr$K
  #plot(dr$s)
  k = ks[j]
  H = dr$H[,1:k]
  
  if(j == 1){
    Theta0 = matrix(0, k, 3)
    diag(Theta0) = 1
    X0 = H %*% Theta0
    rownames(X0) = paste0("loci", index_kb) 
  } else {
    X0 = X0[paste0("loci", index_kb - index_kb %% ress[j-1]), ]
    rownames(X0) = paste0("loci", index_kb)
  }
  
  dfs = c()
  angles = c()
  iters = c()
  conformations = list()
  betas = list()
  spoisms = list(X = X0, beta = log(mean(C)))
  
  for(i in 1:length(loglambdas)){
    lambda = 10^loglambdas[i]
    cat("\nlambda:", lambda)
    #smoothing spline method
    spoisms = PoisMS(C - lambda * n^2/2 * K, H, 
                     beta0 = spoisms$beta, X0 = spoisms$X)
    df = degrees(spoisms$X, dr$s, lambda * n^2)
    cat(" df:", df)
    
    #save info
    dfs = rbind(dfs, data.frame(df, lambda))
    write.csv(dfs, paste0("Results/smoothing/dfs_chr", chr,"_res", res, "_ws.csv"), row.names = F)
    
    #save conformations
    conformations[[i]] = spoisms$X
    betas[[i]] = spoisms$beta
    saveRDS(list(conformations = conformations, betas = betas), paste0("Results/conformations/smoothed_conformations_chr", chr, "_res", res, "_ws.rds"))
  }
  
  X0 = conformations[[1]]
  rownames(X0) = paste0("loci", index_kb)
}

############# Align dfs for different resolutions ####################

dfs = c()
for(res in ress) dfs = rbind(dfs, data.frame(read.csv(paste0("Results/smoothing/dfs_chr", chr,"_res", res,"_ws.csv"), sep = ","), res))

dfs %>% mutate(resolution = paste0(res, "kb")) %>%
  ggplot(aes(log(lambda, 10), df, color = resolution))+
  geom_point()+
  geom_line()+
  xlab(expression(log[10](lambda)))+
  ylab("degrees−of−freedom")+
  labs(color = "resolution")
ggsave(paste0("Plots/resolution/dfs_chr", chr, ".pdf"), width = 4, height = 4)


############# Align conformations and plot them ####################

Lambda_rescale = function(res, num, res0){
  conf = readRDS(paste0("Results/conformations/smoothed_conformations_chr", chr, "_res", res, "_ws.rds"))
  X = conf$conformations[[num]]
  beta = conf$betas[[num]]
  index_kb = readRDS(paste0("Results/info/conformations_chr", chr, "_res", res, ".rds"))$index_kb
  L = Lambda(X, beta)
  Lflat = data.frame(expand_grid(i = index_kb, j = index_kb), count = c(L)) %>% 
    mutate(i = i - i %% res0, j = j - j %% res0) %>%
    group_by(i,j) %>% summarize(count = sum(count)) %>% ungroup()
  data = Lflat %>% mutate(i = (i - min(Lflat$i))/res0 + 1, j = (j - min(Lflat$j))/res0 + 1)
  n = max(c(data$i, data$j))
  L = matrix(0, n, n)
  L[cbind(data$i, data$j)] = data$count
  #filter
  index = which(colSums(L) != 0)
  L = L[index, index]
  return(list(mat = L, flat = Lflat))
}

X_rescale = function(res, num, res0){
  conf = readRDS(paste0("Results/conformations/smoothed_conformations_chr", chr, "_res", res, "_ws.rds"))
  X = conf$conformations[[num]]
  colnames(X) = c("X", "Y", "Z")
  index_kb = readRDS(paste0("Results/info/conformations_chr", chr, "_res", res, ".rds"))$index_kb
  data.frame(i = index_kb, X) %>% mutate(i = i - i %% res0) %>%
    group_by(i) %>% summarize(X = mean(X), Y = mean(Y), Z = mean(Z))
}

plot_3D = function(Xs, color, name){
  nconf = length(Xs)
  ax = list(title = " ", nticks = 8, showticklabels = FALSE)
  ay = list(title = " ", nticks = 8, showticklabels = FALSE)
  az = list(title = " ", nticks = 8, showticklabels = FALSE)
  plt = plot_ly(x = Xs[[1]][,1], y = Xs[[1]][,2], z = Xs[[1]][,3], type = 'scatter3d', mode = 'lines+markers', name = name[1],
                line = list(width = 6, color = color[1]), marker = list(size = 3.5, color = color[1])) %>% 
    layout(title = title, scene = list(xaxis = ax, yaxis = ay, zaxis = az))
  for(i in 2:nconf){
    plt = plt %>% add_trace(x = Xs[[i]][,1], y = Xs[[i]][,2], z = Xs[[i]][,3], mode = 'lines+markers', name = name[i], line = list(width = 6, color = color[i]), 
                            marker = list(size = 3.5, color = color[i]))
  }
  return(plt)
}

adfs = c(15, 25, 50)
for(i in 1:3){
  Ls = list()
  Xs = list()
  Xflats = c()
  name = c()
  for(j in 1:length(ress)){
    Ls[[j]] = Lambda_rescale(ress[j], 3+i-j+1, 100)
    Xs[[j]] = X_rescale(ress[j], 3+i-j+1, 100) %>% dplyr::select(-i) %>% as.matrix()
    if(j == 1){
      Lflats = Ls[[j]]$flat
    } 
    else {
      Lflats = left_join(Lflats, Ls[[j]]$flat, by = c("i" = "i", "j" = "j"))
      Xs[[j]] = align(Xs[[j-1]], Xs[[j]])$Y
    } 
    colnames(Xs[[j]]) = c("X", "Y", "Z")
    rownames(Xs[[j]]) = 1:nrow(Xs[[j]])
    Xflats = rbind(Xflats, data.frame(melt(Xs[[j]]), res = ress[j]))
  }
  Lflats = Lflats %>% filter(i <= j) %>% dplyr::select(-i, -j)
  colnames(Lflats) = paste0("resolution ", ress, "kb")
  colnames(Xflats) = c("index", "coordinate", "value", "resolution")
  plt = ggpairs(Lflats, columns = 1:4, aes(alpha = 0.1),
                upper = list(continuous = wrap("cor", size = 2.5)),
                lower = list(continuous =  wrap("smooth", se = F, color = "black", color = cbPalette[6])))+
    xlab("expected counts")+
    ylab("expected counts")
  ggsave(paste0("Plots/resolution/ecount_correlation_chr", chr, "_res", res, "_df", adfs[i], ".png"), plt, width = 8, height = 8)
  plot_3D(Xs, gg_color_hue(4), paste0(ress, "kb")) %>%
    saveWidget(paste0("Plots/resolution/conformations_aligned_chr", chr, "_res", res, "_df", adfs[i], ".html"), selfcontained = F, libdir = "lib")
  plt = mutate(Xflats, index = as.numeric(index), resolution = paste0(resolution, "kb")) %>%
    ggplot(aes(index, value, color = resolution))+
    geom_line()+
    facet_grid(rows = vars(coordinate))
  ggsave(paste0("Plots/resolution/projections_aligned_chr", chr, "_res", res, "_df", adfs[i], ".pdf"), plt, width = 5, height = 3)
}


