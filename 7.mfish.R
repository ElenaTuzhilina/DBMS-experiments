source("DBMS/WPCMS.R")
source("DBMS/PCMS.R")
source("DBMS/DBMS.R")
source("DBMS/distributions.R")
source("functions.R")
library(dplyr)
library(ggplot2)
library(readr)
library(tidyverse)
library(zoo)
library(reshape2)

# data = read_tsv("Data/chromosome21.tsv") %>% data.frame()
# data = data %>%
#   dplyr::select(c(2,3,1,4,5)) %>%
#   rename(X = 1, Y = 2, Z = 3, coord = 4, copy = 5) %>%
#   mutate(coord = as.numeric(substr(coord, 7, 14)) - 1)
# write.csv(data, "Data/mfish.csv", row.names = F)

############# Load MFISH data ####################

res = 50
start = 15200000
end = 48100000

mfish = read.csv("Data/mfish.csv", header = T)

get_mfish = function(num){
  subfish = mfish %>% filter(copy == num, coord >= start, coord <= end) %>% na.omit()
  index = (subfish$coord - start)/1000/res + 1
  X = matrix(NA, (end - start)/1000/res + 1, 3)
  X[index, ] = subfish %>% dplyr::select(X, Y, Z) %>% as.matrix()
  colnames(X) = c("X", "Y", "Z")
  return(X)
}

############# Find pairwise Procrustes distances between sample of MFISH conformations ####################

set.seed(1)
nrep = 200
sub = sort(sample(1:max(mfish$copy), 200))

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
  #sX = median(sqrt(rowSums(Xsmooth^2)))
  #sY = median(sqrt(rowSums(Ysmooth^2)))
  sX = 1
  sY = 1
  Xsmooth = scale(Xsmooth, center = mX, scale = FALSE)/sX
  Ysmooth = scale(Ysmooth, center = mY, scale = FALSE)/sY
  X = scale(X, center = mX, scale = FALSE)/sX
  Y = scale(Y, center = mY, scale = FALSE)/sY
  flip = TRUE
  while(flip){
    proc = procrustes(Xsmooth, Ysmooth, scale = FALSE)
    R = proc$rotation
    Ysmooth = Ysmooth %*% R
    Y = Y %*% R
    refl = sign(diag(cor(Xsmooth, Ysmooth)))
    Ysmooth = scale(Ysmooth, center = FALSE, scale = refl)
    Y = scale(Y, center = FALSE, scale = refl)
    flip = (sum(refl < 0) > 0)
  }
  return(list(X = X, Y = Y, ss = sum((X - Y)^2, na.rm = T)/length(index), ssmooth = sum((Xsmooth - Ysmooth)^2)/length(index)))
}

SS = matrix(NA, nrep, nrep)
SSmooth = matrix(NA, nrep, nrep)
for(i in 1:(nrep-1)){
  for(j in (i+1):nrep){
    cat("i:", i, "j:", j, "\n")
    X = get_mfish(sub[i])
    Y = get_mfish(sub[j])
    sa = smooth_align(X, Y, 1)
    SS[i,j] = sa$ss
    SSmooth[i,j] = sa$ssmooth
    if((i+j) %% nrep == 0){
      saveRDS(SS, "Results/mfish/mfish_spw_dist_width1.rds")
      #saveRDS(SSmooth, "Results/mfish/mfish_mav_spw_dist.rds")
    } 
  }
}


############# Find Procrustes distances between from each confromation and MFISH conformations ####################

get_conformation = function(method, num){
  res = 50
  chr = 21
  ext = ".RAWobserved"
  if(method == "SPoisMS"){
    conf = readRDS(paste0("Results/conformations/smoothing_conformations_chr", chr, "_res", res, ".rds"))[[method]][[num]]
  } 
  if(method %in% c("PoisMS", "HPoisMS", "ZIPoisMS", "NBMS")) {
    conf = readRDS(paste0("Results/conformations/dbms_conformations_chr", chr, "_res", res, ".csv"))[[method]][[num]]
  }
  if(method == "PM1") conf = read.csv("Data/PM1.csv", header = FALSE, sep = ",") %>% as.matrix()
  index = read.csv("Data/index_hic.csv", header = T)$index
  X = matrix(NA, (end - start)/1000/res + 1, 3)
  X[index, ] = as.matrix(conf)
  colnames(X) = c("X", "Y", "Z")
  return(X)
}

SSconf = c()
for(method in c("SPoisMS", "PoisMS", "HPoisMS","ZIPoisMS", "NBMS")){
  cat(method)
  if(method == "SPoisMS") dfs = c(5, 9, 15, 27, 47, 84) 
  else dfs = seq(10, 100, 10)
  for(num in 1:length(dfs)){
      cat("\ndf:", dfs[num])
      Y = get_conformation(method, num)
      cat("\ni: ")
      for(i in 1:nrep){
      cat(i, " ")
      X = get_mfish(sub[i])
      sa = smooth_align(X, Y, width = 1, smoothY = F)
      SSconf = rbind(SSconf, data.frame(ss = sa$ss, method, df = dfs[num]))
    }
  }
  saveRDS(SSconf, "Results/mfish/mfish_sconf_dist_width1_all.rds")
}

############# Compare the distribution of Procrustes distances ####################

SS = readRDS("Results/mfish/mfish_spw_dist_width1.rds")
SSconf = readRDS("Results/mfish/mfish_sconf_dist_width1.rds")


#rbind(SSconf, data.frame(ss = SS[upper.tri(SS)], method = "MFISH")) %>%
SSconf %>% filter(method == "ZIPoisMS") %>%
  ggplot(aes(x = log(ss,10), y=..density.., fill = factor(df), color = factor(df)))+
  geom_histogram(position = "identity", alpha = 0.3, bins = 30, color = NA)+
  geom_density(alpha = 0.1)+
  facet_wrap(~method)

SSconf %>% group_by(method) %>% summarize(mean(ss))


############# Compare a pair MFISH conformations ####################

# X = get_mfish(500)
# Y = get_mfish(1000)
X = get_conformation("NBMS", 1)
Y = get_conformation("NBMS", 10)

plot_3D(list(X, Y), c("red", "green"), c("X", "Y"))
# plot_flat(X)
# plot_flat(Y)
# cor(X, Y, use  = "complete.obs")
# 
width = 50
Xsmooth = smooth_X(X, width)
Ysmooth = smooth_X(Y, width)
plot(curvature(Xsmooth))
plot(curvature(Ysmooth))
# cor(X, Y, use  = "complete.obs")
# cor(curvature(Xsmooth), curvature(Ysmooth), use  = "complete.obs")
# 
plot_3D(list(Xsmooth, Ysmooth), c("red", "green"), c("X", "Y"))
plot_flat(Xsmooth)
plot_flat(Ysmooth)
# 
sa = smooth_align(X, Y, 1)
plot_3D(list(smooth_X(sa$X, 100), smooth_X(sa$Y, 100)), c("red", "green"), c("X", "Y"))
plot_3D(list(sa$X, sa$Y), c("red", "green"), c("X", "Y"))
# 
# 
# 

X = read.csv("Data/PM1.csv", header = FALSE, sep = ",")
colnames(X) = c("X", "Y", "Z")
rownames(X) = index
plot_3D(list(smooth_X(X, 10)), c("black"), c("PM1"))
plot_flat(X)

melt(X)

SSpastis = c()
for(method in c("PM1")){
  cat(method)
    Y = get_conformation(method, 0)
    cat("\ni: ")
    for(i in 1:nrep){
      cat(i, " ")
      X = get_mfish(sub[i])
      sa = smooth_align(X, Y, width = 1, smoothY = F)
      SSpastis = rbind(SSpastis, data.frame(ss = sa$ss, method, df = NA))
    }
  }
saveRDS(SSpastis, "Results/mfish/mfish_spastis_dist_width1_all.rds")



