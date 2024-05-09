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
library(matlib)
library(htmlwidgets)
library(tidyr)
library(plotly)
library(GGally)

data = read_tsv("Data/chromosome21.tsv") %>% data.frame()
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

SS = matrix(NA, nrep, nrep)
for(i in 1:nrep){
  for(j in 1:nrep){
    cat("i:", i, "j:", j, "\n")
    X = get_mfish(sub[i])
    Y = get_mfish(sub[j])
    sa = smooth_align(X, Y, 0, FALSE, FALSE)
    SS[i,j] = sa$ss
    if((i+j) %% nrep == 0){
      saveRDS(SS, "Results/mfish/mfish_pw_dist_scaled.rds")
    } 
  }
}
saveRDS(SS, "Results/mfish/mfish_pw_dist_scaled.rds")

############# Find Procrustes distances between from each conformation and MFISH conformations ####################

get_conformation = function(method, num = 0){
  res = 50
  chr = 21
  ext = ".RAWobserved"
  if(method == "SPoisMS"){
    conf = readRDS(paste0("Results/conformations/smoothing_conformations_chr", chr, "_res", res, ".rds"))[[method]][[num]]
  } 
  if(method %in% c("PoisMS", "HPoisMS", "ZIPoisMS", "NBMS", "NMS", "EMS")) {
    conf = readRDS(paste0("Results/conformations/dbms_optimal_conformations_chr", chr, "_res", res, ".rds"))[[method]]$X
  }
  if(method %in% c("PM1", "PM2", "poisson")) conf = read.csv(paste0("Data/PASTIS/", method, ".csv"), header = FALSE, sep = ",") %>% as.matrix()
  index = read.csv("Data/index_hic.csv", header = T)$index
  X = matrix(NA, (end - start)/1000/res + 1, 3)
  X[index, ] = conf
  colnames(X) = c("X", "Y", "Z")
  return(X)
}

methods = c("PM1", "PM2", "poisson", "SPoisMS", "PoisMS", "HPoisMS", "ZIPoisMS", "NBMS", "NMS", "EMS")
SSconf = c()
for(method in methods){
  cat(method)
  if(method == "SPoisMS") num = 3
  #if(method %in% c("PoisMS", "SPoisMS")) num = 3
  #if(method %in% c("HPoisMS","ZIPoisMS", "NBMS")) num = 2
  #if(method %in% c("PM1","PM2", "poisson")) num = 0
  Y = get_conformation(method, num)
  cat("\ni: ")
  for(i in 1:nrep){
    cat(i, " ")
    X = get_mfish(sub[i])
    sa = smooth_align(X, Y, width = 0, smoothX = F, smoothY = F)
    SSconf = rbind(SSconf, data.frame(ss = sa$ss, method))
  }
  saveRDS(SSconf, "Results/mfish/mfish_conf_dist_scaled_norm.rds")
}

 ############# Compare the distribution of Procrustes distances to PASTIS and our methods (Figure 10 and 22) ####################

SS = readRDS("Results/mfish/mfish_pw_dist_scaled.rds")
SSconf = readRDS("Results/mfish/mfish_conf_dist_scaled_norm.rds")

methods = c("PM1", "PM2", "poisson", "SPoisMS", "PoisMS", "HPoisMS", "ZIPoisMS", "NBMS")

rbind(data.frame(ss = SS[upper.tri(SS)], method = "MFISH", group = rep(methods, rep(length(SS[upper.tri(SS)]), length(methods)))),
      SSconf %>% filter(method %in% methods) %>%
        mutate(group = method)) %>%
  mutate(group = ifelse(group == "PM1", "PASTIS-PM1", ifelse(group == "PM2", "PASTIS-PM2", ifelse(group == "poisson", "PASTIS-Poisson", group)))) %>%
  mutate(group = factor(group, levels = c("PASTIS-PM1", "PASTIS-PM2", "PASTIS-Poisson", methods[-(1:3)])),
         method = factor(method, levels = c("MFISH", methods))) %>%
  ggplot(aes(x = log(ss, 10), y=..density.., fill = method, color = method))+
  geom_histogram(position = "identity", alpha = 0.3, bins = 25, color = NA)+
  geom_density(alpha = 0.1)+
  facet_wrap(~group)+
  xlab("log(procrustes distance)")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_color_manual(values = c("black", gg_color_hue(length(methods))))+
  scale_fill_manual(values = c("black", gg_color_hue(length(methods))))
ggsave("Plots/mfish/compare_distributions.pdf", height = 5, width = 7)


methods = c("NMS", "EMS")
rbind(data.frame(ss = SS[upper.tri(SS)], method = "MFISH", group = rep(methods, rep(length(SS[upper.tri(SS)]), length(methods)))), 
      SSconf %>% filter(method %in% methods) %>%
        mutate(group = method)) %>%
  mutate(group = factor(group, c(methods)),
         method = factor(method, levels = c("MFISH", methods))) %>%
  ggplot(aes(x = log(ss, 10), y=..density.., fill = method, color = method))+
  geom_histogram(position = "identity", alpha = 0.3, bins = 20, color = NA)+
  geom_density(alpha = 0.1)+
  facet_wrap(~group)+
  xlab("log(procrustes distance)")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_color_manual(values = c("black", cbPalette[c(5, 6)]))+
  scale_fill_manual(values = c("black", cbPalette[c(5, 6)]))
ggsave("Plots/mfish/compare_distributions_norm.pdf", height = 3, width = 5)

for(method in methods){
  cat("\n", method, ":")
  cat(wilcox.test(x = SS[upper.tri(SS)], y = SSconf %>% filter(method == method) %>% pull(ss), alternative = "greater")$p.value)
}

############# Plot PASTIS conformations (Figure 11) ####################

visualize(get_conformation("PM1", 0), "3D") %>%
  saveWidget(paste0("Plots/mfish/conformations_PM1_chr", chr, "_res", res, ".html"), selfcontained = F, libdir = "lib")
visualize(get_conformation("PM2", 0), "3D") %>%
  saveWidget(paste0("Plots/mfish/conformations_PM2_chr", chr, "_res", res, ".html"), selfcontained = F, libdir = "lib")
visualize(get_conformation("poisson", 0), "3D") %>%
  saveWidget(paste0("Plots/mfish/conformations_poisson_chr", chr, "_res", res, ".html"), selfcontained = F, libdir = "lib")

############# Compare a pair MFISH conformations (Figure 18) ####################

X = get_mfish(1)
Y = get_mfish(1000)
Z = get_mfish(1500)
visualize(X, "3D")
visualize(Y, "3D")
visualize(Z, "3D")

Xsmooth = smooth_X(sa$X, 10)
plot_flat(X)
plot_flat(Xsmooth)

plot_3D(list(X, Y), c("red", "green"), c("X", "Y"))
sa = smooth_align(X, Y, 10)
plot_3D(list(smooth_X(sa$X, 10), smooth_X(sa$Y, 10)), c("red", "green"), c("X", "Y"))
plot_3D(list(sa$X, sa$Y), c("red", "green"), c("X", "Y"))

