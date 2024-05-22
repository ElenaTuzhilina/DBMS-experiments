source("DBMS/WPCMS.R")
source("DBMS/PCMS.R")
source("DBMS/DBMS.R")
source("DBMS/distributions.R")
source("functions.R")
library(segmented)
library(htmlwidgets)
library(tictoc)
library(hicrep)
library(HiCdiff)
library(stats)
library(tidyr)
palette = cbPalette[c(4, 2, 3, 7, 5, 6)]

############# Load data ####################

res = 50
chr = 21
start = 15200000
end = 48100000
ext = ".RAWobserved"
filename = paste0("Data/IMR90/", res,"kb_resolution_intrachromosomal/chr", chr, "/MAPQGE30/chr", chr, "_", res, "kb", ext)
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
Ctilde = KRnorm(C)


############# Find train scores for DBMS models  ####################

dfs = seq(10, 100, 10)
methods = c("PoisMS","SPoisMS", "HPoisMS", "ZIPoisMS", "NBMS", "NMS", "EMS")

update_info = function(C, H, method, param, Theta0, info, method_name){
  #get reconstruction
  tic()
  sol = DBMS(C, H, param = param, method = method, Theta = Theta0, eps_wpcms = 1e-5, eps_dbms = 1e-5, verbose_dbms = TRUE)
  save = toc()
  
  if(method_name %in% c("PoisMS", "EMS")) par = NA
  if(method_name %in% c("HPoisMS", "ZIPoisMS")) par =  sol$param$p
  if(method_name == "NBMS") par = sol$param$r
  if(method_name == "NMS") par = sol$param$sigma2
  
  #update info
  return(list(info = rbind(info, data.frame("df" = ncol(H), "beta" = sol$param$beta, "param" = par, "epoch" = sol$epoch, "iter_total" = sol$iter_total, 
                         dbms_evaluate(sol$X, sol$param, C, method, method_name), time = save$toc - save$tic,
                         "method" = method_name)), X = sol$X))
}

info = data.frame()
Xs = list("PoisMS" = list(), "HPoisMS" = list(), "ZIPoisMS" = list(), "NBMS" = list(), "NMS" = list(), "EMS" = list())

for(i in 1:length(dfs)){
    df = dfs[i]
    H = load_H(df, index, orth = TRUE)
    Theta0 = matrix(0, df, 3)
    diag(Theta0) = 1
    
    cat("=======================================\ndf :", df, "method : PoisMS\n=======================================\n")
    param = list(beta = log(mean(C[C > 0])))
    upd = update_info(C, H, pois, param, Theta0, info, "PoisMS")
    info = upd$info
    Xs[["PoisMS"]][[i]] = upd$X

    cat("=======================================\ndf :", df, "method : HPoisMS\n=======================================\n")
    param$p = mean(C == 0)
    upd = update_info(C, H, hpois, param, Theta0, info, "HPoisMS")
    info = upd$info
    Xs[["HPoisMS"]][[i]] = upd$X

    cat("=======================================\ndf :", df, "method : ZIPoisMS\n=======================================\n")
    upd = update_info(C, H, zipois, param, Theta0, info, "ZIPoisMS")
    info = upd$info
    Xs[["ZIPoisMS"]][[i]] = upd$X

    cat("=======================================\ndf :", df, "method : NBMS\n=======================================\n")
    param$p = NULL
    param$r = 1
    upd = update_info(C, H, nbin, param, Theta0, info, "NBMS")
    info = upd$info
    Xs[["NBMS"]][[i]] = upd$X

    cat("=======================================\ndf :", df, "method : NMS\n=======================================\n")
    param = list(beta = mean(Ctilde[Ctilde>0]))
    param$sigma2 = var(c(Ctilde[Ctilde>0]))
    upd = update_info(Ctilde, H, norm, param, Theta0, info, "NMS")
    info = upd$info
    Xs[["NMS"]][[i]] = upd$X
    
    cat("=======================================\ndf :", df, "method : EMS\n=======================================\n")
    param = list(beta = log(mean(Ctilde[Ctilde > 0])))
    upd = update_info(Ctilde, H, expon, param, Theta0, info, "EMS")
    info = upd$info
    Xs[["EMS"]][[i]] = upd$X
    
    write.csv(info, paste0("Results/dbms/bulk_cell_chr", chr, "_res", res, ".csv"), row.names = FALSE)
    saveRDS(Xs, paste0("Results/conformations/dbms_conformations_chr", chr, "_res", res, ".rds"))
}


############# Elbow detection ####################

info = read.csv(paste0("Results/dbms/bulk_cell_chr", chr, "_res", res, ".csv")) %>%
  mutate(method = factor(method, levels = methods))

detect_elbow = function(info, method_name){
  data = info %>% filter(method == method_name)
  LR = lm(loglik ~ df, data = info %>% filter(method == method_name))
  SR = segmented(LR, seg.Z = ~ df, npsi = 1, tol = 1e-6, it.max = 1000)
  elbow = SR$psi[2]
  return(data.frame(df = elbow, method = method_name))
}

elbow = rbind(detect_elbow(info, "PoisMS"), 
              detect_elbow(info, "HPoisMS"), 
              detect_elbow(info, "ZIPoisMS"), 
              detect_elbow(info, "NBMS"), 
              detect_elbow(info, "NMS"),
              detect_elbow(info, "EMS")) 
elbow = elbow %>% mutate(method = factor(method, levels = methods))
write.csv(elbow, paste0("Results/dbms/bulk_cell_elbow_chr", chr, "_res", res, ".csv"), row.names = FALSE)
  
############# Plot train scores (Figure 9 and 20, left) ####################

info %>% filter(!method %in% c("NMS", "EMS")) %>%
ggplot(aes(df, loglik, color = method, linetype = method))+
  geom_line()+
  geom_point()+
  geom_vline(elbow %>% filter(!method %in% c("NMS", "EMS")), 
             mapping = aes(xintercept = df, color = method, linetype = method), size = 0.3)+
  scale_color_manual(values = palette)+
  scale_fill_manual(values = palette)+
  ylab("train loss")+
  xlab("degrees-of-freedom")+
  scale_x_continuous(breaks = dfs)+
  theme_bw()+
  theme(legend.position="none")
ggsave(paste0("Plots/dbms/bulk_cell_chr", chr, "_res", res, ".pdf"), height = 3, width = 3)

info %>% filter(method %in% c("NMS", "EMS")) %>%
  ggplot(aes(df, loglik, color = method, linetype = method))+
  geom_line()+
  geom_point()+
  geom_vline(elbow %>% filter(method %in% c("NMS", "EMS")), mapping = aes(xintercept = df, color = method, linetype = method), size = 0.3)+
  scale_color_manual(values = palette)+
  scale_fill_manual(values = palette)+
  ylab("train loss")+
  xlab("degrees-of-freedom")+
  theme_bw()+
  theme(legend.position="none")+
  scale_x_continuous(breaks = dfs)
ggsave(paste0("Plots/dbms/bulk_cell_chr", chr, "_res", res, "_norm.pdf"), height = 3, width = 3)

############# Plot correlations between E and O counts (Figure 17 and Figure 20, right) ####################

sinfo = info  %>% dplyr::select(-beta, -param, - epoch, -iter_total, -time)

conformations = readRDS(paste0("Results/conformations/smoothing_conformations_chr", chr,"_res", res, ".rds"))
betas = readRDS(paste0("Results/conformations/smoothing_betas_chr", chr,"_res", res, ".rds"))
dfs = c(5, 9, 15, 27, 47, 84)
for(i in 1:length(dfs)){
  sinfo = rbind(sinfo, data.frame(df = dfs[i],  dbms_evaluate(conformations$SPoisMS[[i]], list(beta = betas$SPoisMS[[i]]), C, pois, "PoisMS"), method = "SPoisMS"))
}

sinfo = sinfo %>%
  dplyr::select(df, method,  corE, spcorE, strcorE) %>%
  gather(type, cor, 3:5) %>% 
  mutate(type = ifelse(type == "corE", "Pearson", ifelse(type == "spcorE", "Spearman", "stratified"))) %>%
  mutate(method = factor(method, levels = c("PoisMS","SPoisMS", "HPoisMS", "ZIPoisMS", "NBMS", "NMS", "EMS"))) 

sinfo %>%
  filter(type != "Spearman", !method %in% c("NMS", "EMS")) %>%
  ggplot(aes(df, cor, color = method, linetype = method))+
  geom_line()+
  geom_point()+
  scale_color_manual(values = cbPalette[c(4, 1, 2, 3, 7, 5, 6)])+
  ylab("train correlation")+
  xlab("degrees-of-freedom")+
  scale_x_continuous(breaks = seq(10, 100, 10))+
  facet_wrap(~type, scales = "free")+
  theme_bw()
ggsave(paste0("Plots/dbms/bulk_cell_cors_chr", chr, "_res", res, ".pdf"), height = 3, width = 6)

info  %>% dplyr::select(-beta, -param, - epoch, -iter_total, -time) %>% 
  filter(method %in% c("NMS", "EMS")) %>%
  dplyr::select(df, method, corE, spcorE, strcorE) %>%
  gather(type, cor, 3:5) %>% 
  mutate(type = ifelse(type == "corE", "Pearson", ifelse(type == "spcorE", "Spearman", "stratified"))) %>%
  mutate(method = factor(method, levels = c("NMS", "EMS"))) %>%
  filter(type == "Pearson") %>%
  ggplot(aes(df, cor, color = method, linetype = method))+
  geom_line()+
  geom_point()+
  scale_color_manual(values = palette)+
  ylab("train correlation")+
  xlab("degrees-of-freedom")+
  scale_x_continuous(breaks = seq(10, 100, 10))+
  facet_wrap(~type, scales = "free")+
  theme_bw()
ggsave(paste0("Plots/dbms/bulk_cell_cors_chr", chr, "_res", res, "_norm.pdf"), height = 3, width = 4)

############# Plot convergence performance (Figure 13 and 24) ####################

info %>% filter(!method %in% c("NMS", "EMS")) %>%
  ggplot(aes(df, time, color = method, linetype = method))+
  geom_line()+
  geom_point()+
  xlab("degrees-of-freedom")+
  ylab("time (secs)")+
  scale_color_manual(values = palette) +
  theme_bw()
ggsave(paste0("Plots/dbms/time_chr", chr, "_res", res, ".pdf"), height = 3, width = 4.5)

info %>% filter(!method %in% c("NMS", "EMS")) %>%
  ggplot(aes(df, iter_total, color = method, linetype = method))+
  geom_line()+
  geom_point()+
  xlab("degrees-of-freedom")+
  ylab("total iterations")+
  scale_color_manual(values = palette)+
  theme_bw()
ggsave(paste0("Plots/dbms/iters_chr", chr, "_res", res, ".pdf"), height = 3, width = 4.5)

info %>% filter(!method %in% c("NMS", "EMS")) %>%
  ggplot(aes(df, epoch, color = method, linetype = method))+
  geom_line()+
  geom_point()+
  xlab("degrees-of-freedom")+
  ylab("epochs")+
  scale_color_manual(values = palette)+
  theme_bw()+
  theme(legend.position = "none")
ggsave(paste0("Plots/dbms/epochs_chr", chr, "_res", res, ".pdf"), height = 3, width = 3.5)

############# Find and plot optimal conformations with optimal df (Figure 8 and 21) ####################

save_conformation = function(C, method, param, df, method_name, align_with = NULL, zlim){
  method_name = as.character(method_name)
  H = load_H(df, index, orth = TRUE)
  Theta0 = matrix(0, df, 3)
  diag(Theta0) = 1
  sol = DBMS(C, H, param = param, method = method, Theta = Theta0, eps_wpcms = 1e-5, verbose_wpcms = FALSE, eps_dbms = 1e-5, verbose_dbms = TRUE)
  if(!is.null(align_with)) sol$X = align(align_with, sol$X)$Y
  
  L = Lambda(sol$X, sol$param$beta)
  E = expected(sol$X, sol$param, method_name)

  visualize(sol$X, "3D") %>%
    saveWidget(paste0("Plots/conformations/", method_name,"_chr", chr, "_res", res, ".html"), selfcontained = F, libdir = "lib")
  png(paste0("Plots/conformations/", method_name,"_heatmap_chr", chr, "_res", res, ".png"))
  image.plot(L, axes = F)
  axis(1, at = seq(0, n, by = 100)/n, labels = seq(0, n, by = 100))
  axis(2, at = seq(0, n, by = 100)/n, labels = seq(0, n, by = 100))
  dev.off()

  png(paste0("Plots/conformations/", method_name,"_logEheatmap_chr", chr, "_res", res, ".png"))
  image.plot(log(pmax(E, 0)), axes = F, zlim = zlim)
  axis(1, at = seq(0, n, by = 100)/n, labels = seq(0, n, by = 100))
  axis(2, at = seq(0, n, by = 100)/n, labels = seq(0, n, by = 100))
  dev.off()

  trindex = upper.tri(C, diag = T)
  return(list(X = sol$X, param = sol$param))
}

elbow = read.csv(paste0("Results/dbms/bulk_cell_elbow_chr", chr, "_res", res, ".csv")) %>% mutate(df = round(df))
conf = list()

param = list(beta = log(mean(C[C > 0])), p = mean(C == 0), r = 1)
conf[["PoisMS"]] = save_conformation(C, pois, param, elbow$df[1], "PoisMS", NULL, c(0,9))
conf[["HPoisMS"]] = save_conformation(C, hpois, param, elbow$df[2], "HPoisMS", conf[["PoisMS"]]$X, c(0,9))
conf[["ZIPoisMS"]] = save_conformation(C, zipois, param, elbow$df[3], "ZIPoisMS", conf[["HPoisMS"]]$X, c(0,9))
conf[["NBMS"]] = save_conformation(C, nbin, param, elbow$df[4], "NBMS", conf[["ZIPoisMS"]]$X, c(0,9))

param = list(beta = mean(Ctilde[Ctilde > 0]), sigma2 = var(Ctilde[Ctilde > 0]))
conf[["NMS"]] = save_conformation(Ctilde, norm, param, elbow$df[5], "NMS", NULL, c(-11, 0))
param = list(beta = log(mean(Ctilde[Ctilde > 0])))
conf[["EMS"]] = save_conformation(Ctilde, expon, param, elbow$df[6], "EMS", conf[["NMS"]]$X, c(-11, 0))

saveRDS(conf, paste0("Results/conformations/dbms_optimal_conformations_chr", chr, "_res", res, ".rds"))


