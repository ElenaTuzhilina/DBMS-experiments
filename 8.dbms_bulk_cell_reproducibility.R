source("DBMS/WPCMS.R")
source("DBMS/PCMS.R")
source("DBMS/DBMS.R")
source("DBMS/distributions.R")
source("functions.R")
library(segmented)
library(htmlwidgets)
library(tictoc)
library(hicrep)
library(tidyr)

############# Load data ####################

train = "replicate"
test = "primary"

res = 50
chr = 21
start = 15200000
end = 48100000
ext = ".RAWobserved"
filename = paste0("Data/GM12878/", res,"kb_resolution_intrachromosomal_", train, "/chr", chr, "/MAPQGE30/chr", chr, "_", res, "kb", ext)
data = read.csv(filename, header = FALSE, sep = "\t")
colnames(data) = c("i", "j", "count")
data = data %>% filter(i >= start & j >= start) %>%
  filter(i <= end & j <= end)
data = data %>% mutate(i = (i - start)/res/1000 + 1, j = (j-start)/res/1000 + 1)
n0 = max(c(data$i, data$j))
Ctrain = matrix(0, n0, n0)
Ctrain[cbind(data$i, data$j)] = data$count
Ctrain[cbind(data$j, data$i)] = data$count
#filter
index = which(colSums(Ctrain) != 0)
Ctrain = Ctrain[index, index]
n = ncol(Ctrain)
Ctraintilde = KRnorm(Ctrain)

filename = paste0("Data/GM12878/", res,"kb_resolution_intrachromosomal_", test, "/chr", chr, "/MAPQGE30/chr", chr, "_", res, "kb", ext)
data = read.csv(filename, header = FALSE, sep = "\t")
colnames(data) = c("i", "j", "count")
data = data %>% filter(i >= start & j >= start) %>%
  filter(i <= end & j <= end)
data = data %>% mutate(i = (i - start)/res/1000 + 1, j = (j-start)/res/1000 + 1)
Ctest = matrix(0, n0, n0)
Ctest[cbind(data$i, data$j)] = data$count
Ctest[cbind(data$j, data$i)] = data$count
Ctest = Ctest[index, index]
Ctesttilde = KRnorm(Ctest)

hist(c(Ctrain-Ctest))
hist(c(Ctraintilde-Ctesttilde))

############# Find train scores for the DBMS models ####################

dfs = seq(10, 100, 10)

update_info = function(Ctrain, Ctest, H, method, param, Theta0, info, method_name){
  #get reconstruction
  tic()
  sol = DBMS(Ctrain, H, param = param, method = method, Theta = Theta0, eps_wpcms = 1e-5, eps_dbms = 1e-5, verbose_dbms = TRUE)
  save = toc()
  
  if(method_name %in% c("PoisMS", "EMS")) par = NA
  if(method_name %in% c("HPoisMS", "ZIPoisMS")) par =  sol$param$p
  if(method_name == "NBMS") par = sol$param$r
  if(method_name == "NMS") par = sol$param$sigma2
  
  #update info
  eval = rbind(data.frame(dbms_evaluate(sol$X, sol$param, Ctrain, method, method_name), set = "train"), 
                                                    data.frame(dbms_evaluate(sol$X, sol$param, Ctest, method, method_name), set = "test"))
  if(method_name %in% c("NMS", "EMS")) eval$norm = TRUE
  else{
    eval$norm = F
    eval = rbind(eval,
                 data.frame(dbms_evaluate(sol$X, sol$param, Ctrain, method, method_name, norm = TRUE), set = "train", norm = TRUE), 
                 data.frame(dbms_evaluate(sol$X, sol$param, Ctest, method, method_name, norm = TRUE), set = "test", norm = TRUE))
  } 
  return(list(info = rbind(info, data.frame("df" = ncol(H), "beta" = sol$param$beta, "param" = par, "epoch" = sol$epoch, "iter_total" = sol$iter_total, 
                          eval, time = save$toc - save$tic,
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
    param = list(beta = log(mean(Ctrain[Ctrain > 0])))
    upd = update_info(Ctrain, Ctest, H, pois, param, Theta0, info, "PoisMS")
    info = upd$info
    Xs[["PoisMS"]][[i]] = upd$X
    
    cat("=======================================\ndf :", df, "method : HPoisMS\n=======================================\n")
    param$p = mean(Ctrain == 0)
    upd = update_info(Ctrain, Ctest, H, hpois, param, Theta0, info, "HPoisMS")
    info = upd$info
    Xs[["HPoisMS"]][[i]] = upd$X

    cat("=======================================\ndf :", df, "method : ZIPoisMS\n=======================================\n")
    upd = update_info(Ctrain, Ctest, H, zipois, param, Theta0, info, "ZIPoisMS")
    info = upd$info
    Xs[["ZIPoisMS"]][[i]] = upd$X
    
    cat("=======================================\ndf :", df, "method : NBMS\n=======================================\n")
    param$p = NULL
    param$r = 1
    upd = update_info(Ctrain, Ctest, H, nbin, param, Theta0, info, "NBMS")
    info = upd$info
    Xs[["NBMS"]][[i]] = upd$X
    
    cat("=======================================\ndf :", df, "method : NMS\n=======================================\n")
    param = list(beta = mean(Ctraintilde[Ctraintilde>0]))
    param$sigma2 = var(c(Ctraintilde[Ctraintilde>0]))
    upd = update_info(Ctraintilde, Ctesttilde, H, norm, param, Theta0, info, "NMS")
    info = upd$info
    Xs[["NMS"]][[i]] = upd$X
    
    cat("=======================================\ndf :", df, "method : EMS\n=======================================\n")
    param = list(beta = log(mean(Ctraintilde[Ctraintilde > 0])))
    upd = update_info(Ctraintilde, Ctesttilde, H, expon, param, Theta0, info, "EMS")
    info = upd$info
    Xs[["EMS"]][[i]] = upd$X
    
    write.csv(info, paste0("Results/dbms/bulk_cell_repr_", train, "-", test, "_chr", chr, "_res", res, "_norm.csv"), row.names = FALSE)
    saveRDS(Xs, paste0("Results/conformations/dbms_conformations_repr_", train, "-", test, "_chr", chr, "_res", res, "_norm.csv"))
}


############# Plot correlations between E and O counts (Figure 12) ####################

methods = c("PoisMS", "HPoisMS", "ZIPoisMS", "NBMS")
info = rbind(read.csv(paste0("Results/dbms/bulk_cell_repr_primary-replicate_chr", chr, "_res", res, "_norm.csv")) %>% 
  filter(set == "test", norm = FALSE, method %in% methods) %>%
  mutate(set = "primary"),
  read.csv(paste0("Results/dbms/bulk_cell_repr_replicate-primary_chr", chr, "_res", res, "_norm.csv")) %>% 
    filter(set == "test", norm = FALSE, method %in% methods) %>%
    mutate(set = "replicate"))
  
info =  info %>% filter(method %in% methods) %>%
  mutate(method = factor(method, levels = methods)) %>%
  mutate(set = ifelse(set == "train", "primary", "replicate"))

info %>% dplyr::select(df, method,  corE, spcorE, strcorE, set) %>%
  gather(type, cor, 3:5) %>% 
  mutate(type = ifelse(type == "corE", "Pearson", ifelse(type == "spcorE", "Spearman", "stratified"))) %>%
  mutate(method = factor(method, levels = methodse)) %>%
  filter(type != "Spearman") %>%
  ggplot(aes(df, cor, color = method, linetype = method))+
  geom_line()+
  geom_point()+
  scale_color_manual(values = cbPalette[c(4, 2, 3, 7)])+
  ylab("correlation")+
  xlab("degrees-of-freedom")+
  scale_x_continuous(breaks = seq(10, 100, 10))+
  facet_wrap(set~type, scales = "free_y")+
  theme_bw()
ggsave(paste0("Plots/dbms/bulk_cell_cors_rep_chr", chr, "_res", res, ".pdf"), height = 3.5, width = 6)


############# Plot correlations between E and O counts for different resolutions (Figure 19) ####################

info = rbind(data.frame(read.csv(paste0("Results/dbms/bulk_cell_repr_chr", chr, "_res100.csv")), resolution = "100kb"),
             data.frame(read.csv(paste0("Results/dbms/bulk_cell_repr_chr", chr, "_res50.csv")), resolution = "50kb"),
             data.frame(read.csv(paste0("Results/dbms/bulk_cell_repr_chr", chr, "_res25.csv")), resolution = "25kb"),
             data.frame(read.csv(paste0("Results/dbms/bulk_cell_repr_chr", chr, "_res10.csv")), resolution = "10kb"))

info =  info %>% mutate(method = factor(method, levels = c("PoisMS", "HPoisMS", "ZIPoisMS", "NBMS"))) %>%
  mutate(set = ifelse(set == "train", "primary", "replicate")) %>%
  mutate(resolution = factor(resolution, levels = paste0(c(100,50,25,10), "kb")))

info %>% 
  ggplot(aes(df, strcorE, color = method, linetype = method))+
  geom_line()+
  geom_point()+
  scale_color_manual(values = cbPalette[c(4, 2, 3, 7)])+
  ylab("correlation")+
  xlab("degrees-of-freedom")+
  scale_x_continuous(breaks = seq(10, 100, 10))+
  facet_grid(set~resolution)+ 
  theme_bw()+
  theme(axis.text.x=element_text(size=rel(0.9)))
ggsave(paste0("Plots/dbms/bulk_cell_cors_rep_chr", chr, ".pdf"), height = 5, width = 8)


############# Plot correlations between normalized E and O counts (Figure 23) ####################

methods = c("PoisMS", "HPoisMS", "ZIPoisMS", "NBMS", "NMS", "EMS")
info = rbind(read.csv(paste0("Results/dbms/bulk_cell_repr_primary-replicate_chr", chr, "_res", res, "_norm.csv")) %>% 
               filter(set == "test", norm = TRUE, method %in% methods) %>%
               mutate(set = "primary"),
             read.csv(paste0("Results/dbms/bulk_cell_repr_replicate-primary_chr", chr, "_res", res, "_norm.csv")) %>% 
               filter(set == "test", norm = TRUE, method %in% methods) %>%
               mutate(set = "replicate"))


info =  info %>% filter(method %in% methods, norm == TRUE) %>%
  mutate(method = factor(method, levels = methods))

info %>% dplyr::select(df, method, corE, spcorE, strcorE, set) %>%
  gather(type, cor, 3:5) %>% 
  mutate(type = ifelse(type == "corE", "Pearson", ifelse(type == "spcorE", "Spearman", "stratified"))) %>%
  mutate(method = factor(method, levels = methods)) %>%
  filter(type == "Pearson") %>%
  ggplot(aes(df, cor, color = method, linetype = method))+
  geom_line()+
  geom_point()+
  scale_color_manual(values = cbPalette[c(4, 2, 3, 7, 5, 6)])+
  ylab("correlation")+
  xlab("degrees-of-freedom")+
  scale_x_continuous(breaks = seq(10, 100, 10))+
  facet_grid(set~type)+
  theme_bw()
ggsave(paste0("Plots/dbms/bulk_cell_cors_rep_chr", chr, "_res", res, "_norm.pdf"), height = 3.5, width = 4.5)



