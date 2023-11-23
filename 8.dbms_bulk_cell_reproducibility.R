source("DBMS/WPCMS.R")
source("DBMS/PCMS.R")
source("DBMS/DBMS.R")
source("DBMS/distributions.R")
source("functions.R")
library(segmented)
library(htmlwidgets)
library(tictoc)
library(hicrep)

############# Load data ####################

res = 50
chr = 21
start = 15200000
end = 48100000
ext = ".RAWobserved"
filename = paste0("Data/GM12878/", res,"kb_resolution_intrachromosomal_primary/chr", chr, "/MAPQGE30/chr", chr, "_", res, "kb", ext)
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

filename = paste0("Data/GM12878/", res,"kb_resolution_intrachromosomal_replicate/chr", chr, "/MAPQGE30/chr", chr, "_", res, "kb", ext)
data = read.csv(filename, header = FALSE, sep = "\t")
colnames(data) = c("i", "j", "count")
data = data %>% filter(i >= start & j >= start) %>%
  filter(i <= end & j <= end)
data = data %>% mutate(i = (i - start)/res/1000 + 1, j = (j-start)/res/1000 + 1)
Ctest = matrix(0, n0, n0)
Ctest[cbind(data$i, data$j)] = data$count
Ctest[cbind(data$j, data$i)] = data$count
Ctest = Ctest[index, index]
diff = Ctrain-Ctest
hist(c(diff))

############# Find train scores for PoisMS, HPoisMS, ZIPoisMS, NBMS ####################

dfs = seq(10, 100, 10)

update_info = function(Ctrain, Ctest, H, method, param, Theta0, info, method_name){
  #get reconstruction
  tic()
  sol = DBMS(Ctrain, H, param = param, method = method, Theta = Theta0, eps_wpcms = 1e-5, eps_dbms = 1e-5, verbose_dbms = TRUE)
  save = toc()
  
  if(method_name == "PoisMS") par = NA
  if(method_name %in% c("HPoisMS", "ZIPoisMS")) par =  sol$param$p
  if(method_name == "NBMS") par = sol$param$r
  
  #update info
  eval = rbind(data.frame(dbms_evaluate(sol$X, sol$param, Ctrain, method, method_name), set = "train"), 
        data.frame(dbms_evaluate(sol$X, sol$param, Ctest, method, method_name), set = "test"))
  
  return(list(info = rbind(info, data.frame("df" = ncol(H), "beta" = sol$param$beta, "param" = par, "epoch" = sol$epoch, "iter_total" = sol$iter_total, 
                          eval, time = save$toc - save$tic,
                         "method" = method_name)), X = sol$X))
}

info = data.frame()
Xs = list("PoisMS" = list(), "HPoisMS" = list(), "ZIPoisMS" = list(), "NBMS" = list())

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
    
    write.csv(info, paste0("Results/dbms/bulk_cell_repr_chr", chr, "_res", res, ".csv"), row.names = FALSE)
    saveRDS(Xs, paste0("Results/conformations/dbms_conformations_repr_chr", chr, "_res", res, ".csv"))
}


############# Plot correlations between E and O counts (Figure 12) ####################

info = read.csv(paste0("Results/dbms/bulk_cell_repr_chr", chr, "_res", res, ".csv"))
info =  info %>% mutate(method = factor(method, levels = c("PoisMS", "HPoisMS", "ZIPoisMS", "NBMS"))) %>%
  mutate(set = ifelse(set == "train", "primary", "replicate"))

info %>% dplyr::select(df, method,  corE, spcorE, strcorE, set) %>%
  filter(set == "replicate") %>%
  gather(type, cor, 3:5) %>% 
  mutate(type = ifelse(type == "corE", "Pearson", ifelse(type == "spcorE", "Spearman", "stratified"))) %>%
  mutate(method = factor(method, levels = c("PoisMS","SPoisMS", "HPoisMS", "ZIPoisMS", "NBMS"))) %>%
  ggplot(aes(df, cor, color = method))+
  geom_line()+
  geom_point()+
  scale_color_manual(values = cbPalette[c(4, 2, 3, 7)])+
  ylab("correlation")+
  xlab("degrees-of-freedom")+
  theme(legend.position="none")+
  scale_x_continuous(breaks = seq(10, 100, 10))+
  facet_wrap(~type, scales = "free")
ggsave(paste0("Plots/dbms/bulk_cell_cors_rep_chr", chr, "_res", res, ".pdf"), height = 3, width = 8)


