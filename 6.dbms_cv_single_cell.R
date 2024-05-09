source("DBMS/WPCMS.R")
source("DBMS/PCMS.R")
source("DBMS/DBMS.R")
source("DBMS/distributions.R")
source("functions.R")
library(dplyr)
library(ggplot2)

############# Load data ####################

get_cell = function(num){
  df = as.matrix(read.csv(paste('Data/Single cell/sparse_contact_cell', num, '.csv', sep = ''), header = FALSE))
  n = df[1,1]
  C = matrix(0, n, n)    
  C[df[-1,1:2]] = df[-1,3]
  C = C + t(C)
  diag(C) = diag(C)/2
  return(C)
}

bulk = function(seed){
  set.seed(seed)
  train_index = sample(1:8, 4)
  test_index = (1:8)[-train_index]
  C_train = 0
  C_test = 0
  for(i in 1:4){
    C_train = C_train +  get_cell(train_index[i])
    C_test = C_test +  get_cell(test_index[i])
  }
  index = which(rowSums(C_train) != 0)
  return(list(index = index, 
              train = C_train[index, index], 
              test = C_test[index, index]))
}

############# Find test scores for PoisMS, HPoisMS, ZIPoisMS, NBMS ####################

update_info = function(C, H, method, param, Theta0, seed, info, method_name){
  #get reconstruction
  sol = DBMS(C$train, H, param = param, method = method, Theta = Theta0, eps_wpcms = 1e-5, eps_dbms = 1e-5, verbose_dbms = TRUE)
  
  if(method_name == "PoisMS") par = NA
  if(method_name %in% c("HPoisMS", "ZIPoisMS")) par =  sol$param$p
  if(method_name == "NBMS") par = sol$param$r
  
  #update info
  rbind(info, data.frame("df" = ncol(H), seed, "beta" = sol$param$beta, "param" = par, "epoch" = sol$epoch, "iter_total" = sol$iter_total, 
                         rbind(dbms_evaluate(sol$X, sol$param, C$train, method, method_name),
                               dbms_evaluate(sol$X, sol$param, C$test, method, method_name)), 
                         data = c("train", "test"),
                         "method" = method_name))
}

#bootstrap
dfs = seq(5, 50, 5)
B = 30
info = data.frame()

for(seed in 1:B){
  C = bulk(seed)
  for(df in dfs){
    H = load_H(df, C$index, orth = TRUE)
    set.seed(seed)
    Theta0 = matrix(rnorm(df * 3), df, 3) 
    
    cat("=======================================\nseed :",
        seed, "df :", df, "method : PoisMS\n=======================================\n")
     param = list(beta = log(mean(C$train[C$train>0])))
    info = update_info(C, H, pois, param, Theta0, seed, info, "PoisMS")

    cat("=======================================\nseed :",
        seed, "df :", df, "method : HPoisMS\n=======================================\n")
    param$p = mean(C$train == 0)
    info = update_info(C, H, hpois, param, Theta0, seed, info, "HPoisMS")

    cat("=======================================\nseed :",
        seed, "df :", df, "method : ZIPoisMS\n=======================================\n")
    info = update_info(C, H, zipois, param, Theta0, seed, info, "ZIPoisMS")

    cat("=======================================\nseed :", 
        seed, "df :", df, "method : NBMS\n=======================================\n")
    param$p = NULL
    param$r = 1
    info = update_info(C, H, nbin, param, Theta0, seed, info, "NBMS")
  
    write.csv(info, paste('Results/dbms/cv_single_cell.csv', sep = ''), row.names = FALSE)
  }
}

############# Plot test scores (Figure 9, right) ####################

info = read.csv("Results/dbms/cv_single_cell.csv")
info_sum = info %>% group_by(df, method, data) %>%
  summarize(upper_loglik = quantile(loglik, 0.95),
            upper_cor = quantile(spcorE, 0.95),
            upper_mse = quantile(mseE, 0.95),
            lower_loglik = quantile(loglik, 0.05),
            lower_cor = quantile(spcorE, 0.05),
            lower_mse = quantile(mseE, 0.05),
            loglik = mean(loglik),
            cor = mean(spcorE),
            mse = mean(mseE)) %>%
  mutate(method = factor(method, levels = c("PoisMS", "HPoisMS", "ZIPoisMS", "NBMS")))

info_min = info_sum %>% filter(data == "test") %>%
  group_by(method) %>% 
  slice(which.min(loglik))

info_sum %>% filter(data == "test") %>%
ggplot(aes(df, loglik, color = method, fill = method, linetype = method))+
  geom_ribbon(aes(ymax = upper_loglik, ymin = lower_loglik), alpha = 0.1, color = NA, linetype = "dashed")+
  geom_line()+
  geom_point()+
  geom_point(info_min, mapping = aes(df, loglik, color = method), size = 3, shape = 4)+
  scale_color_manual(values = cbPalette[c(4, 2, 3, 7)])+
  scale_fill_manual(values = cbPalette[c(4, 2, 3, 7)])+
  ylab("test loss")+
  xlab("degrees-of-freedom")+
  theme_bw()
ggsave(paste0("Plots/dbms/cv_single_cell.pdf"), height = 3, width = 4)

