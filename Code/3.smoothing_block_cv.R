source("PoisMS/PCMS.R")
source("PoisMS/WPCMS.R")
source("PoisMS/PoisMS.R")
source("functions.R")
library(matlib)
library(htmlwidgets)
library(tidyr)
library(dplyr)
library(plotly)
library(ggplot2)
library(fields)
library(splines)
cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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

############# Find BCV scores for PoisMS and SPoisMS ####################

dr = DR(index)
K = dr$K
plot(dr$s)
k = 300

train_data = function(C, H, unobs){
  Ctrain = C[-unobs, -unobs]
  nonzero = abs(colSums(H[-unobs,])) > 0
  QR = qr(H[-unobs, nonzero])
  Q = qr.Q(QR)
  R = qr.R(QR)
  Htrain = Q
  A = H[unobs, nonzero] %*% solve(R) %*% t(Q)
  return(list(C = Ctrain, H = Htrain, A = A))
}

train_data_smooth = function(C, H, K, unobs){
  Ctrain = C[-unobs, -unobs]
  Htrain = qr.Q(qr(H[-unobs, ]))
  K11 = K[-unobs, -unobs]
  K12 = K[-unobs, unobs]
  K21 = K[unobs, -unobs]
  K22 = K[unobs, unobs]
  Ktrain = K11 - K12 %*% solve(K22) %*% K21
  A = -solve(K22) %*% K21
  return(list(C = Ctrain, H = Htrain, K = Ktrain, A = A))
}

test_loss = function(C, X, beta, unobs){
  W = matrix(NA, ncol(C), nrow(C))
  W[unobs,] = 1
  W[,unobs] = 1
  D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
  logL = -D^2 + beta
  return(mean(W * (exp(logL) - C * logL), na.rm = T))
}

info = read.csv(paste0("Results/smoothing/dfs_chr", chr, "_res", res, ".csv"))
split = ceiling((1:n) / n*20)
test_scores = c()
Theta0 = matrix(0, k, 3)
Theta0[1:3, 1:3] = matrix(rnorm(9), 3, 3)

for(seed in 1:20){
  cat("\n\nseed:", seed)
  set.seed(0)
  #random initialization of DR
  unobs = which(split == seed)
  for(i in 1:6){
    
    #smoothing spline method
    lambda = info$lambda[i]
    cat("\nlambda:", lambda)
    H = dr$H[,1:k]
    train = train_data_smooth(C, H, K, unobs)
    spoisms = list(Theta = Theta0, beta = log(mean(train$C)))
    spoisms = PoisMS(train$C - lambda/2 * length(train$C) * train$K, train$H,
                     beta0 = spoisms$beta, Theta0 = spoisms$Theta)
    
    #evaluate
    X = matrix(NA, n, 3)
    X[-unobs,] = spoisms$X
    X[unobs,] = train$A %*% spoisms$X
    matplot(X, type = "l")
    test_scores = rbind(test_scores, data.frame(beta = spoisms$beta, score = test_loss(C, X, spoisms$beta, unobs), seed, param = lambda, method = "SPoisMS"))
    
    #B-spline method
    df = round(info$df[i])
    cat(" df:", df)
    H = load_H(df, index, orth = FALSE)
    train = train_data(C, H, unobs)
    poisms = PoisMS(train$C, train$H, beta0 = log(mean(train$C)), Theta0 = Theta0[1:ncol(train$H),])
    
    #evaluate
    X = matrix(NA, n, 3)
    X[-unobs,] = poisms$X
    X[unobs,] = train$A %*% poisms$X
    matplot(X, type = "l")
    test_scores = rbind(test_scores, data.frame(beta = poisms$beta, score = test_loss(C, X, poisms$beta, unobs), seed, param = df, method = "PoisMS"))
    
    write.csv(test_scores, paste0("Results/smoothing/bcv_chr", chr, "_res", res, ".csv"), row.names = F)
  }
}


############# Plot BCV scores for PoisMS and SPoisMS (Figure 5, right) ####################

test_scores = read.csv(paste0("Results/smoothing/bcv_chr", chr, "_res", res, ".csv"))
test_scores_sum = test_scores %>% mutate(lambda = plyr::mapvalues(param, round(info$df), info$lambda)) %>% 
  group_by(method, lambda) %>%
  filter(seed > 1 & seed < 20) %>%
  summarise(mean = mean(score), sd = sd(score)/sqrt(18))
test_scores_min = test_scores_sum %>% group_by(method) %>% 
  slice(which.min(mean))

ggplot(test_scores_sum, aes(x = log(lambda, 10), y = mean, color = method, fill = method, linetype = method))+
  geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd), alpha = 0.2, color = NA, linetype = "dashed")+
  geom_line()+
  geom_point()+
  geom_point(test_scores_min, mapping = aes(log(lambda, 10), mean, color = method), size = 3, shape = 4)+
  ylab("test loss")+
  scale_x_continuous(name = 'degrees-of-freedom',
                     trans='reverse',
                     sec.axis = sec_axis(~ . , name=bquote(log[10](lambda))),
                     breaks = log(info$lambda,10),
                     labels = round(info$df))+
  scale_color_manual(values = cbPalette[c(4, 8)])+
  scale_fill_manual(values = cbPalette[c(4, 8)])+
  coord_cartesian(ylim = c(-55, -40))+
  theme_bw()
ggsave(paste0("Plots/smoothing/bcv_chr", chr, "_res", res, ".pdf"), height = 3, width = 4)
