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
#write.table(C, paste0("Data/hic_chr", chr, "_res", res, ".csv"), row.names = F, col.names = F, sep = ",")

############# Plot contact matrix (Figure 6, right) ####################

png(paste0("Plots/motivation/contact_chr", chr, "_res", res,".png"))
image.plot(ifelse(log(C) == -Inf, NA, log(C)), axes = F, zlim = c(0,9))
axis(1, at = seq(0, n, by = 100)/n, labels = seq(0, n, by = 100))
axis(2, at = seq(0, n, by = 100)/n, labels = seq(0, n, by = 100))
dev.off()

############# Plot correlation of errors (Figure 1, left) ####################

df = 25
H = load_H(df, index)
poisms = PoisMS(C, H)

#correlation of errors
L = Lambda(poisms$X, poisms$beta)
E = C - L
pcorr = cor(E, method = "spearman")
png(paste0("Plots/motivation/error_correlations_chr", chr, "_res", res, ".png"))
image.plot(pcorr, axes = F)
axis(1, at = seq(0, n, by = 100)/n, labels = seq(0, n, by = 100))
axis(2, at = seq(0, n, by = 100)/n, labels = seq(0, n, by = 100))
dev.off()

############# Plot overdispersion (Figure 7) ####################

data.frame(lambda = C[upper.tri(C)], count = L[upper.tri(L)])  %>%
  ggplot(aes(lambda, count))+
  geom_point(alpha = 0.1, color = cbPalette[6])+
  geom_smooth(color = "black", se = F, linewidth = 1)+
  xlab("expected counts")+
  ylab("observed counts")+
  theme(text = element_text(size = 15))+
  geom_abline(aes(slope = 1, intercept = 0), color = "black", linetype = "dashed", size = 1)
ggsave(paste0("Plots/motivation/overdispersion_chr", chr, "_res", res, ".png"), width = 6, height = 3.5)

############# Compute BCV scores ####################

#20-fold cross-validation
train_data = function(C, H, unobs){
  Ctrain = C[-unobs, -unobs]
  nonzero = colSums(H[-unobs,]) != 0
  QR = qr(H[-unobs, nonzero])
  Q = qr.Q(QR)
  R = qr.R(QR)
  Htrain = Q
  A = H[unobs, nonzero] %*% solve(R) %*% t(Q)
  return(list(C = Ctrain, H = Htrain, A = A))
}

test_loss = function(C, X, beta, unobs){
  W = matrix(NA, ncol(C), nrow(C))
  W[unobs,] = 1
  W[,unobs] = 1
  D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
  logL = -D^2 + beta
  return(mean(W * (exp(logL) - C * logL), na.rm = T))
}

split = ceiling((1:n) / n * 20)
dfs = seq(10, 100, 10)

info = c()
for(df in dfs){
  cat("\ndf =", df, "\n")
  for(seed in 1:20){
    cat(seed, " ")
    unobs = which(split == seed)
    H = load_H(df, index, orth = FALSE)
    train = train_data(C, H, unobs)
    k = ncol(train$H)
    poisms = PoisMS(train$C, train$H, 
                  beta0 = log(mean(train$C)),
                  Theta0 = matrix(rnorm(3 * k), k, 3))
    X = matrix(NA, n, 3)
    X[-unobs,] = poisms$X
    X[unobs,] = train$A %*% poisms$X
    matplot(X)
    beta = poisms$beta
    score = test_loss(C, X, beta, unobs)
    info = rbind(info, data.frame(score = score, beta, df, seed))
    write.table(info, paste0("Results/motivation/bcv_poisms_chr", chr, "_res", res, ".csv"), sep = ",", row.names = FALSE)
  }
}


############# Plot variation in BCV scores (Figure 1, right) ####################

info = read.csv(paste0("Results/motivation/bcv_poisms_chr", chr, "_res", res, ".csv"))

ggplot()+
  geom_boxplot(info, mapping = aes(x = df, y = score, group = df), fill = cbPalette[4], alpha = 1, size = 0.2)+
  xlab("degrees-of-freedom")+
  ylab("test loss")+
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = c(-100, 1, 10^seq(2, 14, 2)))+
  scale_x_continuous(breaks = seq(10, 100, 10))+
  geom_hline(aes(yintercept = 0), linetype = "dashed", alpha = 0.5)
ggsave(paste0("Plots/motivation/bcv_poisms_chr", chr, "_res", res, ".pdf"), height = 2.5, width = 4)

