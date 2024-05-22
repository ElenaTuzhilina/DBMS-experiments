source("PoisMS/PCMS.R")
source("PoisMS/WPCMS.R")
source("PoisMS/PoisMS.R")
source("functions.R")
library(matlib)
library(htmlwidgets)
library(tidyr)
library(dplyr)
library(plotly)

############# Load data ####################

res = 10
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

############# Fit SPoisMS and calculate df ####################

dr = DR(index)
K = dr$K
plot(dr$s)
k = 300
H1 = dr$H[,1:k]

degrees = function(X, s, lambda){
  rank = ncol(X)
  dfs = rep(0, rank)
  for(i in 1:rank) dfs[i] = sum(1/(1 + lambda * s / sum(X[,i]^2)))
  mean(dfs)
}

loglambdas = seq(4, -1, -1)
seeds = 1

dfs = c()
angles = c()
iters = c()
conformations = list("PoisMS" = list(), "SPoisMS" = list())
betas = list("PoisMS" = list(), "SPoisMS" = list())

for(seed in 0:(seeds-1)){
  cat("\n\nseed:", seed)
  set.seed(seed)
  #random initialization of DR
  Theta0 = matrix(0, k, 3)
  #diag(Theta0) = 1
  Theta0[1:3, 1:3] = matrix(rnorm(9), 3, 3)
  spoisms = list(Theta = Theta0, beta = log(mean(C)))
  spoisms$X = H1 %*% spoisms$Theta
  for(i in 1:length(loglambdas)){
    lambda = 10^loglambdas[i]
    X0 = spoisms$X
    #smoothing spline method
    spoisms = PoisMS(C - lambda * n^2/2 * K, H1, 
                     beta0 = spoisms$beta, Theta0 = spoisms$Theta)
    df = degrees(spoisms$X, dr$s, lambda * n^2)
    cat("\nlambda:", lambda, "df:", df)
    H2 = load_H(round(df), index)
    #B-spline method
    poisms = PoisMS(C, H2, 
                    beta0 = log(mean(C)), Theta0 = matrix(rnorm(round(df) * 3), round(df), 3))
    
    #save info
    dfs = rbind(dfs, data.frame(df, lambda, seed))
    write.csv(dfs, paste0("Results/smoothing/dfs_chr", chr,"_res", res, ".csv"), row.names = F)
    
    iters = rbind(iters, data.frame(iter = c(spoisms$iter_total, poisms$iter_total), 
                                    method = c("SPoisMS", "PoisMS"), lambda, seed))
    write.csv(iters, paste0("Results/smoothing/iters_chr", chr,"_res", res, ".csv"), row.names = F)
    
    angles = rbind(angles, data.frame(index, angle = c(curvature(spoisms$X), curvature(poisms$X)),
                                      method = c(rep("SPoisMS", n), rep("PoisMS", n)), lambda, seed))
    write.csv(angles, paste0("Results/smoothing/angles_chr", chr,"_res", res, ".csv"), row.names = F)
    
    #save conformations
    conformations[["SPoisMS"]][[seed * length(loglambdas) + i]] = spoisms$X
    conformations[["PoisMS"]][[seed * length(loglambdas) + i]] = poisms$X
    betas[["SPoisMS"]][[seed * length(loglambdas) + i]] = spoisms$beta
    betas[["PoisMS"]][[seed * length(loglambdas) + i]] = poisms$beta
    saveRDS(conformations, paste0("Results/conformations/smoothing_conformations_chr", chr,"_res", res, ".rds"))
    saveRDS(betas, paste0("Results/conformations/smoothing_betas_chr", chr,"_res", res, ".rds"))
  }
}

############# Compare number of iterations required for convergence of SPoisMS and PoisMS (Figure 5, left) ####################

dfs = read.csv(paste0("Results/smoothing/dfs_chr", chr,"_res", res, ".csv"))
iters = read.csv(paste0("Results/smoothing/iters_chr", chr,"_res", res, ".csv"))
ggplot(iters %>% filter(seed == 0), 
       aes(log(lambda, 10), iter, color = method, linetype = method))+
  geom_line()+
  geom_point()+
  ylab("number of iterations")+
  theme(legend.title= element_blank())+
  scale_x_continuous(name = 'degrees-of-freedom',
                     trans='reverse',
                     sec.axis = sec_axis(~ . , name=bquote(log[10](lambda))),
                     breaks = log(dfs$lambda,10),
                                    labels = round(dfs$df))+
  scale_color_manual(values = cbPalette[c(4, 8)])+
  theme_bw()+
  theme(legend.position="none")
ggsave(paste0("Plots/smoothing/iters_chr", chr,"_res", res, ".pdf"), height = 3, width = 3)  

############# Compare reconstruction smoothness measured in angles for SPoisMS and PoisMS (Figure 14) ####################

angles = read.csv(paste0("Results/smoothing/angles_chr", chr,"_res", res, ".csv"))
plt1 = ggplot(angles %>% filter(seed == 0), 
              aes(index, angle, color = method, frame = log(lambda, 10)))+
  geom_line()+
  theme(legend.position = "none")
saveWidget(ggplotly(plt1), paste0("Plots/smoothing/angles_chr", chr,"_res", res, ".html"), selfcontained = F, libdir = "lib")


ggplot(angles %>% group_by(method, lambda, seed) %>% summarise(mean.angle = mean(angle)), 
             aes(log(lambda, 10), mean.angle, color = method, linetype = method))+
  geom_line()+
  geom_point()+
  xlab(bquote(log[10](lambda)))+
  ylab("average angle deficit")+
  theme(legend.title= element_blank())+
  scale_x_continuous(name = 'degrees-of-freedom',
                     trans='reverse',
                     sec.axis = sec_axis(~ . , name=bquote(log[10](lambda))),
                     breaks = log(dfs$lambda,10),
                     labels = round(dfs$df))+
  scale_color_manual(values = cbPalette[c(4, 8)])
  theme_bw()
ggsave(paste0("Plots/smoothing/angles_chr", chr,"_res", res, ".pdf"), height = 3.5, width = 5)  


############# Align and plot SPoisMS and PoisMS conformations (Figure 2-3) ####################

#align conformations
aligned = list("PoisMS" = list(), "SPoisMS" = list())
aligned[["SPoisMS"]][[1]] = conformations[["SPoisMS"]][[1]]
aligned[["PoisMS"]][[1]] =  align(conformations[["SPoisMS"]][[1]], conformations[["PoisMS"]][[1]])$Y
for(i in 2:length(loglambdas)){
  aligned[["SPoisMS"]][[i]] = align(aligned[["SPoisMS"]][[i - 1]], conformations[["SPoisMS"]][[i]])$Y
  aligned[["PoisMS"]][[i]] = align(aligned[["PoisMS"]][[i - 1]], conformations[["PoisMS"]][[i]])$Y
}
saveRDS(aligned, paste0("Results/conformations/smoothing_aligned_chr", chr,"_res", res, ".rds"))

#plot projections
Xs = c()
for(i in 1:length(loglambdas)){
  lambda = 10^loglambdas[i]
  Xs = rbind(Xs, data.frame(aligned[["SPoisMS"]][[i]]) %>% gather("key" = "coordinate") %>%
             mutate(index = rep(index, 3), lambda = lambda, method = "SPoisMS", seed = 0))
  Xs = rbind(Xs, data.frame(aligned[["PoisMS"]][[i]]) %>% gather("key" = "coordinate") %>%
             mutate(index = rep(index, 3), lambda = lambda, method = "PoisMS", seed = 0))
  visualize(aligned[["SPoisMS"]][[i]], "3D") %>%
    saveWidget(paste0("Plots/conformations/SPoisMS", loglambdas[i], "_chr", chr,"_res", res, ".html"), selfcontained = F, libdir = "lib")
  visualize(aligned[["PoisMS"]][[i]], "3D") %>%
    saveWidget(paste0("Plots/conformations/PoisMS", loglambdas[i], "_chr", chr,"_res", res, ".html"), selfcontained = F, libdir = "lib")
}

plt2 = ggplot(Xs, aes(index, value, color = method, frame = log(lambda, 10)))+
  geom_line()+
  facet_grid(coordinate~.)+
  theme(legend.position = "none")+
  ylab(NULL)
saveWidget(ggplotly(plt2), paste0("Plots/smoothing/coordinates_chr", chr,"_res", res, ".html"), selfcontained = F, libdir = "lib")

############# Plot expected counts for SPoisMS (Figure 6, left) ####################

conformations = readRDS(paste0("Results/conformations/smoothing_conformations_chr", chr,"_res", res, ".rds"))
betas = readRDS(paste0("Results/conformations/smoothing_betas_chr", chr,"_res", res, ".rds"))
dfs = c(5, 9, 15, 27, 47, 84)
X = conformations$SPoisMS[[which(dfs == 15)]]
beta = betas$SPoisMS[[which(dfs == 15)]]
E = expected(X, list(beta = beta), "PoisMS")
png(paste0("Plots/conformations/SPoisMS_logEheatmap_chr", chr, "_res", res, ".png"))
image.plot(log(E), axes = F, zlim = c(0,9))
axis(1, at = seq(0, n, by = 100)/n, labels = seq(0, n, by = 100))
axis(2, at = seq(0, n, by = 100)/n, labels = seq(0, n, by = 100))
dev.off()