library(future)
library(doParallel)
library(foreach)
# Test % of K in >0, <0 and ~ 0 groups
###############################################################
# This code is modified from Lebeuf-Taylor et al. 2019 Fig. 1 source code (doi: 10.7554/eLife.45952)
# Code is licensed under CCA 4.0 https://creativecommons.org/licenses/by/4.0/
#
## Testing the shape of the distribution of beneficial mutation ####
#seed <- sample(1:.Machine$integer.max, 1)
# sampled 243945692
seed <- 243945692
set.seed(seed)

# sample fitness effects and return a sorted vector
sampleFitnessEffects <- function(s, n) {
  # Measure fitness relative to the smallest beneficial selection coefficient
  # as recommended by Beisel et al. 2007
  X = s-min(s)
  X = X[X!=0]
  X = sample(X, n)
  X = sort(X)
  return(X)
}

X <- sampleFitnessEffects(mutExp_ben$s, 1000)

LL.GPD <- function(par, d){ # equations from Beisel et al 2007
  # X is the adjusted selection coefficients, i.e. the data
  X = d
  # tau is a parameter of the GPD
  tau = par[1]
  # kappa is a parameter of the GDP
  kappa = par[2]
  # n_1 is the number of observations
  n_1 = length(X)
  LL.a = -(n_1)*log(tau)
  LL.b = NULL
  # for kappa>0
  if(kappa>0) {
    LL.b = -(kappa+1)/kappa * sum (log(1+(kappa*X)/tau))
  }
  # for kappa<0
  if(kappa<0) {
    if(sum(X>-tau/kappa)){LL.b = -1000000}
    else{LL.b = -(kappa+1)/kappa * sum (log(1+(kappa*X)/tau))}
  }
  # for kappa=0
  if(kappa==0) {
    LL.b = -(1/tau)*sum(X)
  }
  LL = LL.a + LL.b
  return(-LL)
}


## Special case kappa=0 (exponential distribution)
LL.Exp <- function(par, d){ # equations from Beisel et al 2007
  # X is the adjusted selection coefficients, i.e. the data
  X = d
  # tau is a parameter of the GPD
  tau = par
  # kappa is a parameter of the GDP
  kappa = 0
  
  # n_1 is the number of observations
  n_1 = length(X)
  LL.a = -(n_1)*log(tau)
  LL.b = -(1/tau)*sum(X)
  LL = LL.a + LL.b
  return(-LL)
}

#1# Optimization ####
require(GenSA)

bootBeisel <- function(d, nullLR) {
  start.par = c(.1, 0) # starting parameter values for the optimization
  ## 1) Find tau, When kappa = 0
  Exp.opt = GenSA(par = start.par[1], fn = LL.Exp, lower = 0.000001, upper = 100, d = d)
  Exp.opt.tau = Exp.opt$par
  LL.Exp.opt = -Exp.opt$value
  ## 2) Find best tau and kappa
  GPD.opt = GenSA(par = start.par, fn = LL.GPD, lower = c(0.000001, -100), upper = c(100, 100), d = d)
  GPD.opt.tau = GPD.opt$par[1]
  GPD.opt.kappa = GPD.opt$par[2]
  LL.GPD.opt = -GPD.opt$value
  
  ## Likelihood ratio ####
  LRT = LL.Exp.opt/LL.GPD.opt
  
  P.value = sum(LRT>nullLR)/length(d)
  return(c(GPD.opt.tau, GPD.opt.kappa, LRT, P.value))
  
}

## Find a null distribution of likelihood ratio ####
calcNullDist <- function(B, X) {
  start.par = c(.1, 0) # starting parameter values for the optimization
  ## 1) Find tau, When kappa = 0 for d = X
  Exp.opt = GenSA(par = start.par[1], fn = LL.Exp, lower = 0.000001, upper = 100, d = X)
  Exp.opt.tau = Exp.opt$par
  LL.Exp.opt = -Exp.opt$value

# Calculate null distribution (gumbel, K = 0)
  LR.null.dist = vector(length = B)
  for(b in 1:B){
    X = rexp(n = length(X), rate = 1/Exp.opt.tau)
    Exp.null.opt = optim(par = start.par[1], fn = LL.Exp, method = "Brent", lower = 0.000001, upper = 100, d = X)
    LL.Exp.null.opt = -Exp.null.opt$value
    ## 2) Find best tau and kappa
    GPD.null.opt = optim(par = c(Exp.opt.tau, 0), fn = LL.GPD, method = "L-BFGS-B", lower = c(0.000001, -100), upper = c(100, 100), d = X)
    # print(b)
    # GPD.null.opt = GenSA(par = start.par, fn = LL.GPD, lower = c(0.000001, -100), upper = c(100, 100))
    LL.GPD.null.opt = -GPD.null.opt$value
    ## Likelihood ratio
    LR.null.dist[b] = LL.Exp.null.opt/LL.GPD.null.opt
  }
  return(LR.null.dist)
}

cl <- makeCluster(future::availableCores())
registerDoParallel(cl)
# Change B to be number of replicates
B = length(unique(mutExp_ben$seed))
bootBeisel_nar_nonpool <- data.frame(tau = numeric(B),
                             kappa = numeric(B),
                             LRT = numeric(B),
                             p.value = numeric(B))

# Subset by replicate
mutExp_ben_seeds <- unique(mutExp_ben$seed)
bootBeisel_nar_nonpool <- foreach(b=1:B, .combine = "rbind") %dopar% {
  require(GenSA)
  # Sample from this replicate's data
  sbst <- mutExp_ben[mutExp_ben$seed == mutExp_ben_seeds[b],]
  samples_n <- min(nrow(sbst), 100)
  d <- sampleFitnessEffects(sbst$s, samples_n)
  # Get the null distribution for this d
  LR.null.dist <- calcNullDist(samples_n, X)
  bootBeisel_nar_nonpool[b, ] <- bootBeisel(d, LR.null.dist)
}

bootBeisel_nar_nonpool <- as.data.frame(bootBeisel_nar_nonpool)
colnames(bootBeisel_nar_nonpool) <- c("tau", "kappa", "LRT", "p.value")
write.csv(bootBeisel_nar_nonpool, "bootBeisel_nar_nonpool.csv", row.names = F)


bootBeisel_add_nonpool <- data.frame(tau = numeric(B),
                             kappa = numeric(B),
                             LRT = numeric(B),
                             p.value = numeric(B))


# additive
X <- sampleFitnessEffects(mutExp_add_ben$s, 1000)
LR.null.dist <- calcNullDist(10000, X)
B = length(unique(mutExp_add_ben$seed))

#Run in parallel
mutExp_add_ben_seeds <- unique(mutExp_add_ben$seed)

bootBeisel_add_nonpool <- foreach(b=1:B, .combine = "rbind") %dopar% {
  require(GenSA)
  sbst <- mutExp_add_ben[mutExp_add_ben$seed == mutExp_add_ben_seeds[b],]
  samples_n <- min(nrow(sbst), 100)
  d <- sampleFitnessEffects(sbst$s, samples_n)
  # Get the null distribution for this d
  LR.null.dist <- calcNullDist(samples_n, X)
  bootBeisel_add_nonpool[b, ] <- bootBeisel(d, LR.null.dist)
}
stopCluster(cl)

bootBeisel_add_nonpool <- as.data.frame(bootBeisel_add_nonpool)
colnames(bootBeisel_add_nonpool) <- c("tau", "kappa", "LRT", "p.value")

write.csv(bootBeisel_add_nonpool, "bootBeisel_add_nonpool.csv", row.names = F)

# Calculate mean and CI of bootstrap
print("bootstrapped GPD fit mean parameters")
print(mean(bootBeisel_nar_nonpool$kappa))
print(CI(bootBeisel_nar_nonpool$kappa))
print(mean(bootBeisel_nar_nonpool$p.value))
print(CI(bootBeisel_nar_nonpool$p.value))
print(mean(bootBeisel_nar_nonpool$LRT))
print(CI(bootBeisel_nar_nonpool$LRT))

print(mean(bootBeisel_add_nonpool$kappa))
print(CI(bootBeisel_add_nonpool$kappa))
print(mean(bootBeisel_add_nonpool$p.value))
print(CI(bootBeisel_add_nonpool$p.value))
print(mean(bootBeisel_add_nonpool$LRT))
print(CI(bootBeisel_add_nonpool$LRT))

hist(bootBeisel_add_nonpool$kappa, breaks = 100)
hist(bootBeisel_add$kappa, breaks = 100)
hist(bootBeisel_nar_nonpool$kappa, breaks = 100)
hist(bootBeisel_nar$kappa, breaks = 100)

##############################################################################
# End Lebeuf-Taylor et al. derived code
##############################################################################
bootBeisel_add_nonpool <- read.csv("bootBeisel_add_nonpool.csv")
bootBeisel_nar_nonpool <- read.csv("bootBeisel_nar_nonpool.csv")

# percent in each bin
print(bootBeisel_add_nonpool %>% 
  summarise(percGumbel = sum(kappa < 1 & kappa > -1)/n(),
            percWeibull = sum(kappa <= -1)/n(),
            percFrechet = sum(kappa >= 1)/n()))

print(bootBeisel_nar_nonpool %>% 
  summarise(percGumbel = sum(kappa < 1 & kappa > -1)/n(),
            percWeibull = sum(kappa <= -1)/n(),
            percFrechet = sum(kappa >= 1)/n()))

print(bootBeisel_add %>% 
  summarise(percGumbel = sum(kappa < 1 & kappa > -1)/n(),
            percWeibull = sum(kappa <= -1)/n(),
            percFrechet = sum(kappa >= 1)/n()))

print(bootBeisel_nar %>% 
  summarise(percGumbel = sum(kappa < 1 & kappa > -1)/n(),
            percWeibull = sum(kappa <= -1)/n(),
            percFrechet = sum(kappa >= 1)/n()))

# Plot each fit overlaid on top of each other
bootBeisel_add$model <- "Additive"
bootBeisel_add$pooled <- "Pooled"
bootBeisel_nar$model <- "Network"
bootBeisel_nar$pooled <- "Pooled"
bootBeisel_add_nonpool$model <- "Additive"
bootBeisel_add_nonpool$pooled <- "Non-pooled"
bootBeisel_nar_nonpool$model <- "Network"
bootBeisel_nar_nonpool$pooled <- "Non-pooled"

d_bootBeisel <- rbind(bootBeisel_add, bootBeisel_nar, 
                    bootBeisel_add_nonpool, bootBeisel_nar_nonpool)

d_gpd <- data.frame(
  tau = rep(d_bootBeisel$tau, each = 1000),
  kappa = rep(d_bootBeisel$kappa, each = 1000),
  p = rep(d_bootBeisel$p.value, each = 1000),
  model = rep(d_bootBeisel$model, each = 1000),
  id = rep(row.names(d_bootBeisel), each = 1000),
  pooled = rep(d_bootBeisel$pooled, each = 1000),
  sample_index = rep(1:1000, times = nrow(d_bootBeisel))
)
d_gpd$val <- devd(rep(seq(0, 0.2, length.out = 1000), times = nrow(d_bootBeisel)), 
                  loc = 0, 
                  scale = d_gpd$tau, 
                  shape = d_gpd$kappa, type = "GP")
d_gpd$density_key <- rep(seq(0, 0.2, length.out = 1000), times = nrow(d_bootBeisel))

# Subset: 5 examples from each group
d_gpd_sbst <- d_gpd %>%
  group_by(model, pooled) %>%
  filter(id %in% sample(unique(d_gpd[d_gpd$model == cur_group()$model & 
                                       d_gpd$pooled == cur_group()$pooled,]$id), 5))

ggplot(d_gpd_sbst,
       aes(x = density_key, y = val, fill = model, group = id)) +
  facet_grid(.~pooled) +
  geom_area(alpha = 0.4, position = "identity") +
  geom_line() +
  scale_fill_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "Network")) +
  labs(x = "s", y = "Density", fill = "Model") + 
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom",
        legend.key = element_rect(colour = "black")) -> plt_beisel_pool
plt_beisel_pool
ggsave("sfig_beisel_pool.png", plt_beisel_pool, width = 8, height = 4, 
       device = png, bg = "white")
