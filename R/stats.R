library(tidyverse)
library(fitdistrplus)
library(multimode)
library(boot)
library(moments)
library(emmeans)
library(future)
library(doParallel)
library(foreach)

# Number of pops adapted
print(d_qg %>% group_by(modelindex) %>%
  filter(gen == 49500) %>%
  summarise(n = n(),
            pAdapted = mean(isAdapted),
            CIAdapted = CI(isAdapted)))

print(d_qg %>% group_by(modelindex, isAdapted) %>%
        filter(gen == 49500) %>%
        summarise(n()))
        

# Fisher test - number of pops adapted
print(fisher.test(table((d_qg %>% distinct(seed, modelindex, .keep_all = T))$modelindex, 
                  (d_qg %>% distinct(seed, modelindex, .keep_all = T))$isAdapted)))


# Fisher test - number of steps taken in walk
print(fisher.test(table((d_fix_ranked_combined %>% filter(rank > 0))$model, 
      (d_fix_ranked_combined %>% filter(rank > 0))$rank)))

# mean number of steps
print(d_fix_ranked_combined %>% 
  group_by(model, seed) %>%
  summarise(nSteps = max(rank)) %>%
  ungroup() %>%
  summarise(meanStepsAcrossBothModels = mean(nSteps),
            CISteps = CI(nSteps)))

print(d_fix_ranked_combined %>% 
  group_by(model, seed) %>%
  summarise(nSteps = max(rank)) %>%
  group_by(nSteps) %>%
  summarise(n_nSteps = n(),
            nRows = nrow(.),
            percentHadnSteps = n_nSteps/nRows))

# Difference in mean generation for an adaptive step
# for populations with more than one step and haven't reached the optimum
d_fix_ranked_combined %>% 
  filter(rank > 1, phenomean < 1.9 | phenomean > 2.1) %>%
  mutate(gen = gen - 50000) -> d_notYetAdapted

print(wilcox.test(gen ~ model, data = d_notYetAdapted, conf.int = T))

# fit gamma distribution to fixed effects
fit_nar <- fitdist((d_fix_ranked %>% filter(rank > 0, s > 0))$s, 
                   distr = "gamma", method = "mle")
fit_add <- fitdist((d_fix_ranked_add %>% filter(rank > 0, s > 0))$s, 
                   distr = "gamma", method = "mle")

# 95% CI
print(fit_nar$sd * qnorm(0.975))
print(fit_add$sd * qnorm(0.975))
print(summary(fit_nar))
print(summary(fit_add))
#plot(fit_add)
#plot(fit_nar)

# Find modes of mutation screen distribution
d <- density(mutExp$s)

modes <- function(d){
  i <- which(diff(sign(diff(d$y))) < 0) + 1
  data.frame(x = d$x[i], y = d$y[i])
}

mutExp_modes <- modes(d)
print(mutExp_modes[order(mutExp_modes$y, decreasing = T),])

d_add <- density(mutExp_add$s)
mutExp_add_modes <- modes(d_add)
print(mutExp_add_modes[order(mutExp_add_modes$y, decreasing = T),])

# Heterozygosity: difference between models
## Average over all generations
print(d_het %>% pivot_longer(cols = c(Ho_l1, Ho_l2), names_to = "locus", values_to = "Ho") %>%
  group_by(modelindex) %>%
  summarise(meanHo = mean(Ho),
            CIHo = CI(Ho)))

# Maximum average
print(d_Ho_sum %>% 
  group_by(model) %>%
  filter(meanHo == max(meanHo)))

# Amount of variation driven by segregating variance
# ratio of fixed/seg effects
print(d_fix_ranked_combined %>% filter(rank > 0) %>%
  mutate(rat = AA_pheno/phenomean) %>%
  group_by(model) %>%
  summarise(meanRatioFixedToSegregating = mean(rat),
            CIRatio = CI(rat)))

# Fit generalised Pareto
mutExp_add_ben <- mutExp_add %>% filter(s > 0)
mutExp_ben <- mutExp %>% filter(s > 0)

###############################################################
# This code is modified from Lebeuf-Taylor et al. 2019 Fig. 1 source code (doi: 10.7554/eLife.45952)
# Code is licensed under CCA 4.0 https://creativecommons.org/licenses/by/4.0/
#
## Testing the shape of the distribution of beneficial mutation ####
# Set seed for waiting time
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

LR.null.dist <- calcNullDist(1000, X)

# Do in parallel
cl <- makeCluster(future::availableCores())
registerDoParallel(cl)

# 10,000 replicates
B = 10000

bootBeisel_nar <- data.frame(tau = numeric(B),
                             kappa = numeric(B),
                             LRT = numeric(B),
                             p.value = numeric(B))

bootBeisel_nar <- foreach(b=1:B, .combine = "rbind") %dopar% {
  require(GenSA)
  d <- sampleFitnessEffects(mutExp_ben$s, 1000)
  bootBeisel_nar[b, ] <- bootBeisel(d, LR.null.dist)
}

bootBeisel_nar <- as.data.frame(bootBeisel_nar)
colnames(bootBeisel_nar) <- c("tau", "kappa", "LRT", "p.value")

write.csv(bootBeisel_nar, "bootBeisel_nar.csv", row.names = F)

bootBeisel_add <- data.frame(tau = numeric(B),
                             kappa = numeric(B),
                             LRT = numeric(B),
                             p.value = numeric(B))

# additive
X <- sampleFitnessEffects(mutExp_add_ben$s, 1000)
LR.null.dist <- calcNullDist(1000, X)

#Run in parallel
bootBeisel_add <- foreach(b=1:B, .combine = "rbind") %dopar% {
  require(GenSA)
  d <- sampleFitnessEffects(mutExp_add_ben$s, 1000)
  bootBeisel_add[b, ] <- bootBeisel(d, LR.null.dist)
}
stopCluster(cl)

bootBeisel_add <- as.data.frame(bootBeisel_add)
colnames(bootBeisel_add) <- c("tau", "kappa", "LRT", "p.value")

write.csv(bootBeisel_add, "bootBeisel_add.csv", row.names = F)

# calculate mean and CI of bootstrap
print("bootstrapped GPD fit mean parameters")
print(mean(bootBeisel_nar$kappa))
print(CI(bootBeisel_nar$kappa))
print(mean(bootBeisel_nar$p.value))
print(CI(bootBeisel_nar$p.value))
print(mean(bootBeisel_nar$LRT))
print(CI(bootBeisel_nar$LRT))

print(mean(bootBeisel_add$kappa))
print(CI(bootBeisel_add$kappa))
print(mean(bootBeisel_add$p.value))
print(CI(bootBeisel_add$p.value))
print(mean(bootBeisel_add$LRT))
print(CI(bootBeisel_add$LRT))

bootBeisel_nar <- read.csv("bootBeisel_nar.csv")
bootBeisel_add <- read.csv("bootBeisel_add.csv")

# Combine p values with Fisher's method
fisherMethod <- function(p, returnStat = F) {
  # calculate chi for vector p
  X <- -2*sum(log(p))
  if (!returnStat)
    return(1-pchisq(X, 2*length(p), lower.tail = T))
  return(paste(X, "df =", 2*length(p)))
}

print(bootBeisel_add %>%
  mutate(p.value = ifelse(p.value == 0, p.value + 1e-4, p.value)) %>% # Since we did 10000 iterations per model, the minimum p-value is actually 0.0001
  summarise(fisherP = fisherMethod(p.value),
            fisherX = fisherMethod(p.value, T)))

print(bootBeisel_nar %>%
  mutate(p.value = ifelse(p.value == 0, p.value + 1e-4, p.value)) %>% # Since we did 10000 iterations per model, the minimum p-value is actually 0.0001
  summarise(fisherP = fisherMethod(p.value),
            fisherX = fisherMethod(p.value, T)))


# Now run a similar analysis but instead of sampling from a joint distribution,
# sample individually within each simulation.
source("./R/beisel_nonpooled.R")
##############################################################################
# End Lebeuf-Taylor et al. derived code
##############################################################################


# Compare means - proportion of beneficial muts, waiting time to beneficial mut
# waiting time to beneficial mut, compare means
mutExp_combined %>% 
  group_by(seed, model, rankFactor) %>%
  summarise(percBeneficial = mean(as.integer(s > 0)),
            waitingTime = 1/(20000 * (9.1528*10^-6) * percBeneficial)) -> mutExp_wt

percBen_lm <- lm(percBeneficial ~ model * rankFactor, 
                   mutExp_wt %>% filter(is.finite(percBeneficial)))
print(summary(percBen_lm))
# plot(percBen_lm)
print(sqrt(diag(vcov(percBen_lm))) * qnorm(0.975))

# Marginal means
percBen_em <- emmeans(percBen_lm, pairwise ~ model | rankFactor)
# emmip(percBen_lm, ~ model | rankFactor, CIs = T)
# percBen_em$emmeans
# percBen_em$contrasts
print(summary(pairs(regrid(percBen_em)), by = NULL) %>% 
        mutate(CI = SE * qnorm(0.975)))

# Bootstrap expected waiting times 
waitingTimeDiff <- function(data, n) {
  data <- data %>% filter(is.finite(waitingTime))
  samples_add <- sample((data %>% filter(model == "Additive"))$waitingTime, n, replace = T)
  samples_net <- sample((data %>% filter(model == "NAR"))$waitingTime, n, replace = T)
  
  result <- data.frame(
    diff = samples_net - samples_add,
    sample_net = samples_net,
    sample_add = samples_add
  )
  return(result)
}
# Set seed for waiting time
#seed <- sample(1:.Machine$integer.max, 1)
# sampled 716691235
seed <- 716691235
set.seed(seed)
d_waitingTime <- waitingTimeDiff(mutExp_wt, 100000)
t.test(d_waitingTime$diff)

print(d_waitingTime %>%
  summarise(meanNARWaitingTime = mean(sample_net),
            meanAddWaitingTime = mean(sample_add),
            CINARWaitingTime = CI(sample_net),
            CIAddWaitingTime = CI(sample_add),
            meanDiffInWaitingTime = mean(diff),
            CIDiffInWaitingTime = CI(diff)))


# alpha/beta regression
aZbZ_lm <- lm(pheno ~ aZbZ, data = plt_aZbZratio$data)
print(summary(aZbZ_lm))

# Optimum phenotype - what bZaZ ratio gives exactly 2: from the above figure,
# it's around 0.8
opt_pheno_ratio <- genRatioLandscapeData(0.8, 0.81)
opt_pheno_ratio$bZaZ <- opt_pheno_ratio$bZ / opt_pheno_ratio$aZ
print(opt_pheno_ratio[order(opt_pheno_ratio$fitness, decreasing = T)[1:2],])


# usage of molecular components (phi)
print(mean(d_molCompDiff$molCompDiff))
print(CI(d_molCompDiff$molCompDiff))

# percentage of models with more than 1 step that used only alpha, only beta, or both
print(d_fix_ranked %>%
  mutate(value_aZ = if_else(mutType == 3, value, 0),
         value_bZ = if_else(mutType == 4, value, 0)) %>%
  group_by(seed) %>%
  filter(n() > 2) %>% # exclude groups with less than 2 steps
  mutate(evoBybZ = all(value_aZ == 0, na.rm = T),
         evoByaZ = all(value_bZ == 0, na.rm = T)) %>%
  ungroup() %>%
  distinct(seed, .keep_all = T) %>% 
  summarise(percEvoByaZ = mean(evoByaZ),
            percEvoBybZ = mean(evoBybZ),
            percEvoByBoth = 1 - (percEvoByaZ + percEvoBybZ),
            countEvoByaZ = sum(evoByaZ),
            countEvoBybZ = sum(evoBybZ),
            countEvoByBoth = n() - (countEvoByaZ + countEvoBybZ)))
