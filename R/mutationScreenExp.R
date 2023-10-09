# From each adaptive step, generate some mutations and get the distribution of 
# fitness effects

.createEffectDataframeAdd <- function(dat, sampled_effects) {
  dat <- dat %>% filter(modelindex == 1)
  
  dat$rowID <- as.integer(rownames(dat))
  
  dat_out <- dat[rep(seq_len(nrow(dat)), each = length(sampled_effects)),]
  rownames(dat_out) <- NULL # reset row names
  
  # calculate fitness with sampled effects added
  Aa <- calcAddFitness(dat_out$fixEffectSum + sampled_effects, 2, 0.05)
  AA <- calcAddFitness(dat_out$fixEffectSum + sampled_effects * 2, 2, 0.05)
  aa <- calcAddFitness(dat_out$fixEffectSum, 2, 0.05)
  
  dat <- dat_out
  dat$avFit <- Aa - aa
  dat$avFit_AA <- AA - aa
  dat$value_AA <- dat$value * 2
  dat$wAA <- AA
  dat$wAa <- Aa
  dat$waa <- aa
  dat$s <- AA - aa
  return(dat)
}

.createEffectDataframe <- function(dat, sampled_effects) {
  dat <- dat %>% filter(modelindex == 2)
  
  # AA = 1+s; Aa = 1+hs; aa = 1
  # AA = d_popfx$fixEffectSum + 2 * value
  # Aa = d_popfx$fixEffectSum + value
  # aa = d_popfx$fixEffectSum
  
  # Calculate base phenotype/fitness with only fixations (reference)
  dat$rowID <- as.integer(rownames(dat))
  
  write.table(dat %>% ungroup() %>% mutate(KZ = 1, KXZ = 1) %>%
                dplyr::select(rowID, fixEffectSum_aZ, fixEffectSum_bZ, KZ, KXZ), 
              "d_grid.csv", sep = ",", col.names = F, row.names = F)
  d_base <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)
  
  # Duplicate rows so that each has the sampled additional effect
  dat_out <- dat[rep(seq_len(nrow(dat)), each = length(sampled_effects)),]
  rownames(dat_out) <- NULL # reset row names
  
  d_popfx <- d_base[rep(seq_len(nrow(d_base)), each = length(sampled_effects)*2),]
  rownames(d_popfx) <- NULL # reset row names
  
  # We need to add to aZ/bZ separately
  d_dat_withFX_aZ <- dat_out
  d_dat_withFX_bZ <- dat_out

  # Now add a random effect: calculate separately for alpha/beta
  d_dat_withFX_aZ$aZ <- exp(log(d_dat_withFX_aZ$fixEffectSum_aZ) + sampled_effects)
  d_dat_withFX_aZ$bZ <- exp(log(d_dat_withFX_aZ$fixEffectSum_bZ))
  d_dat_withFX_bZ$aZ <- exp(log(d_dat_withFX_bZ$fixEffectSum_aZ))
  d_dat_withFX_bZ$bZ <- exp(log(d_dat_withFX_bZ$fixEffectSum_bZ) + sampled_effects)

  # Homozygous estimation - for dominance calculation
  d_dat_withFX_aZ$aZ_AA <- exp(log(d_dat_withFX_aZ$fixEffectSum_aZ) + 2 * sampled_effects)
  d_dat_withFX_aZ$bZ_AA <- exp(log(d_dat_withFX_aZ$fixEffectSum_bZ))
  d_dat_withFX_bZ$aZ_AA <- exp(log(d_dat_withFX_bZ$fixEffectSum_aZ))
  d_dat_withFX_bZ$bZ_AA <- exp(log(d_dat_withFX_bZ$fixEffectSum_bZ) + 2 * sampled_effects)
  
  d_dat_withFX <- rbind(d_dat_withFX_aZ, d_dat_withFX_bZ)
  
  # Add KZ/KXZ values
  d_dat_withFX$KZ <- 1
  d_dat_withFX$KXZ <- 1
  
  write.table(d_dat_withFX %>% ungroup() %>% 
                dplyr::select(rowID, aZ, bZ, KZ, KXZ), 
              "d_grid.csv", sep = ",", col.names = F, row.names = F)
  Aa <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)
  
  write.table(d_dat_withFX %>% ungroup() %>% 
                dplyr::select(rowID, aZ_AA, bZ_AA, KZ, KXZ), 
              "d_grid.csv", sep = ",", col.names = F, row.names = F)
  AA <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)
  
  # Get the effect size by taking away the phenotype missing that fixation
  # Ensure that the tables are aligned by id before we join them
  dat <- d_dat_withFX %>% arrange(rowID)
  Aa <- Aa %>% arrange(id)
  AA <- AA %>% arrange(id)
  
  dat$avFX <- Aa$pheno - d_popfx$pheno
  dat$avFit <- Aa$fitness - d_popfx$fitness
  dat$avFX_AA <- AA$pheno - d_popfx$pheno
  dat$avFit_AA <- AA$fitness - d_popfx$fitness
  dat$wAA <- AA$fitness
  dat$wAa <- Aa$fitness
  dat$waa <- d_popfx$fitness
  dat$s <- dat$wAA - dat$waa
  return(dat)
}

MutationScreenExp <- function(fixed, n, isAdditive = F) {
  # at each step in the walk, sample n mutations for each molecular component
  # and add them to the ancestor (the previous step) - then add that to a dataframe
  
  # sample effects
  fx <- rnorm(n)
  
  # Run mutation sweep
  if (isAdditive) {
    return(.createEffectDataframeAdd(fixed, fx))
  }
  
  .createEffectDataframe(fixed, fx)
}

# seed <- sample(1:.Machine$integer.max, 1)
# sampled 18799215
seed <- 18799215
set.seed(seed)
mutExp <- MutationScreenExp(d_fix_ranked, 1000)

mutExp_add <- MutationScreenExp(d_fix_ranked_add, 1000, T)

mutExp_combined <- rbind(mutExp, mutExp_add)
mutExp_combined$model <- ifelse(mutExp_combined$modelindex == 1, "Additive", "NAR")

# percentage of mutations beneficial, expected waiting time for beneficial mutation,
# mean effect of mutations
mutExp %>% 
  group_by(rankFactor) %>%
  summarise(percBeneficial = mean(as.integer(s > 0)),
            CIperc = CI(as.integer(s > 0)),
            meanEffect = mean(s),
            CIEffect = CI(s)) -> mutExp_sum

mutExp %>%
  group_by(rankFactor, seed) %>%
  summarise(percBeneficial = mean(as.integer(s > 0))) -> mutExp_percBeneficial

mutExp_add %>% 
  group_by(rankFactor) %>%
  summarise(percBeneficial = mean(as.integer(s > 0)),
            CIperc = CI(as.integer(s > 0)),
            meanEffect = mean(s),
            CIEffect = CI(s)) -> mutExp_add_sum

mutExp_add %>%
  group_by(rankFactor, seed) %>%
  summarise(percBeneficial = mean(as.integer(s > 0))) -> mutExp_add_percBeneficial


mutExp_add_sum$model <- "Additive"
mutExp_sum$model <- "NAR"

mutExp_add_percBeneficial$model <- "Additive"
mutExp_percBeneficial$model <- "NAR"


mutExp_sum_combined <- rbind(mutExp_add_sum, mutExp_sum)
mutExp_perc_combined <- rbind(mutExp_add_percBeneficial, mutExp_percBeneficial)
