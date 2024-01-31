# From each adaptive step, generate some mutations and get the distribution of 
# fitness effects

.createEffectDataframeAdd <- function(dat, sampled_effects) {
  dat <- dat %>% filter(modelindex == 1)
  
  dat$rowID <- as.integer(rownames(dat))
  
  dat_out <- dat[rep(seq_len(nrow(dat)), each = length(sampled_effects)),]
  rownames(dat_out) <- NULL # reset row names
  dat_out$sampledEffect <- rep(sampled_effects, length.out = nrow(dat_out))
  
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
  
  dat_out$sampledEffect <- rep(sampled_effects, length.out = nrow(dat_out))
  d_popfx$sampledEffect <- rep(sampled_effects, length.out = nrow(d_popfx))
  
  
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
  
  dat$AA_pheno <- AA$pheno
  dat$Aa_pheno <- Aa$pheno
  dat$aa_pheno <- d_popfx$pheno
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

MutationScreenExp <- function(fixed, n, model = "NAR") {
  # at each step in the walk, sample n mutations for each molecular component
  # and add them to the ancestor (the previous step) - then add that to a dataframe
  
  # sample effects
  fx <- rnorm(n)
  
  # Run mutation sweep
  if (model == "Additive") {
    return(.createEffectDataframeAdd(fixed, fx))
  }
  
  if (model == "Multiplicative") {
    return(.createEffectDataframeMult(fixed, fx))
  }
  
  .createEffectDataframe(fixed, fx)
}


# .createEffectDataframeMult <- function(dat, sampled_effects) {
#   dat <- dat %>% filter(modelindex == 2)
#   
#   dat$rowID <- as.integer(rownames(dat))
#   
#   dat_out <- dat[rep(seq_len(nrow(dat)), each = length(sampled_effects)),]
#   rownames(dat_out) <- NULL # reset row names
#   
#   dat_out$sampledEffect <- rep(sampled_effects, length.out = nrow(dat_out))
#   
#   # calculate fitness with sampled effects added
#   Aa <- calcAddFitness(exp(log(dat_out$fixEffectSum) + sampled_effects), 2, 0.05)
#   AA <- calcAddFitness(exp(log(dat_out$fixEffectSum) + 2 * sampled_effects), 2, 0.05)
#   aa <- calcAddFitness(dat_out$fixEffectSum, 2, 0.05)
#   
#   dat <- dat_out
#   dat$avFit <- Aa - aa
#   dat$avFit_AA <- AA - aa
#   dat$AA_pheno <- exp(log(dat$fixEffectSum) + 2 * sampled_effects)
#   dat$Aa_pheno <- exp(log(dat$fixEffectSum) + sampled_effects)
#   dat$aa_pheno <- dat$fixEffectSum
#   dat$avFX <- dat$AA_pheno - dat$Aa_pheno
#   dat$avFX_AA <- dat$AA_pheno - dat$aa_pheno
#   dat$value_AA <- dat$value * 2
#   dat$wAA <- AA
#   dat$wAa <- Aa
#   dat$waa <- aa
#   dat$s <- AA - aa
#   return(dat)
# }

# MutationScreenExpCompare <- function(fixed, sampledEffects, model = "NAR") {
#   # at each step in the walk, sample n mutations for each molecular component
#   # and add them to the ancestor (the previous step) - then add that to a dataframe
# 
#   # Run mutation sweep
#   if (model == "Additive") {
#     return(.createEffectDataframeAdd(fixed, sampledEffects))
#   }
#   
#   if (model == "Multiplicative") {
#     return(.createEffectDataframeMult(fixed, sampledEffects))
#   }
#   
#   .createEffectDataframe(fixed, sampledEffects)
# }


# seed <- sample(1:.Machine$integer.max, 1)
# sampled 18799215
seed <- 18799215
set.seed(seed)

mutExp_add <- MutationScreenExp(d_fix_ranked_add, 1000, "Additive")

#mutExp <- MutationScreenExp(d_fix_ranked, 1000, "NAR")
#mutExp_mult <- MutationScreenExp(d_fix_ranked_mult, 1000, "Multiplicative")

#mutExp_diff <- CalcNARMultDeviation(mutExp, mutExp_mult)

#mutExp_combined <- rbind(mutExp_diff, mutExp_add)
#mutExp_combined$model <- ifelse(mutExp_combined$modelindex == 1, "Additive", "NAR")

# Plot phenotypes for given effects mult vs NAR
mutExp <- MutationScreenExp(d_fix_ranked, 1000, "NAR")
# sampledFX <- rnorm(1000)
# mutExp <- MutationScreenExpCompare(d_fix_ranked, sampledFX, "NAR")
# mutExp_mult <- MutationScreenExpCompare(d_fix_ranked_mult, sampledFX, "Multiplicative")


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


# mutExp_mult %>% 
#   group_by(rankFactor) %>%
#   summarise(percBeneficial = mean(as.integer(s > 0)),
#             CIperc = CI(as.integer(s > 0)),
#             meanEffect = mean(s),
#             CIEffect = CI(s)) -> mutExp_mult_sum
# 
# mutExp_mult %>%
#   group_by(rankFactor, seed) %>%
#   summarise(percBeneficial = mean(as.integer(s > 0))) -> mutExp_mult_percBeneficial
# 
# mutExp_diff <- CalcNARMultDeviation(mutExp, mutExp_mult)


mutExp_add_sum$model <- "Additive"
mutExp_sum$model <- "NAR"

mutExp_add_percBeneficial$model <- "Additive"
mutExp_percBeneficial$model <- "NAR"

mutExp_sum_combined <- rbind(mutExp_add_sum, mutExp_sum)
mutExp_perc_combined <- rbind(mutExp_add_percBeneficial, mutExp_percBeneficial)

mutExp$model <- "NAR"
mutExp_add$model <- "Additive"
# mutExp_mult$model <- "Multiplicative"
# mutExp_diff$model <- "NAR-multiplicative difference"

# mutExp_combined <- rbind(mutExp, mutExp_mult, mutExp_diff, mutExp_add)
mutExp_combined <- rbind(mutExp, mutExp_add)

# ggplot(mutExp_combined,
#        aes(y = rankFactor, x = s, fill = model)) +
#   facet_grid(.~model) +
#   geom_density_ridges(alpha = 0.4, scale = 1, stat = "binline", bins = 100) +
#   scale_fill_paletteer_d("ggsci::nrc_npg") +
#   scale_y_discrete(labels = parse(text=TeX(step_labs))) +
#   labs(y = "Adaptive step", x = "Fitness effect (s)", fill = "Model") +
#   theme_bw() +
#   theme(text = element_text(size = 12), legend.position = "bottom") -> plt_effectsizerandom_time
# plt_effectsizerandom_time
# 
# ggsave("fig_mutationscreen_diff.png", plt_effectsizerandom_time, 
#        width = 10, height = 6, device = png, bg = "white")
# ggsave("fig_mutationscreen_diff.pdf", plt_effectsizerandom_time, 
#        width = 10, height = 6, bg = "white")

# # Plot regression: NAR pheno vs mult pheno
# mutExp_comp_NAR <-
#   mutExp %>%
#   rename(phenoNAR = AA_pheno) %>%
#   ungroup() %>%
#   select(rankFactor, sampledEffect, phenoNAR)
# 
# mutExp_comp_mult <-
#   mutExp_mult %>%
#   rename(phenoMult = AA_pheno) %>%
#   ungroup() %>%
#   select(rankFactor, sampledEffect, phenoMult)
# 
# mutExp_comp <- mutExp_comp_NAR
# mutExp_comp$phenoMult <- rep(mutExp_comp_mult$phenoMult, times = 2)
# 
# ggplot(mutExp_comp,
#        aes(x = phenoMult, y = phenoNAR)) +
#   facet_grid(rankFactor~.) +
#   geom_scattermost(cbind(mutExp_comp$phenoMult, mutExp_comp$phenoNAR)) +
#   scale_colour_paletteer_d("ggsci::nrc_npg") +
#   scale_y_discrete(labels = parse(text=TeX(step_labs))) +
#   labs(x = "Phenotype (Multiplicative)", y = "Phenotype (NAR)") +
#   theme_bw() +
#   theme(text = element_text(size = 12), legend.position = "bottom")

