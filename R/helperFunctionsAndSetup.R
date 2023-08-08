# Load packages
packageList <- c("dplyr", "ggplot2", "tibble", "tidyr", "cowplot", "ggridges", 
"ggpmisc", "deSolve", "DescTools", "paletteer", "latex2exp", "readr", "RColorBrewer")

lapply(packageList, require, character.only = T)

# Variables
mutType_names <- c(
  TeX("$\\alpha_Z$"),
  TeX("$\\beta_Z$")
)

# Misc functions
se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}

CI <- function(x, quantile = 0.975, na.rm = F) {
  return(qnorm(quantile) * se(x, na.rm))
}

# Fitness effect calculation functions
## Gaussian fitness function, used for both models
calcAddFitness <- function(phenotypes, optimum, width) {
  dists <- (phenotypes - optimum)^2
  return(exp(-(dists * width)))
}

## Calculate the fitness effects in additive populations
CalcAddEffects <- function(dat, isFixed = T, dat_fixed = dat) {
  # If we are calculating fitness for segregating sites, need to evaluate fitness
  # vs the fixed effect background at the timepoint (i.e. all fixations at timepoint gen)
  # Fixation effect is multiplied by 2 because diploid
  dat_fixed <- dat_fixed %>% filter(modelindex == 1)
  dat <- dat %>% 
    group_by(gen, seed) %>%
    mutate(fixEffectSum = 2 * sum(dat_fixed[dat_fixed$gen <= cur_group()$gen &
                                              dat_fixed$seed == cur_group()$seed,]$value))
  
  if (isFixed) {
    # For fixed comparisons:
    # AA = 1+s; Aa = 1+hs; aa = 1
    # AA = fixEffectSum
    # Aa = fixEffectSum - value
    # aa = fixEffectSum - 2 * value
    Aa <- calcAddFitness(dat$fixEffectSum - dat$value, 2, 0.05)
    AA <- calcAddFitness(dat$fixEffectSum, 2, 0.05)
    aa <- calcAddFitness(dat$fixEffectSum - 2 * dat$value, 2, 0.05)
    dat$AA_pheno <- dat$fixEffectSum
    dat$Aa_pheno <- dat$fixEffectSum - dat$value
    dat$aa_pheno <- dat$fixEffectSum - 2 * dat$value
    dat$avFit <- Aa - aa
    dat$avFit_AA <- AA - aa
    dat$value_AA <- dat$value * 2
    dat$wAA <- AA
    dat$wAa <- Aa
    dat$waa <- aa
    dat$s <- AA - aa
    return(dat)
  }
  
  # For segregating comparisons:
  # AA = 1+s; Aa = 1+hs; aa = 1
  # AA = fixEffectSum + 2 * value
  # Aa = fixEffectSum + value
  # aa = fixEffectSum
  # Get effect
  Aa <- calcAddFitness(dat$fixEffectSum + dat$value, 2, 0.05)
  AA <- calcAddFitness(dat$fixEffectSum + dat$value * 2, 2, 0.05)
  aa <- calcAddFitness(dat$fixEffectSum, 2, 0.05)
  
  dat$AA_pheno <- dat$fixEffectSum + 2 * dat$value
  dat$Aa_pheno <- dat$fixEffectSum + dat$value
  dat$aa_pheno <- dat$fixEffectSum
  dat$avFit <- Aa - aa
  dat$avFit_AA <- AA - aa
  dat$value_AA <- dat$value * 2
  dat$wAA <- AA
  dat$wAa <- Aa
  dat$waa <- aa
  dat$s <- AA - aa
  
  return(dat)
}

## Run the ODELandscaper tool to evaluate phenotype and fitness
## for many individuals at once.
runLandscaper <- function(df_path, output, width, optimum, threads, useID = FALSE) {
  command <- "ODELandscaper -i %s -o ./%s -w %f -p %f -t %i"
  if (useID) {
    command <- paste(command, "-I")
  }
  system(sprintf(command,
                 df_path, output, width, optimum, threads))
  result <- read_csv(paste0("./", output), col_names = F, col_types = "d")
  if (useID) {
    names(result) <- c("id", "fitness", "pheno", "aZ", "bZ", "KZ", "KXZ")
  } else {
    names(result) <- c("fitness", "pheno", "aZ", "bZ", "KZ", "KXZ")
  }
  
  return(result)
}

## Calculate fitness in network models
CalcNARPhenotypeEffects <- function(dat, isFixed = T, dat_fixed = dat) {
  dat <- dat %>% filter(modelindex == 2)
  dat_fixed <- dat_fixed %>% filter(modelindex == 2)
  
  # calculate cumulative molecular component values at each step due to only 
  # fixed effects
  # multiply by 2 because diploid
  dat <- dat %>%
    group_by(gen, seed) %>%
    mutate(fixEffectSum_aZ = 2 * sum(dat_fixed[dat_fixed$gen <= cur_group()$gen &
                                                 dat_fixed$mutType == 3 &
                                                 dat_fixed$seed == cur_group()$seed,]$value),
           fixEffectSum_bZ = 2 * sum(dat_fixed[dat_fixed$gen <= cur_group()$gen &
                                                 dat_fixed$mutType == 4 &
                                                 dat_fixed$seed == cur_group()$seed,]$value))
  # Transform to exp scale
  dat$fixEffectSum_aZ <- exp(dat$fixEffectSum_aZ)
  dat$fixEffectSum_bZ <- exp(dat$fixEffectSum_bZ)
  
  dat$rowID <- as.integer(rownames(dat))
  
  # Get phenotypes with the mutation
  write.table(dat %>% ungroup() %>%
                dplyr::select(rowID, fixEffectSum_aZ, fixEffectSum_bZ, KZ, KXZ), 
              "d_grid.csv", sep = ",", col.names = F, row.names = F)
  d_popfx <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)
  
  if (isFixed) {
    # For fixed comparisons:
    # AA = 1+s; Aa = 1+hs; aa = 1
    # AA = d_popfx$fixEffectSum
    # Aa = d_popfx$fixEffectSum - value
    # aa = d_popfx$fixEffectSum - 2 * value
    
    # Now take away the fixed effect: calculate seperately for alpha/beta
    d_dat_withoutFX_aZ <- dat %>% filter(mutType == 3)
    d_dat_withoutFX_bZ <- dat %>% filter(mutType == 4)
    d_dat_withoutFX_aZ$aZ <- exp(log(d_dat_withoutFX_aZ$fixEffectSum_aZ) - d_dat_withoutFX_aZ$value)
    d_dat_withoutFX_aZ$bZ <- exp(log(d_dat_withoutFX_aZ$fixEffectSum_bZ))
    d_dat_withoutFX_bZ$aZ <- exp(log(d_dat_withoutFX_bZ$fixEffectSum_aZ))
    d_dat_withoutFX_bZ$bZ <- exp(log(d_dat_withoutFX_bZ$fixEffectSum_bZ) - d_dat_withoutFX_bZ$value)
    
    # Homozygous estimation
    d_dat_withoutFX_aZ$aZ_aa <- exp(log(d_dat_withoutFX_aZ$fixEffectSum_aZ) - 2 * d_dat_withoutFX_aZ$value)
    d_dat_withoutFX_aZ$bZ_aa <- exp(log(d_dat_withoutFX_aZ$fixEffectSum_bZ))
    d_dat_withoutFX_bZ$aZ_aa <- exp(log(d_dat_withoutFX_bZ$fixEffectSum_aZ))
    d_dat_withoutFX_bZ$bZ_aa <- exp(log(d_dat_withoutFX_bZ$fixEffectSum_bZ) - 2 * d_dat_withoutFX_bZ$value)
    
    d_dat_withoutFX <- rbind(d_dat_withoutFX_aZ, d_dat_withoutFX_bZ)

    write.table(d_dat_withoutFX %>% ungroup() %>% 
                  dplyr::select(rowID, aZ, bZ, KZ, KXZ), 
                "d_grid.csv", sep = ",", col.names = F, row.names = F)
    Aa <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)

    write.table(d_dat_withoutFX %>% ungroup() %>% 
                  dplyr::select(rowID, aZ_aa, bZ_aa, KZ, KXZ), 
                "d_grid.csv", sep = ",", col.names = F, row.names = F)
    aa <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)

    # Get the effect size by taking away the phenotype missing that fixation
    # Ensure that the tables are aligned by id before we join them
    dat <- dat %>% arrange(rowID)
    Aa <- Aa %>% arrange(id)
    aa <- aa %>% arrange(id)
    
    dat$AA_pheno <- d_popfx$pheno
    dat$Aa_pheno <- Aa$pheno
    dat$aa_pheno <- aa$pheno
    dat$avFX <- d_popfx$pheno - Aa$pheno
    dat$avFit <- d_popfx$fitness - Aa$fitness
    dat$avFX_AA <- d_popfx$pheno - aa$pheno
    dat$avFit_AA <- d_popfx$fitness - aa$fitness
    dat$wAA <- d_popfx$fitness
    dat$wAa <- Aa$fitness
    dat$waa <- aa$fitness
    dat$s <- dat$wAA - dat$waa
    return(dat)
  }
  
  # Segregating mutation calculations
  # For segregating comparisons:
  # AA = 1+s; Aa = 1+hs; aa = 1
  # AA = d_popfx$fixEffectSum + 2 * value
  # Aa = d_popfx$fixEffectSum + value
  # aa = d_popfx$fixEffectSum
  
  # Add on the segregating effect to the fixation effects
  d_dat_withFX_aZ <- dat %>% filter(mutType == 3)
  d_dat_withFX_bZ <- dat %>% filter(mutType == 4)
  
  d_dat_withFX_aZ$aZ <- exp(log(d_dat_withFX_aZ$fixEffectSum_aZ) + d_dat_withFX_aZ$value)
  d_dat_withFX_aZ$bZ <- exp(log(d_dat_withFX_aZ$fixEffectSum_bZ))
  d_dat_withFX_bZ$aZ <- exp(log(d_dat_withFX_bZ$fixEffectSum_aZ))
  d_dat_withFX_bZ$bZ <- exp(log(d_dat_withFX_bZ$fixEffectSum_bZ) + d_dat_withFX_bZ$value)
  
  # homozygous effect
  d_dat_withFX_aZ$aZ_AA <- exp(log(d_dat_withFX_aZ$fixEffectSum_aZ) + 2 * d_dat_withFX_aZ$value)
  d_dat_withFX_aZ$bZ_AA <- exp(log(d_dat_withFX_aZ$fixEffectSum_bZ))
  d_dat_withFX_bZ$aZ_AA <- exp(log(d_dat_withFX_bZ$fixEffectSum_aZ))
  d_dat_withFX_bZ$bZ_AA <- exp(log(d_dat_withFX_bZ$fixEffectSum_bZ) + 2 * d_dat_withFX_bZ$value)
  
  d_dat_withFX <- rbind(d_dat_withFX_aZ, d_dat_withFX_bZ)
  
  # Get phenotypes with the mutation
  write.table(d_dat_withFX %>% ungroup() %>% 
                dplyr::select(rowID, aZ, bZ, KZ, KXZ), 
              "d_grid.csv", sep = ",", col.names = F, row.names = F)
  Aa <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)
  
  write.table(d_dat_withFX %>% ungroup() %>% 
                dplyr::select(rowID, aZ_AA, bZ_AA, KZ, KXZ), 
              "d_grid.csv", sep = ",", col.names = F, row.names = F)
  AA <- runLandscaper("d_grid.csv", "data_popfx.csv", 0.05, 2, 8, TRUE)
  
  # Ensure that the tables are aligned by id before we join them
  dat <- dat %>% arrange(rowID)
  Aa <- Aa %>% arrange(id)
  AA <- AA %>% arrange(id)
  
  # Get effect
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


# Rank the fixations in order of adaptive step (first step, second, etc.)
RankFixations <- function(dat, isNAR, dat_burnInFX = dat) {
  index <- as.integer(isNAR) + 1
  
  if (isNAR) {
    # Get fixed effects up to each rank
    dat <- dat %>%
      group_by(gen, seed) %>%
      mutate(fixEffectSum_aZ = exp(2 * sum(dat_burnInFX[dat_burnInFX$gen <= cur_group()$gen &
                                                    dat_burnInFX$mutType == 3 &
                                                    dat_burnInFX$seed == cur_group()$seed,]$value)),
             fixEffectSum_bZ = exp(2 * sum(dat_burnInFX[dat_burnInFX$gen <= cur_group()$gen &
                                                    dat_burnInFX$mutType == 4 &
                                                    dat_burnInFX$seed == cur_group()$seed,]$value)))
    
    d_fix_ranked <- dat %>%
      group_by(seed, modelindex) %>%
      arrange(gen, .by_group = T) %>%
      mutate(rank = row_number()) %>%
      dplyr::select(c(gen, rank, seed, modelindex, mutID, mutType, originGen,
                      value, aZ, bZ, phenomean, w, fixEffectSum_aZ, fixEffectSum_bZ,
                      avFX, avFit, avFit_AA, avFX_AA, AA_pheno, Aa_pheno, aa_pheno,
                      wAA, wAa, waa, s))
    
    step0_pheno <- d_adapted %>% 
      filter(modelindex == index, gen == 49500, interaction(seed, modelindex) %in%
               interaction(d_fix_ranked$seed, d_fix_ranked$modelindex)) %>%
      group_by(gen, seed) %>%
      mutate(rank = 0, value = NA, mutID = NA, avFit = NA,
             fixEffectSum_aZ = exp(2 * sum(dat_burnInFX[dat_burnInFX$gen <= cur_group()$gen &
                                                    dat_burnInFX$mutType == 3 &
                                                    dat_burnInFX$seed == cur_group()$seed,]$value)),
             fixEffectSum_bZ = exp(2 * sum(dat_burnInFX[dat_burnInFX$gen <= cur_group()$gen &
                                                    dat_burnInFX$mutType == 4 &
                                                    dat_burnInFX$seed == cur_group()$seed,]$value))) %>%
      dplyr::select(gen, rank, seed, modelindex, mutID, value, 
                    aZ, bZ, phenomean, w, fixEffectSum_aZ, fixEffectSum_bZ, avFit)
  } else {
    d_fix_ranked <- dat %>%
      group_by(seed, modelindex) %>%
      arrange(gen, .by_group = T) %>%
      mutate(rank = row_number()) %>%
      dplyr::select(c(gen, rank, seed, modelindex, mutID, mutType, originGen,
                      fixEffectSum, value, value_AA, aZ, bZ, phenomean, w, 
                      avFit, avFit_AA, AA_pheno, Aa_pheno, aa_pheno, 
                      wAA, wAa, waa, s))
    
    step0_pheno <- d_adapted %>% 
      filter(modelindex == index, gen == 49500, interaction(seed, modelindex) %in%
               interaction(d_fix_ranked$seed, d_fix_ranked$modelindex)) %>%
      group_by(gen, seed) %>%
      mutate(rank = 0, value = NA, mutID = NA, avFit = NA,
             fixEffectSum = 2 * sum(dat_burnInFX[dat_burnInFX$gen <= cur_group()$gen &
                                                     dat_burnInFX$seed == cur_group()$seed,]$value)) %>%
      dplyr::select(gen, rank, seed, modelindex, mutID, value,
                    aZ, bZ, phenomean, w, fixEffectSum, avFit)
  }
  
  return(rbind(d_fix_ranked, step0_pheno))
}
