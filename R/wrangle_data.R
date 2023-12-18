# Transform/summarise data - 
# calculates fitness effects, ranks fixations in the adaptive walk,
# tidies data for plotting



# load trait evolution data
d_qg <- read.table(paste0(dataPath, "slim_qg.csv"), header = F, 
                   sep = ",", colClasses = c("integer", "factor", "factor", 
                                             rep("numeric", times = 12)), 
                   col.names = c("gen", "seed", "modelindex", "meanH", "VA",
                                 "phenomean", "phenovar", "dist", "w", "deltaPheno",
                                 "deltaw", "aZ", "bZ", "KZ", "KXZ"), 
                   fill = T)

# Filter for models which adapted (within 10% of optimum by the end of the simulation)
d_qg %>%
  distinct() %>%
  group_by(seed, modelindex) %>%
  mutate(isAdapted = any(gen >= 59800 & between(phenomean, 1.9, 2.1))) %>%
  ungroup() -> d_qg
d_adapted <- d_qg %>% filter(isAdapted)

# Check how many there are
d_qg %>% group_by(modelindex) %>%
  filter(gen == 49500) %>%
  summarise(n = n(),
            pAdapted = mean(isAdapted),
            CIAdapted = CI(isAdapted))

# 719 additive models adapted total
# 585 network models adapted total

# load in observed heterozygosity data
d_het <- read.table(paste0(dataPath, "slim_locusHo.csv"), 
                     header = F, sep = ",", 
                     colClasses = c("integer", "factor", "factor", "numeric", "numeric"), 
                                    col.names = c("gen", "seed", "modelindex", "Ho_l1", "Ho_l2"), 
                                    fill = T)

# Pivot to calculate mean heterozygosity across the two loci
d_het %>% pivot_longer(cols = c(Ho_l1, Ho_l2), names_to = "locus", values_to = "Ho") %>%
  group_by(modelindex) %>%
  summarise(meanHo = mean(Ho),
            CIHo = CI(Ho))

d_Ho <- d_het %>% 
  pivot_longer(cols = c(Ho_l1, Ho_l2), names_to = "locus", values_to = "Ho") %>%
  mutate(model = ifelse(modelindex == 1, "Additive", "NAR"))

# Summarise, mean H_O over time
d_Ho %>%
  group_by(gen, model) %>%
  summarise(meanHo = mean(Ho),
            CIHo = CI(Ho)) -> d_Ho_sum


# Load mutation data
d_muts <- read.table(paste0(dataPath, "slim_muts.csv"), header = F, 
                   sep = ",", colClasses = c("integer", "factor", "factor", 
                                             "factor", rep("integer", times = 4),
                                             rep("numeric", times = 3),
                                             rep("integer", times = 2)), 
                   col.names = c("gen", "seed", "modelindex", "mutType", "mutID",
                                 "pos", "constraint", "originGen", "value", "chi",
                                 "Freq", "Count", "fixGen"), 
                   fill = T)

# Filter to include only adapted populations
d_muts_adapted <- d_muts %>% filter(interaction(seed, modelindex) %in% 
                                    interaction(d_adapted$seed, d_adapted$modelindex))

# Combine with d_qg
d_com_adapted <- inner_join(d_adapted, d_muts_adapted, by = c("gen", "seed", "modelindex"))


# Keep only fixations - adaptive walk
d_fix <- d_muts %>%
  filter(Freq == 1) %>%
  group_by(seed, modelindex, mutType) %>%
  distinct(mutID, .keep_all = T)

# Filter for only adapted populations
d_fix_adapted <- d_fix %>% filter(interaction(seed, modelindex) %in% 
                                    interaction(d_adapted$seed, d_adapted$modelindex))

# How long (approximately) did it take to fix?
d_fix_adapted$fixTime <- d_fix_adapted$gen - d_fix_adapted$originGen

# Split into additive/NAR models for separate processing
d_fix_add <- d_fix_adapted %>% filter(modelindex == 1)

# Join with trait data
d_qg_matched_fix <- d_adapted %>% 
  filter(interaction(gen, seed, modelindex) %in% 
           interaction(d_fix_add$gen, d_fix_add$seed, d_fix_add$modelindex)) %>%
  dplyr::select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w) %>% distinct()

d_fix_add <- inner_join(d_fix_add, d_qg_matched_fix, by = c("gen", "seed", "modelindex"))

# Calculate the fitness effects in additive populations (>= 50000 because that's when burn-in ended)
d_fix_add <- CalcAddEffects(d_fix_add) %>% filter(gen >= 50000)

# Filter for NAR populations
d_fix_nar <- d_fix_adapted %>% ungroup() %>% mutate(r = rownames(.)) %>% 
  filter(modelindex == 2)

# First get matched molecular component data for generations where we have fixations
d_qg_matched_fix <- d_adapted %>% 
  filter(interaction(gen, seed, modelindex) %in% 
           interaction(d_fix_nar$gen, d_fix_nar$seed, d_fix_nar$modelindex)) %>%
  dplyr::select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w) %>% distinct()

d_fix_nar <- inner_join(d_fix_nar, d_qg_matched_fix, by = c("gen", "seed", "modelindex"))

# Calculate fitness effects for NAR populations
d_fix_nar <- CalcNARPhenotypeEffects(d_fix_nar) %>% filter(gen >= 50000)

# Get adaptive step order and attach step 0 (phenotype from before the first step in the walk)
d_fix_ranked <- RankFixations(d_fix_nar, T, d_fix %>% filter(modelindex == 2))
d_fix_ranked_add <- RankFixations(d_fix_add, F, d_fix %>% filter(modelindex == 1))

# How many populations took n steps to reach the optimum?
rbind(d_fix_ranked %>% mutate(model = "NAR"), 
      d_fix_ranked_add %>% mutate(model = "Additive")) %>% 
  group_by(model, rank) %>%
  filter(rank != 0) %>%
  summarise(n = n())

# Since there aren't many populations with steps >3, we'll organise the groups
# into 1, 2, >=3
d_fix_ranked %>% 
  mutate(rankFactor = ifelse(rank > 2, "\\geq 3", as.character(rank))) -> d_fix_ranked
d_fix_ranked$rankFactor <- factor(d_fix_ranked$rankFactor, 
                                  levels = c("0", "1", "2", "\\geq 3"))

d_fix_ranked_add %>% 
  mutate(rankFactor = ifelse(rank > 2, "\\geq 3", as.character(rank))) -> d_fix_ranked_add
d_fix_ranked_add$rankFactor <- factor(d_fix_ranked_add$rankFactor, 
                                  levels = c("0", "1", "2", "\\geq 3"))


# get segregating variants
d_muts_adapted %>% 
  filter(modelindex == 2, Freq < 1, interaction(gen, seed, modelindex) %in%
           interaction(d_fix_ranked$gen, d_fix_ranked$seed, 
                       d_fix_ranked$modelindex),
         Freq > 0) -> d_seg_ranked

d_seg_ranked %>%
  group_by(gen, seed) %>%
  distinct() %>%
  mutate(rank = tail(d_fix_ranked[d_fix_ranked$gen == cur_group()$gen &
                               d_fix_ranked$seed == cur_group()$seed,]$rank, 1)) -> d_seg_ranked

d_muts_adapted %>% 
  filter(modelindex == 1, Freq < 1, interaction(gen, seed, modelindex) %in%
           interaction(d_fix_ranked_add$gen, d_fix_ranked_add$seed, 
                       d_fix_ranked_add$modelindex),
         Freq > 0) -> d_seg_ranked_add

d_seg_ranked_add %>%
  distinct() %>%
  group_by(gen, seed) %>%
  mutate(rank = tail(d_fix_ranked_add[d_fix_ranked_add$gen == cur_group()$gen & 
                               d_fix_ranked_add$seed == cur_group()$seed,]$rank, 1)) -> d_seg_ranked_add

# Calculate fitness for the segregating variants
d_qg_matched_seg <- d_qg %>% 
  filter(interaction(gen, seed, modelindex) %in% 
           interaction(d_seg_ranked$gen, d_seg_ranked$seed, d_seg_ranked$modelindex)) %>%
  dplyr::select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w) %>% distinct()

d_seg_ranked <- inner_join(d_seg_ranked, d_qg_matched_seg, 
                           by = c("gen", "seed", "modelindex"))
d_seg_ranked <- CalcNARPhenotypeEffects(d_seg_ranked, F, d_fix_adapted)


d_qg_matched_seg <- d_qg %>% 
  filter(interaction(gen, seed, modelindex) %in% 
           interaction(d_seg_ranked_add$gen, d_seg_ranked_add$seed, 
                       d_seg_ranked_add$modelindex)) %>%
  dplyr::select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w) %>% distinct()

d_seg_ranked_add <- inner_join(d_seg_ranked_add, d_qg_matched_seg, 
                               by = c("gen", "seed", "modelindex"))
d_seg_ranked_add <- CalcAddEffects(d_seg_ranked_add, F, d_fix_adapted)

# contribution to trait seg vs fixed
GetSegFixContributions <- function(seg, fix, isNAR) {
  # weight effects by frequency
  seg %>%
    filter(rank > 0) %>%
    group_by(gen, seed, modelindex, rank) %>%
    summarise(segEffectSum = sum(abs(seg[seg$gen == cur_group()$gen & 
                                   seg$seed == cur_group()$seed,][[ifelse({{ isNAR }}, "avFX_AA", "value_AA")]])*
                                   seg[seg$gen == cur_group()$gen & 
                                         seg$seed == cur_group()$seed,]$Freq),
              weightSumSeg = sum(Freq),
              segWeightedAverage = segEffectSum/weightSumSeg) -> d_segFX

  # Select all fixations before or equal to this one (cumulative sum)
  fix %>%
    filter(rank > 0) %>%
    group_by(gen, seed, modelindex, rank) %>%
    summarise(fixEffectSum = sum(abs(fix[fix$rank <= cur_group()$rank & fix$rank > 0 &
                                           fix$seed == cur_group()$seed,][[ifelse({{ isNAR }}, "avFX_AA", "value_AA")]])),
              fixWeightedAverage = fixEffectSum/rank) -> d_fixFX
  
  return(inner_join(d_fixFX, d_segFX, 
                    by = c("gen", "seed", "modelindex", "rank")) %>%
    group_by(gen, seed, rank) %>%
      # weighted average: rank is the sum of fixation frequencies (1), weightSumFix
    mutate(percFix = (fixEffectSum)/(fixEffectSum + segEffectSum)) %>%
    ungroup(gen, seed) %>%
    summarise(meanPercFix = mean(percFix, na.rm = T),
              CIPercFix = CI(percFix, na.rm = T)))
}

d_segFixRat_add_sum <- GetSegFixContributions(d_seg_ranked_add, d_fix_ranked_add, F)
d_segFixRat_sum <- GetSegFixContributions(d_seg_ranked, d_fix_ranked, T)

# seg effects weighted by frequency
d_seg_ranked$weighteds <- d_seg_ranked$s * d_seg_ranked$Freq
d_seg_ranked_add$weighteds <- d_seg_ranked_add$s * d_seg_ranked_add$Freq


# Some fixations are deleterious in early steps - why?
# Could be they are adaptive in burn-in environment, reached a high freq
# there and then drifted the rest of the way
# so need to calculate fitness effect of the mutation relative to that environment:
# where they originated, and when the mutation first reached >50% frequency
d_fix_ranked_combined <- rbind(d_fix_ranked, d_fix_ranked_add)
d_fix_ranked_combined$model <- if_else(d_fix_ranked_combined$modelindex == 1, "Additive", "NAR")

d_fix_ranked_combined %>% filter(s < 0) -> d_fix_del

# need to find simulations with those muts
d_muts_adapted %>% filter(interaction(seed, modelindex) %in%
                            interaction(d_fix_del$seed, 
                                        d_fix_del$modelindex)) -> d_muts_del

d_muts_del %>% filter(mutID %in% d_fix_del$mutID) %>% distinct() -> d_muts_del

# Calculate fitness effects when they arose, 
# and when they first reached 50% or greater freq - not fixed
d_qg_matched_seg <- d_qg %>% 
  filter(interaction(gen, seed, modelindex) %in% 
           interaction(d_muts_del$gen, d_muts_del$seed, d_muts_del$modelindex)) %>%
  dplyr::select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w) %>% distinct()

d_seg_del <- inner_join(d_muts_del, d_qg_matched_seg, 
                        by = c("gen", "seed", "modelindex"))

# Calculate phenotypes
d_seg_del_add <- CalcAddEffects(d_seg_del %>% filter(modelindex == 1), isFixed = F,
                                dat_fixed = d_fix_adapted %>% filter(modelindex == 1))
d_seg_del <- CalcNARPhenotypeEffects(d_seg_del %>% filter(modelindex == 2), isFixed = F, 
                        dat_fixed = d_fix_adapted %>% filter(modelindex == 2))

# Fitness calculations aren't correct for those before optimum shift: recalculate
d_seg_del[d_seg_del$gen < 50000,]$wAA <- calcAddFitness(d_seg_del[d_seg_del$gen < 50000,]$AA_pheno, 1, 0.05)
d_seg_del[d_seg_del$gen < 50000,]$wAa <- calcAddFitness(d_seg_del[d_seg_del$gen < 50000,]$Aa_pheno, 1, 0.05)
d_seg_del[d_seg_del$gen < 50000,]$waa <- calcAddFitness(d_seg_del[d_seg_del$gen < 50000,]$aa_pheno, 1, 0.05)
d_seg_del$avFit <- d_seg_del$wAa - d_seg_del$waa
d_seg_del$avFit_AA <- d_seg_del$wAA - d_seg_del$waa
d_seg_del$s <- d_seg_del$wAA - d_seg_del$waa

d_seg_del_add[d_seg_del_add$gen < 50000,]$wAA <- calcAddFitness(d_seg_del_add[d_seg_del_add$gen < 50000,]$AA_pheno, 1, 0.05)
d_seg_del_add[d_seg_del_add$gen < 50000,]$wAa <- calcAddFitness(d_seg_del_add[d_seg_del_add$gen < 50000,]$Aa_pheno, 1, 0.05)
d_seg_del_add[d_seg_del_add$gen < 50000,]$waa <- calcAddFitness(d_seg_del_add[d_seg_del_add$gen < 50000,]$aa_pheno, 1, 0.05)
d_seg_del_add$avFit <- d_seg_del_add$wAa - d_seg_del_add$waa
d_seg_del_add$avFit_AA <- d_seg_del_add$wAA - d_seg_del_add$waa
d_seg_del_add$s <- d_seg_del_add$wAA - d_seg_del_add$waa

d_seg_del <- rbind(d_seg_del, d_seg_del_add)
d_seg_del$model <- ifelse(d_seg_del$modelindex == 2, "NAR", "Additive")
rm(d_seg_del_add)

# Get the difference in selection coefficient after optimum shift
# as well as mean frequency just prior to optimum shift
d_seg_del %>%
  filter(gen == 49500 | gen == 50000) %>%
  arrange(gen) %>%
  group_by(seed, model, mutID) %>%
  filter(n() > 1) %>%
  summarise(diff_s = s[2] - s[1],
            Freq = Freq[1]) -> d_del_diffs

# Calculate the contributions of adaptive steps vs segregating mutations on phenotypic variation 
d_fix_ranked %>%
  mutate(value_aZ = if_else(mutType == 3, value, 0),
         value_bZ = if_else(mutType == 4, value, 0)) %>%
  group_by(seed) %>%
  filter(n() > 2) %>% # exclude groups with less than 2 steps
  mutate(molCompDiff = sum(abs(2 * value_aZ), na.rm = T) - sum(abs(2 * value_bZ), na.rm = T),
         molCompCorr = cor(abs(2 * value_aZ), abs(2 * value_bZ))) %>%
  ungroup() %>%
  dplyr::select(seed, molCompDiff, molCompCorr) %>%
  distinct(seed, .keep_all = T) -> d_molCompDiff

# Multiplicative model:
## We need to calculate the deviation between a multiplicative model of phenotypes
## (exp^(sum(alpha, beta)) and the ODE phenotype
## We can take the results from the NAR model, feed in the alpha and beta values
## and use them as a multiplicative model -- same mutations, different phenotype
## function. The difference between phenotypes is the NAR effect on the phenotype
## or the phenotypic variance described by the NAR as opposed to the multiplicative
## scaling

# First get matched molecular component data for generations where we have fixations
d_qg_matched_fix <- d_adapted %>% 
  filter(interaction(gen, seed, modelindex) %in% 
           interaction(d_fix_nar$gen, d_fix_nar$seed, d_fix_nar$modelindex)) %>%
  dplyr::select(gen, seed, modelindex, aZ, bZ, KZ, KXZ, phenomean, w) %>% distinct()

d_fix_mult <- inner_join(d_fix_nar, d_qg_matched_fix, by = c("gen", "seed", "modelindex"))

# Calculate fitness effects for NAR populations
d_fix_mult <- CalcMultEffects(d_fix_mult) %>% filter(gen >= 50000)

# Get adaptive step order and attach step 0 (phenotype from before the first step in the walk)
d_fix_ranked_mult <- RankFixations(d_fix_mult, F, d_fix %>% filter(modelindex == 2))

# How many populations took n steps to reach the optimum?
rbind(d_fix_ranked %>% mutate(model = "NAR"), 
      d_fix_ranked_add %>% mutate(model = "Additive"),
      d_fix_ranked_mult %>% mutate(model = "Multiplicative")) %>% 
  group_by(model, rank) %>%
  filter(rank != 0) %>%
  summarise(n = n())

# Since there aren't many populations with steps >3, we'll organise the groups
# into 1, 2, >=3
d_fix_ranked %>% 
  mutate(rankFactor = ifelse(rank > 2, "\\geq 3", as.character(rank))) -> d_fix_ranked
d_fix_ranked$rankFactor <- factor(d_fix_ranked$rankFactor, 
                                  levels = c("0", "1", "2", "\\geq 3"))

