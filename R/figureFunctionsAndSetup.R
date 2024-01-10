require(igraph)
require(ggraph)
require(extRemes)

# Plot fitness landscape
plotaZbZLandscape <- function(minVal, maxVal) {
  GRID_RES <- 400
  d_grid <- expand.grid(seq(from = minVal, 
                            to = maxVal, length.out = GRID_RES), 
                        seq(from = minVal, 
                            to = maxVal, length.out = GRID_RES), 1, 1)
  write.table(d_grid, "d_pairinput.csv", sep = ",", col.names = F, row.names = F)
  
  d_landscape <- runLandscaper("d_pairinput.csv", "d_pairwiselandscape.csv", 0.05, 2, 8)
  cc <- paletteer_c("viridis::viridis", 3, direction = 1)
  
  minFit <- 0.8 #min(d_landscape$fitness)
  maxFit <- max(d_landscape$fitness)
  # Rescale fitness values so we can change gradient breaks
  wValues <- c(0,
               (0.90-minFit)/(maxFit - minFit),
               (0.99-minFit)/(maxFit - minFit),
               1)
  
  suppressWarnings(
    ggplot(d_landscape %>% mutate(aZbZ = aZ/bZ), 
           aes(x = aZ, y = bZ, colour = fitness)) +
      geom_point() +
      #geom_abline(slope = 1/1.27) +
      # geom_point(data = d_landscape %>% mutate(aZbZ = aZ/bZ) 
      #            %>% filter(aZbZ > 1.25, aZbZ < 1.35), size = 0.1, shape = 4) +
      scale_colour_gradientn(colors = c(cc[1], cc), 
                             limits = c(minFit, 1),
                             values = wValues, na.value = cc[1]) +
      labs(x = TeX("$\\alpha_Z$"), y = TeX("$\\beta_Z$"), 
           colour = "Fitness (w)") +
      theme_bw() + 
      theme(legend.position = "bottom", text = element_text(size = 14)) +
      guides(
        colour = guide_colourbar(barwidth = 10, title.vjust = 0.87)) # i love magic numbers
  )
}


# ODE system for feedback autoregulation:
Hilln <- 8 # constant for Hill function

ODEs_FBA <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- (t > Xstart && t <= Xstop)
    dZ <- bZ * (t > Xstart && t <= Xstop) * (X^Hilln)/(KXZ^Hilln + X^Hilln) * (KZ^Hilln)/(KZ^Hilln + Z^Hilln) - aZ*Z
    dZnoFB <- aZ * (t > Xstart && t <= Xstop) - aZ*ZnoFB
    return(list(c(dZ, dZnoFB)))
  })
}

# Plot feedback dynamics
plotDynamics_FBA <- function(Xstart = 1,
                             Xstop = 6,
                             tmax = 10,
                             dt = 0.1,
                             pars = list(aZ = 1,
                                         bZ = 1,
                                         KXZ = 1,
                                         KZ = 1)) {
  params <- c(Xstart = Xstart, Xstop = Xstop, aZ = pars$aZ, bZ = pars$bZ, KXZ = pars$KXZ, KZ = pars$KZ, Hilln = 8)
  iniState <- c(Z=0, ZnoFB = 0)
  times <- seq(0,tmax,by=dt)
  solution <- ode(iniState, times, ODEs_FBA, params) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(X = ifelse(time > params["Xstart"] & time <= params["Xstop"], 1, 0),
           aZ = as.factor(pars$aZ), bZ = as.factor(pars$bZ), 
           KXZ = as.factor(pars$KXZ), KZ = as.factor(pars$KZ)) %>%
    dplyr::select(time, aZ, bZ, KXZ, KZ, X, Z, ZnoFB)
  return(solution)
}

genRatioLandscapeData <- function(minRatio, maxRatio) {
  GRID_RES <- 1000
  rational <- MASS:::.rat(seq(minRatio, maxRatio, length.out = GRID_RES),
                          max.denominator = 2000)$rat
  d_grid <- data.frame(aZ = rational[,1], 
                       bZ = rational[,2],
                       KZ = 1, 
                       KXZ = 1) %>% distinct()
  
  result <- data.frame(pheno = numeric(nrow(d_grid)), 
                       aZ = numeric(nrow(d_grid)), 
                       bZ = numeric(nrow(d_grid)), 
                       KZ = numeric(nrow(d_grid)),  
                       KXZ = numeric(nrow(d_grid)))
  
  for (i in seq_len(nrow(d_grid))) {
    solution <- plotDynamics_FBA(pars = as.list(d_grid[i,]))
    pheno <- AUC(solution$time, solution$Z, absolutearea = T)
    result[i,] <- c(pheno, d_grid[i,])
  }
  result$fitness <- calcAddFitness(result$pheno, 2, 0.05)
  return(result)
}

plotRatioLandscape <- function(minRatio, maxRatio) {
  d_landscape <- genRatioLandscapeData(minRatio, maxRatio)
  
  cc <- paletteer_c("viridis::viridis", 3, direction = 1)
  
  minFit <- 0.8 #min(d_landscape$fitness)
  maxFit <- max(d_landscape$fitness)
  # Rescale fitness values so we can change gradient breaks
  wValues <- c(0,
               (0.90-minFit)/(maxFit - minFit),
               (0.99-minFit)/(maxFit - minFit),
               1)
  
  suppressWarnings(
    ggplot(d_landscape  %>% mutate(aZbZ = aZ/bZ), 
           aes(x = aZbZ, y = pheno, colour = fitness)) +
      geom_point() +
      scale_colour_gradientn(colors = c(cc[1], cc), 
                             limits = c(ifelse(minFit < 0.8, 0.8, minFit), 1),
                             values = wValues, na.value = cc[1]) +
      labs(x = TeX("$\\alpha_Z / \\beta_Z$"), y = "Phenotype", 
           colour = "Fitness (w)") +
      theme_bw() + 
      theme(legend.position = "bottom", text = element_text(size = 14)) +
      guides(
        colour = guide_colourbar(barwidth = 20, title.vjust = 0.87)) # i love magic numbers
  )
}

plotbZaZvsaZLandscape <- function(minValaZ, maxValaZ, minRatio, maxRatio) {
  GRID_RES <- 400
  aZVals <- seq(minValaZ, maxValaZ, length.out = GRID_RES)
  ratios <- seq(minRatio, maxRatio, length.out = GRID_RES)
  combos <- expand.grid(aZVals, ratios)
  colnames(combos) <- c("aZVals", "ratios")
  
  # calculate beta values from aZ values and ratios
  combos$bZVals <- combos$aZVals * combos$ratios
  
  d_grid <- data.frame(aZ = combos$aZVals,
                       bZ = combos$bZVals,
                       KZ = 1,
                       KXZ = 1) %>% distinct()  
  
  write.table(d_grid, "d_pairinput.csv", sep = ",", col.names = F, row.names = F)

  d_landscape <- runLandscaper("d_pairinput.csv", "d_pairwiselandscape.csv", 0.05, 2, 8)
  cc <- paletteer_c("viridis::viridis", 3, direction = 1)
  
  minFit <- 0.8 #min(d_landscape$fitness)
  maxFit <- max(d_landscape$fitness)
  # Rescale fitness values so we can change gradient breaks
  wValues <- c(0,
               (0.90-minFit)/(maxFit - minFit),
               (0.99-minFit)/(maxFit - minFit),
               1)
  
  suppressWarnings(
    ggplot(d_landscape %>% mutate(bZaZ = bZ/aZ), 
           aes(x = aZ, y = bZaZ, colour = fitness)) +
      geom_point() +
      scale_colour_gradientn(colors = c(cc[1], cc), 
                           limits = c(minFit, 1),
                           values = wValues, na.value = cc[1]) +
      labs(x = TeX("$\\alpha_Z$"), y = TeX("$\\beta_Z/\\alpha_Z$"), 
           colour = "Fitness (w)") +
      theme_bw() + 
      theme(legend.position = "bottom", text = element_text(size = 14)) +
      guides(
        colour = guide_colourbar(barwidth = 10, title.vjust = 0.87)) # i love magic numbers
    )
}

plotbZaZvsbZLandscape <- function(minValbZ, maxValbZ, minRatio, maxRatio) {
  GRID_RES <- 400
  bZVals <- seq(minValbZ, maxValbZ, length.out = GRID_RES)
  ratios <- seq(minRatio, maxRatio, length.out = GRID_RES)
  combos <- expand.grid(bZVals, ratios)
  colnames(combos) <- c("bZVals", "ratios")
  
  # calculate beta values from aZ values and ratios
  combos$aZVals <- combos$bZVals * (1/combos$ratios)
  combos <- drop_na(combos)
  
  d_grid <- data.frame(aZ = combos$aZVals,
                       bZ = combos$bZVals,
                       KZ = 1,
                       KXZ = 1) %>% distinct()  
  
  write.table(d_grid, "d_pairinput.csv", sep = ",", col.names = F, row.names = F)
  
  d_landscape <- runLandscaper("d_pairinput.csv", "d_pairwiselandscape.csv", 0.05, 2, 8)
  cc <- paletteer_c("viridis::viridis", 3, direction = 1)
  
  minFit <- 0.8 #min(d_landscape$fitness)
  maxFit <- max(d_landscape$fitness)
  # Rescale fitness values so we can change gradient breaks
  wValues <- c(0,
               (0.90-minFit)/(maxFit - minFit),
               (0.99-minFit)/(maxFit - minFit),
               1)
  
  suppressWarnings(
    ggplot(d_landscape %>% mutate(bZaZ = bZ/aZ), 
           aes(x = bZ, y = bZaZ, colour = fitness)) +
      geom_point() +
      scale_colour_gradientn(colors = c(cc[1], cc), 
                             limits = c(minFit, 1),
                             values = wValues, na.value = cc[1]) +
      labs(x = TeX("$\\beta_Z$"), y = TeX("$\\beta_Z/\\alpha_Z$"), 
           colour = "Fitness (w)") +
      theme_bw() + 
      theme(legend.position = "bottom", text = element_text(size = 14)) +
      guides(
        colour = guide_colourbar(barwidth = 10, title.vjust = 0.87)) # i love magic numbers
  )
}


genRatioLandscapeData <- function(minRatio, maxRatio) {
  GRID_RES <- 1000
  rational <- MASS:::.rat(seq(minRatio, maxRatio, length.out = GRID_RES),
                          max.denominator = 2000)$rat
  d_grid <- data.frame(aZ = rational[,1], 
                       bZ = rational[,2],
                       KZ = 1, 
                       KXZ = 1) %>% distinct()
  # write.table(d_grid, "d_pairinput.csv", sep = ",", col.names = F, row.names = F)
  # 
  # return(runLandscaper("d_pairinput.csv", "d_pairwiselandscape.csv", 0.05, 2, 8) %>%
  #          filter(pheno > 0))
  
  result <- data.frame(pheno = numeric(nrow(d_grid)), 
                       aZ = numeric(nrow(d_grid)), 
                       bZ = numeric(nrow(d_grid)), 
                       KZ = numeric(nrow(d_grid)),  
                       KXZ = numeric(nrow(d_grid)))
  
  for (i in seq_len(nrow(d_grid))) {
    solution <- plotDynamics_FBA(pars = as.list(d_grid[i,]))
    pheno <- AUC(solution$time, solution$Z, absolutearea = T)
    result[i,] <- c(pheno, d_grid[i,])
  }
  result$fitness <- calcAddFitness(result$pheno, 2, 0.05)
  return(result)
}

plotRatioLandscape <- function(minRatio, maxRatio) {
  d_landscape <- genRatioLandscapeData(minRatio, maxRatio)
  
  cc <- paletteer_c("viridis::viridis", 3, direction = 1)
  
  minFit <- 0.8 #min(d_landscape$fitness)
  maxFit <- max(d_landscape$fitness)
  # Rescale fitness values so we can change gradient breaks
  wValues <- c(0,
               (0.90-minFit)/(maxFit - minFit),
               (0.99-minFit)/(maxFit - minFit),
               1)
  
  suppressWarnings(
    ggplot(d_landscape %>% mutate(aZbZ = aZ/bZ), 
           aes(x = aZbZ, y = pheno, colour = fitness)) +
      geom_point() +
      scale_colour_gradientn(colors = c(cc[1], cc), 
                             limits = c(ifelse(minFit < 0.8, 0.8, minFit), 1),
                             values = wValues, na.value = cc[1]) +
      labs(x = TeX("$\\alpha_Z / \\beta_Z$"), y = "Phenotype", 
           colour = "Fitness (w)") +
      theme_bw() + 
      theme(legend.position = "bottom", text = element_text(size = 14)) +
      guides(
        colour = guide_colourbar(barwidth = 20, title.vjust = 0.87)) # i love magic numbers
  )
}
