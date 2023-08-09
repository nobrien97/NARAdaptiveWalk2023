# Fig 1: NAR network 
colZ   <- brewer.pal(6, "Paired")[4]
colZbg <- brewer.pal(6, "Paired")[3]

colX <- "#0075BD"
colXbg   <- "#CCECFF"

colEg <- brewer.pal(6, "Paired")[6]


genesNAR <- data.frame(name = c("X", "Z"))
activationNAR <- data.frame(from = c("X", "Z"),
                            to =   c("Z", "Z"))
gNAR <- graph_from_data_frame(activationNAR, directed = TRUE, vertices = genesNAR)
plotNARgraph <- ggraph(gNAR, layout = "manual", x = c(0, 0), y=c(0.5, -0.5)) +
  geom_node_point(size = 16, color = c(colX, colZ)) +
  geom_node_point(size = 15, color = c(colXbg, colZbg)) +
  geom_node_text(aes(label = name), size = 10, color = "#5E5E5E") +
  geom_edge_arc(aes(start_cap = circle(8, unit = "mm"),
                    end_cap = circle(8, unit = "mm")),
                arrow = arrow(type = "open", angle = 30, length = unit(3, 'mm')),
                strength = c(0)) +
  geom_edge_loop(aes(start_cap = circle(8, unit = "mm"),
                     end_cap = circle(8, unit = "mm"),
                     span = 90, direction = 0, strength = 0.7),
                 arrow = arrow(type = "open", angle = 90, length = unit(3, 'mm'))) +
  scale_x_continuous(limits = c(-0.7,0.7)) +
  scale_y_continuous(limits = c(-0.7,0.7)) +
  #coord_fixed() +
  theme_graph(plot_margin = unit(rep(1, times = 4), "mm"), background = NA)
plotNARgraph

dat <- plotDynamics_FBA()

plotZ <- ggplot(dat) +
  scale_color_manual(values = colZ) +
  annotate("rect", xmin = 1, xmax = 6, ymin = 0, ymax = 1.05,
           alpha = .2, fill = colX) +
  geom_line(aes(time, Z, color = bZ), linewidth = 1.5) +
  scale_y_continuous(limits = c(0,1.05)) +
  labs(x = "Time", y = "Z expression") +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position = "none")

fig1 <- plot_grid(plotNARgraph, plotZ, rel_widths = c(0.75, 1), labels = "AUTO",
                  label_size = 20)
fig1

ggsave("NARgraph.png", fig1, width = 8, height = 4, device = png, bg = "white")


# phenomean and adaptive walk
# A - phenomean ridgeline plot
d_adapted_walk <- d_adapted %>% filter(gen >= 49000)
breaks <- seq(min(d_adapted_walk$gen - 50000), max(d_adapted_walk$gen - 50000), by = 1000)

d_adapted_walk$gen_group <- breaks[findInterval(d_adapted_walk$gen - 50000, breaks, rightmost.closed = TRUE)]

ggplot(d_adapted_walk,
       aes(y = as.factor(gen_group), x = phenomean, fill = modelindex)) +
  geom_density_ridges(alpha = 0.4) +
  geom_vline(xintercept = 2, linetype = "dashed") +
  scale_fill_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "NAR")) +
  labs(y = "Generations post-optimum shift", x = "Phenotype mean", 
       fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_phenomean_dist
plt_phenomean_dist

# B: phenotype at each step
step_labs <- paste0("$", levels(d_fix_ranked_combined$rankFactor), "$")
ggplot(d_fix_ranked_combined,
       aes(y = rankFactor, x = phenomean, fill = model)) +
  geom_density_ridges(alpha = 0.4) +
  geom_vline(xintercept = 2, linetype = "dashed") +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  scale_y_discrete(labels = parse(text=TeX(step_labs))) +
  labs(y = "Adaptive step", x = "Phenotype mean", 
       fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none") -> plt_adaptivewalk_pheno_dist
plt_adaptivewalk_pheno_dist

leg <- get_legend(plt_phenomean_dist)

plot_grid(plt_phenomean_dist + theme(legend.position = "none"), 
          plt_adaptivewalk_pheno_dist, 
          ncol = 2,
          labels = "AUTO") -> plt_traitevo

plot_grid(plt_traitevo, leg, nrow = 2, rel_heights = c(1, 0.1)) -> plt_traitevo
ggsave("fig_traitevolution.png", plt_traitevo, device = png, bg = "white")


# distribution of fixations and step size over time
d_fix_ranked_combined$model <- if_else(d_fix_ranked_combined$modelindex == 1, "Additive", "NAR")

ggplot(d_fix_ranked_combined %>% filter(rank > 0), aes(x = s, fill = model)) +
  geom_density(alpha = 0.4) + 
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(y = "Density", x = "Fitness effect (s)", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none") -> plt_distfixed
plt_distfixed

ggplot(d_fix_ranked_combined %>% filter(rank > 0), 
       aes(x = s, y = rankFactor, fill = model)) +
  geom_density_ridges(alpha = 0.4) + 
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  scale_y_discrete(labels = parse(text=TeX(step_labs[2:4]))) +
  labs(y = "Adaptive step", x = "Fitness effect (s)", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none") -> plt_distfixed_time
plt_distfixed_time

plot_grid(plt_distfixed, 
          plt_distfixed_time, 
          ncol = 2,
          labels = "AUTO") -> plt_DFEfixations

plot_grid(plt_DFEfixations, leg, nrow = 2, rel_heights = c(1, 0.1)) -> plt_DFEfixations
ggsave("fig_fixations.png", plt_DFEfixations, device = png, bg = "white")


# space of possible mutations (mutation screen)
# A: all possible mutations at each step
ggplot(mutExp_combined, aes(y = rankFactor, x = s, fill = model)) +
  geom_density_ridges(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  scale_y_discrete(labels = parse(text=TeX(step_labs))) +
  labs(y = "Adaptive step", x = "Fitness effect (s)", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "none") -> plt_effectsizerandom_time
plt_effectsizerandom_time

# B: Proportion of mutations that are beneficial
ggplot(mutExp_sum_combined, aes(x = rankFactor, y = percBeneficial, colour = model)) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = percBeneficial - CIperc, 
                              ymax = percBeneficial + CIperc),
                width = 0.2) +
  geom_line(aes(group = model)) +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  scale_x_discrete(labels = parse(text=TeX(step_labs))) +
  labs(x = "Adaptive step", y = "Proportion of\nbeneficial mutations (s > 0)",
       colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "none") -> plt_propbeneficial
plt_propbeneficial

# C: Waiting time to a beneficial mutation
ggplot(mutExp_sum_combined %>% mutate(waitingTime = 1/(10000 * (9.1528*10^-6) * percBeneficial),
                                      CIWaitingTime_lower = 1/(10000 * (9.1528*10^-6) * (percBeneficial - CIperc)),
                                      CIWaitingTime_upper = 1/(10000 * (9.1528*10^-6) * (percBeneficial + CIperc))),
       aes(x = rankFactor, y = waitingTime, colour = model)) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = CIWaitingTime_lower, 
                              ymax = CIWaitingTime_upper),
                width = 0.2) +
  scale_x_discrete(labels = parse(text=TeX(step_labs))) +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  geom_line(aes(group=model)) +
  labs(x = "Adaptive step", y = "Expected waiting time\nto beneficial mutation", colour = "Model") +
  theme_bw() + 
  theme(text = element_text(size = 12), legend.position = "none") -> plt_waitingtime
plt_waitingtime

right_col <- plot_grid(plt_propbeneficial, 
                       plt_waitingtime,
                       ncol = 1,
                       labels = c("B", "C"))

plot_grid(plt_effectsizerandom_time, 
          right_col, 
          ncol = 2,
          labels = "AUTO", rel_heights = c(2, 1)) -> plt_mutationscreen
plot_grid(plt_mutationscreen, leg, nrow = 2, rel_heights = c(1, 0.1)) -> plt_mutationscreen
plt_mutationscreen

ggsave("fig_mutationscreen.png", plt_mutationscreen, width = 10, height = 6, device = png, bg = "white")


# fitness landscape and aZbZ ratio
# A - fitness landscape
plotaZbZLandscape(0, 3) -> plt_aZbZ_landscape
plt_aZbZ_landscape

# B - alpha vs beta/alpha landscape
plotbZaZvsaZLandscape(0, 3, 0, 3) -> plt_bZaZ_aZ_landscape

# C - GPW map of aZbZ
plotRatioLandscape(0.6, 3) -> plt_aZbZratio

plt_aZbZratio +
  stat_poly_line(colour = "#AAAAAA", linetype = "dashed") +
  stat_poly_eq(use_label(c("adj.R2", "p.value"), sep = "*\"; \"*"), 
               label.x = "right", colour = "#000000") +
  geom_hline(yintercept = 2, linetype = "dashed") +
  theme(text = element_text(size = 14)) -> plt_aZbZratio


# D - difference in evolution among alpha and beta
ggplot(d_molCompDiff,
       aes(x = molCompDiff)) +
  geom_density() +
  labs(x = TeX("Molecular component contribution ($\\phi_{\\alpha_Z \\beta_Z}$)"), y = "Density") +
  theme_bw() + 
  theme(text = element_text(size = 14)) -> plt_molCompDiff
plt_molCompDiff

leg <- get_legend(plt_aZbZ_landscape)

plot_grid(plt_aZbZ_landscape + theme(legend.position = "none"), 
          plt_bZaZ_aZ_landscape + theme(legend.position = "none"),
          plt_aZbZratio + theme(legend.position = "none"), 
          plt_molCompDiff,
          nrow = 2,
          labels = "AUTO") -> plt_fitnesslandscape

plot_grid(plt_fitnesslandscape, leg, nrow = 2, rel_heights = c(1, 0.1)) -> plt_fitnesslandscape
plt_fitnesslandscape
ggsave("fig_fitnessLandscape.png", plt_fitnesslandscape, width = 10, height = 8, device = png, bg = "white")


# supp fig - Phenotype fixed vs segregating effects
# A - correlation between phenotype of fixations only with mean phenotype
group_means <- d_fix_ranked_combined %>% filter(rank > 0) %>%
  group_by(model) %>%
  summarise(ratio = mean(AA_pheno/phenomean),
            CIRatio = CI(AA_pheno/phenomean))

# set seed for geom_jitter
set.seed(seed)
ggplot(d_fix_ranked_combined %>% filter(rank > 0), 
       aes(x = as.factor(model), y = AA_pheno/phenomean)) +
  geom_jitter(size = 0.5, shape = 1, alpha = 0.3) +
  geom_point(data = group_means, aes(y = ratio, colour = model), size = 2) +
  geom_errorbar(data = group_means, aes(y = ratio, ymin = ratio - CIRatio,
                                        ymax = ratio + CIRatio,
                                        colour = model), width = 0.1) +
  scale_colour_paletteer_d("ggsci::nrc_npg", guide = NULL) +
  scale_x_discrete(labels = c("Additive", "NAR")) +
  labs(x = "Model", y = "Fixed effect/mean phenotype ratio") +
  theme_bw() + 
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_segfixedcont
plt_segfixedcont

ggsave("sfig_segfixedcont.png", plt_segfixedcont, device = png)


# supp fig - Balancing selection example
ggplot(d_com_adapted %>% filter((modelindex == 1 & seed == 1448101263 & mutID == 4624809) | 
                                  (modelindex == 2 & seed == 2270695859 & mutID == 4607309), 
                                gen >= 49000) %>% 
         distinct() %>%
         mutate(gen = gen - 50000, 
                model = if_else(modelindex == 1, "Additive", "NAR")), 
       aes(x = gen, y = Freq, colour = model)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_x_continuous(labels = scales::comma) +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  labs(x = "Generations post-optimum shift", y = "Allele frequency", 
       colour = "Model") +
  theme(text = element_text(size = 16), legend.position = "bottom")
ggsave("sfig_balsel.png", device = png)

# Supp fig: fitness effect difference in deleterious fixations
ggplot(d_del_diffs, 
       aes(x = diff_s, fill = model)) +
  geom_density(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  theme_bw() +
  labs(x = TeX("Difference in fitness effect after optimum shift $(s_1 - s_0)$"),
       y = "Density", fill = "Model") +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_delfixed_s

ggplot(d_del_diffs, 
       aes(x = Freq, fill = model)) +
  geom_density(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  theme_bw() +
  labs(x = "Frequency at end of burn-in (p)",
       y = "Density", fill = "Model") +
  theme(text = element_text(size = 16), legend.position = "none") -> plt_delfixed_freq

leg <- get_legend(plt_delfixed_s)

plot_grid(plt_delfixed_s + theme(legend.position = "none"),
          plt_delfixed_freq + theme(legend.position = "none"),
          nrow = 2,
          labels = "AUTO") -> plt_delfixed
plt_delfixed <- plot_grid(plt_delfixed, 
                          leg, 
                          nrow = 2,
                          rel_heights = c(1, 0.1))
plt_delfixed
ggsave("sfig_delfixations.png", plt_delfixed, device = png, bg = "white")

# Supp fig: heterozygosity
# should be 0 most of the time if we're under SSWM
ggplot(d_Ho_sum %>% mutate(gen = gen - 50000), 
       aes(x = gen, y = meanHo, colour = model)) +
  geom_line() +
  geom_ribbon(aes(ymin = meanHo - CIHo, ymax = meanHo + CIHo,
                  colour = NA, fill = model),
              alpha = 0.2) +
  scale_x_continuous(labels = scales::comma) +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  scale_fill_paletteer_d("ggsci::nrc_npg", guide = NULL) +
  labs(x = "Generations post-optimum shift", y = "Mean population heterozygosity",
       colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_Ho
plt_Ho
ggsave("sfig_het.png", plt_Ho, device = png)


# supp fig - adaptive step timing for populations not yet at the optimum
ggplot(d_fix_ranked_combined %>% 
         filter(rank > 0, phenomean < 1.9 | phenomean > 2.1) %>%
         mutate(gen = gen - 50000),
       aes(y = rankFactor, x = gen, fill = model)) +
  geom_density_ridges(alpha = 0.4) +
  scale_y_discrete(labels = parse(text=TeX(step_labs[2:4]))) +
  scale_x_continuous(labels = scales::comma) +
  scale_fill_paletteer_d("ggsci::nrc_npg", labels = c("Additive", "NAR")) +
  labs(x = "Fixation generation (post-optimum shift)", y = "Adaptive step", 
       fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_adaptivestepgen_dist
plt_adaptivestepgen_dist
ggsave("sfig_fixgendist.png", plt_adaptivestepgen_dist, device = png)

# Supp fig: dist of beneficial fixations at each step - zoom in of Fig. 2B
ggplot(d_fix_ranked_combined %>% filter(rank > 0), 
       aes(x = s, y = rankFactor, fill = model)) +
  geom_density_ridges(alpha = 0.4) + 
  coord_cartesian(xlim = c(0, 0.1)) +
  scale_y_discrete(labels = parse(text=TeX(step_labs[2:4]))) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(y = "Adaptive step", x = "Fitness effect (s)", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")
ggsave("sfig_driftbarrier.png", device = png)
