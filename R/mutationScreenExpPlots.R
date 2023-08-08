setwd("/mnt/c/GitHub/SLiMTests/tests/fixedK/R")
source("mutationScreenExp.R")


# plot distributions of fitness effects among possible mutations
ggplot(mutExp, aes(x = s, y = as.factor(rank))) +
  geom_density_ridges() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(y = "Adaptive step", x = "Fitness effect (s)") +
  theme_bw()

ggplot(mutExp_add, aes(x = s, y = as.factor(rank))) +
  geom_density_ridges() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(y = "Adaptive step", x = "Fitness effect (s)") +
  theme_bw()

ggplot(mutExp_combined, aes(y = rankFactor, x = s, fill = model)) +
  geom_density_ridges(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(y = "Adaptive step", x = "Fitness effect (s)", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 12)) -> plt_effectsizerandom_time
plt_effectsizerandom_time

# Beneficial mutations only
ggplot(mutExp_combined %>% filter(s > 0), aes(y = rankFactor, x = s, fill = model)) +
  #facet_grid(.~model) +
  geom_density_ridges(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(y = "Adaptive step", x = "Fitness effect (s)", fill = "Model") +
  theme_bw() +
  coord_cartesian(xlim = c(0, 0.25)) +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_effectsizebenrandom_time
plt_effectsizebenrandom_time

# across all steps
ggplot(mutExp_combined, aes(x = s, fill = model)) +
  geom_density(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(x = "Fitness effect (s)", y = "Density", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_effectsizerandom
plt_effectsizerandom

# Beneficial mutations only
ggplot(mutExp_combined %>% filter(s > 0), aes(x = s, fill = model)) +
  geom_density(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(x = "Fitness effect (s)", y = "Density", fill = "Model") +
  theme_bw() +
  coord_cartesian(xlim = c(0, 0.25)) +
  theme(text = element_text(size = 16), legend.position = "bottom") -> plt_beneffectsizerandom
plt_beneffectsizerandom

leg <- get_legend(plt_beneffectsizerandom)
plt_ben_random <- plot_grid(plt_effectsizebenrandom_time + theme(legend.position = "none"),
                            plt_beneffectsizerandom + theme(legend.position = "none"),
                            labels = "AUTO",
                            nrow = 2, align = "v")
plt_ben_random <- plot_grid(plt_ben_random, leg, nrow = 2,
                            rel_heights = c(1, 0.1))
ggsave("fig_s_benrandom.png", plt_ben_random, width = 6, height = 9, 
       device = png, bg = "white")

ggplot(mutExp_sum_combined, aes(x = rankFactor, y = percBeneficial, colour = model)) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = percBeneficial - CIperc, 
                              ymax = percBeneficial + CIperc),
                width = 0.2) +
  geom_line(aes(group = model)) +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  labs(x = "Adaptive step", y = "Proportion of beneficial mutations (s > 0)",
       colour = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16))
ggsave("propbeneficial.png", device = png, bg = "white")

ggplot(mutExp_perc_combined, aes(y = rankFactor, x = percBeneficial, fill = model)) +
  geom_density_ridges(alpha = 0.4) +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  labs(y = "Adaptive step", x = "Beneficial mutation\nprobability (s > 0)", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 12)) -> plt_percBeneficial_time
plt_percBeneficial_time

# across all time points
ggplot(mutExp_perc_combined, aes(x = percBeneficial, fill = model)) +
  geom_density(alpha = 0.4) +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  labs(y = "Density", x = "Beneficial mutation\nprobability (s > 0)", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 12)) -> plt_percBeneficial
plt_percBeneficial


ggplot(mutExp_sum_combined, aes(x = rankFactor, y = meanEffect, colour = model)) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = meanEffect - CIEffect, 
                              ymax = meanEffect + CIEffect),
                width = 0.2) +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  geom_line(aes(group=model)) +
  labs(x = "Adaptive step", y = "Mean effect", colour = "Model") +
  theme_bw() + 
  theme(text = element_text(size = 12))
ggsave("meaneffectstep.png", device = png, bg = "white")

plot_grid(plt_percBeneficial, plt_effectsizerandom,
          nrow = 2, labels = "AUTO")
ggsave("effectstepdist.png", width = 8, height = 6, device = png, bg = "white")

ggplot(mutExp_sum_combined %>% mutate(waitingTime = 1/(10000 * (9.1528*10^-6) * percBeneficial)),
       aes(x = rankFactor, y = waitingTime, colour = model)) +
  geom_point() +
  scale_colour_paletteer_d("ggsci::nrc_npg") +
  geom_line(aes(group=model)) +
  labs(x = "Adaptive step", y = "Expected waiting time to beneficial mutation", colour = "Model") +
  theme_bw() + 
  theme(text = element_text(size = 12))


# average fitness effect among beneficial mutations
mutExp_combined %>%
  filter(s > 0) %>%
  group_by(rankFactor, model) %>%
  summarise(mean_s = mean(s),
            CI_s = CI(s))

# among deleterious muts
mutExp_combined %>%
  filter(s < 0) %>%
  group_by(rankFactor, model) %>%
  summarise(mean_s = mean(s),
            CI_s = CI(s))

# distribution of all beneficial mutations
ggplot(mutExp_combined %>% filter(s > 0), aes(x = s, fill = model)) +
  geom_density(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(x = "Fitness effect (s)", y = "Density", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_evt1
plt_evt1

# distribution of fixed effects
ggplot(d_fix_ranked_combined %>% filter(s > 0),
       aes(x = s, fill = model)) +
  geom_density(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(x = "Fitness effect (s)", y = "Density", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_evt2
plt_evt2

# changing over time
ggplot(mutExp_combined %>% filter(s > 0, rank > 0),
       aes(x = s, y = rankFactor, fill = model)) +
  geom_density_ridges(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(x = "Fitness effect (s)", y = "Adaptive step", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_evt1_time
plt_evt1_time

ggplot(d_fix_ranked_combined %>% filter(s > 0, rank > 0),
       aes(x = s, y = rankFactor, fill = model)) +
  geom_density_ridges(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(x = "Fitness effect (s)", y = "Adaptive step", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_evt2_time
plt_evt2_time

# what about deleterious variants?
ggplot(mutExp_combined %>% filter(s < 0),
       aes(x = s, fill = model)) +
  geom_density(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(x = "Fitness effect (s)", y = "Density", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_del_dist
plt_del_dist

ggplot(mutExp_combined %>% filter(s < 0, rank > 0),
       aes(x = s, y = as.factor(rank), fill = model)) +
  geom_density_ridges(alpha = 0.4) +
  scale_fill_paletteer_d("ggsci::nrc_npg") +
  labs(x = "Fitness effect (s)", y = "Adaptive step", fill = "Model") +
  theme_bw() +
  theme(text = element_text(size = 16)) -> plt_del_dist_time
plt_del_dist_time


