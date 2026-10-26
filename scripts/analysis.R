# Pacheco and Longo (2026)
# Coral bleaching and photochemical efficiency analyses
#
# Expected project structure:
#   data/
#     COLOR.csv
#     exp1_dhw.csv
#     bleaching_prevalence_by_day_R.csv
#     PAM.csv
#   figures/   (created automatically)
#
# Run this script from the root of the GitHub repository.

# 1. Packages -----------------------------------------------------------------

library(tidyverse)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(patchwork)
library(performance)

# 2. General settings ----------------------------------------------------------

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

site_levels <- c("P", "R", "F")
site_labels <- c(
  P = "Tide pool",
  R = "Shallow reef",
  F = "Deep reef"
)
site_colors <- c(
  P = "#e88a25",
  R = "#1ea49d",
  F = "#482344"
)
site_colors_labeled <- c(
  "Tide pool" = "#e88a25",
  "Shallow reef" = "#1ea49d",
  "Deep reef" = "#482344"
)

# 3. Import and prepare data ---------------------------------------------------

# Coral color observations
color_raw <- read.csv2("data/COLOR.csv")

color_data <- color_raw %>%
  mutate(
    color = factor(cor),
    treatment = factor(trat, levels = c("C", "A")),
    site = factor(local, levels = site_levels),
    colony_id = factor(ID),
    experiment_time = factor(tempo),
    condition = case_when(
      cor %in% c("D1", "D2", "rosa") ~ "bleached",
      cor %in% c("D3", "D4", "D5", "D6") ~ "not-bleached",
      cor == "morto" ~ "dead",
      TRUE ~ NA_character_
    ),
    phase = case_when(
      tempo == "0" ~ "pre-stress",
      tempo == "1" ~ "ramp-up",
      tempo %in% c("2", "3") ~ "stress",
      tempo %in% c("4", "5", "6", "7") ~ "post-stress",
      TRUE ~ NA_character_
    ),
    condition = factor(
      condition,
      levels = c("not-bleached", "bleached", "dead")
    ),
    phase = factor(
      phase,
      levels = c("pre-stress", "ramp-up", "stress", "post-stress")
    )
  )

# Daily degree heating weeks (DHW) data
dhw_data <- read.csv("data/exp1_dhw.csv") %>%
  mutate(date = as.Date(date))

# Join DHW data to coral color observations
color_dhw <- color_data %>%
  mutate(date = as.Date(paste(ano, mes, dia, sep = "-"))) %>%
  left_join(dhw_data, by = "date")

# Binary bleaching response; dead colonies are excluded from this analysis
bleaching_data <- color_dhw %>%
  filter(condition != "dead") %>%
  mutate(
    bleaching = case_when(
      condition == "bleached" ~ 1,
      condition == "not-bleached" ~ 0,
      TRUE ~ NA_real_
    ),
    bleaching = as.numeric(bleaching)
  )

# PAM observations joined to DHW data
pam_raw <- read.csv2("data/PAM.csv", dec = ".")

pam_dhw <- pam_raw %>%
  mutate(date = as.Date(paste(ano, mes, dia, sep = "-"))) %>%
  left_join(dhw_data, by = "date") %>%
  mutate(
    site = factor(local, levels = site_levels),
    treatment = factor(trat, levels = c("C", "A")),
    colony_id = factor(ID)
  )

# 4. Bleaching analyses --------------------------------------------------------

# 4.1. Differences among sites in heated colonies during stress and post-stress
bleaching_heated <- bleaching_data %>%
  filter(
    treatment == "A",
    phase %in% c("stress", "post-stress")
  )

model_bleaching_site <- glmmTMB(
  bleaching ~ local + (1 | ID),
  family = binomial(),
  data = bleaching_heated
)

summary(model_bleaching_site)

sim_res_bleaching_site <- DHARMa::simulateResiduals(model_bleaching_site)
plot(sim_res_bleaching_site)
DHARMa::testUniformity(sim_res_bleaching_site)
DHARMa::testDispersion(sim_res_bleaching_site)
DHARMa::testOutliers(sim_res_bleaching_site)

# Estimated bleaching probabilities by site
emm_bleaching_site <- emmeans(
  model_bleaching_site,
  ~ local,
  type = "response"
)

emm_bleaching_site
pairs(emm_bleaching_site, type = "response")

# Secondary contrast: shallow reef (R) versus deep reef (F)
contrast_bleaching_RF <- contrast(
  emm_bleaching_site,
  method = list("R - F" = c(0, 1, -1)),
  adjust = "none"
)

summary(contrast_bleaching_RF, infer = TRUE)

emm_bleaching_site_df <- as.data.frame(emm_bleaching_site)

# Figure 3A: observed responses and estimated bleaching probabilities
fig_bleaching_site <- ggplot() +
  geom_jitter(
    data = bleaching_heated,
    aes(x = local, y = bleaching, color = local),
    alpha = 0.2,
    size = 2,
    position = position_jitter(width = 0.15, height = 0.04)
  ) +
  geom_errorbar(
    data = emm_bleaching_site_df,
    aes(x = local, ymin = asymp.LCL, ymax = asymp.UCL, color = local),
    width = 0.1
  ) +
  geom_point(
    data = emm_bleaching_site_df,
    aes(x = local, y = prob, color = local),
    size = 5
  ) +
  scale_color_manual(values = site_colors) +
  scale_x_discrete(name = NULL, labels = site_labels) +
  labs(y = "Predicted bleaching probability (95% CI)") +
  theme_classic(base_size = 16) +
  theme(legend.position = "none")

fig_bleaching_site

# 4.3. Relationship between DHW and bleaching probability among heated colonies
bleaching_heated_all <- bleaching_data %>%
  filter(treatment == "A")

model_bleaching_dhw <- glmmTMB(
  bleaching ~ local * dhw + (1 | ID),
  family = binomial(),
  data = bleaching_heated_all
)

summary(model_bleaching_dhw)

sim_res_bleaching_dhw <- DHARMa::simulateResiduals(model_bleaching_dhw)
plot(sim_res_bleaching_dhw)
DHARMa::testUniformity(sim_res_bleaching_dhw)
DHARMa::testDispersion(sim_res_bleaching_dhw)
DHARMa::testOutliers(sim_res_bleaching_dhw)

# Estimated DHW slopes for each site
dhw_trends_bleaching <- emtrends(
  model_bleaching_dhw,
  ~ local,
  var = "dhw"
)

dhw_trends_bleaching

# Secondary contrast: DHW slope in R versus F
contrast_dhw_bleaching_RF <- contrast(
  dhw_trends_bleaching,
  method = list("R - F" = c(0, 1, -1)),
  adjust = "none"
)

summary(contrast_dhw_bleaching_RF, infer = TRUE)

# Predicted bleaching probability across the observed DHW range
emm_bleaching_dhw <- emmeans(
  model_bleaching_dhw,
  ~ local | dhw,
  at = list(dhw = seq(0, 8, by = 0.1)),
  type = "response"
)

emm_bleaching_dhw_df <- as.data.frame(emm_bleaching_dhw)

# Figure 3B: predicted bleaching probability across DHW
fig_bleaching_dhw <- ggplot(
  emm_bleaching_dhw_df,
  aes(x = dhw, y = prob, color = local, fill = local)
) +
  geom_ribbon(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    alpha = 0.2,
    color = NA
  ) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) +
  labs(x = "DHW (°C-weeks)", y = NULL) +
  theme_classic(base_size = 16) +
  theme(legend.position = "none")

fig_bleaching_dhw

# Combine Figure 3 panels
fig3 <- (fig_bleaching_site + fig_bleaching_dhw) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 1.01)
  )

fig3
ggsave("figures/Fig3.png", fig3, width = 12, height = 6, dpi = 300)

# 5. Bleaching prevalence and color-class composition --------------------------

# 5.1. Daily bleaching prevalence data
prevalence_data <- read.csv("data/bleaching_prevalence_by_day_R.csv") %>%
  mutate(
    treatment = factor(treatment, levels = c("Heated", "Control")),
    origin_site = factor(
      origin_site,
      levels = c("Tide pool", "Shallow reef", "Deep reef")
    ),
    day = factor(day)
  )

fig_prevalence <- ggplot(
  prevalence_data,
  aes(
    x = day,
    y = bleaching_prevalence_among_alive,
    group = interaction(origin_site, treatment),
    color = interaction(origin_site, treatment),
    fill = interaction(origin_site, treatment)
  )
) +
  geom_ribbon(
    aes(ymin = exact_95pct_CI_lower, ymax = exact_95pct_CI_upper),
    alpha = 0.2,
    color = NA
  ) +
  geom_line(aes(linetype = treatment), linewidth = 0.8) +
  geom_point(aes(shape = treatment)) +
  facet_grid(origin_site ~ .) +
  scale_color_manual(
    values = c(
      "Tide pool.Heated" = "#e88a25",
      "Shallow reef.Heated" = "#1ea49d",
      "Deep reef.Heated" = "#482344",
      "Tide pool.Control" = "#e88a25",
      "Shallow reef.Control" = "#1ea49d",
      "Deep reef.Control" = "#482344"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Tide pool.Heated" = "#e88a25",
      "Shallow reef.Heated" = "#1ea49d",
      "Deep reef.Heated" = "#482344",
      "Tide pool.Control" = "grey70",
      "Shallow reef.Control" = "grey70",
      "Deep reef.Control" = "grey70"
    )
  ) +
  scale_shape_manual(values = c("Heated" = 17, "Control" = 16), name = NULL) +
  scale_linetype_discrete(name = NULL) +
  scale_y_continuous(breaks = c(0, 0.5, 1), name = NULL) +
  labs(
    x = "Experiment day",
    title = "Bleaching prevalence among living colonies"
  ) +
  guides(color = "none", fill = "none") +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 16)
  )

fig_prevalence

# 5.2. Proportion of color classes by site, day, and treatment
color_proportions <- color_dhw %>%
  count(local, condition, dia_exp, trat, .drop = FALSE) %>%
  group_by(local, dia_exp, trat) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    experiment_day = factor(dia_exp),
    condition = factor(
      condition,
      levels = c("not-bleached", "bleached", "dead")
    )
  )

fig_color_composition <- ggplot(
  color_proportions,
  aes(x = experiment_day, y = proportion, fill = condition)
) +
  geom_col() +
  facet_grid(
    local ~ trat,
    labeller = labeller(
      trat = c("C" = "Control treatment", "A" = "Heated treatment"),
      local = function(x) rep("", length(x))
    )
  ) +
  scale_fill_manual(
    values = c(
      "not-bleached" = "lightgrey",
      "bleached" = "#e88a25",
      "dead" = "#0a0310"
    ),
    labels = c(
      "not-bleached" = "Not bleached",
      "bleached" = "Bleached",
      "dead" = "Dead"
    ),
    name = NULL
  ) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1),
    name = "Proportion of colonies"
  ) +
  labs(x = "Experiment day") +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank()
  )

fig_color_composition

# Combine Figure 2 panels
fig2 <- (fig_color_composition + fig_prevalence) +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 0.99)
  )

fig2
ggsave("figures/Fig2.png", fig2, width = 12, height = 8, dpi = 300)

# 6. Photochemical efficiency (YII) over experimental time ---------------------

# Summarize YII by experiment day, treatment, and site
yii_summary <- pam_raw %>%
  group_by(dia_exp, trat, local) %>%
  summarise(
    mean_yii = mean(YIIdec, na.rm = TRUE),
    se_yii = sd(YIIdec, na.rm = TRUE) / sqrt(sum(!is.na(YIIdec))),
    n_colonies = sum(!is.na(YIIdec)),
    .groups = "drop"
  )

# Check expected sample sizes; stop if the data do not match expectations
if (
  any(yii_summary$n_colonies[yii_summary$trat == "A"] != 12) ||
  any(yii_summary$n_colonies[yii_summary$trat == "C"] != 8)
) {
  warning("Some treatment-by-site-by-day groups do not have the expected sample size.")
}

# Mean temperature in heated systems, by experimental day
temperature_by_day <- pam_raw %>%
  group_by(dia_exp) %>%
  summarise(temp = mean(temp, na.rm = TRUE), .groups = "drop")

x_grid <- seq(
  min(temperature_by_day$dia_exp, na.rm = TRUE),
  max(temperature_by_day$dia_exp, na.rm = TRUE),
  length.out = 1200
)

temperature_background <- tibble(
  day = x_grid,
  y = 0.9,
  heated_temp = approx(
    x = temperature_by_day$dia_exp,
    y = temperature_by_day$temp,
    xout = x_grid
  )$y
)

temperature_labels <- temperature_by_day %>%
  filter(dia_exp %in% c(1, 6, 11, 15, 19, 24, 32, 38, 45)) %>%
  mutate(label = sprintf("%.1f", temp))

yii_summary <- yii_summary %>%
  mutate(
    site = factor(local, levels = site_levels, labels = site_labels),
    treatment = factor(trat, levels = c("A", "C"), labels = c("Heated", "Control"))
  )

# Figure 4A: observed YII through the experiment
fig_yii_time <- ggplot() +
  geom_tile(
    data = temperature_background,
    aes(x = day, y = y, fill = heated_temp),
    width = (max(x_grid) - min(x_grid)) / (length(x_grid) - 1),
    height = 0.065
  ) +
  annotate(
    "text",
    x = 8.3,
    y = 0.96,
    label = "Temperature in heated systems (°C)",
    size = 4,
    color = "black"
  ) +
  geom_text(
    data = temperature_labels,
    aes(x = dia_exp, y = 0.9, label = label),
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  geom_vline(
    xintercept = 16,
    linetype = "dotted",
    color = "#eb0258",
    linewidth = 0.5,
    alpha = 0.6
  ) +
  geom_errorbar(
    data = yii_summary,
    aes(
      x = dia_exp,
      ymin = mean_yii - se_yii,
      ymax = mean_yii + se_yii,
      color = site,
      group = interaction(site, treatment)
    ),
    linewidth = 0.5,
    alpha = 0.7,
    width = 0.5
  ) +
  geom_line(
    data = yii_summary,
    aes(
      x = dia_exp,
      y = mean_yii,
      color = site,
      linetype = treatment,
      group = interaction(site, treatment)
    ),
    linewidth = 1,
    alpha = 0.5
  ) +
  geom_point(
    data = yii_summary,
    aes(
      x = dia_exp,
      y = mean_yii,
      color = site,
      shape = treatment
    ),
    size = 5
  ) +
  scale_color_manual(
    name = "Origin",
    values = site_colors_labeled
  ) +
  scale_linetype_manual(
    name = "Treatment",
    values = c("Heated" = "solid", "Control" = "dashed")
  ) +
  scale_shape_manual(
    name = "Treatment",
    values = c("Heated" = 17, "Control" = 16)
  ) +
  scale_fill_gradientn(
    guide = "none",
    colors = c("#F7FAF9", "#F9DFCE", "#f7a29e", "#eb0258"),
    limits = c(26.9, 32),
    breaks = c(27, 29, 32)
  ) +
  scale_x_continuous(
    limits = c(0, 46),
    breaks = c(1, 6, 11, 15, 19, 24, 32, 38, 45)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 0.8, by = 0.2)
  ) +
  labs(
    x = "Experiment day",
    y = "Photochemical efficiency (YII ± SE)"
  ) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2),
    shape = guide_legend(order = 2)
  ) +
  theme_classic(base_size = 16) +
  theme(legend.position = "right")

fig_yii_time
ggsave("figures/Fig4a.png", fig_yii_time, width = 12, height = 7, dpi = 300)

# 7. Photochemical efficiency in relation to DHW -------------------------------

yii_heated <- pam_dhw %>%
  filter(treatment == "A")

# Beta regression with a zero-inflation component
model_yii_dhw <- glmmTMB(
  YIIdec ~ dhw * site + (1 | colony_id),
  ziformula = ~ 1 + (1 | site),
  family = beta_family(link = "logit"),
  data = yii_heated
)

summary(model_yii_dhw)

sim_res_yii_dhw <- DHARMa::simulateResiduals(model_yii_dhw, n = 1000)
plot(sim_res_yii_dhw)
DHARMa::testUniformity(sim_res_yii_dhw)
DHARMa::testDispersion(sim_res_yii_dhw)
DHARMa::testOutliers(sim_res_yii_dhw)

# Estimated YII across DHW values, by site
emm_yii_dhw <- emmeans(
  model_yii_dhw,
  ~ site | dhw,
  at = list(dhw = seq(0, 8, by = 0.1)),
  type = "response"
)

emm_yii_dhw_df <- as.data.frame(emm_yii_dhw)

# Estimated DHW slopes for each site
dhw_trends_yii <- emtrends(
  model_yii_dhw,
  ~ site,
  var = "dhw"
)

dhw_trends_yii

# Secondary contrast: DHW slope in R versus F
contrast_dhw_yii_RF <- contrast(
  dhw_trends_yii,
  method = list("R - F" = c(0, 1, -1)),
  adjust = "none"
)

summary(contrast_dhw_yii_RF, infer = TRUE)

# Figure 4C: model-predicted YII across DHW
fig_yii_dhw_predicted <- ggplot(
  emm_yii_dhw_df,
  aes(x = dhw, y = response, color = site, fill = site)
) +
  geom_ribbon(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    alpha = 0.2,
    color = NA
  ) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = site_colors_labeled, name = NULL) +
  scale_fill_manual(values = site_colors_labeled, name = NULL) +
  labs(
    x = "DHW (°C-weeks)",
    y = "Predicted photochemical efficiency (95% CI)"
  ) +
  theme_classic(base_size = 16) +
  theme(legend.position = "none")

fig_yii_dhw_predicted

# Figure 4B: observed YII and summary trends across DHW
fig_yii_dhw_observed <- ggplot(
  yii_heated,
  aes(x = dhw, y = YIIdec, color = site, group = site)
) +
  geom_point(
    position = position_jitter(width = 0.1, height = 0),
    alpha = 0.2
  ) +
  stat_summary(geom = "pointrange", fun.data = mean_se, size = 1) +
  stat_summary(geom = "line", fun = mean, linewidth = 1) +
  scale_color_manual(values = site_colors_labeled, name = NULL) +
  labs(
    x = "DHW (°C-weeks)",
    y = "Photochemical efficiency (YII ± SE)"
  ) +
  theme_classic(base_size = 16) +
  theme(legend.position = "none")

fig_yii_dhw_observed

# Combine Figure 4BC panels
fig5 <- fig_yii_dhw_observed + fig_yii_dhw_predicted
fig5
ggsave("figures/Fig4BC.png", fig5, width = 12, height = 6, dpi = 300)

# Combine Figures 4A and 4BC as a multi-panel figure
fig6 <- (fig_yii_time / fig5) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 0.99)
  )

fig6
ggsave("figures/Fig4.png", fig6, width = 12, height = 13, dpi = 300)

# End of script ---------------------------------------------------------------
