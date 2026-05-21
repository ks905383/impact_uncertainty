# Standalone figure export for distance_slides.tex.
# Sources distance.R to reuse its fits and objects, then saves PDFs to figs/.

source('impact_uncertainty/tb/distance.R')

fig_dir <- 'impact_uncertainty/tb/figs'
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
save_fig <- function(p, name, w = 7, h = 4) {
  ggsave(file.path(fig_dir, name), p, width = w, height = h, units = 'in')
}

# 1. ME curves with delta-method CIs, one panel per model
p_me <- ggplot(ME, aes(x = temp, y = me)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  facet_wrap(~model, nrow = 1) +
  labs(x = 'Temperature (°C)', y = 'Marginal effect on growth')
save_fig(p_me, 'me_curves.pdf', w = 10, h = 3.2)

# 2. Pairwise distance heatmap
p_D <- ggplot(D_long, aes(x = model_a, y = model_b, fill = distance)) +
  geom_tile() +
  geom_text(aes(label = sprintf('%.4f', distance)), color = 'white', size = 3.2) +
  scale_fill_gradient(low = 'white', high = 'red') +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = 'L1 distance')
save_fig(p_D, 'distance_matrix.pdf', w = 6, h = 5)

# 3-5. LOO bar charts
save_fig(plot_loo(LOO_cont, 'Continent LOO (ERA5-025 vs MERRA2)'),
         'loo_continent.pdf', w = 7, h = 3.2)
save_fig(plot_loo(LOO_dec,  'Decade LOO (ERA5-025 vs MERRA2)'),
         'loo_decade.pdf', w = 7, h = 2.8)
save_fig(plot_loo(LOO_tq,   'Temperature-decile LOO (ERA5-025 vs MERRA2)'),
         'loo_decile.pdf', w = 7, h = 4)

# 6. Temp agreement scatter
p_sc <- ggplot(AB, aes(t_B, t_A)) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  coord_equal() +
  labs(x = 'MERRA2 temp (°C)', y = 'ERA5-025 temp (°C)')
save_fig(p_sc, 'temp_scatter.pdf', w = 5, h = 5)

# 6.1 Temp agreement scatter, residualised
resids <- tibble(
  ERA5=resid(feols(temp_1 ~ 1 | iso + year + iso[time1] + iso[time2], df %>% filter(obs_mod == "ERA5-025"))),
  MERRA=resid(feols(temp_1 ~ 1 | iso + year + iso[time1] + iso[time2], df %>% filter(obs_mod == "MERRA2"))),
       ) 
feols(MERRA ~ ERA5, resids) %>% summary()
resids %>% 
  ggplot(aes(x = MERRA, y = ERA5)) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  coord_equal() +
  labs(x = 'MERRA2 temp (°C)', y = 'ERA5-025 temp (°C)')

# 7. Temp bias by level (loess)
p_bias <- ggplot(AB, aes(t_B, dT)) +
  geom_point(alpha = 0.1, size = 0.4) +
  geom_smooth(method = 'loess', formula = y ~ x, color = 'black') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(x = 'MERRA2 temp (°C)', y = 'ERA5-025 minus MERRA2 (°C)')
save_fig(p_bias, 'temp_bias.pdf', w = 7, h = 4)

# 8. Sierra Leone temp series
p_sle <- df %>%
  filter(iso == 'SLE', obs_mod %in% c(A, B)) %>%
  ggplot(aes(year, temp_1, color = obs_mod)) +
  geom_line() +
  geom_point(size = 1) +
  labs(title = NULL, x = 'Year', y = 'Temperature (°C)', color = NULL) +
  theme(legend.position = 'bottom')
save_fig(p_sle, 'sle_temp.pdf', w = 7, h = 3.5)

cat('\nWrote figures to', fig_dir, '\n')
