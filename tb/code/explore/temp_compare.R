# Diagnostic comparison of the five reanalysis temperature products.
# Self-contained: reads the CSV and applies the balanced filter from
# distance.R, then produces a comprehensive set of diagnostic plots saved
# to figs/ with a `tc_` prefix.

library(tidyverse)
library(fixest)
theme_set(theme_classic())
setwd('/Users/tombearpark/Library/CloudStorage/Dropbox/ks')

fig_dir <- 'impact_uncertainty/tb/figs'
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
save_fig <- function(p, name, w = 7, h = 4)
  ggsave(file.path(fig_dir, name), p, width = w, height = h, units = 'in')

# 1. Load + balance ----------------------------------------------------------
df <- read_csv('reg_inputs/gdppc/reginputs_gdppc_1deg.csv', show_col_types = FALSE) %>%
  mutate(temp_1 = temp_popweight,
         time1  = year - min(year),
         time2  = time1^2)
mods <- sort(unique(df$obs_mod))
keep <- df %>% filter(!is.na(temp_1)) %>%
  distinct(iso, year, obs_mod) %>% count(iso, year) %>%
  filter(n == length(mods)) %>% select(iso, year)
df <- df %>% semi_join(keep, by = c('iso','year'))
message(sprintf('balanced panel: %d (iso, year) pairs × %d datasets',
                nrow(keep), length(mods)))

# Wide table on temp_1, indexed by iso-year ----------------------------------
wide <- df %>%
  select(iso, year, continent, obs_mod, temp_1) %>%
  pivot_wider(names_from = obs_mod, values_from = temp_1)

pairs_df <- combn(mods, 2) %>% t() %>% as.data.frame() %>%
  setNames(c('A','B')) %>% as_tibble()

# 2. Pairwise correlation + mean |diff| ---------------------------------------
pair_summary <- pairs_df %>%
  rowwise() %>%
  mutate(cor = cor(wide[[A]], wide[[B]]),
         mean_abs_diff = mean(abs(wide[[A]] - wide[[B]])),
         mean_signed_diff = mean(wide[[A]] - wide[[B]])) %>%
  ungroup()

cat('\n--- Pairwise temperature agreement (raw) ---\n')
print(pair_summary %>% arrange(desc(mean_abs_diff)))

# Pairwise MAD heatmap
mad_mat <- matrix(0, length(mods), length(mods), dimnames = list(mods, mods))
for (i in seq_len(nrow(pair_summary))) {
  a <- pair_summary$A[i]; b <- pair_summary$B[i]; v <- pair_summary$mean_abs_diff[i]
  mad_mat[a,b] <- v; mad_mat[b,a] <- v
}
mad_long <- as_tibble(mad_mat, rownames = 'A') %>%
  pivot_longer(-A, names_to = 'B', values_to = 'mad')
p_mad <- ggplot(mad_long, aes(A, B, fill = mad)) +
  geom_tile() +
  geom_text(aes(label = sprintf('%.2f', mad)), size = 3, color = 'white') +
  scale_fill_gradient(low = 'white', high = 'red') +
  coord_equal() + labs(x = NULL, y = NULL, fill = 'Mean |diff| (°C)')
save_fig(p_mad, 'tc_pairwise_mad.pdf', w = 5.5, h = 5)

# 3. Country-level disagreement (ERA5 vs MERRA2) ------------------------------
A <- 'ERA5-025'; B <- 'MERRA2'
country_AB <- wide %>%
  mutate(d = .data[[A]] - .data[[B]]) %>%
  group_by(iso, continent) %>%
  summarise(mad = mean(abs(d)), bias = mean(d), sd = sd(d),
            n = n(), .groups = 'drop')

top20 <- country_AB %>% arrange(desc(mad)) %>% head(20)
p_cbars <- ggplot(top20, aes(x = reorder(iso, mad), y = mad, fill = continent)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = 'Mean |ERA5 − MERRA2| (°C)', fill = NULL,
       title = 'Top 20 countries by ERA5 vs MERRA2 disagreement') +
  theme(legend.position = 'bottom')
save_fig(p_cbars, 'tc_country_bars.pdf', w = 8, h = 5.5)

cat('\n--- Top 20 countries by ERA5 vs MERRA2 MAD ---\n')
print(top20)

# 4. Temporal disagreement (all pairs) ---------------------------------------
year_disagr <- pmap_dfr(pair_summary %>% select(A, B),
  function(A, B) {
    wide %>% mutate(d = abs(.data[[A]] - .data[[B]])) %>%
      group_by(year) %>% summarise(mad = mean(d), .groups = 'drop') %>%
      mutate(pair = paste(A, '-', B))
  })

p_years <- ggplot(year_disagr, aes(year, mad, color = pair)) +
  geom_line() +
  labs(x = 'Year', y = 'Mean |diff| across countries (°C)', color = 'Pair',
       title = 'Temporal evolution of pairwise disagreement')
save_fig(p_years, 'tc_year_series.pdf', w = 9, h = 4)

# 5. Residualise each dataset on iso FE + country-specific linear trend ------
resid_df <- map_dfr(mods, function(m) {
  d <- df %>% filter(obs_mod == m) %>%
    drop_na(temp_1)
  r <- resid(feols(
    temp_1 ~ 1 | iso + year + iso[time1] + iso[time2], d, notes = FALSE))
  tibble(iso = d$iso, year = d$year, continent = d$continent,
         obs_mod = m, temp_1 = d$temp_1, resid = r)
})

wide_r <- resid_df %>%
  select(iso, year, continent, obs_mod, resid) %>%
  pivot_wider(names_from = obs_mod, values_from = resid)

# Residualised correlation heatmap — this is what matters for panel ID.
cor_mat_r <- wide_r %>% select(all_of(mods)) %>% cor(use = 'pairwise.complete.obs')
cor_long_r <- as_tibble(cor_mat_r, rownames = 'A') %>%
  pivot_longer(-A, names_to = 'B', values_to = 'cor')
p_cor_r <- ggplot(cor_long_r, aes(A, B, fill = cor)) +
  geom_tile() +
  geom_text(aes(label = sprintf('%.3f', cor)), size = 3) +
  scale_fill_gradient(low = 'white', high = 'darkblue') +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = 'Residual\ncorrelation',
       title = 'Correlation of residualised temperatures\n(iso FE + iso×time trend removed)')
save_fig(p_cor_r, 'tc_cor_heatmap_resid.pdf', w = 6, h = 5.5)

cat('\n--- Residualised pairwise correlation ---\n')
print(round(cor_mat_r, 3))

# Residualised pairwise MAD
pair_summary_r <- pairs_df %>%
  rowwise() %>%
  mutate(cor_resid = cor(wide_r[[A]], wide_r[[B]]),
         mad_resid = mean(abs(wide_r[[A]] - wide_r[[B]]))) %>%
  ungroup()
cat('\n--- Residualised pairwise MAD ---\n')
print(pair_summary_r %>% arrange(desc(mad_resid)))

# Residual ERA5 vs MERRA2 scatter
p_rsc <- wide_r %>%
  ggplot(aes(.data[[B]], .data[[A]])) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  coord_equal() +
  labs(x = paste(B, 'residual (°C)'),
       y = paste(A, 'residual (°C)'),
       title = 'Residualised temp (iso FE + iso×time trend partialled out)')
save_fig(p_rsc, 'tc_resid_scatter.pdf', w = 5.5, h = 5.5)

# Country MAD on residuals (the "weather" disagreement, climate removed)
country_AB_r <- wide_r %>%
  mutate(d = .data[[A]] - .data[[B]]) %>%
  group_by(iso, continent) %>%
  summarise(mad = mean(abs(d)), .groups = 'drop')

top20_r <- country_AB_r %>% arrange(desc(mad)) %>% head(20)
p_cbars_r <- ggplot(top20_r, aes(x = reorder(iso, mad), y = mad, fill = continent)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = 'Mean |residual diff| (°C)', fill = NULL,
       title = 'Top 20 countries by residualised ERA5 vs MERRA2 disagreement') +
  theme(legend.position = 'bottom')
save_fig(p_cbars_r, 'tc_resid_country_bars.pdf', w = 8, h = 5.5)

cat('\n--- Top 20 countries by residualised ERA5 vs MERRA2 MAD ---\n')
print(top20_r)

# 6. Time series of residuals — top 8 vs bottom 8 disagreement countries ----
# Rank by average residualised MAD across ALL 10 pairs per country.
country_ensemble_mad <- pmap_dfr(pair_summary_r %>% select(A, B),
  function(A, B) {
    wide_r %>% transmute(iso, d = abs(.data[[A]] - .data[[B]])) %>%
      group_by(iso) %>% summarise(mad = mean(d), .groups = 'drop')
  }) %>%
  group_by(iso) %>% summarise(avg_mad = mean(mad), .groups = 'drop') %>%
  arrange(desc(avg_mad))

top8_r    <- country_ensemble_mad %>% head(8) %>% pull(iso)
bottom8_r <- country_ensemble_mad %>% tail(8) %>% arrange(avg_mad) %>% pull(iso)

make_rts_panel <- function(isos, title) {
  resid_df %>%
    filter(iso %in% isos) %>%
    mutate(iso = factor(iso, levels = isos)) %>%
    ggplot(aes(year, resid, color = obs_mod)) +
    geom_line(linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40') +
    facet_wrap(~iso, scales = 'free_y', ncol = 2) +
    labs(x = 'Year', y = 'Residualised temp (°C)', color = NULL,
         title = title) +
    theme(legend.position = 'bottom')
}

save_fig(make_rts_panel(top8_r,
           'Top 8 disagreement countries (ensemble residualised MAD)'),
         'tc_resid_ts_top.pdf', w = 7, h = 8)
save_fig(make_rts_panel(bottom8_r,
           'Bottom 8 disagreement countries (ensemble residualised MAD)'),
         'tc_resid_ts_bottom.pdf', w = 7, h = 8)

# 7. Temperature-level dependence of the bias, all pairs ---------------------
bias_by_level <- pmap_dfr(pair_summary %>% select(A, B),
  function(A, B) {
    wide %>% transmute(t_B = .data[[B]], d = .data[[A]] - .data[[B]],
                        pair = paste(A, '-', B))
  })

p_bias <- ggplot(bias_by_level, aes(t_B, d)) +
  geom_point(alpha = 0.05, size = 0.3) +
  geom_smooth(method = 'loess', formula = y ~ x, color = 'black', se = FALSE) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  facet_wrap(~pair, ncol = 5) +
  labs(x = 'Temp (°C) [second dataset]', y = 'Diff A − B (°C)',
       title = 'Pairwise temperature bias vs level')
save_fig(p_bias, 'tc_bias_by_level.pdf', w = 12, h = 6)

# 8. Per-country linear trend, compared across datasets ----------------------
trend_df <- map_dfr(mods, function(m) {
  df %>% filter(obs_mod == m) %>%
    group_by(iso) %>%
    summarise(slope = coef(lm(temp_1 ~ year))['year'], .groups = 'drop') %>%
    mutate(obs_mod = m)
}) %>% pivot_wider(names_from = obs_mod, values_from = slope)

p_trend <- ggplot(trend_df, aes(.data[[B]], .data[[A]])) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  coord_equal() +
  labs(x = sprintf('%s trend (°C/yr)', B),
       y = sprintf('%s trend (°C/yr)', A),
       title = 'Per-country linear trend in temperature: ERA5 vs MERRA2')
save_fig(p_trend, 'tc_trend_scatter.pdf', w = 5.5, h = 5.5)

# 9. Top 20 (iso, year) outliers in the raw ERA5-MERRA2 difference -----------
outliers <- wide %>%
  mutate(d = .data[[A]] - .data[[B]],
         abs_d = abs(d)) %>%
  arrange(desc(abs_d)) %>%
  select(iso, year, continent,
         !!A := !!sym(A), !!B := !!sym(B), d, abs_d) %>%
  head(20)
cat('\n--- Top 20 ERA5 vs MERRA2 (iso, year) outliers ---\n')
print(outliers)

# 9b. 2009-2010 spike: GMFD vs ensemble mean ---------------------------------
# The raw year-series shows a spike in pairs involving GMFD. Isolate it by
# comparing GMFD to the mean of the other four datasets, year by year.
gmfd_dev <- wide %>%
  mutate(others_mean = (`ERA5-025` + `JRA-3Q` + MERRA2 + UDel) / 4,
         dev = GMFD - others_mean,
         abs_dev = abs(dev))

gmfd_year <- gmfd_dev %>%
  group_by(year) %>%
  summarise(mean_abs_dev = mean(abs_dev),
            mean_signed  = mean(dev), .groups = 'drop')

p_gmfd_yr <- ggplot(gmfd_year, aes(year, mean_abs_dev)) +
  geom_line() + geom_point(size = 1) +
  geom_vline(xintercept = c(2009, 2010), linetype = 'dashed', color = 'red') +
  labs(x = 'Year', y = 'Mean |GMFD − others| (°C)',
       title = 'GMFD deviation from the mean of ERA5 / JRA-3Q / MERRA2 / UDel')
save_fig(p_gmfd_yr, 'tc_gmfd_year.pdf', w = 8, h = 3.8)

# Per-pair zoom 2005-2013
pair_zoom <- year_disagr %>% filter(year >= 2005) %>%
  mutate(has_GMFD = str_detect(pair, 'GMFD'))
p_zoom <- ggplot(pair_zoom, aes(year, mad, color = pair, linetype = has_GMFD)) +
  geom_line(linewidth = 0.6) +
  scale_linetype_manual(values = c(`FALSE` = 'dotted', `TRUE` = 'solid'),
                        guide = 'none') +
  geom_vline(xintercept = c(2009, 2010), linetype = 'dashed',
             color = 'grey60') +
  labs(x = 'Year', y = 'Mean |diff| (°C)',
       title = 'Pairwise disagreement 2005–2013 (GMFD pairs solid)',
       color = NULL)
save_fig(p_zoom, 'tc_year_zoom.pdf', w = 9, h = 4)

# Top 15 country deviations in 2010, signed, ordered
top2010 <- gmfd_dev %>%
  filter(year == 2010) %>%
  arrange(desc(abs_dev)) %>% head(15) %>%
  mutate(iso = factor(iso, levels = rev(iso)))

p_2010 <- ggplot(top2010, aes(x = iso, y = dev, fill = dev > 0)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c(`TRUE` = 'firebrick', `FALSE` = 'steelblue'),
                    labels = c(`TRUE` = 'GMFD warmer', `FALSE` = 'GMFD cooler')) +
  geom_hline(yintercept = 0) +
  labs(x = NULL, y = 'GMFD − mean-of-others (°C), 2010',
       title = 'Top 15 countries driving the GMFD 2010 anomaly',
       fill = NULL) +
  theme(legend.position = 'bottom')
save_fig(p_2010, 'tc_2010_countries.pdf', w = 7, h = 5)

# Same for 2009 so the user can compare the two spike years
top2009 <- gmfd_dev %>%
  filter(year == 2009) %>%
  arrange(desc(abs_dev)) %>% head(15) %>%
  mutate(iso = factor(iso, levels = rev(iso)))

p_2009 <- ggplot(top2009, aes(x = iso, y = dev, fill = dev > 0)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c(`TRUE` = 'firebrick', `FALSE` = 'steelblue'),
                    labels = c(`TRUE` = 'GMFD warmer', `FALSE` = 'GMFD cooler')) +
  geom_hline(yintercept = 0) +
  labs(x = NULL, y = 'GMFD − mean-of-others (°C), 2009',
       title = 'Top 15 countries driving the GMFD 2009 anomaly',
       fill = NULL) +
  theme(legend.position = 'bottom')
save_fig(p_2009, 'tc_2009_countries.pdf', w = 7, h = 5)

# 10. Disagreement vs income / climate ---------------------------------------
# Per-country, average residualised MAD across ALL 10 pairs.
country_pair_mad_r <- pmap_dfr(pair_summary %>% select(A, B),
  function(A, B) {
    wide_r %>% transmute(iso, d = abs(.data[[A]] - .data[[B]])) %>%
      group_by(iso) %>% summarise(mad = mean(d), .groups = 'drop') %>%
      mutate(pair = paste(A, '-', B))
  })
country_disagr <- country_pair_mad_r %>%
  group_by(iso) %>%
  summarise(avg_mad_resid = mean(mad), .groups = 'drop')

# Same for raw (bias-inclusive)
country_pair_mad_raw <- pmap_dfr(pair_summary %>% select(A, B),
  function(A, B) {
    wide %>% transmute(iso, d = abs(.data[[A]] - .data[[B]])) %>%
      group_by(iso) %>% summarise(mad = mean(d), .groups = 'drop') %>%
      mutate(pair = paste(A, '-', B))
  })
country_disagr_raw <- country_pair_mad_raw %>%
  group_by(iso) %>%
  summarise(avg_mad_raw = mean(mad), .groups = 'drop')

# Country covariates (one obs_mod, since gdpCAPppp_wdi is obs_mod-invariant).
country_cov <- df %>% filter(obs_mod == mods[1]) %>%
  group_by(iso, continent) %>%
  summarise(log_gdp   = mean(log(gdpCAPppp_wdi), na.rm = TRUE),
            mean_temp = mean(temp_1, na.rm = TRUE),
            .groups = 'drop')

dm <- country_disagr %>%
  left_join(country_disagr_raw, by = 'iso') %>%
  left_join(country_cov,        by = 'iso') %>%
  filter(is.finite(log_gdp), is.finite(mean_temp))

cat('\n--- Disagreement vs covariate correlations ---\n')
print(dm %>% summarise(
  cor_resid_logGDP = cor(avg_mad_resid, log_gdp),
  cor_resid_temp   = cor(avg_mad_resid, mean_temp),
  cor_raw_logGDP   = cor(avg_mad_raw,   log_gdp),
  cor_raw_temp     = cor(avg_mad_raw,   mean_temp)))

label_top <- dm %>% arrange(desc(avg_mad_resid)) %>% head(8) %>% pull(iso)

# Helper: linear fit with R² annotated top-right of the plot.
scatter_lm <- function(data, x, y, xlab, ylab, label_isos = NULL) {
  r2 <- summary(lm(data[[y]] ~ data[[x]]))$r.squared
  p <- ggplot(data, aes(.data[[x]], .data[[y]])) +
    geom_point(aes(color = continent), size = 2, alpha = 0.8) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE,
                color = 'black', linetype = 'dashed') +
    annotate('text', x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = sprintf('R² = %.3f', r2), size = 4) +
    labs(x = xlab, y = ylab, color = NULL) +
    theme(legend.position = 'bottom')
  if (!is.null(label_isos))
    p <- p + ggrepel::geom_text_repel(
      data = data %>% filter(iso %in% label_isos),
      aes(label = iso), size = 3, show.legend = FALSE, max.overlaps = Inf)
  p
}

p_di <- scatter_lm(dm, 'log_gdp', 'avg_mad_resid',
                   'log mean GDP per capita (PPP)',
                   'Average residualised MAD across pairs (°C)',
                   label_isos = label_top)
save_fig(p_di, 'tc_disagr_vs_income.pdf', w = 8, h = 5)

p_dt <- scatter_lm(dm, 'mean_temp', 'avg_mad_resid',
                   'Country mean temperature (°C)',
                   'Average residualised MAD across pairs (°C)',
                   label_isos = label_top)
save_fig(p_dt, 'tc_disagr_vs_temp.pdf', w = 8, h = 5)

# Raw variants too, for contrast — also linear fit + R².
p_di_raw <- scatter_lm(dm, 'log_gdp', 'avg_mad_raw',
                       'log mean GDP per capita (PPP)',
                       'Average RAW MAD across pairs (°C)')
save_fig(p_di_raw, 'tc_disagr_raw_vs_income.pdf', w = 8, h = 5)

p_dt_raw <- scatter_lm(dm, 'mean_temp', 'avg_mad_raw',
                       'Country mean temperature (°C)',
                       'Average RAW MAD across pairs (°C)')
save_fig(p_dt_raw, 'tc_disagr_raw_vs_temp.pdf', w = 8, h = 5)

cat('\nAll figures written to', fig_dir, '\n')
