# Disagreement vs covariate scatters at multiple granularities, plus the
# regression of disagreement on income + climate + interaction, and three
# visualisations of the interaction (heatmap, marginal-effect ribbons, actual
# 2D scatter).
#
# Produces in figs/:
#   tc_disagr_vs_income.pdf      country-level, x = log GDP/cap   (the original)
#   tc_disagr_vs_temp.pdf        country-level, x = mean temp     (the original)
#   tc_disagr_yr_income.pdf      year-level
#   tc_disagr_yr_climate.pdf
#   tc_disagr_cy_income.pdf      country-year level
#   tc_disagr_cy_climate.pdf
#   tc_disagr_heatmap.pdf        interaction predictions surface
#   tc_disagr_marginal.pdf       marginal-effect ribbons
#   tc_disagr_actual_2d.pdf      actual disagreement on the (gdp, temp) plane

library(tidyverse)
library(fixest)
library(patchwork)
theme_set(theme_classic())
setwd('/Users/tombearpark/Library/CloudStorage/Dropbox/ks')
source('impact_uncertainty/tb/code/disagree_helpers.R')

# ============================================================================
# Data load + balance + residualise
# ============================================================================
df <- read_csv('reg_inputs/gdppc/reginputs_gdppc_1deg.csv', show_col_types = FALSE) %>%
  mutate(temp_1 = temp_popweight,
         time1  = year - min(year),
         time2  = time1^2)
mods <- sort(unique(df$obs_mod))
keep <- df %>% filter(!is.na(temp_1)) %>%
  distinct(iso, year, obs_mod) %>% count(iso, year) %>%
  filter(n == length(mods)) %>% select(iso, year)
df <- df %>% semi_join(keep, by = c('iso','year'))

resid_df <- map_dfr(mods, function(m) {
  d <- df %>% filter(obs_mod == m) %>% drop_na(temp_1)
  r <- resid(feols(temp_1 ~ 1 | iso + year + iso[time1] + iso[time2],
                   d, notes = FALSE))
  tibble(iso = d$iso, year = d$year, obs_mod = m, resid = r)
})
wide_r <- resid_df %>% pivot_wider(names_from = obs_mod, values_from = resid)

# ============================================================================
# Project-specific helpers
# ============================================================================

# Aggregate covariates (log GDP, mean temp, continent passthrough) to grouping.
make_country_cov <- function(df, mods, group_cols) {
  df %>% filter(obs_mod == mods[1]) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(log_gdp   = mean(log(gdpCAPppp_wdi), na.rm = TRUE),
              temp      = mean(temp_popweight,    na.rm = TRUE),
              continent = if ('iso' %in% group_cols) first(continent)
                          else NA_character_,
              .groups = 'drop')
}

# Driver: build dm at a given granularity, save the figure, return dm.
build_and_save <- function(group_cols, x, label_var, file, x_lab) {
  dm <- make_disagr_df(wide_r, mods, group_cols, metric = 'mad') %>%
    left_join(make_country_cov(df, mods, group_cols), by = group_cols) %>%
    filter(is.finite(.data[[x]]))
  p <- plot_disagr(dm, x = x, label_var = label_var, x_lab = x_lab,
                   colour_var = 'continent')
  ggsave(file.path('impact_uncertainty/tb/figs', file),
         p, width = 8, height = 5, units = 'in')
  r2 <- summary(lm(reformulate(x, 'disagr'), dm))$r.squared
  cat(sprintf('  %-32s  N = %4d  R² = %.3f\n', file, nrow(dm), r2))
  invisible(dm)
}

# ============================================================================
# Build the six scatter figures
# ============================================================================
cat('\nBuilding disagreement scatters:\n')

build_and_save(c('iso'),        'log_gdp', 'iso',  'tc_disagr_vs_income.pdf',
               'log mean GDP per capita (PPP)')
build_and_save(c('iso'),        'temp',    'iso',  'tc_disagr_vs_temp.pdf',
               'Country mean temperature (°C)')
build_and_save(c('year'),       'log_gdp', 'year', 'tc_disagr_yr_income.pdf',
               'log world-mean GDP per capita (PPP)')
build_and_save(c('year'),       'temp',    'year', 'tc_disagr_yr_climate.pdf',
               'World-mean temperature (°C)')
build_and_save(c('iso','year'), 'log_gdp', NULL,  'tc_disagr_cy_income.pdf',
               'log GDP per capita (PPP)')
build_and_save(c('iso','year'), 'temp',    NULL,  'tc_disagr_cy_climate.pdf',
               'Realised temperature (°C)')

# ============================================================================
# Regressions of disagreement on income + climate
# ============================================================================
dm_iso <- make_disagr_df(wide_r, mods, c('iso'), metric = 'mad') %>%
  left_join(make_country_cov(df, mods, c('iso')), by = 'iso') %>%
  filter(is.finite(log_gdp), is.finite(temp))

m_g  <- feols(disagr ~ log_gdp,        dm_iso, vcov = 'hetero')
m_t  <- feols(disagr ~ temp,           dm_iso, vcov = 'hetero')
m_gt <- feols(disagr ~ log_gdp + temp, dm_iso, vcov = 'hetero')
m_gt_int <- feols(disagr ~ log_gdp + temp + log_gdp * temp,
                  dm_iso, vcov = 'hetero')

cat('\n=== Cross-country regression of disagreement (152 obs, HC1 SE) ===\n')
etable(m_g, m_t, m_gt, m_gt_int,
       headers = c('GDP', 'Temp', 'Both', 'Both + int'),
       fitstat = c('n', 'ar2'),
       digits = 4, digits.stats = 3)

dm_std <- dm_iso %>% mutate(across(c(disagr, log_gdp, temp),
                                    ~ as.vector(scale(.x))))
m_std <- feols(disagr ~ log_gdp + temp, dm_std, vcov = 'hetero')
cat('\nStandardised (z-scored) coefficients:\n')
print(round(coef(m_std)[c('log_gdp','temp')], 3))

# Country-year cross-section.  Cluster on iso for serial correlation only.
dm_cy <- make_disagr_df(wide_r, mods, c('iso','year'), metric = 'mad') %>%
  left_join(make_country_cov(df, mods, c('iso','year')),
            by = c('iso','year')) %>%
  filter(is.finite(log_gdp), is.finite(temp))
m_cy <- feols(disagr ~ log_gdp + temp, dm_cy, cluster = 'iso')
cat('\n=== Country-year regression of disagreement (3096 obs) ===\n')
etable(m_cy, fitstat = c('n', 'r2'), digits = 4, digits.stats = 3)

# ============================================================================
# Heatmap, marginal effects, actual 2D scatter (interaction model)
# ============================================================================
p_heat <- plot_disagr_heatmap(
  m_gt_int, dm_iso, x = 'log_gdp', y = 'temp', label_var = 'iso',
  x_lab = 'log mean GDP per capita (PPP)',
  y_lab = 'Country mean temperature (°C)',
  title = 'Predicted disagreement: log_gdp + temp + log_gdp×temp')
ggsave('impact_uncertainty/tb/figs/tc_disagr_heatmap.pdf', p_heat,
       width = 9, height = 5.5, units = 'in')

p_me <- plot_disagr_marginal(
  m_gt_int, dm_iso, x = 'log_gdp', y = 'temp',
  x_lab = 'log mean GDP per capita (PPP)',
  y_lab = 'Country mean temperature (°C)',
  title_x = 'Marginal effect of income, by climate',
  title_y = 'Marginal effect of climate, by income')
ggsave('impact_uncertainty/tb/figs/tc_disagr_marginal.pdf', p_me,
       width = 12, height = 4.5, units = 'in')

p_act <- plot_disagr_actual_2d(
  dm_iso, x = 'log_gdp', y = 'temp', label_var = 'iso',
  x_lab = 'log mean GDP per capita (PPP)',
  y_lab = 'Country mean temperature (°C)',
  title = paste0('Actual disagreement (', nrow(dm_iso), ' countries)'))
ggsave('impact_uncertainty/tb/figs/tc_disagr_actual_2d.pdf', p_act,
       width = 9, height = 5.5, units = 'in')

cat('\nAll figures written to figs/.\n')
