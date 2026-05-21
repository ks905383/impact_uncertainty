# Mortality version of gdp_disagree.R.
#
# Stage A: replicate the basic mortality regression of regressions_mortality.do
# (block DO_BASIC), per reanalysis, plot the response function as a sanity check.
# Stage B: disagreement-vs-covariate pipeline mirroring gdp_disagree.R, sharing
# helpers in disagree_helpers.R.

library(tidyverse)
library(fixest)
library(patchwork)
theme_set(theme_classic())
setwd('/Users/tombearpark/Library/CloudStorage/Dropbox/ks')
source('impact_uncertainty/tb/code/disagree_helpers.R')

# ============================================================================
# 1. Load + minimal prep
# ============================================================================
df <- read_csv('reg_inputs/mortality/reginputs_mortality_1deg.csv',
               show_col_types = FALSE)
mods <- sort(unique(df$obs_mod))

# ============================================================================
# 2. Prep: replicate prep_data_multclims.do
# ============================================================================
prep_for_reg <- function(df) {
  df %>%
    filter(year <= 2010,
           agegroup != 'total') %>%        # = `drop if agegroup == "0"` in .do
    group_by(iso, agegroup) %>%
    mutate(deathrate_p99 = quantile(deathrate, 0.99, na.rm = TRUE),
           deathrate_w99 = pmin(deathrate, deathrate_p99)) %>%
    ungroup() %>%
    drop_na(deathrate_w99, tavg_poly_1, prcp_poly_1,
            loggdppc_adm1_avg, lr_tavg_1981_2010) %>%
    mutate(adm0_code   = as.integer(factor(iso)),
           adm1_code   = as.integer(factor(paste(iso, adm1_id))),
           adm2_id_eff = if_else(iso == 'JPN', adm1_id, adm2_id),
           adm2_code   = as.integer(factor(paste(iso, adm1_id, adm2_id_eff))),
           agegroup_   = as.integer(factor(agegroup)),
           CHN_ts      = if_else(iso == 'CHN' & year >= 2004, 2L, 1L))
}

df_prep <- prep_for_reg(df)
cat(sprintf('After prep: %s rows, %d adm0, %d adm1, %d adm2, %d agegroups\n',
            format(nrow(df_prep), big.mark = ','),
            n_distinct(df_prep$adm0_code),
            n_distinct(df_prep$adm1_code),
            n_distinct(df_prep$adm2_code),
            n_distinct(df_prep$agegroup_)))

# ============================================================================
# 3. Fit pooled mortality regression per obs_mod
# ============================================================================
# Stata:
# reghdfe deathrate_w99 tavg_poly_1..4 [pw=weight],
#   absorb(adm0_code##c.prcp_poly_1 adm0_code##c.prcp_poly_2
#          i.adm2_code#i.CHN_ts#i.agegroup_  i.adm0_code#i.year#i.agegroup_)
#   cluster(adm1_code)

fit_mortality_pooled <- function(d_one_mod) {
  feols(deathrate_w99 ~ tavg_poly_1 + tavg_poly_2 + tavg_poly_3 + tavg_poly_4 |
          adm0_code[prcp_poly_1] + adm0_code[prcp_poly_2] +
          adm2_code^CHN_ts^agegroup_ + adm0_code^year^agegroup_,
        data    = d_one_mod,
        weights = ~population,
        cluster = ~adm1_code,
        notes   = FALSE)
}

cat('\nFitting per obs_mod...\n')
mortality_fits <- map(mods, function(m) {
  cat('  ', m, '...')
  t0 <- Sys.time()
  fit <- fit_mortality_pooled(df_prep %>% filter(obs_mod == m))
  cat(sprintf(' %.1fs (N = %s)\n',
              as.numeric(difftime(Sys.time(), t0, units = 'secs')),
              format(fit$nobs, big.mark = ',')))
  fit
}) %>% set_names(mods)

# Print coefficients
cat('\n--- tavg_poly_* coefficients per reanalysis ---\n')
etable(mortality_fits)

# ============================================================================
# 4. Plot response function relative to T = 20°C
# ============================================================================
T_grid <- seq(0, 35, by = 0.5); T_ref <- 20
resp_df <- imap_dfr(mortality_fits, function(fit, m) {
  b <- coef(fit)[paste0('tavg_poly_', 1:4)]
  resp <- function(T) b[1]*T + b[2]*T^2 + b[3]*T^3 + b[4]*T^4
  tibble(obs_mod = m, T = T_grid, dy = resp(T_grid) - resp(T_ref))
})

p <- ggplot(resp_df, aes(T, dy, color = obs_mod)) +
  geom_line(linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40') +
  geom_vline(xintercept = T_ref, linetype = 'dotted', color = 'grey60') +
  labs(x = 'Temperature (°C)',
       y = expression(Delta~deathrate ~ vs ~ T == 20*degree*C),
       color = NULL,
       title = 'Mortality response function, by reanalysis (sanity check)') +
  theme(legend.position = 'bottom')

ggsave('impact_uncertainty/tb/figs/mt_response_check.pdf', p,
       width = 9, height = 4.5, units = 'in')
cat('\nWrote mt_response_check.pdf\n')

# ============================================================================
# STAGE B — DISAGREEMENT PIPELINE
# ============================================================================
# Reduce to one age group: temperature is the same across age within
# (adm2, year) so any agegroup gives the same temp series. Add time1, time2.
df_temp <- df_prep %>%
  filter(agegroup == '0-4') %>%
  mutate(time1 = year - min(year),
         time2 = time1^2)

# Balance: keep only (adm2, year) observed in all four reanalyses.
keep_b <- df_temp %>%
  distinct(adm2_code, year, obs_mod) %>% count(adm2_code, year) %>%
  filter(n == length(mods)) %>% select(adm2_code, year)
df_temp <- df_temp %>% semi_join(keep_b, by = c('adm2_code','year'))
cat(sprintf('\nBalanced disagreement panel: %s adm2-years × %d datasets\n',
            format(nrow(keep_b), big.mark = ','), length(mods)))

# Residualise tavg on the FE structure analogous to gdp_disagree.R:
# adm2 + year + adm0[time1] + adm0[time2].
resid_df <- map_dfr(mods, function(m) {
  d <- df_temp %>% filter(obs_mod == m)
  r <- resid(feols(tavg_poly_1 ~ 1 |
                     adm2_code + year + adm0_code[time1] + adm0_code[time2],
                   d, notes = FALSE))
  tibble(adm2_code = d$adm2_code, adm1_code = d$adm1_code,
         adm0_code = d$adm0_code, iso = d$iso, year = d$year,
         obs_mod = m, resid = r)
})
wide_r <- resid_df %>% pivot_wider(names_from = obs_mod, values_from = resid)

# ---------------------------------------------------------------------------
# Aggregate covariates (income, realised temp, iso passthrough) to grouping.
# ---------------------------------------------------------------------------
make_adm_cov <- function(df, mods, group_cols) {
  df %>% filter(obs_mod == mods[1]) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(log_gdp = mean(loggdppc_adm1_avg, na.rm = TRUE),
              temp    = mean(tavg_poly_1,        na.rm = TRUE),
              n_years = n_distinct(year),
              iso     = if (any(c('adm1_code','adm2_code') %in% group_cols))
                        first(iso) else NA_character_,
              .groups = 'drop')
}

# ---------------------------------------------------------------------------
# Driver: build dm at a given granularity, save the figure, return dm.
# ---------------------------------------------------------------------------
build_and_save <- function(group_cols, x, label_var, file, x_lab,
                           min_years = 0) {
  dm <- make_disagr_df(wide_r, mods, group_cols, metric = 'mad') %>%
    left_join(make_adm_cov(df_temp, mods, group_cols), by = group_cols) %>%
    filter(is.finite(.data[[x]]))
  if (min_years > 0 && 'n_years' %in% names(dm))
    dm <- dm %>% filter(n_years >= min_years)
  p <- plot_disagr(dm, x = x, label_var = label_var, x_lab = x_lab,
                   colour_var = 'iso')
  ggsave(file.path('impact_uncertainty/tb/figs', file),
         p, width = 8, height = 5, units = 'in')
  r2 <- summary(lm(reformulate(x, 'disagr'), dm))$r.squared
  cat(sprintf('  %-32s  N = %6d  R² = %.3f\n', file, nrow(dm), r2))
  p
}

# ---------------------------------------------------------------------------
# Build the six scatter figures (mt_disagr_*).
# ---------------------------------------------------------------------------
cat('\nBuilding disagreement scatters:\n')

# Adm1 level (analog of country level — ~152 adm1)
build_and_save(c('adm1_code'),       'log_gdp', 'iso',
               'mt_disagr_vs_income.pdf',  'log mean GDP per capita', min_years = 10)
build_and_save(c('adm1_code'),       'temp',    'iso',
               'mt_disagr_vs_temp.pdf',    'Adm1 mean temperature (°C)')

# Year level
build_and_save(c('year'),            'log_gdp', 'year',
               'mt_disagr_yr_income.pdf',  'log mean GDP per capita (annual avg)')
build_and_save(c('year'),            'temp',    'year',
               'mt_disagr_yr_climate.pdf', 'Sample-mean temperature in year (°C)')

# Adm2-year level (densest)
build_and_save(c('adm2_code','year'), 'log_gdp', NULL,
               'mt_disagr_cy_income.pdf',  'log GDP per capita')
build_and_save(c('adm2_code','year'), 'temp',    NULL,
               'mt_disagr_cy_climate.pdf', 'Realised temperature (°C)')

# ============================================================================
# Regressions of disagreement on income + climate, with interaction
# ============================================================================
min_years <- 10

dm_adm1 <- make_disagr_df(wide_r, mods, c('adm1_code'), metric = 'mad') %>%
  left_join(make_adm_cov(df_temp, mods, c('adm1_code')), by = 'adm1_code') %>%
  filter(is.finite(log_gdp), is.finite(temp)) %>%
  filter(n_years >= !!min_years)

m_g  <- feols(disagr ~ log_gdp,        dm_adm1, vcov = 'hetero')
m_t  <- feols(disagr ~ temp,           dm_adm1, vcov = 'hetero')
m_gt <- feols(disagr ~ log_gdp + temp, dm_adm1, vcov = 'hetero')
m_gt_int <- feols(disagr ~ log_gdp + temp + log_gdp * temp,
                  dm_adm1, vcov = 'hetero')

cat(sprintf('\n=== Cross-adm1 regression of disagreement (%d obs, HC1 SE) ===\n',
            nrow(dm_adm1)))
etable(m_g, m_t, m_gt, m_gt_int,
       headers = c('GDP', 'Temp', 'Both', 'Both + int'),
       fitstat = c('n', 'ar2'),
       digits = 4, digits.stats = 3)

# Standardised version
dm_std <- dm_adm1 %>% mutate(across(c(disagr, log_gdp, temp),
                                    ~ as.vector(scale(.x))))
m_std <- feols(disagr ~ log_gdp + temp, dm_std, vcov = 'hetero')
cat('\nStandardised (z-scored) coefficients:\n')
print(round(coef(m_std)[c('log_gdp','temp')], 3))

# ============================================================================
# Heatmap, marginal effects, actual 2D scatter (interaction model)
# ============================================================================
p_heat <- plot_disagr_heatmap(
  m_gt_int, dm_adm1, x = 'log_gdp', y = 'temp', label_var = 'iso',
  x_lab = 'log mean GDP per capita',
  y_lab = 'Adm1 mean temperature (°C)',
  title = 'Predicted disagreement (mortality, adm1 level)')
ggsave('impact_uncertainty/tb/figs/mt_disagr_heatmap.pdf', p_heat,
       width = 9, height = 5.5, units = 'in')

p_me <- plot_disagr_marginal(
  m_gt_int, dm_adm1, x = 'log_gdp', y = 'temp',
  x_lab = 'log mean GDP per capita',
  y_lab = 'Adm1 mean temperature (°C)',
  title_x = 'Marginal effect of income, by climate',
  title_y = 'Marginal effect of climate, by income')
ggsave('impact_uncertainty/tb/figs/mt_disagr_marginal.pdf', p_me,
       width = 12, height = 4.5, units = 'in')

p_act <- plot_disagr_actual_2d(
  dm_adm1, x = 'log_gdp', y = 'temp', label_var = 'iso',
  x_lab = 'log mean GDP per capita',
  y_lab = 'Adm1 mean temperature (°C)',
  title = paste0('Actual disagreement (', nrow(dm_adm1), ' adm1 units)'))
ggsave('impact_uncertainty/tb/figs/mt_disagr_actual_2d.pdf', p_act,
       width = 9, height = 5.5, units = 'in')

cat('\nAll mortality-disagreement figures written to figs/.\n')
