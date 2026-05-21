# Out-of-sample predictive score per reanalysis, via FWL residualisation +
# fast .lm.fit on the residuals. Self-contained: mirrors the balanced-filter
# logic from distance.R but does not source it.

library(tidyverse)
library(fixest)
theme_set(theme_classic())
setwd('/Users/tombearpark/Library/CloudStorage/Dropbox/ks')

# 1. Load + balance -----------------------------------------------------------
df <- read_csv('reg_inputs/gdppc/reginputs_gdppc_1deg.csv', show_col_types = FALSE) %>%
  mutate(temp_1 = temp_popweight, temp_2 = temp_1^2,
         prec_1 = precip_popweight, prec_2 = prec_1^2,
         time1 = year - min(year), time2 = time1^2)
mods <- sort(unique(df$obs_mod))
keep <- df %>% filter(!is.na(temp_1)) %>%
  distinct(iso, year, obs_mod) %>% count(iso, year) %>%
  filter(n == length(mods)) %>% select(iso, year)
df <- df %>% semi_join(keep, by = c('iso','year'))
message(sprintf('balanced panel: %d (iso, year) pairs × %d models',
                nrow(keep), length(mods)))

# 2. Residualise per dataset (FWL) --------------------------------------------
# Partial out the FE structure and precipitation controls — the nuisance.
residualise_mod <- function(data, mod_name) {
  d <- data %>% filter(obs_mod == mod_name) %>%
    drop_na(growthWDI, temp_1, temp_2, prec_1, prec_2)
  fml <- function(lhs) as.formula(
    paste0(lhs, ' ~ prec_1 + prec_2 | iso + year + iso[time1] + iso[time2]'))
  tibble(iso = d$iso, year = d$year,
         y_t = resid(feols(fml('growthWDI'), d, notes = FALSE)),
         t1  = resid(feols(fml('temp_1'),    d, notes = FALSE)),
         t2  = resid(feols(fml('temp_2'),    d, notes = FALSE)))
}
R <- set_names(mods, mods) %>% map(~residualise_mod(df, .x))

# 3. CV helpers ---------------------------------------------------------------
splits_random <- function(r, k = 10, seed = 1) {
  set.seed(seed)
  folds <- sample(rep(seq_len(k), length.out = nrow(r)))
  split(seq_len(nrow(r)), folds)
}
splits_by_country <- function(r) split(seq_len(nrow(r)), r$iso)
splits_by_year    <- function(r) split(seq_len(nrow(r)), r$year)

score_mod <- function(r, split_fn, label) {
  folds <- split_fn(r)
  X <- cbind(1, r$t1, r$t2); y <- r$y_t
  mses <- map_dbl(folds, function(test) {
    train <- setdiff(seq_len(nrow(r)), test)
    fit <- .lm.fit(X[train, , drop = FALSE], y[train])
    pred <- drop(X[test, , drop = FALSE] %*% fit$coefficients)
    mean((y[test] - pred)^2)
  })
  tibble(split = label,
         mse = mean(mses),
         mse_se = sd(mses) / sqrt(length(mses)),
         oos_r2 = 1 - mean(mses) / var(y),
         n_folds = length(folds))
}

# 4. Score all five datasets under all three splits --------------------------
scores <- map_dfr(mods, function(m) {
  r <- R[[m]]
  bind_rows(
    score_mod(r, function(x) splits_random(x, k = 10), 'random-10fold'),
    score_mod(r, splits_by_country,                    'leave-country-out'),
    score_mod(r, splits_by_year,                       'leave-year-out')
  ) %>% mutate(model = m, .before = 1)
})
print(scores, n = 100)

# 5. Plot ---------------------------------------------------------------------
fig_dir <- 'impact_uncertainty/tb/figs'
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
save_fig <- function(p, name, w = 8, h = 4)
  ggsave(file.path(fig_dir, name), p, width = w, height = h, units = 'in')

scores <- scores %>%
  mutate(split = factor(split, levels = c('random-10fold',
                                          'leave-country-out',
                                          'leave-year-out')))

p_mse <- ggplot(scores, aes(x = reorder(model, -mse), y = mse)) +
  geom_col(fill = 'steelblue') +
  geom_errorbar(aes(ymin = mse - mse_se, ymax = mse + mse_se), width = 0.3) +
  facet_wrap(~split, scales = 'free_x', nrow = 1) +
  coord_flip() +
  labs(x = NULL, y = 'CV-MSE (lower = better)')
save_fig(p_mse, 'oos_mse.pdf', w = 10, h = 3.5)

p_r2 <- ggplot(scores, aes(x = reorder(model, oos_r2), y = oos_r2)) +
  geom_col(fill = 'steelblue') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  facet_wrap(~split, nrow = 1) +
  coord_flip() +
  labs(x = NULL, y = 'OOS R² (higher = better)')
save_fig(p_r2, 'oos_r2.pdf', w = 10, h = 3.5)

cat('\nWrote oos_mse.pdf and oos_r2.pdf to', fig_dir, '\n')
