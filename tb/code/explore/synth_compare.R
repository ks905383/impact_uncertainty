# Synthetic temperature measures: build five constructions and compare each
# against the five raw datasets on ME precision, OOS predictive fit, and
# robustness to leave-one-dataset-out perturbation.
#
# The panel regression throughout uses ensemble-mean precipitation as the
# control, so that temperature is the only thing varying across methods.

library(tidyverse)
library(fixest)
theme_set(theme_classic())
setwd('/Users/tombearpark/Library/CloudStorage/Dropbox/ks')

fig_dir <- 'impact_uncertainty/tb/figs'
save_fig <- function(p, name, w = 8, h = 5)
  ggsave(file.path(fig_dir, name), p, width = w, height = h, units = 'in')

# ============================================================================
# 1. Load + balance + wide format --------------------------------------------
# ============================================================================
df <- read_csv('reg_inputs/gdppc/reginputs_gdppc_1deg.csv', show_col_types = FALSE) %>%
  mutate(temp_1 = temp_popweight,
         prec_1 = precip_popweight,
         time1  = year - min(year),
         time2  = time1^2)
mods <- sort(unique(df$obs_mod))
keep <- df %>% filter(!is.na(temp_1)) %>%
  distinct(iso, year, obs_mod) %>% count(iso, year) %>%
  filter(n == length(mods)) %>% select(iso, year)
df <- df %>% semi_join(keep, by = c('iso','year'))

wide <- df %>%
  select(iso, year, continent, time1, time2, growthWDI, obs_mod, temp_1, prec_1) %>%
  pivot_wider(names_from = obs_mod, values_from = c(temp_1, prec_1))

T_cols <- paste0('temp_1_', mods)
P_cols <- paste0('prec_1_', mods)
T_mat  <- as.matrix(wide[, T_cols]); colnames(T_mat) <- mods
P_mat  <- as.matrix(wide[, P_cols]); colnames(P_mat) <- mods
wide$P_mean    <- rowMeans(P_mat)
wide$P_mean_sq <- wide$P_mean^2

message(sprintf('Balanced panel: %d rows × %d datasets', nrow(wide), length(mods)))

# ============================================================================
# 2. Residualise each raw T on FE structure (shared by PC / Factor diagnostics)
# ============================================================================
nuisance_FE <- ' | iso + year + iso[time1] + iso[time2]'
residualise_T <- function(m) {
  fml <- as.formula(paste0('`temp_1_', m, '` ~ 1', nuisance_FE))
  resid(feols(fml, wide, notes = FALSE))
}
R_mat <- do.call(cbind, map(mods, residualise_T))
colnames(R_mat) <- mods
fit_mat <- T_mat - R_mat   # = fitted FE structure per dataset (country + year + trend)

# ============================================================================
# 3. Build five synthetic measures -------------------------------------------
# ============================================================================

# 3.1 MEAN ----------------------------------------------------------------------
T_MEAN <- rowMeans(T_mat)

# 3.2 MEDIAN (robust variant of MEAN) ------------------------------------------
T_MEDIAN <- apply(T_mat, 1, median)

# 3.3 PC1 (pooled): first PC of residualised temps across all (iso, year) ------
pc_pooled <- prcomp(R_mat, center = TRUE, scale. = FALSE)
F_pooled  <- pc_pooled$x[, 1]
# Rescale to match sd of ensemble-mean residual, and fix sign to align with it.
F_scaled  <- F_pooled / sd(F_pooled) * sd(rowMeans(R_mat))
if (cor(F_scaled, rowMeans(R_mat)) < 0) F_scaled <- -F_scaled
T_PC1 <- rowMeans(fit_mat) + F_scaled   # restore climate + trend in °C units

cat(sprintf('\nPC1 explains %.1f%% of pooled residual variance\n',
            100 * pc_pooled$sdev[1]^2 / sum(pc_pooled$sdev^2)))

# 3.4 FACTOR: per-country factor model via EM ----------------------------------
# Model: R[t, k] = lam_k * F(t) + eps_k,  eps_k ~ N(0, sig2_k)
# Identification: sd(F) = 1; sign chosen to match sign(mean(lam)).
factor_fit_country <- function(R_c, max_iter = 60, tol = 1e-7) {
  K <- ncol(R_c); Tn <- nrow(R_c)
  pc  <- prcomp(R_c, center = FALSE)
  F_h <- pc$x[, 1] / pc$sdev[1]              # unit sd
  lam <- as.vector(crossprod(R_c, F_h)) / Tn
  if (mean(lam) < 0) { lam <- -lam; F_h <- -F_h }
  sig2 <- pmax(colMeans((R_c - F_h %o% lam)^2), 1e-10)
  prev <- -Inf
  for (it in 1:max_iter) {
    wk    <- lam / sig2
    F_h   <- as.vector(R_c %*% wk) / sum(lam^2 / sig2)
    s     <- sd(F_h); if (s > 0) { F_h <- F_h / s; lam <- lam * s }
    lam   <- as.vector(crossprod(R_c, F_h)) / Tn
    sig2  <- pmax(colMeans((R_c - F_h %o% lam)^2), 1e-10)
    ll    <- -sum(log(sig2)) -
             sum((R_c - F_h %o% lam)^2 / matrix(sig2, Tn, K, byrow = TRUE))
    if (abs(ll - prev) < tol) break
    prev <- ll
  }
  list(F = F_h, lam = lam, sig2 = sig2)
}

T_FACTOR <- numeric(nrow(wide))
factor_info <- list()
for (iso_c in unique(wide$iso)) {
  idx <- which(wide$iso == iso_c)
  fit_c <- factor_fit_country(R_mat[idx, , drop = FALSE])
  factor_info[[iso_c]] <- fit_c
  # Back to °C: common-factor path scaled by mean lambda, added to fitted trend.
  T_FACTOR[idx] <- rowMeans(fit_mat[idx, , drop = FALSE]) +
                   fit_c$F * mean(fit_c$lam)
}

nuis_pfml <- ' ~ P_mean + P_mean_sq | iso + year + iso[time1] + iso[time2]'
y_r <- resid(feols(as.formula(paste0('growthWDI', nuis_pfml)), wide, notes = FALSE))

# 3.5 OOS-OPT attempt: characterise the objective surface --------------------
# We *wanted* weights minimising CV-MSE of the growth regression. But MSE is
# essentially flat over the simplex here (we already knew from oos_score.R that
# the five raw datasets are indistinguishable on OOS MSE). Quantify the flatness
# by random Dirichlet sampling over the 5-simplex.
oos_mse_at_w <- function(w, folds) {
  Tw <- as.vector(T_mat %*% w)
  wloc <- wide; wloc$Tw <- Tw; wloc$Tw_sq <- Tw^2
  Tr  <- resid(feols(as.formula(paste0('Tw',    nuis_pfml)), wloc, notes = FALSE))
  Tsr <- resid(feols(as.formula(paste0('Tw_sq', nuis_pfml)), wloc, notes = FALSE))
  X <- cbind(Tr, Tsr)
  mean(map_dbl(folds, function(test) {
    train <- setdiff(seq_along(y_r), test)
    fit <- .lm.fit(X[train, , drop = FALSE], y_r[train])
    if (any(is.na(fit$coefficients))) return(NA_real_)
    pred <- drop(X[test, , drop = FALSE] %*% fit$coefficients)
    mean((y_r[test] - pred)^2)
  }), na.rm = TRUE)
}
set.seed(1)
folds10 <- split(seq_along(y_r), sample(rep(1:10, length.out = length(y_r))))

set.seed(2)
W_rand <- matrix(rexp(500 * length(mods)), ncol = length(mods))
W_rand <- W_rand / rowSums(W_rand)   # Dirichlet(1,...,1) samples on 5-simplex
mse_rand <- apply(W_rand, 1, oos_mse_at_w, folds = folds10)
cat(sprintf(
  '\nOOS-OPT sensitivity: across 500 random weight vectors, MSE range = [%.6f, %.6f]\n',
  min(mse_rand), max(mse_rand)))
cat(sprintf('  relative spread = %.3f%% of mean\n',
            100 * (max(mse_rand) - min(mse_rand)) / mean(mse_rand)))
cat('  => objective is flat => OOS-OPT is ill-identified.\n')
cat('  We will NOT report OOS-OPT as a standalone method.  Dropped.\n')

# 3.5′ MIN-VAR: inverse-variance weighted combination ------------------------
# Factor model already gave us sigma_k(c). Under the model
#   T_k(c,t) = T*(c,t) + eta_k(c,t), eta_k ⊥ eta_j, Var(eta_k(c,.)) = sigma_k^2(c),
# the efficient (min-variance, unbiased-for-T*) weighted mean is
#   w_k(c) = (1/sigma_k^2(c)) / sum_k' 1/sigma_k'^2(c).
# This differs from ANCHORED in that the noise variance is estimated from the
# factor model, not as RMSE-against-UDel.
T_MINVAR <- numeric(nrow(wide))
for (iso_c in unique(wide$iso)) {
  ii   <- which(wide$iso == iso_c)
  sig2 <- pmax(factor_info[[iso_c]]$sig2, 1e-10)
  w_c  <- (1 / sig2) / sum(1 / sig2)
  T_MINVAR[ii] <- as.vector(T_mat[ii, , drop = FALSE] %*% w_c)
}

# 3.6 ANCHORED: UDel as station-proxy, inverse-MSE weights per country --------
reanalyses <- setdiff(mods, 'UDel')
T_ANCHORED <- numeric(nrow(wide))
for (iso_c in unique(wide$iso)) {
  idx <- which(wide$iso == iso_c)
  t_u <- wide$`temp_1_UDel`[idx]
  mses_c <- sapply(reanalyses, function(k)
    mean((wide[[paste0('temp_1_', k)]][idx] - t_u)^2))
  w_c <- (1 / mses_c) / sum(1 / mses_c)
  T_combo <- sapply(reanalyses, function(k) wide[[paste0('temp_1_', k)]][idx]) %*% w_c
  T_ANCHORED[idx] <- as.vector(T_combo)
}

# ============================================================================
# 4. Evaluate every method on ME + CI + OOS ---------------------------------
# ============================================================================
method_temps <- list(
  `ERA5-025` = T_mat[, 'ERA5-025'],
  GMFD       = T_mat[, 'GMFD'],
  `JRA-3Q`   = T_mat[, 'JRA-3Q'],
  MERRA2     = T_mat[, 'MERRA2'],
  UDel       = T_mat[, 'UDel'],
  MEAN       = T_MEAN,
  MEDIAN     = T_MEDIAN,
  PC1        = T_PC1,
  FACTOR     = T_FACTOR,
  `MIN-VAR`  = T_MINVAR,
  ANCHORED   = T_ANCHORED
)

# 4.1 Panel regression and ME(T) + delta-method SE -----------------------------
evaluate_method <- function(T_vec, method_name) {
  dat <- wide; dat$T_ <- T_vec; dat$T_sq <- T_vec^2
  fit <- feols(growthWDI ~ T_ + T_sq + P_mean + P_mean_sq |
                 iso + year + iso[time1] + iso[time2],
               dat, cluster = 'iso', notes = FALSE)
  b  <- coef(fit)[c('T_','T_sq')]
  V  <- vcov(fit)[c('T_','T_sq'), c('T_','T_sq')]
  T_grid <- c(5, 14, 20, 25, 28, 30)
  me     <- b[1] + 2 * b[2] * T_grid
  var_me <- V[1,1] + 4 * T_grid^2 * V[2,2] + 4 * T_grid * V[1,2]
  tibble(method = method_name, T_eval = T_grid,
         me = as.numeric(me), se_me = sqrt(as.numeric(var_me)),
         beta1 = b[1], beta2 = b[2], se_beta2 = sqrt(V[2,2]))
}
me_results <- map2_dfr(method_temps, names(method_temps), evaluate_method)

me_at_25 <- me_results %>% filter(T_eval == 25) %>%
  mutate(lo = me - 1.96*se_me, hi = me + 1.96*se_me, ci_width = hi - lo,
         t_beta2 = beta2 / se_beta2,
         kind = if_else(method %in% mods, 'raw', 'synthetic')) %>%
  arrange(ci_width)
cat('\nME at 25°C, sorted by CI width (also showing t-statistic on β₂):\n')
print(me_at_25 %>% select(method, beta2, se_beta2, t_beta2, me, ci_width, kind), n = Inf)

cat('\n(Same) sorted by |t-statistic on β₂| — the econometric precision measure:\n')
print(me_at_25 %>% arrange(desc(abs(t_beta2))) %>%
        select(method, beta2, se_beta2, t_beta2, me, ci_width, kind), n = Inf)

# 4.2 OOS MSE on three CV splits ----------------------------------------------
oos_method <- function(T_vec, method_name) {
  dat <- wide; dat$T_ <- T_vec; dat$T_sq <- T_vec^2
  Tr  <- resid(feols(as.formula(paste0('T_',    nuis_pfml)), dat, notes = FALSE))
  Tsr <- resid(feols(as.formula(paste0('T_sq', nuis_pfml)), dat, notes = FALSE))
  X <- cbind(Tr, Tsr)
  score <- function(folds) {
    mses <- map_dbl(folds, function(test) {
      train <- setdiff(seq_along(y_r), test)
      fit <- .lm.fit(X[train, , drop = FALSE], y_r[train])
      if (any(is.na(fit$coefficients))) return(NA_real_)
      pred <- drop(X[test, , drop = FALSE] %*% fit$coefficients)
      mean((y_r[test] - pred)^2)
    })
    ok <- !is.na(mses)
    c(mean(mses[ok]), sd(mses[ok]) / sqrt(sum(ok)))
  }
  set.seed(42)
  f_r <- split(seq_along(y_r), sample(rep(1:10, length.out = length(y_r))))
  f_c <- split(seq_along(y_r), wide$iso)
  f_y <- split(seq_along(y_r), wide$year)
  mr <- score(f_r); mc <- score(f_c); my <- score(f_y)
  tibble(method = method_name,
         mse_random  = mr[1], se_random  = mr[2],
         mse_country = mc[1], se_country = mc[2],
         mse_year    = my[1], se_year    = my[2])
}
oos_results <- map2_dfr(method_temps, names(method_temps), oos_method)
oos_results <- oos_results %>%
  mutate(kind = if_else(method %in% mods, 'raw', 'synthetic'))
cat('\nOOS MSE by split (sorted by random-10-fold MSE), 6 sig figs:\n')
oos_results %>% arrange(mse_random) %>%
  mutate(across(where(is.numeric), ~signif(.x, 6))) %>%
  print(n = Inf)

# ============================================================================
# 5. Jackknife: rebuild each synthetic leaving out one dataset ---------------
# ============================================================================
jackknife_synth <- function(method_name, leave_out) {
  idx  <- setdiff(mods, leave_out)
  Ts   <- T_mat[, idx, drop = FALSE]
  Rs   <- R_mat[, idx, drop = FALSE]
  Fs   <- fit_mat[, idx, drop = FALSE]
  switch(method_name,
    MEAN   = rowMeans(Ts),
    MEDIAN = apply(Ts, 1, median),
    PC1    = {
      pc <- prcomp(Rs, center = TRUE, scale. = FALSE)
      Fp <- pc$x[,1] / sd(pc$x[,1]) * sd(rowMeans(Rs))
      if (cor(Fp, rowMeans(Rs)) < 0) Fp <- -Fp
      rowMeans(Fs) + Fp
    },
    FACTOR = {
      out <- numeric(nrow(wide))
      for (iso_c in unique(wide$iso)) {
        ii <- which(wide$iso == iso_c)
        fc <- factor_fit_country(Rs[ii, , drop = FALSE])
        out[ii] <- rowMeans(Fs[ii, , drop = FALSE]) + fc$F * mean(fc$lam)
      }
      out
    },
    `MIN-VAR` = {
      out <- numeric(nrow(wide))
      for (iso_c in unique(wide$iso)) {
        ii <- which(wide$iso == iso_c)
        fc <- factor_fit_country(Rs[ii, , drop = FALSE])
        sig2 <- pmax(fc$sig2, 1e-10)
        w_c <- (1 / sig2) / sum(1 / sig2)
        out[ii] <- as.vector(Ts[ii, , drop = FALSE] %*% w_c)
      }
      out
    },
    ANCHORED = {
      if (leave_out == 'UDel') return(rowMeans(Ts))
      reans <- setdiff(idx, 'UDel')
      out <- numeric(nrow(wide))
      for (iso_c in unique(wide$iso)) {
        ii <- which(wide$iso == iso_c)
        t_u <- wide$`temp_1_UDel`[ii]
        ms <- sapply(reans, function(k)
          mean((wide[[paste0('temp_1_', k)]][ii] - t_u)^2))
        w <- (1 / ms) / sum(1 / ms)
        tc <- sapply(reans, function(k) wide[[paste0('temp_1_', k)]][ii]) %*% w
        out[ii] <- as.vector(tc)
      }
      out
    })
}

synth_names <- c('MEAN', 'MEDIAN', 'PC1', 'FACTOR', 'MIN-VAR', 'ANCHORED')
jk_df <- map_dfr(synth_names, function(meth) {
  map_dfr(mods, function(lo) {
    T_vec <- jackknife_synth(meth, lo)
    evaluate_method(T_vec, paste0(meth, '_no_', lo)) %>% filter(T_eval == 25) %>%
      mutate(synthetic = meth, leave_out = lo)
  })
})
jk_summary <- jk_df %>%
  group_by(synthetic) %>%
  summarise(me_mean = mean(me), me_sd = sd(me),
            me_range = max(me) - min(me))
cat('\nJackknife stability of ME(25°C):\n')
print(jk_summary)

# ============================================================================
# 6. Factor-model diagnostics: measurement-error variance per dataset --------
# ============================================================================
sigma_by_dataset <- tibble(
  iso = names(factor_info),
  !!!set_names(
    map(mods, function(k) map_dbl(factor_info, ~sqrt(.x$sig2[match(k, mods)]))),
    paste0('sigma_', mods))
) %>%
  pivot_longer(-iso, names_to = 'dataset', values_to = 'sigma') %>%
  mutate(dataset = str_remove(dataset, 'sigma_'))

cat('\nFactor-model sigma_k summary (mean across countries, °C):\n')
sigma_summary <- sigma_by_dataset %>%
  group_by(dataset) %>%
  summarise(sigma_mean = mean(sigma), sigma_sd = sd(sigma))
print(sigma_summary)

# ============================================================================
# 7. Plots -------------------------------------------------------------------
# ============================================================================
# 7.1 ME curves: raw vs synthetic
me_curve_df <- map2_dfr(method_temps, names(method_temps), function(T_vec, m) {
  dat <- wide; dat$T_ <- T_vec; dat$T_sq <- T_vec^2
  fit <- feols(growthWDI ~ T_ + T_sq + P_mean + P_mean_sq |
                 iso + year + iso[time1] + iso[time2],
               dat, cluster = 'iso', notes = FALSE)
  b <- coef(fit)[c('T_','T_sq')]; V <- vcov(fit)[c('T_','T_sq'),c('T_','T_sq')]
  Tg <- seq(0, 32, by = 0.5)
  me <- b[1] + 2 * b[2] * Tg
  se <- sqrt(V[1,1] + 4 * Tg^2 * V[2,2] + 4 * Tg * V[1,2])
  tibble(method = m, temp = Tg, me = me, lo = me - 1.96*se, hi = me + 1.96*se,
         kind = if_else(m %in% mods, 'raw', 'synthetic'))
})
p_me <- ggplot(me_curve_df, aes(temp, me, color = method, fill = method)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40') +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.1, color = NA) +
  geom_line() +
  facet_wrap(~kind, nrow = 1) +
  labs(x = 'Temperature (°C)', y = 'Marginal effect on growth',
       color = NULL, fill = NULL)
save_fig(p_me, 'sc_me_curves.pdf', w = 12, h = 4.5)

# 7.2 CI width at 25°C
p_ci <- ggplot(me_at_25, aes(x = reorder(method, -ci_width), y = ci_width, fill = kind)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c(raw = 'grey60', synthetic = 'steelblue')) +
  labs(x = NULL, y = '95% CI width of ME at 25°C', fill = NULL)
save_fig(p_ci, 'sc_ci_width.pdf', w = 7, h = 5)

# 7.2b Coefficient plot for beta_2 with ±1.96·SE bars
p_b2 <- me_at_25 %>%
  mutate(b2_lo = beta2 - 1.96*se_beta2, b2_hi = beta2 + 1.96*se_beta2) %>%
  ggplot(aes(x = reorder(method, -beta2), y = beta2, color = kind)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = b2_lo, ymax = b2_hi), width = 0.25) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40') +
  scale_color_manual(values = c(raw = 'grey40', synthetic = 'steelblue')) +
  coord_flip() +
  labs(x = NULL, y = expression(hat(beta)[2]~" with 95% CI"), color = NULL)
save_fig(p_b2, 'sc_coef_beta2.pdf', w = 8, h = 5)

# 7.2c Coefficient plot for ME(25°C) with 95% CI
p_me25 <- me_at_25 %>%
  ggplot(aes(x = reorder(method, -me), y = me, color = kind)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40') +
  scale_color_manual(values = c(raw = 'grey40', synthetic = 'steelblue')) +
  coord_flip() +
  labs(x = NULL, y = 'ME at 25°C with 95% CI', color = NULL)
save_fig(p_me25, 'sc_coef_me25.pdf', w = 8, h = 5)

# 7.3 OOS MSE by split
oos_long <- oos_results %>%
  pivot_longer(-c(method, kind), names_to = c('.value', 'split'),
               names_pattern = '(mse|se)_(.+)') %>%
  mutate(split = factor(split, levels = c('random','country','year')))
p_oos <- ggplot(oos_long, aes(x = reorder(method, -mse), y = mse, fill = kind)) +
  geom_col() +
  geom_errorbar(aes(ymin = mse - se, ymax = mse + se), width = 0.3) +
  facet_wrap(~split, scales = 'free_x', nrow = 1) +
  coord_flip() +
  scale_fill_manual(values = c(raw = 'grey60', synthetic = 'steelblue')) +
  labs(x = NULL, y = 'CV-MSE', fill = NULL)
save_fig(p_oos, 'sc_oos_mse.pdf', w = 12, h = 5)

# 7.4 Jackknife stability
p_jk <- ggplot(jk_df, aes(x = leave_out, y = me)) +
  geom_point(size = 2, color = 'steelblue') +
  geom_hline(aes(yintercept = me_mean),
             data = jk_df %>% group_by(synthetic) %>% summarise(me_mean = mean(me)),
             linetype = 'dashed') +
  facet_wrap(~synthetic, nrow = 1) +
  labs(x = 'Dataset left out', y = 'ME at 25°C under jackknife') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_jk, 'sc_jackknife.pdf', w = 12, h = 4)

# 7.5 Factor-model sigma_k by dataset (boxplot across countries)
p_sig <- ggplot(sigma_by_dataset, aes(reorder(dataset, sigma, median), sigma)) +
  geom_boxplot(fill = 'steelblue', alpha = 0.3) +
  labs(x = NULL, y = 'sigma_k (°C) — factor-model measurement-error SD')
save_fig(p_sig, 'sc_factor_sigma.pdf', w = 7, h = 4.5)

# 7.6 Spaghetti plots of raw + synthetic series for top-disagreement countries
# ----------------------------------------------------------------------------
# Rank countries by ensemble-level residualised MAD (average over 10 pairs).
pairs_idx <- combn(length(mods), 2)
country_ensemble_mad <- tibble(
  iso = wide$iso,
  avg_mad = apply(pairs_idx, 2, function(p)
      abs(R_mat[, p[1]] - R_mat[, p[2]])) %>% rowMeans()
) %>% group_by(iso) %>% summarise(avg_mad = mean(avg_mad), .groups = 'drop') %>%
  arrange(desc(avg_mad))

top8 <- country_ensemble_mad %>% head(8) %>% pull(iso)

# Residualise each synthetic on the same FE structure (for the resid plot).
residualise_synth <- function(T_vec) {
  wloc <- wide; wloc$.T <- T_vec
  resid(feols(.T ~ 1 | iso + year + iso[time1] + iso[time2],
              wloc, notes = FALSE))
}
synth_names_plot <- c('MEAN','MEDIAN','PC1','FACTOR','MIN-VAR','ANCHORED')
synth_vecs       <- list(T_MEAN, T_MEDIAN, T_PC1, T_FACTOR, T_MINVAR, T_ANCHORED)
R_synth <- map(synth_vecs, residualise_synth) %>% set_names(synth_names_plot)

# Long-format data frames
raw_long <- map_dfr(mods, function(m) {
  tibble(iso = wide$iso, year = wide$year,
         dataset = m, temp = T_mat[, m], resid = R_mat[, m])
})
synth_long <- map_dfr(synth_names_plot, function(m) {
  tibble(iso = wide$iso, year = wide$year,
         dataset = m,
         temp = synth_vecs[[match(m, synth_names_plot)]],
         resid = R_synth[[m]])
})

colours <- c(MEAN='#1f77b4', MEDIAN='#ff7f0e', PC1='#2ca02c',
             FACTOR='#d62728', `MIN-VAR`='#9467bd', ANCHORED='#8c564b')

make_spaghetti <- function(val_col, ylab, title) {
  ggplot() +
    geom_line(data = raw_long %>% filter(iso %in% top8),
              aes(year, .data[[val_col]], group = dataset),
              color = 'grey40', alpha = 0.35, linewidth = 0.4) +
    geom_line(data = synth_long %>% filter(iso %in% top8),
              aes(year, .data[[val_col]], color = dataset),
              linewidth = 0.6) +
    scale_color_manual(values = colours) +
    facet_wrap(~ factor(iso, levels = top8), scales = 'free_y', ncol = 2) +
    labs(x = 'Year', y = ylab, color = 'Synthetic', title = title) +
    theme(legend.position = 'bottom')
}

p_sp_raw <- make_spaghetti('temp', 'Temperature (°C)',
                           'Raw + synthetic temperatures, top-8 disagreement countries')
save_fig(p_sp_raw, 'sc_spaghetti_raw.pdf', w = 10, h = 9)

p_sp_resid <- make_spaghetti('resid', 'Residualised temp (°C)',
                             'Residualised raw + synthetic, top-8 disagreement countries')
save_fig(p_sp_resid, 'sc_spaghetti_resid.pdf', w = 10, h = 9)

cat('\nAll figures written to', fig_dir, '\n')
