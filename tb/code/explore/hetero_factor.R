# Temperature-dependent measurement-error variance (sigma_k(T)):
# a step-by-step build. This first step is a diagnostic — does the constant-
# sigma factor fit from synth_compare.R leave residuals whose variance depends
# on temperature level?
#
# Self-contained: mirrors the data load + residualisation from synth_compare.R.

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
T_mat  <- as.matrix(wide[, T_cols]); colnames(T_mat) <- mods
message(sprintf('Balanced panel: %d rows × %d datasets', nrow(wide), length(mods)))

# ============================================================================
# 2. Residualise each raw T on FE structure (same as synth_compare.R) --------
# ============================================================================
residualise_T <- function(m) {
  fml <- as.formula(paste0('`temp_1_', m, '` ~ 1 | iso + year + iso[time1] + iso[time2]'))
  resid(feols(fml, wide, notes = FALSE))
}
R_mat <- do.call(cbind, map(mods, residualise_T))
colnames(R_mat) <- mods

# ============================================================================
# 3. Constant-sigma per-country factor model (EM) ----------------------------
# ============================================================================
factor_fit_country <- function(R_c, max_iter = 60, tol = 1e-7) {
  K <- ncol(R_c); Tn <- nrow(R_c)
  pc  <- prcomp(R_c, center = FALSE)
  F_h <- pc$x[, 1] / pc$sdev[1]
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

F_hat      <- numeric(nrow(wide))
lam_mat    <- matrix(0, nrow(wide), length(mods),
                     dimnames = list(NULL, mods))
factor_info_c <- list()
for (iso_c in unique(wide$iso)) {
  ii <- which(wide$iso == iso_c)
  fc <- factor_fit_country(R_mat[ii, , drop = FALSE])
  factor_info_c[[iso_c]] <- fc
  F_hat[ii] <- fc$F
  lam_mat[ii, ] <- matrix(fc$lam, nrow = length(ii), ncol = length(mods),
                          byrow = TRUE)
}

# Residuals from the constant-sigma fit: epsilon_k(c,t) = R_k - lam_k(c) * F(c,t)
eps_mat <- R_mat - lam_mat * F_hat
colnames(eps_mat) <- mods

# ============================================================================
# 4. STEP 1 — DIAGNOSTIC: sigma_k(decile) ------------------------------------
# ============================================================================
# Bin by pooled decile of country-mean-centred ensemble temperature. This is
# "hot for this country" rather than "hot" absolutely — separates climate
# from weather.
T_bar   <- rowMeans(T_mat)
wide$T_bar        <- T_bar
wide$T_bar_dev    <- T_bar - ave(T_bar, wide$iso)              # within-country
wide$decile_abs   <- ntile(wide$T_bar,     10)                  # absolute-T pool
wide$decile_dev   <- ntile(wide$T_bar_dev, 10)                  # weather-T pool

# sigma_k(d) for each binning. cut = "abs" (climate axis) vs "dev" (weather axis)
sigma_by_bin <- function(bin_col, bin_label) {
  map_dfr(mods, function(k) {
    tibble(dataset = k, bin_type = bin_label,
           decile  = wide[[bin_col]],
           eps_sq  = eps_mat[, k]^2) %>%
      group_by(dataset, bin_type, decile) %>%
      summarise(sigma_hat = sqrt(mean(eps_sq)), n = n(), .groups = 'drop')
  })
}
diag_sigma <- bind_rows(
  sigma_by_bin('decile_abs', 'Absolute T (climate axis)'),
  sigma_by_bin('decile_dev', 'Within-country deviation (weather axis)'))

cat('\n--- sigma_k by absolute-T decile (d1=coldest country-years, d10=hottest) ---\n')
print(diag_sigma %>% filter(bin_type == 'Absolute T (climate axis)') %>%
        pivot_wider(id_cols = dataset, names_from = decile,
                    values_from = sigma_hat, names_prefix = 'd') %>%
        mutate(across(-dataset, ~round(.x, 3))), n = Inf)

cat('\n--- sigma_k by within-country-deviation decile ---\n')
print(diag_sigma %>% filter(bin_type == 'Within-country deviation (weather axis)') %>%
        pivot_wider(id_cols = dataset, names_from = decile,
                    values_from = sigma_hat, names_prefix = 'd') %>%
        mutate(across(-dataset, ~round(.x, 3))), n = Inf)

ratio_df <- diag_sigma %>% filter(decile %in% c(1, 10)) %>%
  pivot_wider(id_cols = c(dataset, bin_type),
              names_from = decile, values_from = sigma_hat,
              names_prefix = 'd') %>%
  mutate(ratio_hot_over_cold = d10 / d1)
cat('\n--- Ratio sigma(d10) / sigma(d1) per dataset and binning ---\n')
print(ratio_df)

# --- Plot A: sigma_k(decile), faceted by binning ---------------------------
p_dec <- ggplot(diag_sigma, aes(decile, sigma_hat, color = dataset)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:10) +
  facet_wrap(~bin_type, scales = 'free_y') +
  labs(x = 'Decile (1 = coldest, 10 = hottest)',
       y = expression(hat(sigma)[k] ~ '(°C)'),
       color = NULL,
       title = 'Residual SD by dataset and temperature decile')
save_fig(p_dec, 'hf_diag_sigma_decile.pdf', w = 11, h = 5)

# --- Plot B: sigma-loess vs raw T level ------------------------------------
# Here we plot |epsilon_k| against absolute temperature level (ensemble mean),
# with a loess smoother per dataset.
loess_long <- map_dfr(mods, function(k) {
  tibble(dataset = k,
         T_level = T_bar,
         abs_eps = abs(eps_mat[, k]))
})
p_loess <- ggplot(loess_long, aes(T_level, abs_eps, color = dataset)) +
  geom_smooth(method = 'loess', formula = y ~ x, se = FALSE, linewidth = 0.6) +
  labs(x = 'Ensemble-mean temperature level (°C)',
       y = expression(group('|', hat(epsilon)[k], '|') ~ '(°C)'),
       color = NULL,
       title = 'Absolute residual vs temperature level, loess smoother per dataset')
save_fig(p_loess, 'hf_diag_sigma_vs_T.pdf', w = 9, h = 5)

cat('\nStep 1 figures written to', fig_dir, '\n')

# ============================================================================
# 5. STEP 2 — σ_k(d) AS INPUT FOR T-DEPENDENT SYNTHESIS ----------------------
# ============================================================================
# We originally planned iterative GLS: alternately update F, λ, and σ_k(d)
# under the heteroskedastic factor model R_k = λ_k F + ε_k(d). That fit is
# pathological: the lowest-noise dataset (ERA5) has σ → 0 under GLS because
# the weight 1/σ² feeds back to concentrate F onto it, shrinking residuals,
# shrinking σ further. Even with a σ-floor the fit just rests on the floor
# for ERA5 at every decile — no useful heteroskedasticity estimate.
#
# Rather than report a degenerate MLE, we use the **descriptive** σ_k(d) from
# Step 1 (residuals of the constant-σ factor fit): these are biased (MLE
# under the wrong model) but they are the correct first-order estimate of
# residual SD per (dataset, decile) and are suitable for constructing
# inverse-variance weights in Step 3.

bin_vec <- wide$decile_abs
D       <- max(bin_vec)

sigma_hat_d <- diag_sigma %>%
  filter(bin_type == 'Absolute T (climate axis)') %>%
  select(dataset, decile, sigma_hat) %>%
  pivot_wider(names_from = dataset, values_from = sigma_hat) %>%
  arrange(decile) %>%
  select(-decile) %>%
  as.matrix()                                     # D × K
colnames(sigma_hat_d) <- mods
cat('\nσ_k(d) from Step-1 residuals (°C), by absolute-T decile:\n')
print(round(sigma_hat_d, 3))

# Implied GLS weights (row-normalised 1/σ²):
sigma2_d <- sigma_hat_d^2
weights_d <- t(apply(1 / sigma2_d, 1, function(x) x / sum(x)))
colnames(weights_d) <- mods
cat('\nImplied inverse-variance weights w_k(d), row sums to 1:\n')
print(round(weights_d, 3))

# Extremes of dataset weight across deciles — shows which dataset dominates
# and where.
cat('\nPer-dataset weight range across deciles (min / max):\n')
print(round(apply(weights_d, 2, range), 3))

# --- Plot A: σ_k(d) per dataset, big clean version --------------------------
p_sig <- as_tibble(sigma_hat_d) %>%
  mutate(decile = row_number()) %>%
  pivot_longer(-decile, names_to = 'dataset', values_to = 'sigma_hat') %>%
  ggplot(aes(decile, sigma_hat, color = dataset)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:10) +
  labs(x = 'Absolute-T decile (1 = coldest country-years, 10 = hottest)',
       y = expression(hat(sigma)[k](d) ~ '(°C)'),
       color = NULL,
       title = 'Residual SD per dataset and absolute-T decile')
save_fig(p_sig, 'hf_sigma_curves.pdf', w = 9, h = 4.5)

# --- Plot B: GLS weight heatmap ---------------------------------------------
weights_long <- as_tibble(weights_d) %>%
  mutate(decile = row_number()) %>%
  pivot_longer(-decile, names_to = 'dataset', values_to = 'weight')
p_heat <- ggplot(weights_long, aes(factor(decile), dataset, fill = weight)) +
  geom_tile() +
  geom_text(aes(label = sprintf('%.2f', weight)), size = 3) +
  scale_fill_gradient(low = 'white', high = 'darkblue') +
  labs(x = 'Absolute-T decile', y = NULL, fill = 'GLS\nweight',
       title = 'Optimal inverse-variance weights (row normalised)')
save_fig(p_heat, 'hf_weights_heatmap.pdf', w = 10, h = 3.5)

cat('\nStep 2 figures written to', fig_dir, '\n')

# ============================================================================
# 6. STEP 3 — T-DEPENDENT MIN-VAR SYNTHESIS ----------------------------------
# ============================================================================

# Ensemble mean precipitation, shared precip control across all methods.
P_cols <- paste0('prec_1_', mods)
P_mat  <- as.matrix(wide[, P_cols])
wide$P_mean    <- rowMeans(P_mat)
wide$P_mean_sq <- wide$P_mean^2

# 6.1 Build the synthetic temperatures ---------------------------------------
T_MEAN <- rowMeans(T_mat)

# Constant-σ MIN-VAR (per-country) — matches the method in synth_compare.R
T_MINVAR_const <- numeric(nrow(wide))
for (iso_c in unique(wide$iso)) {
  ii   <- which(wide$iso == iso_c)
  sig2 <- pmax(factor_info_c[[iso_c]]$sig2, 1e-10)
  w_c  <- (1 / sig2) / sum(1 / sig2)
  T_MINVAR_const[ii] <- as.vector(T_mat[ii, , drop = FALSE] %*% w_c)
}

# T-dependent MIN-VAR: weights depend on absolute-T decile.
# weights_d is 10 × 5; bin_vec gives each obs its decile; apply row-wise.
W_obs <- weights_d[bin_vec, ]   # N × K, weights per observation
T_MINVAR_T <- rowSums(T_mat * W_obs)

# FACTOR synthetic (per-country F from Step 1, rescaled) — as in synth_compare.
fit_mat <- T_mat - R_mat
T_FACTOR <- numeric(nrow(wide))
for (iso_c in unique(wide$iso)) {
  ii <- which(wide$iso == iso_c)
  fc <- factor_info_c[[iso_c]]
  T_FACTOR[ii] <- rowMeans(fit_mat[ii, , drop = FALSE]) + fc$F * mean(fc$lam)
}

# 6.2 Panel regression evaluator (mirrors synth_compare.R) ------------------
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

method_temps <- list(
  `ERA5-025`    = T_mat[, 'ERA5-025'],
  MERRA2        = T_mat[, 'MERRA2'],
  MEAN          = T_MEAN,
  FACTOR        = T_FACTOR,
  `MIN-VAR`     = T_MINVAR_const,
  `MIN-VAR-T`   = T_MINVAR_T
)

me_results <- map2_dfr(method_temps, names(method_temps), evaluate_method)
me_at_25 <- me_results %>% filter(T_eval == 25) %>%
  mutate(lo = me - 1.96*se_me, hi = me + 1.96*se_me, ci_width = hi - lo,
         t_beta2 = beta2 / se_beta2,
         kind = if_else(method %in% mods, 'raw', 'synthetic'))

cat('\nME at 25°C comparison (sorted by |t| on β₂):\n')
print(me_at_25 %>% arrange(desc(abs(t_beta2))) %>%
        select(method, beta2, se_beta2, t_beta2, me, ci_width, kind), n = Inf)

# 6.3 Forest plot with MIN-VAR-T added --------------------------------------
p_forest <- me_at_25 %>%
  mutate(b2_lo = beta2 - 1.96*se_beta2, b2_hi = beta2 + 1.96*se_beta2,
         highlight = method == 'MIN-VAR-T') %>%
  ggplot(aes(x = reorder(method, -beta2), y = beta2,
             color = kind, size = highlight)) +
  geom_point() +
  geom_errorbar(aes(ymin = b2_lo, ymax = b2_hi), width = 0.25, size = 0.6) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40') +
  scale_color_manual(values = c(raw = 'grey40', synthetic = 'steelblue')) +
  scale_size_manual(values = c(`FALSE` = 2.5, `TRUE` = 4),
                    guide = 'none') +
  coord_flip() +
  labs(x = NULL, y = expression(hat(beta)[2]~' with 95% CI'),
       color = NULL,
       title = 'T-dependent MIN-VAR joins the synthesis comparison')
save_fig(p_forest, 'hf_forest_update.pdf', w = 9, h = 4.5)

# 6.4 Applied-weight distribution for MIN-VAR-T -----------------------------
weight_long <- as_tibble(W_obs) %>%
  mutate(iso = wide$iso, decile = wide$decile_abs) %>%
  pivot_longer(all_of(mods), names_to = 'dataset', values_to = 'weight')

p_weights <- ggplot(weight_long, aes(weight, fill = dataset)) +
  geom_histogram(bins = 30, color = 'white', linewidth = 0.1) +
  facet_wrap(~dataset, scales = 'free_y', nrow = 1) +
  labs(x = 'Applied MIN-VAR-T weight across observations',
       y = 'N observations', fill = NULL,
       title = 'Distribution of per-observation dataset weights') +
  theme(legend.position = 'none')
save_fig(p_weights, 'hf_weight_distribution.pdf', w = 12, h = 3.5)

# 6.5 Jackknife: leave-one-dataset-out stability of MIN-VAR-T ----------------
# Drop dataset k*, refit Step-1 per-country factor on the remaining 4,
# recompute σ_{k'}(d), rebuild weights, rebuild synthetic, evaluate ME(25).
jackknife_mvt <- function(leave_out) {
  idx  <- setdiff(mods, leave_out)
  R_s  <- R_mat[, idx, drop = FALSE]
  T_s  <- T_mat[, idx, drop = FALSE]
  # 1. per-country factor fit on 4 datasets
  F_s <- numeric(nrow(wide)); lam_s <- matrix(0, nrow(wide), length(idx))
  for (iso_c in unique(wide$iso)) {
    ii <- which(wide$iso == iso_c)
    fc <- factor_fit_country(R_s[ii, , drop = FALSE])
    F_s[ii] <- fc$F
    lam_s[ii, ] <- matrix(fc$lam, nrow = length(ii), ncol = length(idx), byrow = TRUE)
  }
  # 2. σ_{k'}(d) pooled
  eps_s <- R_s - lam_s * F_s
  sig2_s <- matrix(0, D, length(idx))
  for (d in 1:D) {
    sig2_s[d, ] <- pmax(colMeans(eps_s[bin_vec == d, , drop = FALSE]^2), 1e-10)
  }
  # 3. weights, 4. synthetic
  w_s  <- t(apply(1 / sig2_s, 1, function(x) x / sum(x)))
  W_s  <- w_s[bin_vec, ]
  T_s_synth <- rowSums(T_s * W_s)
  evaluate_method(T_s_synth, paste0('MVT_no_', leave_out)) %>%
    filter(T_eval == 25) %>% mutate(leave_out = leave_out)
}
jk_mvt <- map_dfr(mods, jackknife_mvt)
cat('\nJackknife stability of MIN-VAR-T (leave one dataset out):\n')
print(jk_mvt %>% select(leave_out, me, beta2, se_beta2))
cat(sprintf('  range of ME(25°C) = %.5f\n', max(jk_mvt$me) - min(jk_mvt$me)))

cat('\nStep 3 figures written to', fig_dir, '\n')

# ============================================================================
# 7. STEP 4 — T-DEPENDENT ATTENUATION CORRECTION ----------------------------
# ============================================================================
# Classical EIV correction.  For the MIN-VAR-T synthesis, the effective
# measurement-error variance at observation (c,t) with decile d is
#   sigma2_synth(d) = 1 / sum_k (1 / sigma_k^2(d)).
# The residualised regressors are (tilde T, tilde T^2).  Their
# measurement-error components are (tilde eta, 2 bar_T * tilde eta + eta^2).
# We apply a scalar correction separately to beta_1 and beta_2 using the
# classical formula beta_corr = beta_obs / kappa where
#   kappa = (Var(regressor observed) - Var(noise)) / Var(regressor observed).

# sigma^2 of the MIN-VAR-T synthesis per decile
sigma2_synth <- 1 / rowSums(1 / sigma_hat_d^2)   # length D
cat('\nSynthetic σ per decile of absolute T (°C):\n')
print(round(sqrt(sigma2_synth), 4))

# Per-observation noise variance
sigma2_synth_obs <- sigma2_synth[bin_vec]         # length N

# Residualise MIN-VAR-T and its square on FE + precip
nuis_pfml <- ' ~ P_mean + P_mean_sq | iso + year + iso[time1] + iso[time2]'
wide$.Tmv    <- T_MINVAR_T
wide$.Tmv_sq <- T_MINVAR_T^2
Tmv_r  <- resid(feols(as.formula(paste0('.Tmv',    nuis_pfml)), wide, notes = FALSE))
Tmv_sq_r <- resid(feols(as.formula(paste0('.Tmv_sq', nuis_pfml)), wide, notes = FALSE))

# Per-observation noise variances (classical MEM + leading-order chain rule)
var_noise_T     <- mean(sigma2_synth_obs)
var_noise_Tsq   <- mean(4 * T_MINVAR_T^2 * sigma2_synth_obs)

# Observed regressor variances (after residualisation on FE + precip)
var_obs_T   <- var(Tmv_r)
var_obs_Tsq <- var(Tmv_sq_r)

kappa_1 <- max((var_obs_T   - var_noise_T)   / var_obs_T,   1e-6)
kappa_2 <- max((var_obs_Tsq - var_noise_Tsq) / var_obs_Tsq, 1e-6)

cat(sprintf('\nAttenuation factors for MIN-VAR-T:\n'))
cat(sprintf('  kappa_1 (beta_1) = %.4f   (1-kappa_1 = %.2f%% noise)\n',
            kappa_1, 100 * (1 - kappa_1)))
cat(sprintf('  kappa_2 (beta_2) = %.4f   (1-kappa_2 = %.2f%% noise)\n',
            kappa_2, 100 * (1 - kappa_2)))

# Corrected coefficients
mv_fit <- evaluate_method(T_MINVAR_T, 'MIN-VAR-T')
b1_obs  <- unique(mv_fit$beta1)
b2_obs  <- unique(mv_fit$beta2)
se_b1   <- sqrt(mean(c((mv_fit$se_me - 2 * mv_fit$T_eval * mv_fit$se_beta2)^2)))  # unused; keep SEs scaled below
se_b2   <- unique(mv_fit$se_beta2)

b1_corr <- b1_obs / kappa_1
b2_corr <- b2_obs / kappa_2
# SE scales inversely with kappa (since beta_corr = beta/kappa)
se_b1_corr <- se_b1 / kappa_1
se_b2_corr <- se_b2 / kappa_2

cat(sprintf('\nObserved  vs  corrected coefficients:\n'))
cat(sprintf('  beta_1:  %+.6f  →  %+.6f  (×%.3f)\n', b1_obs, b1_corr, 1/kappa_1))
cat(sprintf('  beta_2:  %+.6f  →  %+.6f  (×%.3f)\n', b2_obs, b2_corr, 1/kappa_2))

# --- Plot A: kappa(d) per decile (diagnostic of local identification) -------
kappa_d <- tibble(
  decile      = 1:D,
  var_F_d     = sapply(1:D, function(d) var(F_hat[bin_vec == d])),
  sigma2_d    = sigma2_synth,
  kappa       = var_F_d / (var_F_d + sigma2_d)
)
p_kappa <- ggplot(kappa_d, aes(decile, kappa)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey40') +
  scale_x_continuous(breaks = 1:10) +
  ylim(0, 1.05) +
  labs(x = 'Absolute-T decile',
       y = expression(kappa(d) == frac(Var(F | d), Var(F|d) + bar(sigma)^2(d))),
       title = 'Local identification strength under MIN-VAR-T synthesis')
save_fig(p_kappa, 'hf_attenuation_factor.pdf', w = 8, h = 4)

# --- Plot B: corrected vs uncorrected ME curve ------------------------------
T_grid <- seq(0, 32, by = 0.5)
me_obs  <- b1_obs + 2 * b2_obs * T_grid
me_corr <- b1_corr + 2 * b2_corr * T_grid
# Combined SE for ME via delta method.  For observed, from fit.
fit_me <- feols(growthWDI ~ T_ + T_sq + P_mean + P_mean_sq |
                  iso + year + iso[time1] + iso[time2],
                wide %>% mutate(T_ = T_MINVAR_T, T_sq = T_MINVAR_T^2),
                cluster = 'iso', notes = FALSE)
V <- vcov(fit_me)[c('T_','T_sq'), c('T_','T_sq')]
se_obs  <- sqrt(V[1,1] + 4*T_grid^2*V[2,2] + 4*T_grid*V[1,2])
# Corrected SE: inflate β's by 1/kappa's.
# Var(ME_corr) = V11/κ₁² + 4T²V22/κ₂² + 4T*Cov(β̂1/κ₁, β̂2/κ₂)
#             = V11/κ₁² + 4T²V22/κ₂² + 4T*V12/(κ₁κ₂)
se_corr <- sqrt(V[1,1]/kappa_1^2 + 4*T_grid^2*V[2,2]/kappa_2^2 +
                4*T_grid*V[1,2]/(kappa_1*kappa_2))

me_df <- bind_rows(
  tibble(temp = T_grid, me = me_obs,  lo = me_obs  - 1.96*se_obs,  hi = me_obs  + 1.96*se_obs,  stage = 'Uncorrected (MIN-VAR-T)'),
  tibble(temp = T_grid, me = me_corr, lo = me_corr - 1.96*se_corr, hi = me_corr + 1.96*se_corr, stage = 'EIV-corrected')
)
p_me_corr <- ggplot(me_df, aes(temp, me, color = stage, fill = stage)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40') +
  labs(x = 'Temperature (°C)', y = 'Marginal effect on growth',
       color = NULL, fill = NULL,
       title = 'EIV-corrected ME curve (MIN-VAR-T synthesis)')
save_fig(p_me_corr, 'hf_me_corrected.pdf', w = 9, h = 4.5)

# Report ME(25°C) with and without correction
me_25_obs  <- b1_obs  + 2 * b2_obs  * 25
me_25_corr <- b1_corr + 2 * b2_corr * 25
cat(sprintf('\nME(25°C)  observed:  %+.5f\n',  me_25_obs))
cat(sprintf('ME(25°C)  corrected: %+.5f\n',   me_25_corr))
cat(sprintf('  amplification factor = %.3f\n', me_25_corr / me_25_obs))

# --- Comparison: attenuation correction applied to each raw dataset ---------
# Treat each raw dataset on its own: noise variance = its own sigma_k(d).
compare_kappa <- map_dfr(c(mods, 'MIN-VAR-T'), function(m) {
  T_vec <- if (m == 'MIN-VAR-T') T_MINVAR_T else T_mat[, m]
  sig2_vec <- if (m == 'MIN-VAR-T') sigma2_synth_obs
              else (sigma_hat_d[bin_vec, m])^2
  wloc <- wide; wloc$.T <- T_vec; wloc$.Tsq <- T_vec^2
  Tr   <- resid(feols(as.formula(paste0('.T',    nuis_pfml)), wloc, notes = FALSE))
  Tsr  <- resid(feols(as.formula(paste0('.Tsq', nuis_pfml)), wloc, notes = FALSE))
  vn_T  <- mean(sig2_vec)
  vn_T2 <- mean(4 * T_vec^2 * sig2_vec)
  k1    <- max((var(Tr)  - vn_T ) / var(Tr ), 1e-6)
  k2    <- max((var(Tsr) - vn_T2) / var(Tsr), 1e-6)
  mv    <- evaluate_method(T_vec, m)
  b1 <- unique(mv$beta1); b2 <- unique(mv$beta2)
  tibble(method = m,
         kappa_1 = k1, kappa_2 = k2,
         beta2_obs = b2, beta2_corr = b2 / k2,
         me25_obs = b1 + 2*b2*25,
         me25_corr = (b1/k1) + 2*(b2/k2)*25)
})
cat('\nAttenuation correction per method (raw datasets vs MIN-VAR-T):\n')
print(compare_kappa %>% mutate(across(where(is.numeric), ~signif(.x, 4))))

# --- Plot C: cross-method ME(25°C) with and without EIV correction ----------
ck_long <- compare_kappa %>%
  select(method, me25_obs, me25_corr) %>%
  pivot_longer(-method, names_to = 'stage', values_to = 'me25') %>%
  mutate(stage = recode(stage, me25_obs = 'Uncorrected',
                                me25_corr = 'EIV-corrected'),
         stage = factor(stage, levels = c('Uncorrected','EIV-corrected')))

p_cross <- ggplot(ck_long, aes(x = reorder(method, me25), y = me25,
                                color = stage, shape = stage)) +
  geom_point(size = 3) +
  geom_line(aes(group = method), color = 'grey60', linewidth = 0.3,
            position = position_nudge(x = 0)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40') +
  scale_color_manual(values = c(Uncorrected = 'grey40',
                                 `EIV-corrected` = 'firebrick')) +
  coord_flip() +
  labs(x = NULL, y = 'ME at 25°C',
       color = NULL, shape = NULL,
       title = 'EIV correction collapses the cross-dataset ME spread')
save_fig(p_cross, 'hf_me_corrected_compare.pdf', w = 9, h = 4.5)

cat('\nStep 4 figures written to', fig_dir, '\n')
