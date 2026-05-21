library(tidyverse)
library(useful)
library(fixest)
theme_set(theme_classic())

setwd('/Users/tombearpark/Library/CloudStorage/Dropbox/ks')


# LOAD DATA ---------------------------------------------------------------

balanced_sample <- TRUE  # keep only (iso, year) with non-NA temp in every obs_mod

df <- read_csv('reg_inputs/gdppc/reginputs_gdppc_1deg.csv') %>%
  mutate(temp_1 = temp_popweight,
         temp_2 = temp_1^2,
         prec_1 = precip_popweight,
         prec_2 = prec_1^2,
         time1 = year - min(year),
         time2 = time1^2,
         decade = 10L * (year %/% 10L))

mods <- sort(unique(df$obs_mod))

if (balanced_sample) {
  keep <- df %>%
    filter(!is.na(temp_1)) %>%
    distinct(iso, year, obs_mod) %>%
    count(iso, year) %>%
    filter(n == length(mods)) %>%
    select(iso, year)
  df <- df %>% semi_join(keep, by = c('iso', 'year'))
  message(sprintf('balanced_sample: kept %d (iso, year) pairs across all %d models',
                  nrow(keep), length(mods)))
}

# Temperature decile (pooled across obs_mod × iso × year).
df <- df %>%
  mutate(temp_decile = ntile(temp_1, 10))

df %>% 
  ggplot() + 
  geom_histogram(aes( x= temp_1), bins = 30) + 
  facet_wrap(~obs_mod)

df %>%
  filter(iso == "USA") %>%
  ggplot() +
  geom_line(aes(x = year, y = temp_1)) +
  facet_wrap(~obs_mod, nrow=1)

# RUN MODELS --------------------------------------------------------------

MMs <- lapply(mods, function(m) {
  feols(
    growthWDI ~ temp_1 + temp_2 + prec_1 + prec_2 | iso + year + iso[time1] + iso[time2],
    df %>% filter(obs_mod == m),
    cluster = 'iso',
    lean = TRUE
    )
  }
  ) %>% 
  set_names(mods)

Preds <- map_dfr(mods, function(m)
  predict_poly(MMs[[m]], var = 'temp', min = 0, max = 30, ref = 14) %>%
    mutate(model = m)
)

get_me <- function(m, xlims = 0:30, model){
  b  <- coef(m)[c('temp_1', 'temp_2')]
  V  <- vcov(m)[c('temp_1', 'temp_2'), c('temp_1', 'temp_2')]
  me <- b['temp_1'] + 2 * b['temp_2'] * xlims
  se <- sqrt(V['temp_1','temp_1'] +
             4 * xlims^2 * V['temp_2','temp_2'] +
             4 * xlims      * V['temp_1','temp_2'])
  tibble(temp = xlims,
         me = as.numeric(me),
         se = as.numeric(se),
         lo = me - 1.96 * se,
         hi = me + 1.96 * se,
         model = model)
}

ME <- map_dfr(mods, function(m) get_me(MMs[[m]], model = m))

## Plot ME and predictions 
plot_rf_poly(Preds) + facet_wrap(~model, nrow=1)
ggplot(Preds) + geom_line(aes(x = temp, y = response, color = model)) 

ggplot(ME, aes(x = temp, y = me)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  facet_wrap(~model, nrow=1)

ggplot(ME %>% filter(temp %in% c(5, 25))) + 
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_errorbar(aes(x = model, ymin = lo, ymax = hi, color = model), width = 0.5) + 
  facet_wrap(~temp)

# distance between models -------------------------------------------------
# Scalar distance between two ME curves, integrated over the empirical
# distribution of sample temperatures: E_T[ |ME_a(T) - ME_b(T)| ].

me_at <- function(m, t){
  cf <- coef(m)
  cf['temp_1'] + 2 * cf['temp_2'] * t
}

me_distance <- function(m1, m2, t, type = c('L1','L2')){
  type <- match.arg(type)
  d <- me_at(m1, t) - me_at(m2, t)
  if (type == 'L1') mean(abs(d)) else sqrt(mean(d^2))
}

# Per-model empirical temp distributions (df is long in obs_mod).
temps_by_mod <- split(df$temp_1, df$obs_mod) %>% lapply(na.omit)

D <- outer(mods, mods, Vectorize(function(a, b)
  me_distance(MMs[[a]], MMs[[b]],
              t = c(temps_by_mod[[a]], temps_by_mod[[b]]),
              type = 'L1')))
dimnames(D) <- list(mods, mods)
print(round(D, 4))

D_long <- as_tibble(D, rownames = 'model_a') %>%
  pivot_longer(-model_a, names_to = 'model_b', values_to = 'distance')

ggplot(D_long, aes(x = model_a, y = model_b, fill = distance)) +
  geom_tile() +
  geom_text(aes(label = sprintf('%.4f', distance)), color = 'white', size = 3) +
  labs(x = NULL, y = NULL, fill = 'L1 distance') +
  coord_equal() + 
  scale_fill_gradient(low = 'white', high = 'red') 

# explain -----------------------------------------------------------------
# "ERA5-025" vs "MERRA2"
# Leave-one-(continent|decade)-out: for each group g, drop all rows with that
# group, refit BOTH models on the remainder, recompute the distance.
# reduction = d_full - d_LOO.  Large positive => that group drives the gap.

fit_mod <- function(data, mod_name) {
  feols(
    growthWDI ~ temp_1 + temp_2 + prec_1 + prec_2 |
      iso + year + iso[time1] + iso[time2],
    data %>% filter(obs_mod == mod_name),
    cluster = 'iso', lean = TRUE
  )
}

pair_distance <- function(data, mod_a, mod_b, type = 'L1') {
  ma <- fit_mod(data, mod_a)
  mb <- fit_mod(data, mod_b)
  ta <- na.omit(data$temp_1[data$obs_mod == mod_a])
  tb <- na.omit(data$temp_1[data$obs_mod == mod_b])
  me_distance(ma, mb, t = c(ta, tb), type = type)
}

loo_scan <- function(data, group_var, mod_a, mod_b, type = 'L1') {
  groups <- sort(unique(data[[group_var]]))
  d_full <- pair_distance(data, mod_a, mod_b, type = type)
  map_dfr(groups, function(g) {
    d <- pair_distance(data %>% filter(.data[[group_var]] != g),
                       mod_a, mod_b, type = type)
    tibble(group_var = group_var,
           group = as.character(g),
           distance = d,
           full_distance = d_full,
           reduction = d_full - d)
  })
}

A <- 'ERA5-025'; B <- 'MERRA2'

LOO_cont <- loo_scan(df, 'continent',   A, B)
LOO_dec  <- loo_scan(df, 'decade',      A, B)
LOO_tq   <- loo_scan(df, 'temp_decile', A, B)
LOO <- bind_rows(LOO_cont, LOO_dec, LOO_tq)
print(LOO)

plot_loo <- function(loo_df, title) {
  ggplot(loo_df, aes(x = reorder(group, -distance), y = distance)) +
    geom_col(fill = 'steelblue') +
    geom_hline(aes(yintercept = full_distance), linetype = 'dashed') +
    coord_flip() +
    labs(x = NULL, y = 'LOO distance', title = title)
}

plot_loo(LOO_cont, paste('Continent LOO:',    A, 'vs', B))
plot_loo(LOO_dec,  paste('Decade LOO:',       A, 'vs', B))
plot_loo(LOO_tq,   paste('Temp-decile LOO:',  A, 'vs', B))

# investigate -------------------------------------------------------------
# Hypotheses for the ERA5 vs MERRA2 gap, concentrated in the hot tail:
#   H1: the two datasets disagree on temperature itself — bias or hot-tail divergence.
#   H2: a handful of countries drive the gap (leave-one-country-out).
#   H3: residual growth signal differs in the hot tail.

AB <- df %>%
  filter(obs_mod %in% c(A, B)) %>%
  select(iso, year, continent, obs_mod, temp_1, growthWDI) %>%
  pivot_wider(names_from = obs_mod, values_from = temp_1) %>%
  rename(t_A = !!A, t_B = !!B) %>%
  mutate(dT = t_A - t_B)

# H1a: scatter of the two temps + 45-degree line
ggplot(AB, aes(t_B, t_A)) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  labs(x = paste(B, 'temp'), y = paste(A, 'temp'),
       title = 'Temperature agreement: ERA5-025 vs MERRA2')

# H1b: how does the bias (dT) depend on temperature level?
ggplot(AB, aes(t_B, dT)) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_smooth(method = 'loess', formula = y ~ x) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(x = paste(B, 'temp'), y = paste(A, '-', B, 'temp'),
       title = 'Temperature bias vs. level')

cat('\n--- Temperature agreement (ERA5 vs MERRA2) ---\n')
AB %>%
  summarise(cor = cor(t_A, t_B),
            mean_bias = mean(dT),
            sd_bias   = sd(dT),
            mad_bias  = mean(abs(dT))) %>%
  print()

cat('\n--- Bias by pooled temp decile ---\n')
AB %>%
  mutate(dec = ntile(t_B, 10)) %>%
  group_by(dec) %>%
  summarise(mean_t_B = mean(t_B), mean_dT = mean(dT),
            sd_dT = sd(dT), n = n()) %>%
  print()

# H2: leave-one-country-out (top 15 most-influential)
LOO_iso <- loo_scan(df, 'iso', A, B)
cat('\n--- Top 15 countries by distance reduction ---\n')
LOO_iso %>% arrange(desc(reduction)) %>% head(15) %>% print()
cat('\n--- Bottom 10 (removing them widens the gap) ---\n')
LOO_iso %>% arrange(reduction) %>% head(10) %>% print()

# H3: compare mean growth in the hot tail across the two data frames.
# (df is long in obs_mod but growthWDI per iso-year is the same outcome —
#  so this is really asking: which obs enter the regression in the hot tail?)
cat('\n--- Hot-tail (top decile) composition ---\n')
df %>%
  filter(obs_mod %in% c(A, B), temp_decile == 10) %>%
  group_by(obs_mod) %>%
  summarise(n = n(),
            mean_temp = mean(temp_1, na.rm = TRUE),
            mean_growth = mean(growthWDI, na.rm = TRUE)) %>%
  print()

# Which iso-years are in the top decile for one dataset but not the other?
# (i.e. the reanalyses disagree about which country-years are "hot")
hot_disagree <- df %>%
  filter(obs_mod %in% c(A, B)) %>%
  select(iso, year, obs_mod, temp_decile) %>%
  pivot_wider(names_from = obs_mod, values_from = temp_decile) %>%
  rename(d_A = !!A, d_B = !!B) %>%
  filter((d_A == 10) != (d_B == 10))
cat(sprintf('\n%d (iso, year) pairs are top-decile in one dataset but not the other\n',
            nrow(hot_disagree)))

# check 1: drop the 140 hot-tail-disagreement pairs, recompute distance
df_no_disagree <- df %>% anti_join(hot_disagree %>% select(iso, year), by = c('iso','year'))
d_full      <- pair_distance(df,              A, B)
d_no_disagr <- pair_distance(df_no_disagree,  A, B)
cat(sprintf('\n--- Drop hot-tail-disagreement pairs ---\nfull distance: %.5f   after drop: %.5f   reduction: %.5f\n',
            d_full, d_no_disagr, d_full - d_no_disagr))

# check 2: Sierra Leone temperature series, both datasets
df %>%
  filter(iso == 'SLE', obs_mod %in% c(A, B)) %>%
  ggplot(aes(year, temp_1, color = obs_mod)) +
  geom_line() +
  geom_point(size = 1) +
  labs(title = 'Sierra Leone: ERA5-025 vs MERRA2 temperature', y = 'Temp (°C)')

df %>%
  filter(iso == 'SLE', obs_mod %in% c(A, B)) %>%
  select(iso, year, obs_mod, temp_1) %>%
  pivot_wider(names_from = obs_mod, values_from = temp_1) %>%
  mutate(diff = !!sym(A) - !!sym(B)) %>%
  ggplot(aes(year, diff)) +
  geom_col() +
  geom_hline(yintercept = 0) +
  labs(title = 'Sierra Leone: ERA5-025 minus MERRA2 temperature', y = 'diff (°C)')
