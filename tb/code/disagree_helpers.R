# Shared helpers for the disagreement-vs-covariate analysis.
# Used by gdp_disagree.R and mortality_disagree.R.
#
# Expects: tidyverse, fixest, patchwork, ggrepel loaded by the caller.

# ---------------------------------------------------------------------------
# make_disagr_df: pairwise residualised disagreement aggregated to a grouping.
#
# wide_r:     tibble with one row per observation and one column per dataset
#             holding residualised values (after partialling out the FE
#             structure of the main regression).
# mods:       character vector of dataset names (column names in wide_r).
# group_cols: any subset of columns to group by (e.g. c('iso'), c('year'),
#             c('iso','year'), c('adm1_id'), ...).
# metric:     'mad' = mean of pairwise |A - B| across all pairs (default),
#             'sd'  = mean of cross-dataset SD per row.
#
# Returns:    tibble with group_cols + a single `disagr` column.
# ---------------------------------------------------------------------------
make_disagr_df <- function(wide_r, mods, group_cols = c('iso'),
                           metric = c('mad', 'sd')) {
  metric <- match.arg(metric)
  if (metric == 'mad') {
    pairs_df <- combn(mods, 2) %>% t() %>% as.data.frame() %>%
      setNames(c('A', 'B')) %>% as_tibble()
    pmap_dfr(pairs_df, function(A, B) {
      wide_r %>%
        mutate(d = abs(.data[[A]] - .data[[B]])) %>%
        group_by(across(all_of(group_cols))) %>%
        summarise(mad = mean(d), .groups = 'drop')
    }) %>%
      group_by(across(all_of(group_cols))) %>%
      summarise(disagr = mean(mad), .groups = 'drop')
  } else {
    wide_r %>%
      mutate(sd_row = apply(across(all_of(mods)), 1, sd)) %>%
      group_by(across(all_of(group_cols))) %>%
      summarise(disagr = mean(sd_row), .groups = 'drop')
  }
}

# ---------------------------------------------------------------------------
# plot_disagr: scatter of disagr vs an x-covariate, with OLS line + R².
# colour_var: column to map to point colour (default 'continent').  Auto-
# suppresses if the column is missing or all-NA.
# ---------------------------------------------------------------------------
plot_disagr <- function(dm, x = 'log_gdp', label_var = NULL, label_top = 8,
                        x_lab = NULL, y_lab = 'Disagreement (°C)',
                        colour_var = 'continent',
                        point_alpha = NULL) {
  r2 <- summary(lm(reformulate(x, response = 'disagr'), dm))$r.squared
  has_colour <- colour_var %in% names(dm) && !all(is.na(dm[[colour_var]]))
  if (is.null(point_alpha))
    point_alpha <- if (nrow(dm) > 500) 0.3 else 0.8

  p <- ggplot(dm, aes(.data[[x]], disagr))
  p <- if (has_colour)
         p + geom_point(aes(color = .data[[colour_var]]),
                        size = 2, alpha = point_alpha)
       else
         p + geom_point(size = 2, alpha = point_alpha)
  p <- p +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE,
                color = 'black', linetype = 'dashed') +
    annotate('text', x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = sprintf('R² = %.3f', r2), size = 4) +
    labs(x = x_lab %||% x, y = y_lab, color = NULL) +
    theme(legend.position = if (has_colour) 'bottom' else 'none')

  if (!is.null(label_var)) {
    p <- p + ggrepel::geom_text_repel(
      data = dm %>% arrange(desc(disagr)) %>% head(label_top),
      aes(label = .data[[label_var]]),
      size = 3, show.legend = FALSE, max.overlaps = Inf)
  }
  p
}

# ---------------------------------------------------------------------------
# .extract_int: pull (intercept, β_x, β_y, β_xy) from a fitted model whose
# formula was `disagr ~ x + y + x:y` (regardless of fixest's coefficient
# ordering).  Returns named numeric vector and matching covariance matrix.
# ---------------------------------------------------------------------------
.extract_int <- function(model, x, y) {
  b   <- coef(model)
  V   <- vcov(model)
  ix  <- which(names(b) == '(Intercept)')
  ig  <- which(names(b) == x)
  it  <- which(names(b) == y)
  iigt <- which(names(b) %in% c(paste(x, y, sep = ':'),
                                 paste(y, x, sep = ':')))
  list(b = c(int = unname(b[ix]),
             x   = unname(b[ig]),
             y   = unname(b[it]),
             xy  = unname(b[iigt])),
       V = V[c(ix, ig, it, iigt), c(ix, ig, it, iigt)])
}

# ---------------------------------------------------------------------------
# plot_disagr_heatmap: viridis raster of predicted ŷ from an interaction
# model, with country dots and 5 contour lines.
# ---------------------------------------------------------------------------
plot_disagr_heatmap <- function(model, dm, x = 'log_gdp', y = 'temp',
                                label_var = NULL, label_top = 6,
                                fill_lab = 'predicted\ndisagr (°C)',
                                x_lab = NULL, y_lab = NULL,
                                title = NULL, grid_n = 80) {
  e <- .extract_int(model, x, y)
  pred <- function(xv, yv) e$b['int'] + e$b['x']*xv + e$b['y']*yv +
                            e$b['xy']*xv*yv

  x_seq <- seq(min(dm[[x]]), max(dm[[x]]), length.out = grid_n)
  y_seq <- seq(min(dm[[y]]), max(dm[[y]]), length.out = grid_n)
  grid <- expand_grid(!!x := x_seq, !!y := y_seq) %>%
    mutate(disagr_hat = pred(.data[[x]], .data[[y]]))

  p <- ggplot() +
    geom_raster(data = grid,
                aes(.data[[x]], .data[[y]], fill = disagr_hat)) +
    geom_contour(data = grid,
                 aes(.data[[x]], .data[[y]], z = disagr_hat),
                 bins = 5, color = 'white', linewidth = 0.3, alpha = 0.6) +
    geom_point(data = dm, aes(.data[[x]], .data[[y]]),
               shape = 21, size = 2, fill = 'white', color = 'black',
               stroke = 0.4) +
    scale_fill_viridis_c(name = fill_lab) +
    labs(x = x_lab %||% x, y = y_lab %||% y, title = title)

  if (!is.null(label_var)) {
    top <- dm %>% arrange(desc(disagr)) %>% head(label_top)
    p <- p + ggrepel::geom_text_repel(
      data = top,
      aes(.data[[x]], .data[[y]], label = .data[[label_var]]),
      size = 3, color = 'white', max.overlaps = Inf,
      segment.color = 'white')
  }
  p
}

# ---------------------------------------------------------------------------
# plot_disagr_marginal: two side-by-side panels (∂y/∂x by y, ∂y/∂y by x)
# with delta-method 95% CI ribbons and rug plots of observed values.
# ---------------------------------------------------------------------------
plot_disagr_marginal <- function(model, dm, x = 'log_gdp', y = 'temp',
                                 x_lab = NULL, y_lab = NULL,
                                 me_x_lab = NULL, me_y_lab = NULL,
                                 title_x = NULL, title_y = NULL) {
  e <- .extract_int(model, x, y)
  b <- e$b; V <- e$V

  me_x_fn <- function(yv) b['x'] + b['xy'] * yv
  se_x    <- function(yv) sqrt(V[2,2] + yv^2*V[4,4] + 2*yv*V[2,4])
  me_y_fn <- function(xv) b['y'] + b['xy'] * xv
  se_y    <- function(xv) sqrt(V[3,3] + xv^2*V[4,4] + 2*xv*V[3,4])

  y_grid <- seq(min(dm[[y]]), max(dm[[y]]), length.out = 200)
  x_grid <- seq(min(dm[[x]]), max(dm[[x]]), length.out = 200)
  me_x_df <- tibble(yv = y_grid, me = me_x_fn(y_grid),
                    se = se_x(y_grid)) %>%
    mutate(lo = me - 1.96*se, hi = me + 1.96*se)
  me_y_df <- tibble(xv = x_grid, me = me_y_fn(x_grid),
                    se = se_y(x_grid)) %>%
    mutate(lo = me - 1.96*se, hi = me + 1.96*se)

  p1 <- ggplot(me_x_df, aes(yv, me)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
    geom_line(linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40') +
    geom_rug(data = dm, aes(.data[[y]]), inherit.aes = FALSE,
             alpha = 0.3, length = unit(0.02, 'npc')) +
    labs(x = y_lab %||% y, y = me_x_lab %||% paste0('d disagr / d ', x),
         title = title_x)

  p2 <- ggplot(me_y_df, aes(xv, me)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
    geom_line(linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40') +
    geom_rug(data = dm, aes(.data[[x]]), inherit.aes = FALSE,
             alpha = 0.3, length = unit(0.02, 'npc')) +
    labs(x = x_lab %||% x, y = me_y_lab %||% paste0('d disagr / d ', y),
         title = title_y)

  p1 + p2
}

# ---------------------------------------------------------------------------
# plot_disagr_actual_2d: scatter of (x, y) coloured by actual disagr.
# Same axes/colour scale as plot_disagr_heatmap so they can be compared.
# ---------------------------------------------------------------------------
plot_disagr_actual_2d <- function(dm, x = 'log_gdp', y = 'temp',
                                  label_var = NULL, label_top = 6,
                                  fill_lab = 'actual\ndisagr (°C)',
                                  x_lab = NULL, y_lab = NULL, title = NULL) {
  p <- ggplot(dm, aes(.data[[x]], .data[[y]], fill = disagr)) +
    geom_point(shape = 21, size = 3, color = 'black', stroke = 0.3) +
    scale_fill_viridis_c(name = fill_lab) +
    labs(x = x_lab %||% x, y = y_lab %||% y, title = title)
  if (!is.null(label_var)) {
    top <- dm %>% arrange(desc(disagr)) %>% head(label_top)
    p <- p + ggrepel::geom_text_repel(
      data = top, aes(label = .data[[label_var]]),
      size = 3, max.overlaps = Inf)
  }
  p
}
