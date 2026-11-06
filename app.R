#==============================================================================
# BASIN: Bayesian Sediment Source Apportionment Tool
# Edition with Ghost Source Detector (Semi-Blind G+1 Unmixing)
# Part 1 of 2: Packages, Helpers, Stan Code Generator, UI
#==============================================================================

suppressPackageStartupMessages({
  suppressWarnings({
    if (!require("pacman", quietly = TRUE)) install.packages("pacman", quiet = TRUE)
    options(repos = c(CRAN = "https://cloud.r-project.org"))

    pacman::p_load(
      shiny, rstan, tidyverse, caret, reshape2, dplyr, ggplot2,
      DT, shinythemes, data.table, MASS, zip, compositions,
      scoringRules, loo
    )
    
    # Устанавливаем groundhog только один раз, если его нет
    #if (!requireNamespace("groundhog", quietly = TRUE)) {
    #  install.packages("groundhog")
    #}
    
    #library(groundhog)
    
    # Указываем дату, на которую всё работало стабильно.
    # Например, 1 июня 2025 года.
    #groundhog.library(
    #  c("shiny", "rstan", "tidyverse", "caret", "reshape2", "dplyr", 
    #    "ggplot2", "DT", "shinythemes", "data.table", "MASS", "zip", 
    #    "compositions", "scoringRules", "loo"),
    #  date = "2026-06-01"
    #)
  })
})

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


#==============================================================================
# 1. HELPER FUNCTIONS
#==============================================================================

# Безопасное масштабирование матрицы с защитой от нулевой дисперсии
safe_scale_matrix <- function(x) {
  x <- as.matrix(x)

  ctr <- colMeans(x, na.rm = TRUE)
  scl <- apply(x, 2, sd, na.rm = TRUE)

  ctr[is.na(ctr)] <- 0
  scl[is.na(scl)] <- 1
  scl[scl < 1e-12] <- 1

  list(
    y = scale(x, center = ctr, scale = scl),
    center = as.numeric(ctr),
    scale = as.numeric(scl)
  )
}


# Приводим выборки пропорций к 3D: [iterations, N_mixtures, N_groups]
ensure_3d_samples <- function(x, n_groups) {
  if (is.null(dim(x))) {
    return(matrix(x, ncol = n_groups))
  }

  if (length(dim(x)) == 2) {
    return(array(x, dim = c(dim(x)[1], 1, dim(x)[2])))
  }

  x
}


# log(mean(exp(x))) устойчиво
log_mean_exp <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)

  m <- max(x)
  m + log(mean(exp(x - m)))
}


# Безопасное создание group_id
make_group_id <- function(labels, groups) {
  match(as.character(labels), as.character(groups))
}


# Геометрическое расстояние от точки p до отрезка [a, b] на плоскости
dist_point_to_segment_2d <- function(p, a, b) {
  ab <- b - a
  ap <- p - a
  ab_len_sq <- sum(ab^2)
  if (ab_len_sq < 1e-12) return(sqrt(sum(ap^2)))
  t_proj <- sum(ap * ab) / ab_len_sq
  t_proj <- max(0, min(1, t_proj))
  closest <- a + t_proj * ab
  sqrt(sum((p - closest)^2))
}


# Проверка нахождения точки внутри выпуклого многоугольника (2D) и расстояние до его границы
dist_to_convex_hull_2d <- function(pts_matrix, target_pt) {
  pts_matrix <- as.matrix(pts_matrix[, 1:2, drop = FALSE])
  target_pt <- as.numeric(target_pt[1:2])

  if (nrow(pts_matrix) < 3) {
    return(sqrt(sum((target_pt - colMeans(pts_matrix))^2)))
  }

  hull_idx <- grDevices::chull(pts_matrix[, 1], pts_matrix[, 2])
  hull_pts <- pts_matrix[hull_idx, , drop = FALSE]
  n_h <- nrow(hull_pts)

  if (n_h < 3) {
    return(sqrt(sum((target_pt - colMeans(pts_matrix))^2)))
  }

  cross_signs <- numeric(n_h)
  min_edge_dist <- Inf

  for (i in seq_len(n_h)) {
    j <- if (i == n_h) 1L else i + 1L
    a <- hull_pts[i, ]
    b <- hull_pts[j, ]
    edge <- b - a
    to_p <- target_pt - a
    cross_signs[i] <- edge[1] * to_p[2] - edge[2] * to_p[1]
    d_edge <- dist_point_to_segment_2d(target_pt, a, b)
    if (d_edge < min_edge_dist) min_edge_dist <- d_edge
  }

  is_inside <- all(cross_signs <= 1e-9) || all(cross_signs >= -1e-9)
  if (is_inside) return(0.0)
  min_edge_dist
}


# Проекция вектора смеси y на симплекс средних известных источников M [G x T]
# с помощью быстрого алгоритма Франка-Вульфа и расчет якоря для Ghost Source
# Геометрический якорь Ghost Source с принципом внешнего отталкивания (Outward Repulsion)
calc_ghost_anchor <- function(y_src, group_id, n_groups, y_mix) {
  y_src <- as.matrix(y_src)
  y_mix <- as.matrix(y_mix)
  n_tr <- ncol(y_src)

  # Матрица средних известных источников [n_groups x n_tr]
  M <- matrix(0, nrow = n_groups, ncol = n_tr)
  for (g in seq_len(n_groups)) {
    idx <- which(group_id == g)
    if (length(idx) > 0) {
      M[g, ] <- colMeans(y_src[idx, , drop = FALSE], na.rm = TRUE)
    }
  }
  M[is.na(M)] <- 0

  # Проецируем каждую смесь на выпуклую оболочку известных источников M
  res_mat   <- matrix(0, nrow = nrow(y_mix), ncol = n_tr)
  proj_mat  <- matrix(0, nrow = nrow(y_mix), ncol = n_tr)
  res_norms <- numeric(nrow(y_mix))

  for (i in seq_len(nrow(y_mix))) {
    yi <- as.numeric(y_mix[i, ])
    w <- rep(1 / n_groups, n_groups)
    pred <- as.numeric(crossprod(w, M))

    for (iter in 1:80) {
      grad <- as.numeric(M %*% (pred - yi))
      k_min <- which.min(grad)
      v_dir <- -M[k_min, ] + pred
      denom <- sum(v_dir^2)
      if (denom < 1e-12) break
      gamma <- sum((pred - yi) * v_dir) / denom
      gamma <- max(0, min(1, gamma))
      if (gamma < 1e-7) break
      w <- (1 - gamma) * w
      w[k_min] <- w[k_min] + gamma
      pred <- as.numeric(crossprod(w, M))
    }

    r_i <- yi - pred
    proj_mat[i, ] <- pred
    res_mat[i, ]  <- r_i
    res_norms[i]  <- sqrt(sum(r_i^2))
  }

  max_norm <- max(res_norms, na.rm = TRUE)

  if (!is.finite(max_norm) || max_norm < 1e-6) {
    anchor <- colMeans(M) + rep(2.5 / sqrt(n_tr), n_tr)
  } else {
    top_idx <- which(res_norms >= 0.80 * max_norm)
    w_top   <- res_norms[top_idx]^2 / sum(res_norms[top_idx]^2)

    extreme_proj <- as.numeric(crossprod(w_top, proj_mat[top_idx, , drop = FALSE]))
    extreme_res  <- as.numeric(crossprod(w_top, res_mat[top_idx, , drop = FALSE]))
    res_len      <- sqrt(sum(extreme_res^2))

    if (res_len < 1e-6) {
      unit_dir <- rep(1 / sqrt(n_tr), n_tr)
    } else {
      unit_dir <- extreme_res / res_len
    }

    # Внешнее отталкивание: якорь всегда снаружи симплекса известных источников (>= 2.5 сигмы)
    dist_out <- max(2.5, 1.15 * res_len)
    anchor   <- extreme_proj + dist_out * unit_dir
  }

  anchor <- pmax(-4.5, pmin(4.5, anchor))
  anchor[!is.finite(anchor)] <- 0
  as.numeric(anchor)
}


# Сводная таблица параметров для UI и автосохранения
build_summary_df <- function(fit, groups, p_param, add_phi = FALSE) {
  res_p <- as.data.frame(summary(fit, pars = p_param)$summary)
  samples_p <- rstan::extract(fit, pars = p_param)[[1]]
  if (p_param %in% c("P_ind", "P_level")) {
    samples_p <- ensure_3d_samples(samples_p, length(groups))
    samples_p <- apply(samples_p, c(1, 3), mean)
  }
  if (is.null(dim(samples_p)) || length(dim(samples_p)) < 2) {
    samples_p <- matrix(samples_p, ncol = length(groups))
  }
  q25 <- apply(samples_p, 2, quantile, probs = 0.25)
  q75 <- apply(samples_p, 2, quantile, probs = 0.75)
  df_out <- data.frame(
    Source = groups,
    mean = colMeans(samples_p),
    sd = apply(samples_p, 2, sd),
    q025 = apply(samples_p, 2, quantile, probs = 0.025),
    q25 = q25,
    median = apply(samples_p, 2, median),
    q75 = q75,
    q975 = apply(samples_p, 2, quantile, probs = 0.975),
    W50 = q75 - q25,
    check.names = FALSE
  )
  if ("n_eff" %in% colnames(res_p)) {
    if (length(res_p$Rhat) > length(groups)) {
      rhat_mat <- matrix(res_p$Rhat, ncol = length(groups), byrow = TRUE)
      neff_mat <- matrix(res_p$n_eff, ncol = length(groups), byrow = TRUE)
      df_out$min_neff <- apply(neff_mat, 2, min, na.rm = TRUE)
      df_out$max_Rhat <- apply(rhat_mat, 2, max, na.rm = TRUE)
    } else {
      df_out$min_neff <- res_p$n_eff[seq_len(length(groups))]
      df_out$max_Rhat <- res_p$Rhat[seq_len(length(groups))]
    }
  }
  if (add_phi) {
    res_phi <- as.data.frame(summary(fit, pars = "phi_ind")$summary)
    s_phi <- rstan::extract(fit, pars = "phi_ind")[[1]]
    phi_row <- data.frame(
      Source = "phi_ind", mean = mean(s_phi), sd = sd(s_phi),
      q025 = quantile(s_phi, 0.025), q25 = quantile(s_phi, 0.25),
      median = median(s_phi), q75 = quantile(s_phi, 0.75),
      q975 = quantile(s_phi, 0.975),
      W50 = quantile(s_phi, 0.75) - quantile(s_phi, 0.25),
      check.names = FALSE
    )
    if ("n_eff" %in% colnames(res_phi)) {
      phi_row$min_neff <- res_phi$n_eff
      phi_row$max_Rhat <- res_phi$Rhat
    }
    df_out <- rbind(df_out, phi_row)
  }
  num_cols <- sapply(df_out, is.numeric)
  df_out[num_cols] <- lapply(df_out[num_cols], function(x) round(x, 3))
  df_out
}


#==============================================================================
# 2. STAN CODE GENERATOR - WITH GHOST SOURCE (G+1) & ORTHOGONALIZED BIAS
#==============================================================================

#==============================================================================
# 2. STAN CODE GENERATOR - WITH GHOST DETECTOR, ORTHO-BIAS & HOTSPOT MODES
#==============================================================================

generate_stan_code <- function(
    pool_mode = "no_pool",
    b = FALSE,
    b_src = FALSE,
    cov_mode = "mix",
    coda = FALSE,
    has_groups = FALSE,
    bias = FALSE,
    ghost = FALSE,
    hotspot_mode = "mean"
) {
  cov_mix <- cov_mode %in% c("mix", "full")
  cov_src <- cov_mode == "full"
  use_levels  <- isTRUE(has_groups && pool_mode != "no_pool")
  use_b_src     <- isTRUE(b && b_src && !coda)
  use_k_scale   <- isTRUE(hotspot_mode == "k_scale")
  use_proc_only <- isTRUE(hotspot_mode == "proc_only")
  use_nested    <- isTRUE(hotspot_mode == "nested")

  n_p_str <- if (ghost) "N_groups + 1" else "N_groups"
  w_str   <- if (ghost) "conn_weights_ext" else "conn_weights"

  #--------------------------------------------------------------------------
  # DATA BLOCK
  #--------------------------------------------------------------------------
  data_blk <- "data {
  int<lower=1> N_mix;
  int<lower=1> N_sources;
  int<lower=1> N_tracers;
  int<lower=1> N_groups;
  matrix[N_mix, N_tracers] y_mix;
  matrix[N_sources, N_tracers] y_source;
  int<lower=1, upper=N_groups> group_id[N_sources];
  vector[N_groups] conn_weights;
  vector[N_mix] ssa_mix;
  vector[N_sources] ssa_source;
  vector[N_groups] ssa_group;
  vector[N_tracers] y_center;
  vector<lower=0>[N_tracers] y_scale;
  int<lower=0, upper=N_mix> holdout_id;
"

  if (ghost) {
    data_blk <- paste0(
      data_blk,
      "  vector[N_tracers] ghost_anchor;\n",
      "  real<lower=0> alpha_ghost;\n"
    )
  }

  if (use_k_scale) {
    data_blk <- paste0(
      data_blk,
      "  real<lower=0.1> k_hotspot;\n"
    )
  }

  if (use_nested) {
    data_blk <- paste0(
      data_blk,
      "  real<lower=0.1> sigma_hotspot;\n"
    )
  }

  if (use_levels) {
    data_blk <- paste0(
      data_blk,
      "  int<lower=1> N_mix_levels;\n",
      "  int<lower=1, upper=N_mix_levels> mix_level_id[N_mix];\n"
    )
  }

  if (coda) {
    data_blk <- paste0(
      data_blk,
      "  int<lower=1> N_tracers_raw;\n",
      "  matrix[N_tracers_raw, N_tracers] V;\n",
      "  vector<lower=0>[N_groups] src_totals;\n"
    )
  }

  data_blk <- paste0(data_blk, "}\n")

  #--------------------------------------------------------------------------
  # TRANSFORMED DATA BLOCK
  #--------------------------------------------------------------------------
  tdata_blk <- ""
  if (ghost) {
    tdata_blk <- "transformed data {
  vector[N_groups + 1] conn_weights_ext;
  for (g in 1:N_groups) conn_weights_ext[g] = conn_weights[g];
  conn_weights_ext[N_groups + 1] = alpha_ghost;
}
"
  }

  #--------------------------------------------------------------------------
  # PARAMETERS BLOCK
  #--------------------------------------------------------------------------
  param_blk <- "parameters {\n"

  if (pool_mode == "hier") {
    param_blk <- paste0(
      param_blk,
      sprintf("  simplex[%s] P_global;\n", n_p_str),
      sprintf("  simplex[%s] P_ind[N_mix];\n", n_p_str),
      "  real<lower=0> phi_ind;\n"
    )

    if (use_levels) {
      param_blk <- paste0(
        param_blk,
        sprintf("  simplex[%s] P_level[N_mix_levels];\n", n_p_str),
        "  real<lower=0> phi_level;\n"
      )
    }

  } else if (pool_mode == "no_pool") {
    param_blk <- paste0(
      param_blk,
      sprintf("  simplex[%s] P_ind[N_mix];\n", n_p_str)
    )

  } else {
    if (use_levels) {
      param_blk <- paste0(
        param_blk,
        sprintf("  simplex[%s] P_level[N_mix_levels];\n", n_p_str)
      )
    } else {
      param_blk <- paste0(
        param_blk,
        sprintf("  simplex[%s] P;\n", n_p_str)
      )
    }
  }

  param_blk <- paste0(param_blk, "  vector[N_tracers] mu_source[N_groups];\n")

  if (ghost) {
    param_blk <- paste0(param_blk, "  vector[N_tracers] mu_ghost;\n")
  }

  if (use_nested) {
    param_blk <- paste0(param_blk, "  matrix[N_mix, N_sources] eta_sub;\n")
  }

  if (cov_src) {
    param_blk <- paste0(
      param_blk,
      "  cholesky_factor_corr[N_tracers] L_Omega_src[N_groups];\n",
      "  vector<lower=0>[N_tracers] L_sigma_src[N_groups];\n"
    )
  } else {
    param_blk <- paste0(
      param_blk,
      "  vector<lower=0>[N_tracers] sigma_source[N_groups];\n"
    )
  }

  if (b) {
    if (coda) {
      param_blk <- paste0(param_blk, "  vector[N_tracers_raw] beta_sorting;\n")
    } else {
      param_blk <- paste0(param_blk, "  vector[N_tracers] beta_sorting;\n")
    }
  }

  if (cov_mix) {
    param_blk <- paste0(
      param_blk,
      "  cholesky_factor_corr[N_tracers] L_Omega;\n"
    )
    if (!use_proc_only) {
      param_blk <- paste0(
        param_blk,
        "  vector<lower=0>[N_tracers] L_sigma;\n"
      )
    }
  } else if (!use_proc_only) {
    param_blk <- paste0(param_blk, "  vector<lower=0>[N_tracers] sigma_mix;\n")
  }

  if (bias) {
    param_blk <- paste0(param_blk, "  vector[N_tracers] delta_bias;\n")
  }

  param_blk <- paste0(param_blk, "}\n")


  #--------------------------------------------------------------------------
  # MODEL BLOCK
  #--------------------------------------------------------------------------
  model_blk <- "model {\n"

  # 0. Объявления локальных переменных в начале блока model
  if (cov_src) {
    model_blk <- paste0(
      model_blk,
      "  matrix[N_tracers, N_tracers] L_cov_src[N_groups];\n"
    )
  }
  if (cov_mix && !use_k_scale && !use_proc_only) {
    model_blk <- paste0(
      model_blk,
      "  matrix[N_tracers, N_tracers] L_cov_mix;\n"
    )
  }
  if (coda) {
    if (!use_nested) {
      model_blk <- paste0(
        model_blk,
        "  vector[N_tracers_raw] src_prop_matrix[N_groups];\n"
      )
    }
    if (ghost) {
      model_blk <- paste0(
        model_blk,
        "  vector[N_tracers_raw] ghost_prop_vec;\n",
        "  real ghost_total = mean(src_totals);\n"
      )
    }
  }
  if (bias) {
    model_blk <- paste0(
      model_blk,
      "  vector[N_tracers] delta_eff;\n"
    )
  }

  #--- 1. Proportion pooling priors ------------------------------------------
  if (pool_mode == "hier") {
    model_blk <- paste0(
      model_blk,
      "  phi_ind ~ gamma(2, 0.1);\n",
      sprintf("  P_global ~ dirichlet(%s);\n", w_str)
    )

    if (use_levels) {
      model_blk <- paste0(
        model_blk,
        "  phi_level ~ gamma(2, 0.1);\n",
        "  for (L in 1:N_mix_levels) P_level[L] ~ dirichlet(P_global * phi_level);\n"
      )
    }

  } else if (pool_mode == "no_pool") {
    model_blk <- paste0(
      model_blk,
      sprintf("  for (i in 1:N_mix) P_ind[i] ~ dirichlet(%s);\n", w_str)
    )

  } else {
    if (use_levels) {
      model_blk <- paste0(
        model_blk,
        sprintf("  for (L in 1:N_mix_levels) P_level[L] ~ dirichlet(%s);\n", w_str)
      )
    } else {
      model_blk <- paste0(model_blk, sprintf("  P ~ dirichlet(%s);\n", w_str))
    }
  }

  #--- 2. Other priors -------------------------------------------------------
  if (ghost) {
    model_blk <- paste0(
      model_blk,
      "  mu_ghost ~ normal(ghost_anchor, 0.35);\n"
    )
  }

  if (use_nested) {
    model_blk <- paste0(
      model_blk,
      "  to_vector(eta_sub) ~ normal(0, sigma_hotspot);\n"
    )
  }

  if (b) {
    model_blk <- paste0(model_blk, "  beta_sorting ~ normal(0, 0.2);\n")
  }

  if (cov_mix) {
    model_blk <- paste0(
      model_blk,
      "  L_Omega ~ lkj_corr_cholesky(4);\n"
    )
    if (!use_proc_only) {
      model_blk <- paste0(
        model_blk,
        "  L_sigma ~ student_t(3, 0, 1);\n"
      )
    }
    if (!use_k_scale && !use_proc_only) {
      model_blk <- paste0(
        model_blk,
        "  L_cov_mix = diag_pre_multiply(L_sigma, L_Omega);\n"
      )
    }
  } else if (!use_proc_only) {
    model_blk <- paste0(model_blk, "  sigma_mix ~ student_t(3, 0, 1);\n")
  }

  if (bias) {
    model_blk <- paste0(model_blk, "  delta_bias ~ student_t(3, 0, 0.1);\n")
  }

  #--- 3. Source model (Known sources 1..N_groups) ---------------------------
  if (cov_src) {
    model_blk <- paste0(
      model_blk,
      "  for (g in 1:N_groups) {\n",
      "    L_Omega_src[g] ~ lkj_corr_cholesky(2);\n",
      "    L_sigma_src[g] ~ student_t(3, 0, 1);\n",
      "    mu_source[g] ~ normal(0, 5);\n",
      "    L_cov_src[g] = diag_pre_multiply(L_sigma_src[g], L_Omega_src[g]);\n",
      "  }\n"
    )

    if (use_b_src) {
      model_blk <- paste0(
        model_blk,
        "  for (i in 1:N_sources) {\n",
        "    int g_id = group_id[i];\n",
        "    vector[N_tracers] mu_src_i;\n",
        "    for (t in 1:N_tracers) {\n",
        "      real raw_gt = fmax(1e-6, y_center[t] + y_scale[t] * mu_source[g_id, t]);\n",
        "      real corr_gt = raw_gt * exp(beta_sorting[t] * (ssa_source[i] - ssa_group[g_id]));\n",
        "      mu_src_i[t] = (corr_gt - y_center[t]) / y_scale[t];\n",
        "    }\n",
        "    y_source[i]' ~ multi_normal_cholesky(mu_src_i, L_cov_src[g_id]);\n",
        "  }\n"
      )
    } else {
      model_blk <- paste0(
        model_blk,
        "  for (i in 1:N_sources) {\n",
        "    y_source[i]' ~ multi_normal_cholesky(mu_source[group_id[i]], L_cov_src[group_id[i]]);\n",
        "  }\n"
      )
    }

  } else {
    model_blk <- paste0(
      model_blk,
      "  for (g in 1:N_groups) {\n",
      "    mu_source[g] ~ normal(0, 5);\n",
      "    sigma_source[g] ~ student_t(3, 0, 1);\n",
      "  }\n"
    )

    if (use_b_src) {
      model_blk <- paste0(
        model_blk,
        "  for (i in 1:N_sources) {\n",
        "    int g_id = group_id[i];\n",
        "    vector[N_tracers] mu_src_i;\n",
        "    for (t in 1:N_tracers) {\n",
        "      real raw_gt = fmax(1e-6, y_center[t] + y_scale[t] * mu_source[g_id, t]);\n",
        "      real corr_gt = raw_gt * exp(beta_sorting[t] * (ssa_source[i] - ssa_group[g_id]));\n",
        "      mu_src_i[t] = (corr_gt - y_center[t]) / y_scale[t];\n",
        "    }\n",
        "    y_source[i]' ~ normal(mu_src_i, sigma_source[g_id]);\n",
        "  }\n"
      )
    } else {
      model_blk <- paste0(
        model_blk,
        "  for (i in 1:N_sources) {\n",
        "    y_source[i]' ~ normal(mu_source[group_id[i]], sigma_source[group_id[i]]);\n",
        "  }\n"
      )
    }
  }

  if (coda) {
    if (!use_nested) {
      model_blk <- paste0(
        model_blk,
        "  for (g in 1:N_groups) src_prop_matrix[g] = softmax(V * mu_source[g]);\n"
      )
    }
    if (ghost) {
      model_blk <- paste0(
        model_blk,
        "  ghost_prop_vec = softmax(V * mu_ghost);\n"
      )
    }
  }

  if (bias) {
    if (ghost) {
      model_blk <- paste0(
        model_blk,
        "  {\n",
        "    vector[N_tracers] mu_bar = rep_vector(0.0, N_tracers);\n",
        "    vector[N_tracers] v_gh;\n",
        "    real v_gh_sq;\n",
        "    for (g in 1:N_groups) mu_bar += mu_source[g] / N_groups;\n",
        "    v_gh = mu_ghost - mu_bar;\n",
        "    v_gh_sq = dot_self(v_gh) + 1e-9;\n",
        "    delta_eff = delta_bias - (dot_product(delta_bias, v_gh) / v_gh_sq) * v_gh;\n",
        "  }\n"
      )
    } else {
      model_blk <- paste0(
        model_blk,
        "  delta_eff = delta_bias;\n"
      )
    }
  }

  #--- 4. Mixture likelihood loop --------------------------------------------
  model_blk <- paste0(model_blk, "\n  for (i in 1:N_mix) {\n")

  if (pool_mode %in% c("hier", "no_pool")) {
    current_p_g     <- "P_ind[i][g]"
    current_p_ghost <- "P_ind[i][N_groups + 1]"
  } else if (use_levels) {
    current_p_g     <- "P_level[mix_level_id[i]][g]"
    current_p_ghost <- "P_level[mix_level_id[i]][N_groups + 1]"
  } else {
    current_p_g     <- "P[g]"
    current_p_ghost <- "P[N_groups + 1]"
  }

  # 1. ВСЕ ОБЪЯВЛЕНИЯ ЛОКАЛЬНЫХ ПЕРЕМЕННЫХ В НАЧАЛЕ ЦИКЛА
  model_blk <- paste0(
    model_blk,
    "    vector[N_tracers] mu_f;\n",
    "    vector[N_tracers] y_f = y_mix[i]';\n"
  )

  if (use_nested) {
    model_blk <- paste0(
      model_blk,
      "    matrix[N_groups, N_tracers] mu_eff_i = rep_matrix(0.0, N_groups, N_tracers);\n"
    )
  }

  if (coda) {
    model_blk <- paste0(
      model_blk,
      "    vector[N_tracers_raw] mu_pred_prop;\n"
    )
    if (use_nested) {
      model_blk <- paste0(
        model_blk,
        "    vector[N_tracers_raw] src_prop_matrix[N_groups];\n"
      )
    }
  }

  # 2. ИСПОЛНЯЕМЫЕ ОПЕРАТОРЫ (СЭМПЛИРОВАНИЕ И ВЫЧИСЛЕНИЯ)
  if (pool_mode == "hier") {
    if (use_levels) {
      model_blk <- paste0(
        model_blk,
        "    P_ind[i] ~ dirichlet(P_level[mix_level_id[i]] * phi_ind);\n"
      )
    } else {
      model_blk <- paste0(
        model_blk,
        "    P_ind[i] ~ dirichlet(P_global * phi_ind);\n"
      )
    }
  }

  if (use_nested) {
    model_blk <- paste0(
      model_blk,
      "    {\n",
      "      vector[N_groups] max_eta_g = rep_vector(-1e9, N_groups);\n",
      "      vector[N_groups] sum_exp_g = rep_vector(0.0, N_groups);\n",
      "      vector[N_sources] w_sub;\n",
      "      for (s in 1:N_sources) {\n",
      "        int g_s = group_id[s];\n",
      "        if (eta_sub[i, s] > max_eta_g[g_s]) max_eta_g[g_s] = eta_sub[i, s];\n",
      "      }\n",
      "      for (s in 1:N_sources) {\n",
      "        int g_s = group_id[s];\n",
      "        w_sub[s] = exp(eta_sub[i, s] - max_eta_g[g_s]);\n",
      "        sum_exp_g[g_s] += w_sub[s];\n",
      "      }\n",
      "      for (s in 1:N_sources) {\n",
      "        int g_s = group_id[s];\n",
      "        real w_norm = w_sub[s] / sum_exp_g[g_s];\n",
      "        for (t in 1:N_tracers) mu_eff_i[g_s, t] += w_norm * y_source[s, t];\n",
      "      }\n",
      "    }\n"
    )
  }

  active_mu_gt <- if (use_nested) "mu_eff_i[g, t]" else "mu_source[g, t]"

  if (coda) {
    if (use_nested) {
      model_blk <- paste0(
        model_blk,
        "    for (g in 1:N_groups) src_prop_matrix[g] = softmax(V * to_vector(mu_eff_i[g]));\n"
      )
    }
    model_blk <- paste0(
      model_blk,
      "    for (t in 1:N_tracers_raw) {\n",
      "      real val = 0;\n",
      "      for (g in 1:N_groups) {\n"
    )

    if (b) {
      model_blk <- paste0(
        model_blk,
        sprintf(
          "        val += %s * src_totals[g] * src_prop_matrix[g, t] * exp(beta_sorting[t] * (ssa_mix[i] - ssa_group[g]));\n",
          current_p_g
        )
      )
    } else {
      model_blk <- paste0(
        model_blk,
        sprintf(
          "        val += %s * src_totals[g] * src_prop_matrix[g, t];\n",
          current_p_g
        )
      )
    }

    model_blk <- paste0(model_blk, "      }\n")

    if (ghost) {
      if (b) {
        model_blk <- paste0(
          model_blk,
          sprintf(
            "      val += %s * ghost_total * ghost_prop_vec[t] * exp(beta_sorting[t] * ssa_mix[i]);\n",
            current_p_ghost
          )
        )
      } else {
        model_blk <- paste0(
          model_blk,
          sprintf(
            "      val += %s * ghost_total * ghost_prop_vec[t];\n",
            current_p_ghost
          )
        )
      }
    }

    model_blk <- paste0(
      model_blk,
      "      mu_pred_prop[t] = val + 1e-9;\n",
      "    }\n"
    )

    if (bias) {
      model_blk <- paste0(model_blk, "    mu_f = (log(mu_pred_prop)' * V)' + delta_eff;\n")
    } else {
      model_blk <- paste0(model_blk, "    mu_f = (log(mu_pred_prop)' * V)';\n")
    }

  } else {
    model_blk <- paste0(
      model_blk,
      "    for (t in 1:N_tracers) {\n",
      "      real val = 0;\n"
    )

    if (b) {
      model_blk <- paste0(
        model_blk,
        "      for (g in 1:N_groups) {\n",
        sprintf("        real raw_gt = fmax(1e-6, y_center[t] + y_scale[t] * %s);\n", active_mu_gt),
        sprintf(
          "        val += %s * raw_gt * exp(beta_sorting[t] * (ssa_mix[i] - ssa_group[g]));\n",
          current_p_g
        ),
        "      }\n"
      )
      if (ghost) {
        model_blk <- paste0(
          model_blk,
          "      {\n",
          "        real raw_gh = fmax(1e-6, y_center[t] + y_scale[t] * mu_ghost[t]);\n",
          sprintf(
            "        val += %s * raw_gh * exp(beta_sorting[t] * ssa_mix[i]);\n",
            current_p_ghost
          ),
          "      }\n"
        )
      }
      model_blk <- paste0(model_blk, "      val = (val - y_center[t]) / y_scale[t];\n")
    } else {
      model_blk <- paste0(
        model_blk,
        sprintf(
          "      for (g in 1:N_groups) val += %s * %s;\n",
          current_p_g, active_mu_gt
        )
      )
      if (ghost) {
        model_blk <- paste0(
          model_blk,
          sprintf(
            "      val += %s * mu_ghost[t];\n",
            current_p_ghost
          )
        )
      }
    }

    if (bias) {
      model_blk <- paste0(model_blk, "      mu_f[t] = val + delta_eff[t];\n    }\n")
    } else {
      model_blk <- paste0(model_blk, "      mu_f[t] = val;\n    }\n")
    }
  }

  # Правдоподобие смеси с учетом (или без учета) дисперсии процесса / маргинализации k-хотспотов
  if (use_k_scale || use_proc_only) {
    base_sig_str <- if (cov_mix) "L_sigma[t]" else "sigma_mix[t]"
    src_sig_str  <- if (cov_src) "L_sigma_src[g, t]" else "sigma_source[g, t]"
    init_var_str <- if (use_proc_only) "1e-6" else sprintf("square(%s)", base_sig_str)
    div_k_str    <- if (use_proc_only) "" else " / k_hotspot"

    model_blk <- paste0(
      model_blk,
      "    {\n",
      "      vector[N_tracers] sigma_tot_i;\n",
      "      for (t in 1:N_tracers) {\n",
      sprintf("        real var_t = %s;\n", init_var_str),
      sprintf("        for (g in 1:N_groups) var_t += square(%s) * square(%s)%s;\n", current_p_g, src_sig_str, div_k_str),
      "        sigma_tot_i[t] = sqrt(var_t);\n",
      "      }\n"
    )
    if (cov_mix) {
      model_blk <- paste0(
        model_blk,
        "      if (i != holdout_id) y_f ~ multi_student_t_cholesky(4, mu_f, diag_pre_multiply(sigma_tot_i, L_Omega));\n"
      )
    } else {
      model_blk <- paste0(
        model_blk,
        "      if (i != holdout_id) y_f ~ student_t(4, mu_f, sigma_tot_i);\n"
      )
    }
    model_blk <- paste0(model_blk, "    }\n")
  } else {
    if (cov_mix) {
      model_blk <- paste0(
        model_blk,
        "    if (i != holdout_id) y_f ~ multi_student_t_cholesky(4, mu_f, L_cov_mix);\n"
      )
    } else {
      model_blk <- paste0(
        model_blk,
        "    if (i != holdout_id) y_f ~ student_t(4, mu_f, sigma_mix);\n"
      )
    }
  }

  model_blk <- paste0(model_blk, "  }\n}\n")


  #--------------------------------------------------------------------------
  # GENERATED QUANTITIES BLOCK
  #--------------------------------------------------------------------------
  gen_blk <- "generated quantities {
  vector[N_mix] log_lik;
  {
"
  if (cov_mix && !use_k_scale && !use_proc_only) {
    gen_blk <- paste0(
      gen_blk,
      "    matrix[N_tracers, N_tracers] L_cov_mix = diag_pre_multiply(L_sigma, L_Omega);\n"
    )
  }
  if (coda) {
    if (!use_nested) {
      gen_blk <- paste0(
        gen_blk,
        "    vector[N_tracers_raw] src_prop_matrix[N_groups];\n"
      )
    }
    if (ghost) {
      gen_blk <- paste0(
        gen_blk,
        "    vector[N_tracers_raw] ghost_prop_vec = softmax(V * mu_ghost);\n",
        "    real ghost_total = mean(src_totals);\n"
      )
    }
  }
  if (bias) {
    gen_blk <- paste0(
      gen_blk,
      "    vector[N_tracers] delta_eff;\n"
    )
  }

  if (coda && !use_nested) {
    gen_blk <- paste0(
      gen_blk,
      "    for (g in 1:N_groups) src_prop_matrix[g] = softmax(V * mu_source[g]);\n"
    )
  }

  if (bias) {
    if (ghost) {
      gen_blk <- paste0(
        gen_blk,
        "    {\n",
        "      vector[N_tracers] mu_bar = rep_vector(0.0, N_tracers);\n",
        "      vector[N_tracers] v_gh;\n",
        "      real v_gh_sq;\n",
        "      for (g in 1:N_groups) mu_bar += mu_source[g] / N_groups;\n",
        "      v_gh = mu_ghost - mu_bar;\n",
        "      v_gh_sq = dot_self(v_gh) + 1e-9;\n",
        "      delta_eff = delta_bias - (dot_product(delta_bias, v_gh) / v_gh_sq) * v_gh;\n",
        "    }\n"
      )
    } else {
      gen_blk <- paste0(
        gen_blk,
        "    delta_eff = delta_bias;\n"
      )
    }
  }

  gen_blk <- paste0(
    gen_blk,
    "    for (i in 1:N_mix) {\n",
    "      vector[N_tracers] mu_f;\n",
    "      vector[N_tracers] y_f = y_mix[i]';\n"
  )

  if (use_nested) {
    gen_blk <- paste0(
      gen_blk,
      "      matrix[N_groups, N_tracers] mu_eff_i = rep_matrix(0.0, N_groups, N_tracers);\n"
    )
  }

  if (coda) {
    gen_blk <- paste0(
      gen_blk,
      "      vector[N_tracers_raw] mu_pred_prop;\n"
    )
    if (use_nested) {
      gen_blk <- paste0(
        gen_blk,
        "      vector[N_tracers_raw] src_prop_matrix[N_groups];\n"
      )
    }
  }

  if (use_nested) {
    gen_blk <- paste0(
      gen_blk,
      "      {\n",
      "        vector[N_groups] max_eta_g = rep_vector(-1e9, N_groups);\n",
      "        vector[N_groups] sum_exp_g = rep_vector(0.0, N_groups);\n",
      "        vector[N_sources] w_sub;\n",
      "        for (s in 1:N_sources) {\n",
      "          int g_s = group_id[s];\n",
      "          if (eta_sub[i, s] > max_eta_g[g_s]) max_eta_g[g_s] = eta_sub[i, s];\n",
      "        }\n",
      "        for (s in 1:N_sources) {\n",
      "          int g_s = group_id[s];\n",
      "          w_sub[s] = exp(eta_sub[i, s] - max_eta_g[g_s]);\n",
      "          sum_exp_g[g_s] += w_sub[s];\n",
      "        }\n",
      "        for (s in 1:N_sources) {\n",
      "          int g_s = group_id[s];\n",
      "          real w_norm = w_sub[s] / sum_exp_g[g_s];\n",
      "          for (t in 1:N_tracers) mu_eff_i[g_s, t] += w_norm * y_source[s, t];\n",
      "        }\n",
      "      }\n"
    )
  }

  if (coda) {
    if (use_nested) {
      gen_blk <- paste0(
        gen_blk,
        "      for (g in 1:N_groups) src_prop_matrix[g] = softmax(V * to_vector(mu_eff_i[g]));\n"
      )
    }
    gen_blk <- paste0(
      gen_blk,
      "      for (t in 1:N_tracers_raw) {\n",
      "        real val = 0;\n",
      "        for (g in 1:N_groups) {\n"
    )

    if (b) {
      gen_blk <- paste0(
        gen_blk,
        sprintf(
          "          val += %s * src_totals[g] * src_prop_matrix[g, t] * exp(beta_sorting[t] * (ssa_mix[i] - ssa_group[g]));\n",
          current_p_g
        )
      )
    } else {
      gen_blk <- paste0(
        gen_blk,
        sprintf(
          "          val += %s * src_totals[g] * src_prop_matrix[g, t];\n",
          current_p_g
        )
      )
    }

    gen_blk <- paste0(gen_blk, "        }\n")

    if (ghost) {
      if (b) {
        gen_blk <- paste0(
          gen_blk,
          sprintf(
            "        val += %s * ghost_total * ghost_prop_vec[t] * exp(beta_sorting[t] * ssa_mix[i]);\n",
            current_p_ghost
          )
        )
      } else {
        gen_blk <- paste0(
          gen_blk,
          sprintf(
            "        val += %s * ghost_total * ghost_prop_vec[t];\n",
            current_p_ghost
          )
        )
      }
    }

    gen_blk <- paste0(
      gen_blk,
      "        mu_pred_prop[t] = val + 1e-9;\n",
      "      }\n"
    )

    if (bias) {
      gen_blk <- paste0(gen_blk, "      mu_f = (log(mu_pred_prop)' * V)' + delta_eff;\n")
    } else {
      gen_blk <- paste0(gen_blk, "      mu_f = (log(mu_pred_prop)' * V)';\n")
    }

  } else {
    gen_blk <- paste0(
      gen_blk,
      "      for (t in 1:N_tracers) {\n",
      "        real val = 0;\n"
    )

    if (b) {
      gen_blk <- paste0(
        gen_blk,
        "        for (g in 1:N_groups) {\n",
        sprintf("          real raw_gt = fmax(1e-6, y_center[t] + y_scale[t] * %s);\n", active_mu_gt),
        sprintf(
          "          val += %s * raw_gt * exp(beta_sorting[t] * (ssa_mix[i] - ssa_group[g]));\n",
          current_p_g
        ),
        "        }\n"
      )
      if (ghost) {
        gen_blk <- paste0(
          gen_blk,
          "        {\n",
          "          real raw_gh = fmax(1e-6, y_center[t] + y_scale[t] * mu_ghost[t]);\n",
          sprintf(
            "          val += %s * raw_gh * exp(beta_sorting[t] * ssa_mix[i]);\n",
            current_p_ghost
          ),
          "        }\n"
        )
      }
      gen_blk <- paste0(gen_blk, "        val = (val - y_center[t]) / y_scale[t];\n")
    } else {
      gen_blk <- paste0(
        gen_blk,
        sprintf(
          "        for (g in 1:N_groups) val += %s * %s;\n",
          current_p_g, active_mu_gt
        )
      )
      if (ghost) {
        gen_blk <- paste0(
          gen_blk,
          sprintf(
            "        val += %s * mu_ghost[t];\n",
            current_p_ghost
          )
        )
      }
    }

    if (bias) {
      gen_blk <- paste0(gen_blk, "        mu_f[t] = val + delta_eff[t];\n      }\n")
    } else {
      gen_blk <- paste0(gen_blk, "        mu_f[t] = val;\n      }\n")
    }
  }

  if (use_k_scale || use_proc_only) {
    base_sig_str <- if (cov_mix) "L_sigma[t]" else "sigma_mix[t]"
    src_sig_str  <- if (cov_src) "L_sigma_src[g, t]" else "sigma_source[g, t]"
    init_var_str <- if (use_proc_only) "1e-6" else sprintf("square(%s)", base_sig_str)
    div_k_str    <- if (use_proc_only) "" else " / k_hotspot"

    gen_blk <- paste0(
      gen_blk,
      "      {\n",
      "        vector[N_tracers] sigma_tot_i;\n",
      "        for (t in 1:N_tracers) {\n",
      sprintf("          real var_t = %s;\n", init_var_str),
      sprintf("          for (g in 1:N_groups) var_t += square(%s) * square(%s)%s;\n", current_p_g, src_sig_str, div_k_str),
      "          sigma_tot_i[t] = sqrt(var_t);\n",
      "        }\n"
    )
    if (cov_mix) {
      gen_blk <- paste0(
        gen_blk,
        "        log_lik[i] = multi_student_t_cholesky_lpdf(y_f | 4, mu_f, diag_pre_multiply(sigma_tot_i, L_Omega));\n"
      )
    } else {
      gen_blk <- paste0(
        gen_blk,
        "        log_lik[i] = student_t_lpdf(y_f | 4, mu_f, sigma_tot_i);\n"
      )
    }
    gen_blk <- paste0(gen_blk, "      }\n")
  } else {
    if (cov_mix) {
      gen_blk <- paste0(
        gen_blk,
        "      log_lik[i] = multi_student_t_cholesky_lpdf(y_f | 4, mu_f, L_cov_mix);\n"
      )
    } else {
      gen_blk <- paste0(
        gen_blk,
        "      log_lik[i] = student_t_lpdf(y_f | 4, mu_f, sigma_mix);\n"
      )
    }
  }

  gen_blk <- paste0(gen_blk, "    }\n  }\n}\n")

  return(paste0(data_blk, tdata_blk, param_blk, model_blk, gen_blk))
}


#==============================================================================
# 3. USER INTERFACE (UI)
#==============================================================================

ui <- fluidPage(
  theme = shinytheme("flatly"),

  # --- ВСТАВИТЬ СРАЗУ ПОСЛЕ theme = shinytheme("flatly"), (строка 618 app_49.R) ---
  tags$head(
    tags$style(HTML("
      .exp-locked {
        opacity: 0.50 !important;
        pointer-events: none !important;
        cursor: not-allowed !important;
      }
      .exp-badge {
        font-size: 0.78em;
        color: #7f8c8d;
        font-style: italic;
        margin-left: 4px;
      }
    ")),
    tags$script(HTML("
      $(document).on('shiny:connected', function() {
        // 1. Блокируем чекбокс Ghost Source (G+1) и чекбокс opt_b_src
        $('#opt_ghost, #opt_b_src').prop('disabled', true)
          .closest('.checkbox').addClass('exp-locked')
          .find('span').append(' <span class=\"exp-badge\">[Experimental — disabled]</span>');

        // 2. Блокируем только пункты 3 (k_scale) и 4 (nested) в переключателе hotspot_mode,
        //    оставляя пункты 1 (mean) и 2 (proc_only) полностью активными!
        $(\"input[name='src_rep_mode'][value='k_scale'], input[name='src_rep_mode'][value='nested']\")
          .prop('disabled', true)
          .closest('.radio').addClass('exp-locked')
          .find('span').append(' <span class=\"exp-badge\">[Experimental — disabled]</span>');
      });
    "))
  ),

  titlePanel(
    title = HTML("<b>BASIN:</b> Bayesian Sediment Source Apportionment Tool"),
    windowTitle = "BASIN: Bayesian Sediment Source Apportionment Tool"
  ),

  sidebarLayout(
    sidebarPanel(

      #--------------------------------------------------
      # 1. Data input
      #--------------------------------------------------
      h4("1. Data & Role Assignment"),
      fileInput("file", "Upload dataset (CSV):", accept = ".csv"),

      uiOutput("ui_group_col"),
      uiOutput("ui_mix_label"),
      uiOutput("ui_mix_cov"),
      uiOutput("ui_size_col"),
      uiOutput("ui_tracer_select"),

      hr(),

      #--------------------------------------------------
      # 2. Tracer Pre-Filtering
      #--------------------------------------------------
      h4("2. Tracer Pre-Filtering"),
      radioButtons(
        "filter_mode",
        "Range-Test Normalization (Sorting Proxy):",
        choices = c(
          "Raw only" = "raw",
          "Proxy-normalized" = "corr",
          "Strict (Both)" = "both",
          "Lenient (Either)" = "any"
        ),
        selected = "any",
        inline = TRUE
      ),

      checkboxInput(
        "opt_buddy",
        "1. Inter-Tracer Association Preservation Test (APT / GAPT)",
        FALSE
      ),

      checkboxInput(
        "opt_hull",
        "2. Multivariate PCA Convex-Hull Leverage Filter",
        FALSE
      ),

      actionButton(
        "run_filter",
        "Run Automated Tracer Pre-Filtering",
        class = "btn-warning",
        style = "width:100%"
      ),

      hr(),

      #--------------------------------------------------
      # 3. Model Constructor (4 Design Axes)
      #--------------------------------------------------
      h4("3. Model Constructor"),

      # Axis 3.1: Information Sharing
      selectInput(
        "pool_mode",
        "3.1. Cross-Sample Information Sharing (Pooling):",
        choices = c(
          "Hierarchical Partial Pooling (Ensemble Regularization)" = "hier",
          "Non-Hierarchical / Isolated Proportions (No Pooling)" = "no_pool",
          "Complete Pooling (Technical Replicates / Single Mean)" = "complete"
        ),
        selected = "hier"
      ),
      helpText(
        style = "font-size:0.85em;",
        "Use 'Hierarchical' for stationary ensembles; use 'No Pooling' for isolated samples (N_mix = 1) or event-dominated sediment cores."
      ),

      # Axis 3.2: Error-Scale & Covariance
      conditionalPanel(
        condition = "input.src_rep_mode == 'centroid'",
        radioButtons(
          "err_scale_mode",
          "3.2a. Residual Error-Scale Specification:",
          choices = c(
            "Free Regularized Scales (Adaptive M-Estimator, Default)" = "mean",
            "Process-Error Only (Source-Variance Propagated)" = "proc_only"
          ),
          selected = "mean"
        )
      ),

      radioButtons(
        "cov_mode",
        "3.2b. Inter-Tracer Covariance Architecture:",
        choices = c(
          "Independent (Diagonal Error Matrix)" = "none",
          "Target Mixture Covariance (Recommended for Ensembles)" = "mix",
          "Full Bayes (Source + Mixture Covariance, Ideal for N_mix = 1)" = "full"
        ),
        selected = "mix"
      ),

      # Axis 3.3: Source End-Member Representation
      radioButtons(
        "src_rep_mode",
        "3.3. Source End-Member Representation:",
        choices = c(
          "Parametric Group Centroids (Uniform Source Mobilization)" = "centroid",
          "Sub-Sampled Erosion Hotspots (Hybrid: Free Scale + SD^2 / k)" = "k_scale",
          "Discrete Sample Simplex (Non-Parametric Nested Convex Hull)" = "nested"
        ),
        selected = "centroid"
      ),

      conditionalPanel(
        condition = "input.src_rep_mode == 'k_scale'",
        numericInput(
          "k_hotspot",
          "Active hotspots / sub-samples per source (k):",
          value = 3, min = 1, max = 30, step = 1
        ),
        helpText(
          style = "margin-bottom:0; font-size:0.85em;",
          "Propagates intra-source variance SD^2 / k into the mixture likelihood."
        )
      ),

      conditionalPanel(
        condition = "input.src_rep_mode == 'nested'",
        sliderInput(
          "sigma_hotspot",
          "Hotspot Selectivity (Logit SD):",
          min = 0.2, max = 3.0, value = 1.5, step = 0.1
        ),
        helpText(
          style = "margin-bottom:0; font-size:0.85em;",
          "Models each mixture as a convex combination of individual source samples within each group."
        )
      ),

      # Скрытый переключатель hotspot_mode для 100% совместимости со всем серверным кодом
      conditionalPanel(
        condition = "false",
        radioButtons(
          "hotspot_mode", NULL,
          choices = c("mean", "proc_only", "k_scale", "nested"),
          selected = "mean"
        )
      ),

      # Axis 3.4: Non-Conservatism & Unsampled-Source Buffers
      tags$label("3.4. Non-Conservatism & Compositional Modules:"),
      fluidRow(
        column(6, checkboxInput("opt_b", "Sorting Correction (Beta)", FALSE)),
        column(6, checkboxInput("opt_bias", "Systematic Bias Absorption", FALSE))
      ),

      conditionalPanel(
        condition = "input.opt_b == true && input.opt_coda == false",
        checkboxInput(
          "opt_b_src",
          "Anchor Beta on intra-source proxy variance",
          FALSE
        )
      ),

      fluidRow(
        column(12, checkboxInput("opt_coda", "Compositional Transform (CoDA / ILR)", FALSE))
      ),

      tags$div(
        style = "background-color: #eef6ff; padding: 10px; border-left: 4px solid #2c3e50; border-radius: 4px; margin-bottom: 12px;",
        checkboxInput(
          "opt_ghost",
          HTML("<b>Unsampled 'Ghost' Source Detector (G+1)</b>"),
          FALSE
        ),
        conditionalPanel(
          condition = "input.opt_ghost == true",
          sliderInput(
            "ghost_alpha",
            "Ghost Sparsity Prior (Dirichlet α_ghost):",
            min = 0.02,
            max = 0.50,
            value = 0.12,
            step = 0.02
          ),
          helpText(
            style = "margin-bottom:0; font-size:0.85em;",
            "Adds a latent (G+1)-th end-member anchored to out-of-simplex residuals."
          )
        )
      ),

      hr(),

      #--------------------------------------------------
      # 4. Inference Engine Settings
      #--------------------------------------------------
      h4("4. Bayesian Inference Engine"),
      selectInput(
        "calc_method",
        "Sampling Algorithm:",
        choices = c(
          "Hamiltonian Monte Carlo (Full HMC / NUTS)" = "mcmc",
          "Variational Inference (Fast ADVI Approximation)" = "vb"
        )
      ),

      fluidRow(
        column(6, numericInput("iter", "Iterations:", 2000)),
        column(6, numericInput("chains", "Chains:", 4))
      ),

      hr(),

      textOutput("count_status"),
      br(),

      actionButton(
        "run_model",
        "5. Compile & Run Unmixing Model",
        class = "btn-primary",
        style = "width:100%"
      ),

      hr(),

      #--------------------------------------------------
      # Geographic / connectivity weights
      #--------------------------------------------------
      uiOutput("conn_inputs")
    ),

    mainPanel(
      tabsetPanel(

        #------------------------------------------------
        # Tab 1: Raw data
        #------------------------------------------------
        tabPanel(
          "Data table",
          DTOutput("raw_table")
        ),

        #------------------------------------------------
        # Tab 2: Tracer selection and plots
        #------------------------------------------------
        tabPanel(
          "Tracer selection",
          fluidRow(
            column(
              5,
              h4("1. Final choice"),
              helpText("Check or uncheck variables before starting the model."),
              uiOutput("ui_final_tracers_chk"),
              hr(),
              h4("2. Filtering report"),
              verbatimTextOutput("filter_log", placeholder = TRUE),
              tags$style(
                type = "text/css",
                "#filter_log { max-height: 300px; overflow-y: scroll; white-space: pre-wrap; }"
              )
            ),
            column(
              7,
              h4("3. Multi-proxy Discriminant Visualization"),
              tabsetPanel(
                tabPanel(
                  "PCA (Variance Explorer)",
                  plotOutput("pca_plot", height = "500px")
                ),
                tabPanel(
                  "LDA (Source Separation)",
                  plotOutput("lda_plot", height = "500px"),
                  helpText(
                    "LDA maximizes the distance between source groups. ",
                    "Ellipses represent 95% confidence intervals."
                  )
                )
              )
            )
          )
        ),

        #------------------------------------------------
        # Tab 3: Model results
        #------------------------------------------------
        tabPanel(
          "Results",
          plotOutput("density_plot", height = "500px"),
          hr(),

          conditionalPanel(
            condition = "input.opt_b == true",
            plotOutput("beta_plot", height = "400px")
          ),

          hr(),
          tableOutput("summary_table"),
          hr(),

          fluidRow(
            column(
              6,
              downloadButton(
                "download_ind_res",
                "Download Individual Results (CSV)",
                class = "btn-info",
                style = "width:100%"
              )
            ),
            column(
              6,
              downloadButton(
                "download_chains",
                "Download Raw MCMC Chains (CSV)",
                class = "btn-warning",
                style = "width:100%"
              )
            )
          ),

          hr(),

          h4("Individual Sample Posterior Distributions"),

          fluidRow(
            column(4, uiOutput("ui_show_sample_id")),
            column(8, plotOutput("ind_density_plot", height = "350px"))
          ),

          br(),

          downloadButton(
            "download_all_plots",
            "Download All Density Plots (ZIP)",
            class = "btn-primary",
            style = "width:100%"
          )
        ),

        #------------------------------------------------
        # Tab 4: Reconstruction analysis
        #------------------------------------------------
        tabPanel(
          "Reconstruction Analysis",
          uiOutput("reconstruction_ui")
        ),

        #------------------------------------------------
        # Tab 5: NEW - Ghost Source Detector Diagnostics
        #------------------------------------------------
        tabPanel(
          "Ghost Source Detector",
          uiOutput("ghost_detector_ui")
        ),

        #------------------------------------------------
        # Tab 6: Validation
        #------------------------------------------------
        tabPanel(
          "Validation (virtual mixtures)",
          br(),

          helpText(
            "Validation generates virtual mixtures from the source data. ",
            "If 'Detect Unsampled Ghost Source' is active, you can either test the False-Alarm Rate (all sources known) ",
            "or hide a real source to stress-test whether BASIN recovers the missing source proportion!"
          ),

          fluidRow(
            column(
              3,
              radioButtons(
                "val_gen_mode",
                "VM Generation Logic:",
                choices = c(
                  "1. Full Stoch. (Dirichlet Props + MVN Tracers)" = "stoch_rand",
                  "2. Semi-Stoch. (Simplex Grid + MVN Tracers)" = "stoch_simp",
                  "3. Semi-Det. (Dirichlet Props + Mean Tracers)" = "det_rand",
                  "4. Full Det. (Simplex Grid + Mean Tracers)" = "det_simp"
                ),
                selected = "stoch_rand"
              )
            ),

            column(
              3,
              conditionalPanel(
                condition = "input.val_gen_mode == 'stoch_rand' || input.val_gen_mode == 'det_rand'",
                numericInput("val_n", "Number of VMs:", 20, min = 5, max = 500)
              ),

              conditionalPanel(
                condition = "input.val_gen_mode == 'det_simp' || input.val_gen_mode == 'stoch_simp'",
                numericInput("val_step", "Simplex Step (e.g. 0.05):", 0.05, min = 0.01, max = 0.50, step = 0.01)
              ),

              uiOutput("ui_val_hide_source")
            ),

            column(
              3,
              selectInput(
                "val_noise",
                "Data noise level:",
                choices = c(
                  "Not (0%)" = 0,
                  "Low (5%)" = 0.05,
                  "Medium (15%)" = 0.15,
                  "High (25%)" = 0.25
                ),
                selected = 0.05
              )
            ),

            column(
              3,
              actionButton(
                "run_val",
                "1. Generate and estimate",
                class = "btn-success",
                style = "width:100%"
              ),
              br(), br(),
              downloadButton(
                "download_val_data",
                "2. Download results (CSV)",
                style = "width:100%"
              )
            )
          ),

          hr(),

          h4("The model accuracy (True vs Predicted):"),
          plotOutput("val_plot", height = "600px"),

          h4("Errors & Performance Metrics:"),
          tableOutput("val_metrics")
        ),

        #------------------------------------------------
        # Tab 7: Model comparison
        #------------------------------------------------
        tabPanel(
          "Model Comparison",
          br(),

          h4("Bayesian Model Comparison & Diagnostics"),

          helpText(
            "Select multiple models to rank their predictive performance. ",
            "The model with elpd_diff = 0 is the best. ",
            "Compare only models fitted to the same data, same tracer set, and same transformation."
          ),

          hr(),

          fluidRow(
            column(
              4,
              uiOutput("ui_multi_compare_selector"),
              hr(),
              actionButton(
                "run_loo_compare",
                "Run LOO-CV",
                class = "btn-success",
                style = "width:100%"
              ),
              hr(),
              actionButton(
                "run_exact_loo",
                "Run Exact LOO-CV (SLOW!)",
                class = "btn-warning",
                style = "width:100%"
              )
            ),
            column(
              8,
              verbatimTextOutput("loo_results_text"),
              tags$style(
                type = "text/css",
                "#loo_results_text { font-weight: bold; background-color: #f8f9fa; }"
              )
            )
          ),

          hr(),

          h4("Pareto k Diagnostics (Influential Observations)"),

          helpText(
            "Points above the red dashed line (k > 0.7) strongly influence the model ",
            "and might be outliers or poorly fitted mixtures."
          ),

          fluidRow(
            column(4, uiOutput("ui_pareto_selector")),
            column(8, plotOutput("pareto_plot", height = "400px"))
          )
        ),

        #------------------------------------------------
        # Tab 8: Help
        #------------------------------------------------
        tabPanel(
          "User Guide & Math",
          if (file.exists("basin_help.html")) {
            includeHTML("basin_help.html")
          } else {
            tags$div(
              style = "padding:30px;",
              h4("Help file not found"),
              p("The file basin_help.html was not found in the application directory.")
            )
          }
        )
      )
    )
  ),

  hr(),

  tags$div(
    style = "text-align: center; color: #7f8c8d; padding-bottom: 20px; font-size: 0.9em;",
    tags$p(
      HTML("Developed by <b>Sergey Kharchenko</b>. Supported by <b>RSF</b>, grant No. 25-17-00143.")
    )
  )
)

#==============================================================================
# 4. SERVER LOGIC - WITH GHOST SOURCE DETECTOR (Part 2 of 2)
#==============================================================================

server <- function(input, output, session) {

  values <- reactiveValues(
    df = NULL,
    raw_df = NULL,
    groups = NULL,
    tracers = NULL,
    fit = NULL,
    stan_out = "",
    compiled_list = list(),
    level_names = NULL,
    saved_models = list(),
    loo_objects = list(),
    exact_loo_results = list(),
    suggested_tracers = NULL,
    filter_msg = "Click Tracers filtering before running the model",
    y_center = NULL,
    y_scale = NULL,
    ssa_grp_z = NULL,
    ghost_anchor = NULL,
    stan_data = NULL,
    compiled_model = NULL,
    last_model_name = NULL,
    current_run_folder = NULL,
    val_data = NULL,
    val_export = NULL,
    val_export_full = NULL
  )


  #==========================================================================
  # Helpers inside server
  #==========================================================================

  # Возвращает список имен источников с учетом латентного Ghost_Source (G+1)
  get_effective_groups <- function() {
    req(values$groups)
    if (isTRUE(input$opt_ghost)) {
      c(values$groups, "Ghost_Source")
    } else {
      values$groups
    }
  }

  reset_current_model <- function() {
    values$fit <- NULL
    values$stan_data <- NULL
    values$compiled_model <- NULL
    values$last_model_name <- NULL
    values$stan_out <- "Settings have been changed. Re-run the model."
  }

  reset_comparison <- function() {
    values$saved_models <- list()
    values$loo_objects <- list()
    values$exact_loo_results <- list()
    values$last_model_name <- NULL
  }

  reset_validation <- function() {
    values$val_data <- NULL
    values$val_export <- NULL
    values$val_export_full <- NULL
  }

  get_use_levels <- function(df_mix_rows = NULL) {
    if (is.null(input$mix_cov) || input$mix_cov == "None") return(FALSE)
    if (is.null(input$pool_mode) || input$pool_mode == "no_pool") return(FALSE)

    if (!is.null(df_mix_rows)) {
      return(input$mix_cov %in% colnames(df_mix_rows))
    }

    if (!is.null(values$df)) {
      return(input$mix_cov %in% colnames(values$df))
    }

    FALSE
  }

  sanitize_filename <- function(x) {
    x <- as.character(x)
    x <- gsub("[^A-Za-z0-9_-]+", "_", x)
    x <- gsub("_+", "_", x)
    ifelse(x == "", "sample", x)
  }


  #==========================================================================
  # Reset logic
  #==========================================================================

  observeEvent(input$file, {
    req(input$file)

    values$raw_df <- data.table::fread(
      input$file$datapath,
      data.table = FALSE,
      check.names = FALSE
    )
    values$df <- values$raw_df

    values$groups <- NULL
    values$suggested_tracers <- NULL
    values$filter_msg <- "Click Tracers filtering before running the model"

    reset_current_model()
    reset_comparison()
    reset_validation()
  })

  observeEvent(list(input$group_col, input$mix_label), {
    req(values$raw_df, input$group_col, input$mix_label)

    df <- values$raw_df

    if (input$group_col %in% colnames(df)) {
      if (input$mix_label %in% as.character(df[[input$group_col]])) {
        df_src <- df[as.character(df[[input$group_col]]) != input$mix_label, ]
        df_mix <- df[as.character(df[[input$group_col]]) == input$mix_label, ]

        df_src <- df_src[order(as.character(df_src[[input$group_col]])), ]
        values$df <- rbind(df_src, df_mix)
      } else {
        values$df <- df
      }
    }

    reset_current_model()
    reset_comparison()
    reset_validation()
  }, ignoreInit = TRUE)

  observeEvent(
    list(
      input$pool_mode,
      input$opt_b,
      input$opt_b_src,
      input$cov_mode,
      input$opt_coda,
      input$opt_ghost,
      input$ghost_alpha,
      input$hotspot_mode,
      input$k_hotspot,
      input$sigma_hotspot,
      input$calc_method,
      input$opt_buddy,
      input$opt_hull,
      input$opt_bias,
      input$size_col,
      input$mix_cov
    ),
    {
      reset_current_model()
    },
    ignoreInit = TRUE
  )

  observeEvent(
    list(input$final_selected_tracers, input$opt_coda),
    {
      reset_current_model()
      reset_comparison()
      reset_validation()
    },
    ignoreInit = TRUE
  )

  observeEvent(input$pool_mode, {
    if (!is.null(input$pool_mode) && input$pool_mode == "no_pool") {
      updateSelectInput(session, "mix_cov", selected = "None")
    }
  }, ignoreInit = TRUE)

  # Синхронизация разделенных осей (3.2a Error-Scale и 3.3 Source Representation) с внутренним hotspot_mode
  observeEvent(list(input$src_rep_mode, input$err_scale_mode), {
    req(input$src_rep_mode)
    effective_mode <- if (input$src_rep_mode == "centroid") {
      if (!is.null(input$err_scale_mode)) input$err_scale_mode else "mean"
    } else {
      input$src_rep_mode
    }
    updateRadioButtons(session, "hotspot_mode", selected = effective_mode)
  }, ignoreInit = TRUE)


  #==========================================================================
  # Dynamic selectors
  #==========================================================================

  output$ui_group_col <- renderUI({
    req(values$raw_df)

    selectInput(
      "group_col",
      "Sources column:",
      choices = colnames(values$raw_df),
      selected = colnames(values$raw_df)[min(2, ncol(values$raw_df))]
    )
  })

  output$ui_mix_label <- renderUI({
    req(values$raw_df, input$group_col)

    choices <- unique(as.character(values$raw_df[[input$group_col]]))

    selectInput(
      "mix_label",
      "Mixture label:",
      choices = choices,
      selected = tail(choices, 1)
    )
  })

  output$ui_mix_cov <- renderUI({
    req(values$df)

    cols <- colnames(values$df)

    selectInput(
      "mix_cov",
      "Mixture grouping covariate (e.g. Season/Trap):",
      choices = c("None", cols),
      selected = "None"
    )
  })

  output$ui_size_col <- renderUI({
    req(values$df)

    num_cols <- colnames(values$df)[sapply(values$df, is.numeric)]

    selectInput(
      "size_col",
      "Sorting / Enrichment Proxy Column (e.g., SSA, d50, TOC; optional):",
      choices = c("None", num_cols),
      selected = "None"
    )
  })

  output$ui_tracer_select <- renderUI({
    req(values$df)

    all_nums <- colnames(values$df)[sapply(values$df, is.numeric)]
    exclude <- c("ID", "id", "Id")

    if (!is.null(input$size_col) && input$size_col != "None") {
      exclude <- c(exclude, input$size_col)
    }

    if (!is.null(input$mix_cov) && input$mix_cov != "None") {
      exclude <- c(exclude, input$mix_cov)
    }

    cand <- setdiff(all_nums, exclude)

    selectizeInput(
      "tracers_to_use",
      "Choose tracers:",
      choices = cand,
      selected = cand,
      multiple = TRUE
    )
  })

  output$ui_show_sample_id <- renderUI({
    req(values$fit, values$df, input$mix_label, input$group_col)

    df_mix <- values$df[as.character(values$df[[input$group_col]]) == input$mix_label, ]

    mix_ids <- if ("ID" %in% colnames(df_mix)) {
      as.character(df_mix$ID)
    } else {
      paste0("Mix_", seq_len(nrow(df_mix)))
    }

    selectInput(
      "show_sample_id",
      "Select Sample:",
      choices = mix_ids
    )
  })

  output$ui_val_hide_source <- renderUI({
    req(values$groups)
    if (!isTRUE(input$opt_ghost)) return(NULL)

    choices_list <- c(
      "1. Hide nothing (Include Ghost in VMs: G+1)" = "include_ghost",
      "2. Hide Ghost Source (False-Alarm Test: Ghost = 0)" = "hide_ghost"
    )

    if (length(values$groups) >= 3) {
      real_src_choices <- setNames(
        values$groups,
        paste0(seq_along(values$groups) + 2, ". Hide ", values$groups, " (Blind Test)")
      )
      choices_list <- c(choices_list, real_src_choices)
    }

    selectInput(
      "val_hide_source",
      "Ghost Validation Mode:",
      choices = choices_list,
      selected = "include_ghost"
    )
  })

  output$raw_table <- renderDT({
    req(values$df)
    datatable(values$df, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$count_status <- renderText({
    paste("Selected for analysis:", length(input$final_selected_tracers), "tracers")
  })

  output$conn_inputs <- renderUI({
    req(values$df, input$group_col, input$mix_label)

    grs <- unique(as.character(values$df[[input$group_col]]))
    grs <- grs[grs != input$mix_label]
    grs <- sort(grs)

    values$groups <- grs

    tagList(
      h4("Geographical weights:"),
      p(tags$small("e.g. Area * Erosion * IC")),
      lapply(
        grs,
        function(g) {
          numericInput(
            paste0("w_", g),
            paste("Source:", g),
            value = 1.0,
            step = 0.1
          )
        }
      )
    )
  })


  #==========================================================================
  # Filtering (Range + Robust GAPT + True Convex Hull Penalty + KW + Corr)
  #==========================================================================

  observeEvent(input$run_filter, {
    req(values$df, input$tracers_to_use)

    withProgress(message = "Evaluating tracer conservatism and source-space geometry...", {

      df <- values$df
      gr_col <- input$group_col
      m_lab <- input$mix_label
      ssa_col <- input$size_col
      candidates <- input$tracers_to_use

      src <- df[as.character(df[[gr_col]]) != m_lab, ]
      mix <- df[as.character(df[[gr_col]]) == m_lab, ]

      if (nrow(src) == 0 || nrow(mix) == 0) {
        showNotification("Error: Source or Mixture group is empty.", type = "error")
        return()
      }

      ssa_factor <- 1.0

      if (!is.null(ssa_col) && ssa_col != "None" && ssa_col %in% colnames(df)) {
        mean_ssa_mix <- mean(mix[[ssa_col]], na.rm = TRUE)
        mean_ssa_src <- mean(src[[ssa_col]], na.rm = TRUE)

        if (!is.na(mean_ssa_mix) && !is.na(mean_ssa_src) && mean_ssa_src > 0 && mean_ssa_mix > 0) {
          ssa_factor <- mean_ssa_mix / mean_ssa_src
        }
      }

      log_text <- c(paste("--- PROXY-AWARE RANGE FILTERING (Mode:", input$filter_mode, ") ---"))
      log_text <- c(log_text, paste("[INFO] Sorting proxy enrichment factor:", round(ssa_factor, 2)))

      pass_range <- c()

      for (t in candidates) {
        if (!is.numeric(df[[t]])) next

        mv_raw <- mean(mix[[t]], na.rm = TRUE)
        mv_norm <- mv_raw / ssa_factor

        s_min <- min(src[[t]], na.rm = TRUE)
        s_max <- max(src[[t]], na.rm = TRUE)
        s_span <- s_max - s_min
        low_b  <- s_min - 0.05 * max(abs(s_min), s_span)
        high_b <- s_max + 0.05 * max(abs(s_max), s_span)

        check_raw <- (mv_raw >= low_b && mv_raw <= high_b)
        check_norm <- (mv_norm >= low_b && mv_norm <= high_b)

        res_check <- switch(
          input$filter_mode,
          "raw" = check_raw,
          "corr" = check_norm,
          "both" = (check_raw && check_norm),
          "any" = (check_raw || check_norm)
        )

        if (is.na(res_check)) res_check <- FALSE

        if (res_check) {
          pass_range <- c(pass_range, t)
          if (!check_raw && check_norm) {
            log_text <- c(log_text, paste0("[RANGE] ", t, ": SALVAGED by SSA correction"))
          }
        } else {
          log_text <- c(log_text, paste0("[RANGE] ", t, ": REMOVED (Outside source polygon)"))
        }
      }

      pass_buddy <- pass_range

      if (input$opt_buddy && length(pass_range) > 2) {
        log_text <- c(log_text, "--- INTER-TRACER ASSOCIATION PRESERVATION TEST (APT) ---")

        src_cor <- suppressWarnings(
          cor(src[, pass_range, drop = FALSE], use = "pairwise.complete.obs", method = "spearman")
        )
        src_cor[is.na(src_cor)] <- 0
        diag(src_cor) <- 0

        bad_buddies <- c()

        for (t in pass_range) {
          buddies <- pass_range[abs(src_cor[t, pass_range]) > 0.85]
          if (length(buddies) < 2) next

          broken_links <- 0

          for (b in buddies) {
            coeffs <- tryCatch(
              suppressWarnings(coef(MASS::rlm(src[[t]] ~ src[[b]], maxit = 50))),
              error = function(e) {
                tryCatch(coef(lm(src[[t]] ~ src[[b]])), error = function(e2) NULL)
              }
            )

            if (is.null(coeffs) || any(is.na(coeffs))) next

            pred_raw <- coeffs[1] + coeffs[2] * mean(mix[[b]], na.rm = TRUE)
            pred_norm <- coeffs[1] + coeffs[2] * (mean(mix[[b]], na.rm = TRUE) / ssa_factor)

            err_raw <- abs(mean(mix[[t]], na.rm = TRUE) - pred_raw) / (abs(pred_raw) + 1e-9)
            err_norm <- abs((mean(mix[[t]], na.rm = TRUE) / ssa_factor) - pred_norm) / (abs(pred_norm) + 1e-9)

            active_error <- switch(
              input$filter_mode,
              "raw" = err_raw,
              "corr" = err_norm,
              "both" = max(err_raw, err_norm),
              "any" = min(err_raw, err_norm)
            )

            if (is.na(active_error) || active_error > 0.10) {
              broken_links <- broken_links + 1
            }
          }

          if (broken_links / length(buddies) >= 0.5) {
            bad_buddies <- c(bad_buddies, t)
            log_text <- c(
              log_text,
              paste0("[APT] ", t, ": REMOVED (Broken ", broken_links, "/", length(buddies), " links)")
            )
          }
        }

        pass_buddy <- setdiff(pass_range, bad_buddies)
      }

      pass_hull <- pass_buddy

      if (input$opt_hull && length(pass_buddy) > 4) {
        log_text <- c(log_text, "--- PCA CONVEX HULL PENALTY ---")

        mix_norm <- mix
        for (tr_name in pass_buddy) {
          mix_norm[[tr_name]] <- mix[[tr_name]] / ssa_factor
        }

        calc_pca_hull_leverage <- function(tr_list, target_mix_df) {
          pca <- tryCatch(
            prcomp(src[, tr_list, drop = FALSE], scale. = TRUE),
            error = function(e) NULL
          )

          if (is.null(pca) || ncol(pca$x) < 2) return(NA_real_)

          s_pcs <- predict(pca, src[, tr_list, drop = FALSE])[, 1:2, drop = FALSE]
          m_pcs <- predict(pca, target_mix_df[, tr_list, drop = FALSE])[, 1:2, drop = FALSE]

          r_cloud <- mean(sqrt(rowSums(sweep(s_pcs, 2, colMeans(s_pcs), "-")^2)))
          if (is.na(r_cloud) || r_cloud < 1e-9) return(NA_real_)

          d_hull_vec <- apply(m_pcs, 1, function(pt) dist_to_convex_hull_2d(s_pcs, pt))
          d_hull <- mean(d_hull_vec, na.rm = TRUE)

          d_cent <- sqrt(sum((colMeans(m_pcs) - colMeans(s_pcs))^2))

          if (d_hull > 0) {
            return((d_cent + d_hull) / (r_cloud + 1e-9))
          } else {
            return(0.0)
          }
        }

        get_active_dist <- function(tr_list) {
          d_raw <- calc_pca_hull_leverage(tr_list, mix)
          d_norm <- calc_pca_hull_leverage(tr_list, mix_norm)
          vals <- c(d_raw, d_norm)
          if (all(is.na(vals))) return(NA_real_)
        
          switch(
            input$filter_mode,
            "raw" = d_raw,
            "corr" = d_norm,
            "both" = max(vals, na.rm = TRUE),
            "any" = min(vals, na.rm = TRUE)
          )
        }

        base_dist <- get_active_dist(pass_buddy)
        bad_hull <- c()

        if (!is.na(base_dist) && base_dist > 1.2) {
          for (t in pass_buddy) {
            new_dist <- get_active_dist(setdiff(pass_buddy, t))

            if (!is.na(new_dist) && new_dist < (base_dist * 0.95)) {
              bad_hull <- c(bad_hull, t)
              log_text <- c(log_text, paste0("[HULL] ", t, ": REMOVED (High PCA hull leverage)"))
            }
          }
        }

        pass_hull <- setdiff(pass_buddy, bad_hull)
      }

      pass_kw <- c()

      for (t in pass_hull) {
        pval <- tryCatch({
          kruskal.test(
            as.formula(paste0("`", t, "` ~ `", gr_col, "`")),
            data = src
          )$p.value
        }, error = function(e) NA)

        if (!is.na(pval) && pval < 0.05) {
          pass_kw <- c(pass_kw, t)
        } else {
          log_text <- c(log_text, paste0("[KW-TEST] ", t, ": REMOVED (No discrimination)"))
        }
      }

      final_pass <- pass_kw

      if (length(final_pass) > 1) {
        cor_mat <- cor(src[, final_pass, drop = FALSE], use = "pairwise.complete.obs")
        cor_mat[is.na(cor_mat)] <- 0

        high_c_idx <- caret::findCorrelation(cor_mat, cutoff = 0.95)

        if (length(high_c_idx) > 0) {
          for (rt in final_pass[high_c_idx]) {
            log_text <- c(log_text, paste0("[CORR] ", rt, ": REMOVED (Redundant)"))
          }
          final_pass <- final_pass[-high_c_idx]
        }
      }

      log_text <- c(
        log_text,
        "------------------------",
        paste("PASSED:", length(final_pass), "tracers.")
      )

      values$suggested_tracers <- final_pass
      values$filter_msg <- paste(log_text, collapse = "\n")

      showNotification(
        paste("Complete. Suggested:", length(final_pass), "tracers."),
        type = "message"
      )
    })
  })

  output$filter_log <- renderText({
    values$filter_msg
  })

  output$ui_final_tracers_chk <- renderUI({
    req(input$tracers_to_use)

    choices <- input$tracers_to_use
    selected <- if (!is.null(values$suggested_tracers)) values$suggested_tracers else choices

    checkboxGroupInput(
      "final_selected_tracers",
      label = NULL,
      choices = choices,
      selected = selected,
      inline = TRUE
    )
  })

  outputOptions(output, "ui_final_tracers_chk", suspendWhenHidden = FALSE)


  #==========================================================================
  # PCA / LDA plots
  #==========================================================================

  output$pca_plot <- renderPlot({
    req(input$final_selected_tracers, input$group_col, input$mix_label, values$df)

    tr <- input$final_selected_tracers
    df <- values$df
    m_lab <- input$mix_label
    gr_col <- input$group_col

    src <- df[as.character(df[[gr_col]]) != m_lab, ]
    mix <- df[as.character(df[[gr_col]]) == m_lab, ]

    pca <- prcomp(src[, tr, drop = FALSE], scale. = TRUE)
    var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

    plot_df <- rbind(
      data.frame(predict(pca, src[, tr, drop = FALSE]), Group = as.character(src[[gr_col]])),
      data.frame(predict(pca, mix[, tr, drop = FALSE]), Group = m_lab)
    )

    hull <- plot_df %>%
      filter(Group != m_lab) %>%
      group_by(Group) %>%
      slice(chull(PC1, PC2))

    loadings <- as.data.frame(pca$rotation[, 1:2])
    loadings$Tracer <- rownames(loadings)

    scale_factor <- min(max(abs(plot_df$PC1)), max(abs(plot_df$PC2))) /
      (max(abs(pca$rotation[, 1:2])) + 1e-9) * 0.8

    loadings$PC1 <- loadings$PC1 * scale_factor
    loadings$PC2 <- loadings$PC2 * scale_factor

    ggplot() +
      geom_polygon(data = hull, aes(x = PC1, y = PC2, fill = Group), alpha = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
      geom_point(
        data = filter(plot_df, Group != m_lab),
        aes(x = PC1, y = PC2, color = Group),
        size = 5,
        alpha = 0.8
      ) +
      geom_segment(
        data = loadings,
        aes(x = 0, y = 0, xend = PC1, yend = PC2),
        arrow = arrow(length = unit(0.2, "cm")),
        color = "black",
        linewidth = 0.7,
        alpha = 0.6
      ) +
      geom_text(
        data = loadings,
        aes(x = PC1 * 1.15, y = PC2 * 1.15, label = Tracer),
        color = "black",
        fontface = "bold",
        size = 5
      ) +
      geom_point(
        data = filter(plot_df, Group == m_lab),
        aes(x = PC1, y = PC2),
        shape = 23,
        fill = "red",
        size = 8,
        stroke = 2
      ) +
      scale_color_viridis_d(option = "turbo") +
      scale_fill_viridis_d(option = "turbo") +
      theme_bw(base_size = 18) +
      theme(
        axis.line = element_line(linewidth = 1.2),
        axis.text = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      ) +
      labs(
        title = paste("PCA Biplot (Mixing Polygon):", m_lab),
        x = paste0("PC1 (", var_exp[1], "%)"),
        y = paste0("PC2 (", var_exp[2], "%)")
      )
  })

  output$lda_plot <- renderPlot({
    req(input$final_selected_tracers, input$group_col, input$mix_label, values$df)

    tr <- input$final_selected_tracers
    df <- values$df
    m_lab <- input$mix_label
    gr_col <- input$group_col

    src_clean <- df[as.character(df[[gr_col]]) != m_lab, ] %>%
      drop_na(all_of(tr))

    if (length(unique(src_clean[[gr_col]])) < 2 || nrow(src_clean) <= length(tr)) {
      return(NULL)
    }

    lda_mod <- tryCatch(
      MASS::lda(src_clean[, tr, drop = FALSE], grouping = src_clean[[gr_col]]),
      error = function(e) NULL
    )

    if (is.null(lda_mod)) return(NULL)

    lda_values <- predict(lda_mod)
    lda_plot_df <- data.frame(lda_values$x, Group = as.character(src_clean[[gr_col]]))
    prop_lda <- lda_mod$svd^2 / sum(lda_mod$svd^2)

    is_1d <- ncol(lda_values$x) == 1

    if (is_1d) {
      lda_plot_df$LD2 <- rnorm(nrow(lda_plot_df), 0, 0.05)
      colnames(lda_plot_df)[1] <- "LD1"
      x_label <- "LD1 (100%)"
      y_label <- "Jitter"

      ld_loadings <- data.frame(
        LD1 = lda_mod$scaling[, 1],
        LD2 = 0,
        Tracer = rownames(lda_mod$scaling)
      )
    } else {
      x_label <- paste0("LD1 (", round(prop_lda[1] * 100, 1), "%)")
      y_label <- paste0("LD2 (", round(prop_lda[2] * 100, 1), "%)")

      ld_loadings <- data.frame(
        LD1 = lda_mod$scaling[, 1],
        LD2 = lda_mod$scaling[, 2],
        Tracer = rownames(lda_mod$scaling)
      )
    }

    lda_plot_df$SampleID <- if ("ID" %in% colnames(src_clean)) {
      src_clean$ID
    } else {
      seq_len(nrow(lda_plot_df))
    }

    max_val_points <- max(abs(lda_plot_df$LD1), if (!is_1d) abs(lda_plot_df$LD2) else 0)
    max_val_loadings <- max(abs(ld_loadings$LD1), abs(ld_loadings$LD2))

    scale_factor <- (max_val_points / (max_val_loadings + 1e-9)) * 0.8

    ld_loadings$LD1 <- ld_loadings$LD1 * scale_factor
    ld_loadings$LD2 <- ld_loadings$LD2 * scale_factor

    ggplot() +
      stat_ellipse(
        data = lda_plot_df,
        aes(x = LD1, y = LD2, color = Group, fill = Group),
        geom = "polygon",
        alpha = 0.2,
        show.legend = FALSE
      ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
      geom_point(
        data = lda_plot_df,
        aes(x = LD1, y = LD2, color = Group, fill = Group),
        size = 4,
        alpha = 0.8
      ) +
      geom_text(
        data = lda_plot_df,
        aes(x = LD1, y = LD2, label = SampleID, color = Group),
        vjust = -1.2,
        size = 4,
        show.legend = FALSE
      ) +
      geom_segment(
        data = ld_loadings,
        aes(x = 0, y = 0, xend = LD1, yend = LD2),
        arrow = arrow(length = unit(0.2, "cm")),
        color = "black",
        linewidth = 0.7,
        alpha = 0.6
      ) +
      geom_text(
        data = ld_loadings,
        aes(x = LD1 * 1.15, y = LD2 * 1.15, label = Tracer),
        color = "black",
        fontface = "bold",
        size = 5
      ) +
      theme_bw(base_size = 18) +
      scale_color_viridis_d(option = "turbo") +
      scale_fill_viridis_d(option = "turbo") +
      labs(title = "Linear Discriminant Analysis", x = x_label, y = y_label) +
      theme(
        legend.position = "bottom",
        axis.line = element_line(linewidth = 1.2),
        panel.grid.minor = element_blank()
      )
  })


  #==========================================================================
  # MAIN MODEL EXECUTION (CONTINUATION OF PART 2)
  #==========================================================================

  observeEvent(input$run_model, {
    req(input$final_selected_tracers, values$df, input$group_col, input$mix_label)

    if (length(input$final_selected_tracers) < 2) {
      showNotification("Error: Select at least 2 tracers!", type = "error")
      return()
    }

    size_active <- !is.null(input$size_col) &&
      input$size_col != "None" &&
      input$size_col %in% colnames(values$df)

    if (input$opt_b && !size_active) {
      showNotification(
        "Beta correction requires an SSA / SSA-proxy column. Select a column or disable Beta correction.",
        type = "error"
      )
      return()
    }

    withProgress(message = "Preparing and calculating...", value = 0, {

      tr <- input$final_selected_tracers
      m_lab <- input$mix_label
      gr_col <- input$group_col

      df_mix_rows <- values$df[as.character(values$df[[gr_col]]) == m_lab, ]
      df_src_rows <- values$df[as.character(values$df[[gr_col]]) != m_lab, ]

      if (nrow(df_mix_rows) == 0 || nrow(df_src_rows) == 0) {
        showNotification("Error: Source or Mixture group is empty.", type = "error")
        return()
      }

      if (is.null(values$groups) || length(values$groups) == 0) {
        showNotification("Error: Source groups are not initialized yet.", type = "error")
        return()
      }

      y_src_raw <- as.matrix(df_src_rows[, tr, drop = FALSE])
      y_mix_raw <- as.matrix(df_mix_rows[, tr, drop = FALSE])

      if (any(is.na(y_src_raw)) || any(is.na(y_mix_raw))) {
        showNotification("Error: Selected tracers contain NA values.", type = "error")
        return()
      }

      has_nonpos <- any(y_src_raw <= 0) || any(y_mix_raw <= 0)
      if (has_nonpos && input$opt_coda) {
        showNotification(
          "Warning: Non-positive tracer values detected with CoDA (ILR) active! Clamping to 1e-6.",
          type = "warning",
          duration = 8
        )
        y_src_raw[y_src_raw <= 0] <- 1e-6
        y_mix_raw[y_mix_raw <= 0] <- 1e-6
      } else if (has_nonpos && input$opt_b) {
        for (col_i in seq_len(ncol(y_src_raw))) {
          min_v <- min(c(y_src_raw[, col_i], y_mix_raw[, col_i]), na.rm = TRUE)
          if (min_v <= 0) {
            shift_v <- abs(min_v) + 1.0
            y_src_raw[, col_i] <- y_src_raw[, col_i] + shift_v
            y_mix_raw[, col_i] <- y_mix_raw[, col_i] + shift_v
          }
        }
      }

      if (!input$opt_coda) {
        sc_src <- safe_scale_matrix(y_src_raw)
        y_src <- sc_src$y
        y_mix <- scale(y_mix_raw, center = sc_src$center, scale = sc_src$scale)

        values$y_center <- sc_src$center
        values$y_scale <- sc_src$scale

        V_matrix <- matrix(0, nrow = 0, ncol = 0)
        N_tr_raw <- ncol(y_mix)
      } else {
        ys_prop <- y_src_raw / rowSums(y_src_raw)
        ym_prop <- y_mix_raw / rowSums(y_mix_raw)

        N_tr_raw <- ncol(ys_prop)
        V_matrix <- compositions::ilrBase(D = N_tr_raw)

        y_src <- as.matrix(log(ys_prop) %*% V_matrix)
        y_mix <- as.matrix(log(ym_prop) %*% V_matrix)

        values$y_center <- NULL
        values$y_scale <- NULL
      }

      y_src <- as.matrix(y_src)
      y_mix <- as.matrix(y_mix)

      if (any(!is.finite(y_src)) || any(!is.finite(y_mix))) {
        showNotification("Error: Non-finite values after data transformation.", type = "error")
        return()
      }

      #--------------------------------------------------------------
      # SSA preparation, standardized
      #--------------------------------------------------------------
      ssa_m <- rep(0, nrow(y_mix))
      ssa_s <- rep(0, nrow(y_src))
      ssa_grp_vec <- rep(0, length(values$groups))

      if (size_active) {
        ssa_m_raw <- df_mix_rows[[input$size_col]]
        ssa_s_raw <- df_src_rows[[input$size_col]]

        if (all(is.na(ssa_s_raw))) {
          ssa_s_raw <- rep(0, length(ssa_s_raw))
        } else {
          ssa_s_raw[is.na(ssa_s_raw)] <- mean(ssa_s_raw, na.rm = TRUE)
        }

        if (all(is.na(ssa_m_raw))) {
          ssa_m_raw <- rep(mean(ssa_s_raw, na.rm = TRUE), length(ssa_m_raw))
        } else {
          ssa_m_raw[is.na(ssa_m_raw)] <- mean(ssa_m_raw, na.rm = TRUE)
        }

        mean_s <- mean(ssa_s_raw, na.rm = TRUE)
        sd_s <- sd(ssa_s_raw, na.rm = TRUE)

        if (!is.na(sd_s) && sd_s >= 1e-8) {
          ssa_m <- (ssa_m_raw - mean_s) / sd_s
          ssa_s <- (ssa_s_raw - mean_s) / sd_s
          src_labels <- as.character(df_src_rows[[gr_col]])
          ssa_grp_vec <- sapply(values$groups, function(g) {
            x <- ssa_s[src_labels == g]
            if (length(x) == 0 || all(is.na(x))) 0 else mean(x, na.rm = TRUE)
          })
        }
      }

      values$ssa_grp_z <- ssa_grp_vec

      #--------------------------------------------------------------
      # Levels & group_id
      #--------------------------------------------------------------
      use_levels <- get_use_levels(df_mix_rows)
      N_mix_levels <- NULL
      mix_level_id <- NULL

      if (use_levels) {
        cov_vec <- as.character(df_mix_rows[[input$mix_cov]])
        cov_vec[is.na(cov_vec)] <- "NA"
        cov_factor <- factor(cov_vec)
        mix_level_id <- as.integer(cov_factor)
        N_mix_levels <- nlevels(cov_factor)
        values$level_names <- levels(cov_factor)
      } else {
        values$level_names <- "Overall"
      }

      src_labels <- as.character(df_src_rows[[gr_col]])
      group_id <- make_group_id(src_labels, values$groups)

      if (any(is.na(group_id))) {
        showNotification("Error: Source labels do not match source groups.", type = "error")
        return()
      }

      conn_weights <- vapply(
        values$groups,
        function(g) {
          x <- input[[paste0("w_", g)]]
          if (is.null(x) || is.na(x)) 1 else as.numeric(x)
        },
        numeric(1)
      )

      if (any(is.na(conn_weights)) || any(conn_weights <= 0)) {
        showNotification("Error: Geographical weights must be positive numbers.", type = "error")
        return()
      }

      #--------------------------------------------------------------
      # Compile model (including Ghost Detector option)
      #--------------------------------------------------------------
      use_b_src_flag <- isTRUE(input$opt_b && input$opt_b_src && !input$opt_coda)
      use_ghost_flag <- isTRUE(input$opt_ghost)
      h_mode <- if (!is.null(input$hotspot_mode)) input$hotspot_mode else "mean"

      model_id <- paste0(
        input$pool_mode,
        as.integer(input$opt_b),
        as.integer(use_b_src_flag),
        input$cov_mode,
        as.integer(input$opt_coda),
        as.integer(use_levels),
        as.integer(input$opt_bias),
        as.integer(use_ghost_flag),
        "_", h_mode
      )

      if (is.null(values$compiled_list[[model_id]])) {
        showNotification("Compiling Stan model (this may take time)...", type = "message")

        code <- generate_stan_code(
          pool_mode = input$pool_mode,
          b = input$opt_b,
          b_src = use_b_src_flag,
          cov_mode = input$cov_mode,
          coda = input$opt_coda,
          has_groups = use_levels,
          bias = input$opt_bias,
          ghost = use_ghost_flag,
          hotspot_mode = h_mode
        )

        values$compiled_list[[model_id]] <- stan_model(model_code = code)
      }

      compiled_mod <- values$compiled_list[[model_id]]
      values$compiled_model <- compiled_mod

      #--------------------------------------------------------------
      # Prepare Stan data
      #--------------------------------------------------------------
      y_ctr_stan <- if (!is.null(values$y_center)) values$y_center else rep(0, ncol(y_mix))
      y_scl_stan <- if (!is.null(values$y_scale))  values$y_scale  else rep(1, ncol(y_mix))

      s_data <- list(
        N_mix = nrow(y_mix),
        N_sources = nrow(y_src),
        N_tracers = ncol(y_mix),
        N_groups = length(values$groups),
        y_mix = y_mix,
        y_source = y_src,
        group_id = group_id,
        conn_weights = conn_weights,
        ssa_mix = as.array(ssa_m),
        ssa_source = as.array(ssa_s),
        ssa_group = as.array(ssa_grp_vec),
        y_center = as.array(as.numeric(y_ctr_stan)),
        y_scale = as.array(as.numeric(y_scl_stan)),
        holdout_id = 0L
      )

      if (use_ghost_flag) {
        gh_anchor <- calc_ghost_anchor(y_src, group_id, length(values$groups), y_mix)
        values$ghost_anchor <- gh_anchor
        s_data$ghost_anchor <- as.array(as.numeric(gh_anchor))
        s_data$alpha_ghost  <- if (!is.null(input$ghost_alpha)) as.numeric(input$ghost_alpha) else 0.12
      } else {
        values$ghost_anchor <- NULL
      }

      if (h_mode == "k_scale") {
        s_data$k_hotspot <- if (!is.null(input$k_hotspot) && !is.na(input$k_hotspot)) max(0.1, as.numeric(input$k_hotspot)) else 3.0
      }
      if (h_mode == "nested") {
        s_data$sigma_hotspot <- if (!is.null(input$sigma_hotspot) && !is.na(input$sigma_hotspot)) max(0.1, as.numeric(input$sigma_hotspot)) else 1.5
      }

      if (use_levels) {
        s_data$N_mix_levels <- N_mix_levels
        s_data$mix_level_id <- as.array(mix_level_id)
      }

      if (input$opt_coda) {
        s_data$N_tracers_raw <- N_tr_raw
        s_data$V <- V_matrix
        s_data$src_totals <- as.array(
          sapply(values$groups, function(g) {
            mean(rowSums(df_src_rows[df_src_rows[[gr_col]] == g, tr, drop = FALSE], na.rm = TRUE))
          })
        )
      }

      values$stan_data <- s_data

      #--------------------------------------------------------------
      # Fit model
      #--------------------------------------------------------------
      t_start <- Sys.time()
      logs <- capture.output({
        tryCatch({
          if (input$calc_method == "vb") {
            values$fit <- vb(
              compiled_mod,
              data = s_data,
              init = "random",
              output_samples = 2000,
              seed = 123,
              tol_rel_obj = 0.005
            )
          } else {
            values$fit <- sampling(
              compiled_mod,
              data = s_data,
              iter = input$iter,
              chains = input$chains,
              init = 0.5,
              control = list(adapt_delta = 0.98, max_treedepth = 14)
            )
          }
          suppressWarnings(gc())
        }, error = function(e) {
          showNotification(paste("Stan Error:", e$message), type = "error")
          values$stan_out <- paste("ERROR:", e$message)
          values$fit <- NULL
        })
      }, type = "output")
      t_end <- Sys.time()

      values$stan_out <- paste(logs, collapse = "\n")

      #--------------------------------------------------------------
      # Auto-save model & workspace
      #--------------------------------------------------------------
      if (!is.null(values$fit) && isS4(values$fit)) {
        eff_groups <- get_effective_groups()

        pool_label <- switch(
          input$pool_mode,
          "no_pool" = "NoPool",
          "hier" = "Hier",
          "complete" = "Complete"
        )
        cov_label <- switch(
          input$cov_mode,
          "none" = "Cov:No",
          "mix" = "Cov:Mix",
          "full" = "Cov:Full"
        )

        opts <- c(pool_label, cov_label)
        if (h_mode == "proc_only") opts <- c(opts, "ProcErr")
        if (h_mode == "k_scale") opts <- c(opts, paste0("HotK(", s_data$k_hotspot, ")"))
        if (h_mode == "nested")  opts <- c(opts, paste0("Nested(", s_data$sigma_hotspot, ")"))
        if (input$opt_coda) opts <- c(opts, "ILR") else opts <- c(opts, "Scale")
        if (input$opt_b) opts <- c(opts, "+Beta")
        if (use_b_src_flag) opts <- c(opts, "+BetaSrc")
        if (input$opt_bias) opts <- c(opts, "+Bias")
        if (use_ghost_flag) opts <- c(opts, paste0("+Ghost(", s_data$alpha_ghost, ")"))

        config_str <- paste(opts, collapse = "-")
        tracer_tag <- paste(substr(sort(tr), 1, 4), collapse = "+")
        base_name <- paste0("[", toupper(input$calc_method), "] ", length(tr), "Tr (", tracer_tag, ") | ", config_str)

        model_name <- base_name
        counter <- 1
        while (model_name %in% names(values$saved_models)) {
          model_name <- paste0(base_name, " (", counter, ")")
          counter <- counter + 1
        }

        values$saved_models[[model_name]] <- values$fit
        values$last_model_name <- model_name

        safe_folder_name <- gsub("[^A-Za-z0-9_-]", "_", model_name)
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        run_folder <- file.path("BASIN_Workspaces", paste0(timestamp, "_", safe_folder_name))
        dir.create(run_folder, recursive = TRUE, showWarnings = FALSE)
        values$current_run_folder <- run_folder

        tryCatch({
          write.csv(as.data.frame(values$fit), file.path(run_folder, "01_MCMC_Raw_Chains.csv"), row.names = FALSE)
          saveRDS(values$fit, file.path(run_folder, "02_Stanfit_Object.rds"))

          p_name <- if (input$pool_mode == "hier") {
            "P_global"
          } else if (input$pool_mode == "no_pool") {
            "P_ind"
          } else if (use_levels) {
            "P_level"
          } else {
            "P"
          }
          samples <- rstan::extract(values$fit, pars = p_name)[[1]]
          if (p_name %in% c("P_ind", "P_level")) {
            samples <- ensure_3d_samples(samples, length(eff_groups))
          }
          if (length(dim(samples)) == 3) samples <- apply(samples, c(1, 3), mean)
          if (is.null(dim(samples)) || length(dim(samples)) < 2) {
            samples <- matrix(samples, ncol = length(eff_groups))
          }
          colnames(samples) <- eff_groups
          p_plot <- ggplot(reshape2::melt(samples), aes(x = value, fill = Var2)) +
            geom_density(aes(y = after_stat(scaled)), alpha = 0.6, linewidth = 1.2) +
            scale_fill_viridis_d(option = "mako") +
            theme_classic() +
            labs(title = "Global Posterior")
          ggsave(file.path(run_folder, "03_Global_Posterior_Density.png"), plot = p_plot, width = 8, height = 5)

          if (input$opt_b && "beta_sorting" %in% names(rstan::extract(values$fit))) {
            betas <- rstan::extract(values$fit)$beta_sorting
            if (is.null(dim(betas)) || length(dim(betas)) < 2) betas <- matrix(betas, ncol = 1)
            if (ncol(betas) == length(tr)) colnames(betas) <- tr
            else colnames(betas) <- paste0("V", seq_len(ncol(betas)))
            b_plot <- ggplot(reshape2::melt(betas), aes(x = Var2, y = value, fill = Var2)) +
              geom_boxplot(alpha = 0.7, show.legend = FALSE) +
              geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
              theme_minimal()
            ggsave(file.path(run_folder, "04_Beta_Coefficients.png"), plot = b_plot, width = 8, height = 5)
          }

          wall_s <- as.numeric(difftime(t_end, t_start, units = "secs"))
          chain_times <- tryCatch({
            el <- as.data.frame(rstan::get_elapsed_time(values$fit))
            el$chain <- seq_len(nrow(el))
            el$total_s <- round(el$warmup + el$sample, 3)
            el$warmup  <- round(el$warmup, 3)
            el$sample  <- round(el$sample, 3)
            el
          }, error = function(e) NULL)
          if (!is.null(chain_times)) {
            write.csv(chain_times, file.path(run_folder, "05_Chain_Times.csv"), row.names = FALSE)
          }

          agg_df <- build_summary_df(values$fit, eff_groups, p_name, add_phi = (input$pool_mode == "hier"))
          write.csv(agg_df, file.path(run_folder, "06_Aggregated_Summary.csv"), row.names = FALSE)
          src_agg <- agg_df[agg_df$Source %in% eff_groups, , drop = FALSE]
          max_rhat <- tryCatch({x <- max(src_agg$max_Rhat, na.rm = TRUE); if (is.finite(x)) x else NA_real_}, error = function(e) NA_real_)
          min_neff <- tryCatch({x <- min(src_agg$min_neff, na.rm = TRUE); if (is.finite(x)) x else NA_real_}, error = function(e) NA_real_)

          smape_raw <- NA_real_; smape_corr <- NA_real_
          tryCatch({
            rc <- isolate(recon_calc())
            smape_raw  <- mean(rc$sMAPE_Raw,  na.rm = TRUE)
            smape_corr <- mean(rc$sMAPE_Corr, na.rm = TRUE)
          }, error = function(e) NULL)

          meta <- data.frame(
            Model = model_name,
            Dataset = if (!is.null(input$file)) input$file$name else NA_character_,
            Dataset_KB = if (!is.null(input$file)) round(input$file$size / 1024, 1) else NA_real_,
            Pooling = input$pool_mode,
            Covariance = input$cov_mode,
            CoDA = input$opt_coda,
            Beta = input$opt_b,
            Beta_Src = use_b_src_flag,
            Bias = input$opt_bias,
            Ghost = use_ghost_flag,
            Tracers = paste(tr, collapse = ", "),
            n_tracers = length(tr),
            wall_time_s = wall_s,
            mean_chain_s = if (!is.null(chain_times)) mean(chain_times$total_s) else NA_real_,
            min_chain_s  = if (!is.null(chain_times)) min(chain_times$total_s)  else NA_real_,
            max_chain_s  = if (!is.null(chain_times)) max(chain_times$total_s)  else NA_real_,
            max_Rhat = max_rhat,
            min_neff = min_neff,
            mean_sMAPE_raw = smape_raw,
            mean_sMAPE_corr = smape_corr,
            Date = Sys.time()
          )
          num_meta <- sapply(meta, is.numeric)
          meta[num_meta] <- lapply(meta[num_meta], function(x) round(x, 3))
          write.csv(meta, file.path(run_folder, "00_Run_Metadata.csv"), row.names = FALSE)

          runtime_log <- file.path("BASIN_Workspaces", "Runtime_Log.csv")
          if (file.exists(runtime_log)) {
            old_log <- tryCatch(read.csv(runtime_log, check.names = FALSE), error = function(e) NULL)
            if (!is.null(old_log)) {
              common <- intersect(colnames(old_log), colnames(meta))
              meta <- rbind(old_log[, common, drop = FALSE], meta[, common, drop = FALSE])
            }
          }
          write.csv(meta, runtime_log, row.names = FALSE)
        }, error = function(e) {
          showNotification(paste("Auto-export failed:", e$message), type = "warning")
        })

        tryCatch({
          output_dir <- "BASIN_Saved_Results"
          if (!dir.exists(output_dir)) dir.create(output_dir)
          timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
          mix_ids <- if ("ID" %in% colnames(df_mix_rows)) as.character(df_mix_rows$ID) else paste0("Mix_", seq_len(nrow(df_mix_rows)))
          N_mix <- nrow(df_mix_rows)
          auto_df <- data.frame(Sample_ID = mix_ids, check.names = FALSE)

          if (input$pool_mode %in% c("hier", "no_pool")) {
            raw_samples <- rstan::extract(values$fit, pars = "P_ind")$P_ind
            raw_samples <- ensure_3d_samples(raw_samples, length(eff_groups))
            for (i in seq_along(eff_groups)) {
              src_data <- raw_samples[, , i, drop = FALSE]
              auto_df[[paste0(eff_groups[i], " (Mean)")]]   <- as.numeric(colMeans(src_data))
              auto_df[[paste0(eff_groups[i], " (Median)")]] <- as.numeric(apply(src_data, 2, median))
            }
          } else if (use_levels) {
            raw_samples <- rstan::extract(values$fit, pars = "P_level")$P_level
            raw_samples <- ensure_3d_samples(raw_samples, length(eff_groups))
            for (i in seq_along(eff_groups)) {
              src_data <- raw_samples[, , i, drop = FALSE]
              auto_df[[paste0(eff_groups[i], " (Mean)")]]   <- colMeans(src_data)[mix_level_id]
              auto_df[[paste0(eff_groups[i], " (Median)")]] <- apply(src_data, 2, median)[mix_level_id]
            }
          } else {
            raw_samples <- rstan::extract(values$fit, pars = "P")$P
            if (is.vector(raw_samples)) raw_samples <- matrix(raw_samples, ncol = length(eff_groups))
            for (i in seq_along(eff_groups)) {
              auto_df[[paste0(eff_groups[i], " (Mean)")]]   <- rep(mean(raw_samples[, i]), N_mix)
              auto_df[[paste0(eff_groups[i], " (Median)")]] <- rep(median(raw_samples[, i]), N_mix)
            }
          }
          write.csv(auto_df, file.path(output_dir, paste0("Results_", timestamp, ".csv")), row.names = FALSE)
        }, error = function(e) {
          showNotification(paste("Auto-save failed:", e$message), type = "warning")
        })
      }
    })
  })


  #==========================================================================
  # Results plots / tables (Supports G+1 Ghost_Source)
  #==========================================================================

  output$density_plot <- renderPlot({
    req(values$fit)
    if (!isS4(values$fit)) return(NULL)

    eff_groups <- get_effective_groups()
    use_levels <- get_use_levels()

    p_name <- if (input$pool_mode == "hier") {
      "P_global"
    } else if (input$pool_mode == "no_pool") {
      "P_ind"
    } else if (use_levels) {
      "P_level"
    } else {
      "P"
    }

    samples_list <- tryCatch(rstan::extract(values$fit), error = function(e) NULL)
    if (is.null(samples_list) || !(p_name %in% names(samples_list))) return(NULL)

    samples <- samples_list[[p_name]]
    if (p_name %in% c("P_ind", "P_level")) {
      samples <- ensure_3d_samples(samples, length(eff_groups))
    }
    if (length(dim(samples)) == 3) {
      samples <- apply(samples, c(1, 3), mean)
    }
    if (is.null(dim(samples)) || length(dim(samples)) < 2) {
      samples <- matrix(samples, ncol = length(eff_groups))
    }

    colnames(samples) <- eff_groups
    df_p <- reshape2::melt(samples)

    ggplot(df_p, aes(x = value, fill = Var2)) +
      geom_density(aes(y = after_stat(scaled)), alpha = 0.6, linewidth = 1.2) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
      scale_fill_viridis_d(option = "mako") +
      theme_classic(base_size = 20) +
      theme(axis.line = element_line(linewidth = 1.5)) +
      labs(x = "Proportion", y = "Normalized Density (Max=1)", fill = "Source")
  })

  output$ind_density_plot <- renderPlot({
    req(values$fit, input$show_sample_id, values$df, input$group_col, input$mix_label)

    eff_groups <- get_effective_groups()
    df_mix <- values$df[as.character(values$df[[input$group_col]]) == input$mix_label, ]
    mix_ids <- if ("ID" %in% colnames(df_mix)) as.character(df_mix$ID) else paste0("Mix_", seq_len(nrow(df_mix)))

    sample_idx <- which(mix_ids == input$show_sample_id)
    if (length(sample_idx) == 0) return(NULL)

    all_samples <- rstan::extract(values$fit)
    ind_data <- NULL

    if ("P_ind" %in% names(all_samples)) {
      arr <- ensure_3d_samples(all_samples$P_ind, length(eff_groups))
      tmp <- arr[, sample_idx, , drop = FALSE]
      ind_data <- matrix(tmp, nrow = dim(tmp)[1], ncol = dim(tmp)[3])
    } else if ("P_level" %in% names(all_samples) && get_use_levels()) {
      arr <- ensure_3d_samples(all_samples$P_level, length(eff_groups))
      lvl_names <- if (!is.null(values$level_names)) values$level_names else sort(unique(as.character(df_mix[[input$mix_cov]])))
      level_idx <- match(as.character(df_mix[[input$mix_cov]])[sample_idx], lvl_names)
      if (!is.na(level_idx)) {
        tmp <- arr[, level_idx, , drop = FALSE]
        ind_data <- matrix(tmp, nrow = dim(tmp)[1], ncol = dim(tmp)[3])
      }
    } else if ("P" %in% names(all_samples)) {
      ind_data <- all_samples$P
    }

    if (is.null(ind_data)) return(NULL)
    if (is.null(dim(ind_data))) ind_data <- matrix(ind_data, ncol = length(eff_groups))

    colnames(ind_data) <- eff_groups
    df_p <- reshape2::melt(ind_data)
    df_p$Var2 <- factor(df_p$Var2, levels = eff_groups)

    ggplot(df_p, aes(x = value, fill = Var2)) +
      geom_density(aes(y = after_stat(scaled)), alpha = 0.6, linewidth = 0.8) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
      scale_fill_viridis_d(option = "turbo") +
      theme_classic(base_size = 16) +
      labs(
        title = paste("Posterior Density for Sample:", input$show_sample_id),
        x = "Proportion",
        y = "Normalized Density (Max=1)",
        fill = "Source"
      )
  })

  output$download_all_plots <- downloadHandler(
    filename = function() {
      paste0("All_Sample_Plots_", Sys.Date(), ".zip")
    },
    content = function(file) {
      req(values$fit, values$df, input$group_col, input$mix_label)

      eff_groups <- get_effective_groups()
      temp_dir <- tempdir()
      plot_folder <- file.path(temp_dir, "plots")
      unlink(plot_folder, recursive = TRUE, force = TRUE)
      dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

      df_mix <- values$df[as.character(values$df[[input$group_col]]) == input$mix_label, ]
      mix_ids <- if ("ID" %in% colnames(df_mix)) as.character(df_mix$ID) else paste0("Mix_", seq_len(nrow(df_mix)))

      all_samples <- rstan::extract(values$fit)
      use_levels <- get_use_levels()
      cov_factor <- NULL

      if (use_levels) {
        lvl_names <- if (!is.null(values$level_names)) values$level_names else sort(unique(as.character(df_mix[[input$mix_cov]])))
        cov_factor <- match(as.character(df_mix[[input$mix_cov]]), lvl_names)
      }

      withProgress(message = "Generating plots...", value = 0, {
        for (i in seq_along(mix_ids)) {
          ind_data <- NULL
          if ("P_ind" %in% names(all_samples)) {
            arr <- ensure_3d_samples(all_samples$P_ind, length(eff_groups))
            tmp <- arr[, i, , drop = FALSE]
            ind_data <- matrix(tmp, nrow = dim(tmp)[1], ncol = dim(tmp)[3])
          } else if ("P_level" %in% names(all_samples) && use_levels && !is.null(cov_factor[i])) {
            arr <- ensure_3d_samples(all_samples$P_level, length(eff_groups))
            tmp <- arr[, cov_factor[i], , drop = FALSE]
            ind_data <- matrix(tmp, nrow = dim(tmp)[1], ncol = dim(tmp)[3])
          } else if ("P" %in% names(all_samples)) {
            ind_data <- all_samples$P
          }

          if (is.null(dim(ind_data))) ind_data <- matrix(ind_data, ncol = length(eff_groups))
          colnames(ind_data) <- eff_groups
          df_p <- reshape2::melt(ind_data)
          df_p$Var2 <- factor(df_p$Var2, levels = eff_groups)

          p <- ggplot(df_p, aes(x = value, fill = Var2)) +
            geom_density(aes(y = after_stat(scaled)), alpha = 0.6) +
            scale_x_continuous(limits = c(0, 1)) +
            scale_fill_viridis_d(option = "turbo") +
            theme_classic() +
            labs(title = paste("Sample:", mix_ids[i]), x = "Proportion", y = "Density")

          safe_id <- sanitize_filename(mix_ids[i])
          ggsave(file.path(plot_folder, paste0(safe_id, ".png")), plot = p, width = 6, height = 4, dpi = 150)
          incProgress(1 / length(mix_ids))
        }
      })

      zip::zip(zipfile = file, files = list.files(plot_folder, full.names = FALSE), root = plot_folder)
    }
  )

  output$beta_plot <- renderPlot({
    req(values$fit, input$opt_b)
    samples_all <- rstan::extract(values$fit)
    if (!"beta_sorting" %in% names(samples_all)) return(NULL)

    betas <- samples_all$beta_sorting
    if (is.null(dim(betas)) || length(dim(betas)) < 2) betas <- matrix(betas, ncol = 1)
    if (ncol(betas) == length(input$final_selected_tracers)) {
      colnames(betas) <- input$final_selected_tracers
    } else {
      colnames(betas) <- paste0("V", seq_len(ncol(betas)))
    }

    df_b <- reshape2::melt(betas)
    ggplot(df_b, aes(x = Var2, y = value, fill = Var2)) +
      geom_boxplot(alpha = 0.7, show.legend = FALSE) +
      geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
      theme_minimal(base_size = 18) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      labs(title = "Beta Coefficients (Particle Size Correction)", x = "Tracer", y = "Beta Value")
  })

  output$summary_table <- renderTable({
    req(values$fit)
    eff_groups <- get_effective_groups()
    use_levels <- get_use_levels()

    p_param <- if (input$pool_mode == "hier") {
      "P_global"
    } else if (input$pool_mode == "no_pool") {
      "P_ind"
    } else if (use_levels) {
      "P_level"
    } else {
      "P"
    }

    res_p <- as.data.frame(summary(values$fit, pars = p_param)$summary)
    samples_p <- rstan::extract(values$fit, pars = p_param)[[1]]

    if (p_param %in% c("P_ind", "P_level")) {
      samples_p <- ensure_3d_samples(samples_p, length(eff_groups))
    }
    if (length(dim(samples_p)) == 3) {
      samples_p <- apply(samples_p, c(1, 3), mean)
    }
    if (is.null(dim(samples_p))) {
      samples_p <- matrix(samples_p, ncol = length(eff_groups))
    }

    q025 <- apply(samples_p, 2, quantile, probs = 0.025)
    q25  <- apply(samples_p, 2, quantile, probs = 0.25)
    q75  <- apply(samples_p, 2, quantile, probs = 0.75)
    q975 <- apply(samples_p, 2, quantile, probs = 0.975)

    medians_p <- apply(samples_p, 2, median)
    median_p_corr <- medians_p / sum(medians_p)

    modes_p <- apply(samples_p, 2, function(x) {
      d <- density(x)
      d$x[which.max(d$y)]
    })
    modes_p_corr <- modes_p / sum(modes_p)

    df_final <- data.frame(
      Source = eff_groups,
      mean = colMeans(samples_p),
      sd = apply(samples_p, 2, sd),
      `2.5%` = q025,
      `25% (Q1)` = q25,
      median = medians_p,
      `median corr` = median_p_corr,
      mode = modes_p,
      `mode corr` = modes_p_corr,
      `75% (Q3)` = q75,
      `97.5%` = q975,
      `W50 (IQR)` = q75 - q25,
      check.names = FALSE
    )

    if ("n_eff" %in% colnames(res_p)) {
      if (length(res_p$Rhat) > length(eff_groups)) {
        rhat_mat <- matrix(res_p$Rhat, ncol = length(eff_groups), byrow = TRUE)
        neff_mat <- matrix(res_p$n_eff, ncol = length(eff_groups), byrow = TRUE)
        df_final$n_eff <- apply(neff_mat, 2, min, na.rm = TRUE)
        df_final$Rhat  <- apply(rhat_mat, 2, max, na.rm = TRUE)
      } else {
        df_final$n_eff <- res_p$n_eff[seq_len(length(eff_groups))]
        df_final$Rhat  <- res_p$Rhat[seq_len(length(eff_groups))]
      }
    }

    if (input$pool_mode == "hier") {
      res_phi <- as.data.frame(summary(values$fit, pars = "phi_ind")$summary)
      samples_phi <- rstan::extract(values$fit, pars = "phi_ind")[[1]]
      q25phi <- quantile(samples_phi, 0.25)
      q75phi <- quantile(samples_phi, 0.75)

      df_phi <- data.frame(
        Source = "Phi (Precision)",
        mean = mean(samples_phi),
        sd = sd(samples_phi),
        `2.5%` = quantile(samples_phi, 0.025),
        `25% (Q1)` = q25phi,
        median = median(samples_phi),
        `median corr` = NA,
        mode = density(samples_phi)$x[which.max(density(samples_phi)$y)],
        `mode corr` = NA,
        `75% (Q3)` = q75phi,
        `97.5%` = quantile(samples_phi, 0.975),
        `W50 (IQR)` = q75phi - q25phi,
        check.names = FALSE
      )
      if ("n_eff" %in% colnames(res_phi)) {
        df_phi$n_eff <- res_phi$n_eff
        df_phi$Rhat  <- res_phi$Rhat
      }
      df_final <- rbind(df_final, df_phi)
    }

    return(df_final)
  }, digits = 3, na = "-")

  output$download_ind_res <- downloadHandler(
    filename = function() {
      paste0("Individual_Results_Detailed_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$fit, values$df)
      if (!isS4(values$fit)) return(NULL)

      eff_groups <- get_effective_groups()
      df_mix_rows <- values$df[as.character(values$df[[input$group_col]]) == input$mix_label, ]
      mix_ids <- if ("ID" %in% colnames(df_mix_rows)) as.character(df_mix_rows$ID) else paste0("Mix_", seq_len(nrow(df_mix_rows)))
      N_mix <- nrow(df_mix_rows)
      final_df <- data.frame(Sample_ID = mix_ids, check.names = FALSE)

      tryCatch({
        use_levels <- get_use_levels(df_mix_rows)

        if (input$pool_mode %in% c("hier", "no_pool")) {
          raw_samples <- rstan::extract(values$fit, pars = "P_ind")$P_ind
          raw_samples <- ensure_3d_samples(raw_samples, length(eff_groups))

          for (i in seq_along(eff_groups)) {
            src_matrix <- matrix(raw_samples[, , i, drop = FALSE], nrow = dim(raw_samples)[1], ncol = N_mix)
            q25 <- as.numeric(apply(src_matrix, 2, quantile, probs = 0.25))
            q75 <- as.numeric(apply(src_matrix, 2, quantile, probs = 0.75))

            final_df[[paste0(eff_groups[i], " (Mean)")]]   <- as.numeric(colMeans(src_matrix))
            final_df[[paste0(eff_groups[i], " (Median)")]] <- as.numeric(apply(src_matrix, 2, median))
            final_df[[paste0(eff_groups[i], " (2.5%)")]]   <- as.numeric(apply(src_matrix, 2, quantile, probs = 0.025))
            final_df[[paste0(eff_groups[i], " (97.5%)")]]  <- as.numeric(apply(src_matrix, 2, quantile, probs = 0.975))
            final_df[[paste0(eff_groups[i], " (25%)")]]    <- q25
            final_df[[paste0(eff_groups[i], " (75%)")]]    <- q75
            final_df[[paste0(eff_groups[i], " (W50)")]]    <- q75 - q25
          }
        } else if (use_levels) {
          raw_samples <- rstan::extract(values$fit, pars = "P_level")$P_level
          raw_samples <- ensure_3d_samples(raw_samples, length(eff_groups))
          lvl_names <- if (!is.null(values$level_names)) values$level_names else sort(unique(as.character(df_mix_rows[[input$mix_cov]])))
          cov_factor <- match(as.character(df_mix_rows[[input$mix_cov]]), lvl_names)

          for (i in seq_along(eff_groups)) {
            src_data <- raw_samples[, , i, drop = FALSE]
            final_df[[paste0(eff_groups[i], " (Mean)")]]   <- colMeans(src_data)[cov_factor]
            final_df[[paste0(eff_groups[i], " (Median)")]] <- apply(src_data, 2, median)[cov_factor]
          }
          final_df$Note <- "Level-Specific Model"
        } else {
          raw_samples <- rstan::extract(values$fit, pars = "P")$P
          if (is.vector(raw_samples)) raw_samples <- matrix(raw_samples, ncol = length(eff_groups))

          for (i in seq_along(eff_groups)) {
            src_vec <- raw_samples[, i]
            q25 <- quantile(src_vec, 0.25)
            q75 <- quantile(src_vec, 0.75)

            final_df[[paste0(eff_groups[i], " (Mean)")]]   <- rep(mean(src_vec), N_mix)
            final_df[[paste0(eff_groups[i], " (Median)")]] <- rep(median(src_vec), N_mix)
            final_df[[paste0(eff_groups[i], " (2.5%)")]]   <- rep(quantile(src_vec, 0.025), N_mix)
            final_df[[paste0(eff_groups[i], " (97.5%)")]]  <- rep(quantile(src_vec, 0.975), N_mix)
            final_df[[paste0(eff_groups[i], " (25%)")]]    <- rep(q25, N_mix)
            final_df[[paste0(eff_groups[i], " (75%)")]]    <- rep(q75, N_mix)
            final_df[[paste0(eff_groups[i], " (W50)")]]    <- rep(q75 - q25, N_mix)
          }
          final_df$Note <- "Static Model"
        }
        write.csv(final_df, file, row.names = FALSE)
      }, error = function(e) {
        showNotification(paste("Download Error:", e$message), type = "error")
      })
    }
  )

  output$download_chains <- downloadHandler(
    filename = function() {
      paste0("MCMC_Raw_Chains_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$fit)
      write.csv(as.data.frame(values$fit), file, row.names = FALSE)
    }
  )

  #==========================================================================
  # Reconstruction Analysis (Supports Ghost_Source G+1, Beta & Nested Hotspots)
  #==========================================================================

  output$reconstruction_ui <- renderUI({
    req(values$fit)
    dynamic_height <- max(450, ceiling(length(input$final_selected_tracers) / 3) * 280)
    size_ok <- !is.null(input$size_col) && input$size_col != "None"

    fit_names <- tryCatch(names(rstan::extract(values$fit)), error = function(e) character(0))
    has_beta   <- isTRUE(input$opt_b && size_ok && ("beta_sorting" %in% fit_names))
    has_nested <- isTRUE(input$hotspot_mode == "nested" && ("eta_sub" %in% fit_names))

    tagList(
      br(),
      h4("Multi-Tracer Mass-Balance Reconstruction"),
      hr(),
      fluidRow(
        column(
          4,
          if (has_beta) {
            radioButtons(
              "recon_mode",
              "1. Sorting Correction Mode:",
              choices = c(
                "Conservative (No Beta)" = "raw",
                "Corrected (+Beta SSA)" = "corr"
              ),
              selected = "corr",
              inline = FALSE
            )
          } else {
            p(tags$b("1. Size-Sorting Mode: "), "Conservative (Beta disabled)")
          }
        ),
        column(
          5,
          if (has_nested) {
            radioButtons(
              "recon_src_mode",
              "2. Source End-Member Basis (Nested Mode):",
              choices = c(
                "Group Centroids (Uniform Group Mean)" = "centroid",
                "Active Hotspots (eta_sub Sample Simplex)" = "nested"
              ),
              selected = "nested",
              inline = FALSE
            )
          } else {
            p(tags$b("2. Source End-Member Basis: "), "Group Centroids (Uniform Group Mean)")
          }
        ),
        column(
          3,
          div(
            style = "text-align:right; margin-top: 20px;",
            downloadButton("download_recon_csv", "Download Full Details", class = "btn-info")
          )
        )
      ),
      br(),
      fluidRow(
        column(8, plotOutput("recon_plot", height = paste0(dynamic_height, "px"))),
        column(4, h4("Error per Tracer (%)"), tableOutput("recon_metrics_table"))
      ),
      hr(),
      h4("Sample-specific reconstruction performance"),
      plotOutput("recon_error_boxplot", height = paste0(dynamic_height, "px"))
    )
  })

  recon_calc <- reactive({
    req(values$fit, values$df, input$final_selected_tracers, values$groups)

    tr <- input$final_selected_tracers
    gr_col <- input$group_col
    m_lab <- input$mix_label
    groups <- values$groups
    eff_groups <- get_effective_groups()

    df_mix_rows <- values$df[as.character(values$df[[gr_col]]) == m_lab, ]
    df_src <- values$df[as.character(values$df[[gr_col]]) != m_lab, ]
    N_mix <- nrow(df_mix_rows)
    N_src <- nrow(df_src)
    mix_ids <- if ("ID" %in% colnames(df_mix_rows)) as.character(df_mix_rows$ID) else paste0("Mix_", seq_len(N_mix))

    use_levels <- get_use_levels(df_mix_rows)
    p_param <- if (input$pool_mode %in% c("hier", "no_pool")) "P_ind" else if (use_levels) "P_level" else "P"

    samples <- rstan::extract(values$fit)
    p_samples <- samples[[p_param]]

    if (p_param %in% c("P_ind", "P_level")) {
      p_samples <- ensure_3d_samples(p_samples, length(eff_groups))
    }

    if (p_param == "P_ind") {
      medians <- matrix(apply(p_samples, c(2, 3), median), ncol = length(eff_groups))
    } else if (p_param == "P_level") {
      level_med <- matrix(apply(p_samples, c(2, 3), median), ncol = length(eff_groups))
      cov_vec <- as.character(df_mix_rows[[input$mix_cov]])
      cov_vec[is.na(cov_vec)] <- "NA"
      lvl_names <- if (!is.null(values$level_names)) values$level_names else sort(unique(cov_vec))
      idx <- match(cov_vec, lvl_names)
      idx[is.na(idx)] <- 1
      medians <- level_med[idx, , drop = FALSE]
    } else {
      if (is.null(dim(p_samples))) p_samples <- matrix(p_samples, ncol = length(eff_groups))
      med <- apply(p_samples, 2, median)
      medians <- matrix(rep(med, N_mix), ncol = length(eff_groups), byrow = TRUE)
    }

    row_sums <- rowSums(medians)
    row_sums[row_sums == 0] <- 1
    medians_corr <- medians / row_sums

    # 1. Классическая матрица групповых средних (Group Centroids) [N_eff_groups x N_tr]
    src_means_mat <- matrix(0, nrow = length(eff_groups), ncol = length(tr))
    src_range <- numeric(length(tr))

    for (i in seq_along(tr)) {
      t_vals <- df_src[[tr[i]]]
      rng <- if (all(is.na(t_vals))) 1 else max(t_vals, na.rm = TRUE) - min(t_vals, na.rm = TRUE)
      if (is.na(rng) || rng <= 0) rng <- 1
      src_range[i] <- rng

      for (g in seq_along(groups)) {
        src_means_mat[g, i] <- mean(df_src[as.character(df_src[[gr_col]]) == groups[g], tr[i]], na.rm = TRUE)
      }
    }

    # Если активен Ghost Source, добавляем его реконструированные физические концентрации в (G+1)-ю строку
    if (isTRUE(input$opt_ghost) && "mu_ghost" %in% names(samples)) {
      mu_gh_med <- apply(samples$mu_ghost, 2, median)
      if (!input$opt_coda && !is.null(values$y_center) && !is.null(values$y_scale)) {
        src_means_mat[length(eff_groups), ] <- values$y_center + values$y_scale * mu_gh_med
      } else if (input$opt_coda && !is.null(values$stan_data$V)) {
        V_mat <- values$stan_data$V
        gh_tot <- mean(values$stan_data$src_totals)
        gh_raw_vec <- as.numeric(V_mat %*% mu_gh_med)
        gh_prop <- exp(gh_raw_vec - max(gh_raw_vec))
        gh_prop <- gh_prop / sum(gh_prop)
        src_means_mat[length(eff_groups), ] <- gh_prop * gh_tot
      }
    }
    src_means_mat[is.na(src_means_mat)] <- 0

    # 2. 3D-массив эффективных концентраций для режима Nested Hotspots [N_mix x N_eff_groups x N_tr]
    src_nested_arr <- array(0, dim = c(N_mix, length(eff_groups), length(tr)))
    for (i in seq_len(N_mix)) {
      src_nested_arr[i, , ] <- src_means_mat
    }

    has_nested <- isTRUE(input$hotspot_mode == "nested" && ("eta_sub" %in% names(samples)))
    if (has_nested) {
      eta_draws <- samples$eta_sub # [iterations, N_mix, N_src]
      if (length(dim(eta_draws)) == 3 && dim(eta_draws)[2] == N_mix && dim(eta_draws)[3] == N_src) {
        eta_med <- apply(eta_draws, c(2, 3), median) # [N_mix x N_src]
      } else {
        eta_med <- matrix(0, nrow = N_mix, ncol = N_src)
      }

      src_raw_mat <- as.matrix(df_src[, tr, drop = FALSE])
      src_raw_mat[is.na(src_raw_mat)] <- 0
      src_labels_chr <- as.character(df_src[[gr_col]])

      for (i in seq_len(N_mix)) {
        for (g in seq_along(groups)) {
          s_idx <- which(src_labels_chr == groups[g])
          if (length(s_idx) > 0) {
            e_vec <- eta_med[i, s_idx]
            w_sub <- exp(e_vec - max(e_vec))
            w_sub <- w_sub / sum(w_sub)
            src_nested_arr[i, g, ] <- as.numeric(crossprod(w_sub, src_raw_mat[s_idx, , drop = FALSE]))
          }
        }
      }
    }

    obs_conc <- as.matrix(df_mix_rows[, tr, drop = FALSE])
    obs_conc[is.na(obs_conc)] <- 0

    size_active <- !is.null(input$size_col) &&
      input$size_col != "None" &&
      input$size_col %in% colnames(df_mix_rows) &&
      input$size_col %in% colnames(df_src)

    has_beta <- input$opt_b && size_active && "beta_sorting" %in% names(samples)
    beta_med <- rep(0, length(tr))
    if (has_beta) {
      beta_tmp <- apply(samples$beta_sorting, 2, median)
      beta_med <- if (length(beta_tmp) >= length(tr)) beta_tmp[seq_along(tr)] else c(beta_tmp, rep(0, length(tr) - length(beta_tmp)))
    }

    delta_med <- rep(0, length(tr))
    if (input$opt_bias && "delta_bias" %in% names(samples)) {
      delta_tmp <- apply(samples$delta_bias, 2, median)
      delta_med <- if (length(delta_tmp) >= length(tr)) delta_tmp[seq_along(tr)] else c(delta_tmp, rep(0, length(tr) - length(delta_tmp)))

      if (isTRUE(input$opt_ghost) && "mu_ghost" %in% names(samples) && !input$opt_coda) {
        mu_src_med <- apply(samples$mu_source, c(2, 3), median)
        mu_bar_vec <- colMeans(matrix(mu_src_med, ncol = length(tr)))
        v_gh_vec   <- apply(samples$mu_ghost, 2, median) - mu_bar_vec
        delta_med  <- delta_med - (sum(delta_med * v_gh_vec) / (sum(v_gh_vec^2) + 1e-9)) * v_gh_vec
      }
    }

    beta_ok <- FALSE
    ssa_m_z <- rep(0, N_mix)
    ssa_grp_z <- rep(0, length(eff_groups))

    if (has_beta) {
      ssa_m_raw <- df_mix_rows[[input$size_col]]
      ssa_s_raw <- df_src[[input$size_col]]
      ssa_s_raw[is.na(ssa_s_raw)] <- mean(ssa_s_raw, na.rm = TRUE)
      ssa_m_raw[is.na(ssa_m_raw)] <- mean(ssa_s_raw, na.rm = TRUE)
      mean_s <- mean(ssa_s_raw, na.rm = TRUE)
      sd_s   <- sd(ssa_s_raw, na.rm = TRUE)

      if (!is.na(sd_s) && sd_s >= 1e-8) {
        ssa_m_z <- (ssa_m_raw - mean_s) / sd_s
        ssa_s_z <- (ssa_s_raw - mean_s) / sd_s
        src_labels <- as.character(df_src[[gr_col]])
        ssa_grp_z_known <- sapply(groups, function(g) {
          x <- ssa_s_z[src_labels == g]
          if (length(x) == 0 || all(is.na(x))) 0 else mean(x, na.rm = TRUE)
        })
        ssa_grp_z <- if (isTRUE(input$opt_ghost)) c(ssa_grp_z_known, 0) else ssa_grp_z_known
        beta_ok <- TRUE
      }
    }

    # Вспомогательная функция расчета предсказаний для выбранного базиса источников (Centroid или Nested)
    compute_predictions <- function(use_nested_basis = FALSE) {
      p_raw_mat  <- matrix(0, nrow = N_mix, ncol = length(tr))
      p_corr_mat <- matrix(0, nrow = N_mix, ncol = length(tr))

      for (i in seq_len(N_mix)) {
        M_i <- if (use_nested_basis) src_nested_arr[i, , ] else src_means_mat
        if (is.null(dim(M_i))) M_i <- matrix(M_i, nrow = length(eff_groups), ncol = length(tr))

        base_i <- as.numeric(medians_corr[i, ] %*% M_i)
        if (!input$opt_coda && !is.null(values$y_scale) && length(values$y_scale) == length(tr)) {
          p_raw_mat[i, ] <- base_i + values$y_scale * delta_med
        } else {
          p_raw_mat[i, ] <- base_i
        }

        if (has_beta && beta_ok) {
          for (t_idx in seq_along(tr)) {
            raw_gt_pos <- pmax(1e-6, M_i[, t_idx])
            corr_val <- sum(medians_corr[i, ] * raw_gt_pos * exp(beta_med[t_idx] * (ssa_m_z[i] - ssa_grp_z)))
            if (!input$opt_coda && !is.null(values$y_scale) && length(values$y_scale) == length(tr)) {
              p_corr_mat[i, t_idx] <- corr_val + values$y_scale[t_idx] * delta_med[t_idx]
            } else {
              p_corr_mat[i, t_idx] <- corr_val
            }
          }
        } else {
          p_corr_mat[i, ] <- p_raw_mat[i, ]
        }
      }
      list(raw = p_raw_mat, corr = p_corr_mat)
    }

    res_cent   <- compute_predictions(use_nested_basis = FALSE)
    res_nested <- compute_predictions(use_nested_basis = has_nested)

    res_list <- list()
    for (t_idx in seq_along(tr)) {
      o    <- obs_conc[, t_idx]
      pr_c <- res_cent$raw[, t_idx]
      pc_c <- res_cent$corr[, t_idx]
      pr_n <- res_nested$raw[, t_idx]
      pc_n <- res_nested$corr[, t_idx]

      res_list[[t_idx]] <- data.frame(
        Sample = mix_ids,
        Tracer = tr[t_idx],
        Observed = o,
        # 1. Групповые средние (Group Centroid)
        Pred_Raw = pr_c,
        Pred_Corr = pc_c,
        Bias_Raw_pct = ((pr_c - o) / (abs(o) + 1e-9)) * 100,
        Bias_Corr_pct = ((pc_c - o) / (abs(o) + 1e-9)) * 100,
        sMAPE_Raw = (abs(pr_c - o) / ((abs(pr_c) + abs(o) + 1e-9) / 2)) * 100,
        sMAPE_Corr = (abs(pc_c - o) / ((abs(pc_c) + abs(o) + 1e-9) / 2)) * 100,
        NRMSE_src_Raw = (abs(pr_c - o) / (src_range[t_idx] + 1e-9)) * 100,
        NRMSE_src_Corr = (abs(pc_c - o) / (src_range[t_idx] + 1e-9)) * 100,
        # 2. Вложенные хотспоты (Nested Sample Simplex)
        Pred_Nested_Raw = pr_n,
        Pred_Nested_Corr = pc_n,
        Bias_Nested_Raw_pct = ((pr_n - o) / (abs(o) + 1e-9)) * 100,
        Bias_Nested_Corr_pct = ((pc_n - o) / (abs(o) + 1e-9)) * 100,
        sMAPE_Nested_Raw = (abs(pr_n - o) / ((abs(pr_n) + abs(o) + 1e-9) / 2)) * 100,
        sMAPE_Nested_Corr = (abs(pc_n - o) / ((abs(pc_n) + abs(o) + 1e-9) / 2)) * 100,
        NRMSE_src_Nested_Raw = (abs(pr_n - o) / (src_range[t_idx] + 1e-9)) * 100,
        NRMSE_src_Nested_Corr = (abs(pc_n - o) / (src_range[t_idx] + 1e-9)) * 100,
        # Сдвиги за счет сортировки (Beta) и за счет селективности хотспотов (Nested)
        Beta_Shift_pct = if (has_beta && beta_ok) ((pc_c - pr_c) / (abs(pr_c) + 1e-9)) * 100 else 0,
        Hotspot_Shift_pct = if (has_nested) ((pc_n - pc_c) / (abs(pc_c) + 1e-9)) * 100 else 0,
        check.names = FALSE
      )
    }

    df_final <- do.call(rbind, res_list)
    clamp_cols <- c(
      "Bias_Raw_pct", "Bias_Corr_pct",
      "Bias_Nested_Raw_pct", "Bias_Nested_Corr_pct",
      "Beta_Shift_pct", "Hotspot_Shift_pct"
    )
    for (cc in clamp_cols) {
      if (cc %in% colnames(df_final)) {
        df_final[[cc]][df_final[[cc]] > 500]  <- 500
        df_final[[cc]][df_final[[cc]] < -500] <- -500
      }
    }
    return(df_final)
  })

  active_recon_data <- reactive({
    df <- recon_calc()
    mode <- if (is.null(input$recon_mode)) "corr" else input$recon_mode
    src_basis <- if (!is.null(input$recon_src_mode) && input$hotspot_mode == "nested") {
      input$recon_src_mode
    } else {
      "centroid"
    }

    res <- data.frame(Sample = df$Sample, Tracer = df$Tracer, Observed = df$Observed)

    if (src_basis == "nested") {
      if (mode == "corr") {
        res$Predicted <- df$Pred_Nested_Corr; res$Bias_pct <- df$Bias_Nested_Corr_pct
        res$sMAPE <- df$sMAPE_Nested_Corr;     res$NRMSE_src <- df$NRMSE_src_Nested_Corr
      } else {
        res$Predicted <- df$Pred_Nested_Raw;  res$Bias_pct <- df$Bias_Nested_Raw_pct
        res$sMAPE <- df$sMAPE_Nested_Raw;     res$NRMSE_src <- df$NRMSE_src_Nested_Raw
      }
    } else {
      if (mode == "corr") {
        res$Predicted <- df$Pred_Corr; res$Bias_pct <- df$Bias_Corr_pct
        res$sMAPE <- df$sMAPE_Corr;     res$NRMSE_src <- df$NRMSE_src_Corr
      } else {
        res$Predicted <- df$Pred_Raw;  res$Bias_pct <- df$Bias_Raw_pct
        res$sMAPE <- df$sMAPE_Raw;     res$NRMSE_src <- df$NRMSE_src_Raw
      }
    }
    return(res)
  })

  output$recon_error_boxplot <- renderPlot({
    req(active_recon_data())
    ggplot(active_recon_data(), aes(x = Tracer, y = Bias_pct, fill = Tracer)) +
      geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +
      geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.fill = "white") +
      theme_minimal(base_size = 16) +
      scale_fill_viridis_d(option = "mako", guide = "none") +
      labs(title = "Reconstruction Bias per Tracer", y = "Bias (%)", x = "") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x = element_blank()) +
      coord_cartesian(ylim = c(-150, 150))
  })

  output$recon_metrics_table <- renderTable({
    active_recon_data() %>%
      group_by(Tracer) %>%
      summarise(
        `Mean Bias (%)` = mean(Bias_pct, na.rm = TRUE),
        `sMAPE (%)` = mean(sMAPE, na.rm = TRUE),
        `Error vs Src Range (%)` = mean(NRMSE_src, na.rm = TRUE),
        .groups = "drop"
      )
  }, digits = 1)

  output$recon_plot <- renderPlot({
    req(active_recon_data())
    ggplot(active_recon_data(), aes(x = Observed, y = Predicted, color = Tracer)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
      geom_point(size = 3, alpha = 0.7) +
      facet_wrap(~Tracer, scales = "free") +
      theme_bw(base_size = 14) +
      labs(title = "Mass-Balance Tracer Reconstruction", x = "Observed Tracer Value", y = "Reconstructed Tracer Value")
  })

  output$download_recon_csv <- downloadHandler(
    filename = function() paste0("Reconstruction_Details_", Sys.Date(), ".csv"),
    content = function(file) { req(recon_calc()); write.csv(recon_calc(), file, row.names = FALSE) }
  )


  #==========================================================================
  # NEW TAB: GHOST SOURCE DETECTOR (G+1 SEMI-BLIND UNMIXING DIAGNOSTICS)
  #==========================================================================

  ghost_calc <- reactive({
    req(values$fit, values$df, input$final_selected_tracers, values$groups)
    if (!isTRUE(input$opt_ghost)) return(NULL)

    samples <- tryCatch(rstan::extract(values$fit), error = function(e) NULL)
    if (is.null(samples) || !("mu_ghost" %in% names(samples))) return(NULL)

    tr <- input$final_selected_tracers
    groups <- values$groups
    eff_groups <- get_effective_groups()
    g_ghost_idx <- length(eff_groups)

    df_mix_rows <- values$df[as.character(values$df[[input$group_col]]) == input$mix_label, ]
    df_src <- values$df[as.character(values$df[[input$group_col]]) != input$mix_label, ]
    N_mix <- nrow(df_mix_rows)
    mix_ids <- if ("ID" %in% colnames(df_mix_rows)) as.character(df_mix_rows$ID) else paste0("Mix_", seq_len(N_mix))

    # 1. Извлекаем выборки пропорций Ghost_Source по каждой пробе и глобально
    use_levels <- get_use_levels(df_mix_rows)
    if ("P_ind" %in% names(samples)) {
      p_arr <- ensure_3d_samples(samples$P_ind, length(eff_groups))
      p_ghost_mat <- matrix(p_arr[, , g_ghost_idx], nrow = dim(p_arr)[1], ncol = N_mix)
    } else if ("P_level" %in% names(samples) && use_levels) {
      p_arr <- ensure_3d_samples(samples$P_level, length(eff_groups))
      cov_vec <- as.character(df_mix_rows[[input$mix_cov]])
      cov_vec[is.na(cov_vec)] <- "NA"
      lvl_names <- if (!is.null(values$level_names)) values$level_names else sort(unique(cov_vec))
      idx <- match(cov_vec, lvl_names); idx[is.na(idx)] <- 1
      p_ghost_mat <- matrix(p_arr[, idx, g_ghost_idx], nrow = dim(p_arr)[1], ncol = N_mix)
    } else {
      p_raw <- samples$P
      if (is.null(dim(p_raw))) p_raw <- matrix(p_raw, ncol = length(eff_groups))
      p_ghost_mat <- matrix(rep(p_raw[, g_ghost_idx], N_mix), nrow = nrow(p_raw), ncol = N_mix)
    }

    p_ghost_global_draws <- rowMeans(p_ghost_mat)

    sample_df <- data.frame(
      Sample_ID     = factor(mix_ids, levels = unique(mix_ids)),
      Mean          = colMeans(p_ghost_mat),
      Median        = apply(p_ghost_mat, 2, median),
      q025          = apply(p_ghost_mat, 2, quantile, probs = 0.025),
      q25           = apply(p_ghost_mat, 2, quantile, probs = 0.25),
      q75           = apply(p_ghost_mat, 2, quantile, probs = 0.75),
      q975          = apply(p_ghost_mat, 2, quantile, probs = 0.975),
      Prob_gt_5pct  = colMeans(p_ghost_mat > 0.05),
      Prob_gt_10pct = colMeans(p_ghost_mat > 0.10),
      check.names   = FALSE
    )

    # 2. Реконструкция геохимического профиля Ghost_Source (в Z-шкале и физических единицах)
    mu_gh_draws <- samples$mu_ghost
    if (is.null(dim(mu_gh_draws)) || length(dim(mu_gh_draws)) < 2) {
      mu_gh_draws <- matrix(mu_gh_draws, ncol = length(tr))
    }

    src_raw_mat <- as.matrix(df_src[, tr, drop = FALSE])
    sc_ref <- safe_scale_matrix(src_raw_mat)

    if (!input$opt_coda && !is.null(values$y_center) && !is.null(values$y_scale)) {
      z_draws <- as.matrix(mu_gh_draws)
      phys_draws <- sweep(sweep(z_draws, 2, values$y_scale, "*"), 2, values$y_center, "+")
      phys_draws <- as.matrix(phys_draws)
      # Обнуляем отрицательные значения только для неотрицательных трассеров, сохраняя dim матрицы
      for (j in seq_len(ncol(phys_draws))) {
        if (min(src_raw_mat[, j], na.rm = TRUE) >= 0) {
          phys_draws[phys_draws[, j] < 0, j] <- 0
        }
      }
    } else if (input$opt_coda && !is.null(values$stan_data$V)) {
      V_mat <- values$stan_data$V
      gh_tot <- mean(values$stan_data$src_totals)
      phys_draws <- matrix(0, nrow = nrow(mu_gh_draws), ncol = length(tr))
      for (it in seq_len(nrow(mu_gh_draws))) {
        v_raw <- as.numeric(V_mat %*% mu_gh_draws[it, ])
        p_tr <- exp(v_raw - max(v_raw))
        p_tr <- p_tr / sum(p_tr)
        phys_draws[it, ] <- p_tr * gh_tot
      }
      z_draws <- as.matrix(scale(phys_draws, center = sc_ref$center, scale = sc_ref$scale))
    } else {
      z_draws <- as.matrix(mu_gh_draws)
      phys_draws <- as.matrix(mu_gh_draws)
    }

    sig_df <- data.frame(
      Tracer          = tr,
      Z_Median        = apply(z_draws, 2, median),
      Z_q025          = apply(z_draws, 2, quantile, probs = 0.025),
      Z_q25           = apply(z_draws, 2, quantile, probs = 0.25),
      Z_q75           = apply(z_draws, 2, quantile, probs = 0.75),
      Z_q975          = apply(z_draws, 2, quantile, probs = 0.975),
      Ghost_Phys_Med  = apply(phys_draws, 2, median),
      Ghost_Phys_2.5  = apply(phys_draws, 2, quantile, probs = 0.025),
      Ghost_Phys_97.5 = apply(phys_draws, 2, quantile, probs = 0.975),
      check.names     = FALSE
    )

    for (g in groups) {
      sig_df[[paste0("Mean_", g)]] <- colMeans(df_src[as.character(df_src[[input$group_col]]) == g, tr, drop = FALSE], na.rm = TRUE)
    }

    list(
      global_med     = median(p_ghost_global_draws),
      global_mean    = mean(p_ghost_global_draws),
      global_q025    = as.numeric(quantile(p_ghost_global_draws, 0.025)),
      global_q975    = as.numeric(quantile(p_ghost_global_draws, 0.975)),
      prob_gt_5pct   = mean(p_ghost_global_draws > 0.05),
      prob_gt_10pct  = mean(p_ghost_global_draws > 0.10),
      max_sample_med = max(sample_df$Median, na.rm = TRUE),
      max_sample_id  = as.character(sample_df$Sample_ID[which.max(sample_df$Median)]),
      sample_df      = sample_df,
      sig_df         = sig_df
    )
  })

  output$ghost_detector_ui <- renderUI({
    if (!isTRUE(input$opt_ghost)) {
      return(
        tags$div(
          style = "padding: 25px; background-color: #f8f9fa; border-left: 5px solid #3498db; margin-top: 20px;",
          h4("Ghost Source Detector is currently disabled"),
          p("To screen for unsampled erosion hotspots (semi-blind G+1 unmixing):"),
          tags$ol(
            tags$li("Check ", tags$b("Detect Unsampled 'Ghost' Source (G+1)"), " in Section 3 of the left sidebar."),
            tags$li("Adjust the Dirichlet sparsity parameter ", tags$b("alpha"), " (default 0.12 strongly penalizes false alarms)."),
            tags$li("Click ", tags$b("5. Compile and run the model"), ".")
          )
        )
      )
    }

    req(values$fit)
    gc <- ghost_calc()
    req(gc)

    # Формируем автоматический геоморфологический диагноз
    med_pct  <- round(gc$global_med * 100, 1)
    max_pct  <- round(gc$max_sample_med * 100, 1)
    p5_pct   <- round(gc$prob_gt_5pct * 100, 1)
    p10_pct  <- round(gc$prob_gt_10pct * 100, 1)

    enriched_tr <- gc$sig_df$Tracer[gc$sig_df$Z_Median > 0.5]
    depleted_tr <- gc$sig_df$Tracer[gc$sig_df$Z_Median < -0.5]

    enr_txt <- if (length(enriched_tr) > 0) paste(enriched_tr, collapse = ", ") else "None (near catchment average)"
    dep_txt <- if (length(depleted_tr) > 0) paste(depleted_tr, collapse = ", ") else "None (near catchment average)"

    banner_style <- if (gc$global_med >= 0.15 || gc$max_sample_med >= 0.25) {
      "background-color: #fdecea; border-left: 6px solid #e74c3c; padding: 15px; border-radius: 4px;"
    } else if (gc$global_med >= 0.05 || gc$max_sample_med >= 0.10) {
      "background-color: #fef9e7; border-left: 6px solid #f39c12; padding: 15px; border-radius: 4px;"
    } else {
      "background-color: #eafaf1; border-left: 6px solid #2ecc71; padding: 15px; border-radius: 4px;"
    }

    verdict_title <- if (gc$global_med >= 0.15 || gc$max_sample_med >= 0.25) {
      "WARNING: Major Unsampled 'Ghost' Erosion Hotspot Detected!"
    } else if (gc$global_med >= 0.05 || gc$max_sample_med >= 0.10) {
      "NOTICE: Moderate / Episodic Unsampled Source Signal Detected"
    } else {
      "PASS: Known Sources Sufficiently Explain the Mixtures (No Ghost Hotspot)"
    }

    tagList(
      br(),
      tags$div(
        style = banner_style,
        h4(style = "margin-top:0; font-weight:bold;", verdict_title),
        p(
          HTML(paste0(
            "<b>Global Ghost Contribution (Median):</b> ", med_pct, "% ",
            "(95% CI: ", round(gc$global_q025 * 100, 1), "% – ", round(gc$global_q975 * 100, 1), "%) | ",
            "<b>Posterior Prob(Ghost > 5%):</b> ", p5_pct, "% | ",
            "<b>Prob(Ghost > 10%):</b> ", p10_pct, "%<br/>",
            "<b>Peak Sample Contribution:</b> ", max_pct, "% (in sample <code>", gc$max_sample_id, "</code>)<br/>",
            "<b>Reconstructed Ghost Tracer Anomaly:</b> Enriched in [<b>", enr_txt,
            "</b>]; Depleted in [<b>", dep_txt, "</b>]."
          ))
        )
      ),
      hr(),
      fluidRow(
        column(
          6,
          h4("1. Reconstructed Ghost Source Fingerprint (Z-Score)"),
          helpText("Values > 0 indicate tracers enriched in the missing source relative to sampled sources; values < 0 indicate depletion."),
          plotOutput("ghost_profile_plot", height = "420px")
        ),
        column(
          6,
          h4("2. Ghost Source Contribution Across Mixture Samples"),
          helpText("Identifies specific flood events, traps, or core horizons affected by unsampled erosion hotspots."),
          plotOutput("ghost_sample_plot", height = "420px")
        )
      ),
      hr(),
      fluidRow(
        column(8, h4("3. Reconstructed Tracer Values (Original Units) of the Ghost Source vs. Known Sources")),
        column(4, div(style = "text-align:right;", downloadButton("download_ghost_csv", "Download Ghost Report (CSV)", class = "btn-info")))
      ),
      tableOutput("ghost_signature_table")
    )
  })

  output$ghost_profile_plot <- renderPlot({
    gc <- ghost_calc()
    req(gc)
    df_s <- gc$sig_df
    df_s$Tracer <- factor(df_s$Tracer, levels = df_s$Tracer[order(df_s$Z_Median)])
    df_s$Role <- ifelse(df_s$Z_Median >= 0, "Enriched (+)", "Depleted (-)")

    ggplot(df_s, aes(x = Tracer, y = Z_Median, color = Role)) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
      geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray60") +
      geom_linerange(aes(ymin = Z_q025, ymax = Z_q975), linewidth = 0.8, alpha = 0.5) +
      geom_linerange(aes(ymin = Z_q25, ymax = Z_q75), linewidth = 2.2, alpha = 0.9) +
      geom_point(size = 4) +
      coord_flip() +
      scale_color_manual(values = c("Enriched (+)" = "#d35400", "Depleted (-)" = "#2980b9")) +
      theme_bw(base_size = 15) +
      theme(legend.position = "bottom") +
      labs(
        x = "Tracer",
        y = "Standardized Z-Score (relative to known sources)",
        color = "Ghost Anomaly"
      )
  })

  output$ghost_sample_plot <- renderPlot({
    gc <- ghost_calc()
    req(gc)
    df_m <- gc$sample_df

    ggplot(df_m, aes(x = Sample_ID, y = Median)) +
      geom_hline(yintercept = 0.05, linetype = "dotted", color = "#f39c12", linewidth = 0.9) +
      geom_hline(yintercept = 0.15, linetype = "dashed", color = "#e74c3c", linewidth = 0.9) +
      geom_linerange(aes(ymin = q025, ymax = q975), color = "#7f8c8d", linewidth = 0.7) +
      geom_linerange(aes(ymin = q25, ymax = q75), color = "#2c3e50", linewidth = 2.0) +
      geom_point(aes(color = Median > 0.10), size = 3.5) +
      scale_color_manual(values = c("FALSE" = "#2c3e50", "TRUE" = "#e74c3c"), guide = "none") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      theme_bw(base_size = 15) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
      labs(
        x = "Mixture Sample",
        y = "Ghost Source Proportion (Median & 50%/95% CI)"
      )
  })

  output$ghost_signature_table <- renderTable({
    gc <- ghost_calc()
    req(gc)
    gc$sig_df
  }, digits = 2)

  output$download_ghost_csv <- downloadHandler(
    filename = function() paste0("Ghost_Source_Diagnostics_", Sys.Date(), ".csv"),
    content = function(file) {
      gc <- ghost_calc()
      req(gc)
      write.csv(gc$sig_df, file, row.names = FALSE)
    }
  )


  #==========================================================================
  # Validation: virtual mixtures (Supports Blind Leave-One-Source-Out Test)
  #==========================================================================

  observeEvent(input$run_val, {
    req(values$df, input$final_selected_tracers, values$groups)

    if (length(input$final_selected_tracers) < 2) {
      showNotification("Error: Select at least 2 tracers for validation.", type = "error")
      return()
    }

    withProgress(message = "Validation (Virtual Mixtures)...", value = 0, {

      tr <- input$final_selected_tracers
      gr_col <- input$group_col
      m_lab <- input$mix_label

      all_groups <- values$groups
      df_src_all <- values$df[as.character(values$df[[gr_col]]) != m_lab, ]
      df_mix_real <- values$df[as.character(values$df[[gr_col]]) == m_lab, ]

      if (nrow(df_src_all) == 0) {
        showNotification("Error: Source group is empty.", type = "error")
        return(NULL)
      }

      # 1. Определяем выбранный режим валидации при активном Ghost Source Detector:
      #    - "include_ghost" : 1) Не скрывать ничего (G известных + Ghost_Source с долями от 0 до 1)
      #    - "hide_ghost"    : 2) Скрыть призрак (False-Alarm Test: смеси только из G известных, True Ghost = 0)
      #    - <имя источника> : 3+) Скрыть реальный источник (Blind Test: его истинная доля идет в Ghost)
      use_ghost_val <- isTRUE(input$opt_ghost)
      val_mode_sel  <- if (use_ghost_val && !is.null(input$val_hide_source)) input$val_hide_source else "hide_ghost"

      is_blind_hide  <- use_ghost_val && (val_mode_sel %in% all_groups) && (length(all_groups) >= 3)
      is_synth_ghost <- use_ghost_val && (val_mode_sel == "include_ghost")
      is_false_alarm <- use_ghost_val && (val_mode_sel == "hide_ghost")

      # 2. Параметры всех реальных источников
      src_params <- list()
      cov_accum <- matrix(0, nrow = length(tr), ncol = length(tr))

      for (g in all_groups) {
        sub_data <- as.matrix(df_src_all[as.character(df_src_all[[gr_col]]) == g, tr, drop = FALSE])
        sub_data[is.na(sub_data)] <- 0
        mu <- colMeans(sub_data)
        if (nrow(sub_data) < 2) {
          cov_mat <- diag(1e-6, length(tr))
        } else {
          cov_mat <- cov(sub_data)
          cov_mat[is.na(cov_mat)] <- 0
          diag(cov_mat) <- pmax(diag(cov_mat), 0) + 1e-6
        }
        src_params[[g]] <- list(mu = mu, sigma = cov_mat)
        cov_accum <- cov_accum + (cov_mat / length(all_groups))
      }

      # Для Режима 1 ("Не скрывать ничего"): добавляем внешний Ghost_Source в пул генерации ВМ
      if (is_synth_ghost) {
        y_s_tmp  <- as.matrix(df_src_all[, tr, drop = FALSE])
        y_m_tmp  <- as.matrix(df_mix_real[, tr, drop = FALSE])
        sc_tmp   <- safe_scale_matrix(y_s_tmp)
        y_m_sc   <- scale(y_m_tmp, center = sc_tmp$center, scale = sc_tmp$scale)
        g_id_tmp <- make_group_id(as.character(df_src_all[[gr_col]]), all_groups)

        # Вычисляем внешний якорь призрака (с отталкиванием >= 2.5 сигмы наружу от полигона)
        anch_z <- calc_ghost_anchor(sc_tmp$y, g_id_tmp, length(all_groups), y_m_sc)
        ghost_mu_phys <- sc_tmp$center + sc_tmp$scale * anch_z

        for (j in seq_along(tr)) {
          if (min(y_s_tmp[, j], na.rm = TRUE) >= 0 && ghost_mu_phys[j] < 1e-6) {
            ghost_mu_phys[j] <- 1e-6
          }
        }
        names(ghost_mu_phys) <- tr

        src_params[["Ghost_Source"]] <- list(mu = ghost_mu_phys, sigma = cov_accum)
        gen_groups <- c(all_groups, "Ghost_Source")
      } else {
        gen_groups <- all_groups
      }

      n_gen_groups <- length(gen_groups)

      # 3. Генерация истинных пропорций виртуальных смесей
      is_simplex <- input$val_gen_mode %in% c("det_simp", "stoch_simp")
      is_stoch_tracers <- input$val_gen_mode %in% c("stoch_rand", "stoch_simp")

      if (n_gen_groups == 1) {
        N_val <- 1
        true_props <- matrix(1, nrow = 1, ncol = 1)
        colnames(true_props) <- gen_groups
      } else if (is_simplex) {
        step <- if (is.null(input$val_step) || input$val_step <= 0) 0.05 else input$val_step
        seq_vals <- seq(0, 1, by = step)
        grid <- expand.grid(rep(list(seq_vals), n_gen_groups - 1))
        grid$last <- 1 - rowSums(grid)
        grid <- grid[grid$last >= -1e-9 & grid$last <= 1 + 1e-9, , drop = FALSE]

        if (nrow(grid) == 0) {
          true_props <- matrix(1 / n_gen_groups, nrow = 1, ncol = n_gen_groups)
          colnames(true_props) <- gen_groups
          N_val <- 1
        } else {
          grid[grid < 0] <- 0
          true_props <- as.matrix(t(apply(grid, 1, function(x) x / sum(x))))
          colnames(true_props) <- gen_groups
          N_val <- nrow(true_props)
          if (N_val > 300) {
            set.seed(123)
            true_props <- true_props[sample(seq_len(N_val), 300), , drop = FALSE]
            N_val <- 300
          }
        }
      } else {
        N_val <- input$val_n
        true_props <- matrix(0, nrow = N_val, ncol = n_gen_groups)
        colnames(true_props) <- gen_groups
        for (i in seq_len(N_val)) {
          p_rand <- rexp(n_gen_groups, rate = 1)
          true_props[i, ] <- p_rand / sum(p_rand)
        }
      }

      # 4. Генерация концентраций трассеров в виртуальных смесях
      y_val_raw <- matrix(0, nrow = N_val, ncol = length(tr))
      colnames(y_val_raw) <- tr
      noise <- as.numeric(input$val_noise)

      for (i in seq_len(N_val)) {
        mix_vec <- numeric(length(tr))
        p_current <- true_props[i, ]

        for (k in seq_len(n_gen_groups)) {
          g_name <- gen_groups[k]
          if (is_stoch_tracers) {
            fallback_draw <- function() {
              s2 <- diag(src_params[[g_name]]$sigma)
              s2[is.na(s2) | s2 <= 0] <- 1e-6
              rnorm(length(tr), src_params[[g_name]]$mu, sqrt(s2))
            }
            sim <- tryCatch(
              as.numeric(MASS::mvrnorm(1, src_params[[g_name]]$mu, src_params[[g_name]]$sigma)),
              error = function(e) fallback_draw()
            )
            if (length(sim) != length(tr) || any(!is.finite(sim))) sim <- src_params[[g_name]]$mu
          } else {
            sim <- src_params[[g_name]]$mu
          }
          mix_vec <- mix_vec + (sim * p_current[k])
        }

        if (!is.na(noise) && noise > 0) {
          mix_vec <- mix_vec * runif(length(mix_vec), 1 - noise, 1 + noise)
        }
        mix_vec[mix_vec <= 1e-6] <- 1e-6
        y_val_raw[i, ] <- mix_vec
      }

      # Микро-стабилизация при 0% шума против вырождения матрицы Холецкого
      for (j in seq_len(ncol(y_val_raw))) {
        src_sd_j <- sd(df_src_all[[tr[j]]], na.rm = TRUE)
        if (is.na(src_sd_j) || src_sd_j < 1e-6) src_sd_j <- 1.0
        y_val_raw[, j] <- y_val_raw[, j] + rnorm(N_val, 0, 0.002 * src_sd_j)
      }

      # 5. Формирование обучающей выборки источников и эталонных долей для каждого из 3 режимов
      if (is_blind_hide) {
        # Режим 3+: Скрыт один реальный источник
        train_groups <- setdiff(all_groups, val_mode_sel)
        df_src <- df_src_all[as.character(df_src_all[[gr_col]]) != val_mode_sel, ]
        true_props_eval <- cbind(true_props[, train_groups, drop = FALSE], Ghost_Source = true_props[, val_mode_sel])
        colnames(true_props_eval) <- c(train_groups, paste0("Ghost (", val_mode_sel, ")"))
      } else if (is_false_alarm) {
        # Режим 2: Скрыт призрак (все реальные источники в обучении, истинная доля призрака = 0)
        train_groups <- all_groups
        df_src <- df_src_all
        true_props_eval <- cbind(true_props, Ghost_Source = rep(0, N_val))
      } else {
        # Режим 1 (Не скрывать ничего, G+1) или обычный режим без Ghost
        train_groups <- all_groups
        df_src <- df_src_all
        true_props_eval <- true_props
      }

      n_train_groups  <- length(train_groups)
      n_out_groups    <- ncol(true_props_eval)
      out_group_names <- colnames(true_props_eval)

      export_df <- data.frame(ID = paste0("Virtual_Mix_", seq_len(N_val)))
      export_df <- cbind(export_df, as.data.frame(y_val_raw))
      props_df <- as.data.frame(true_props_eval)
      colnames(props_df) <- paste0("TRUE_", colnames(props_df))
      values$val_export <- cbind(export_df, props_df)

      # 6. Масштабирование или CoDA относительно обучающих источников
      y_src_real <- as.matrix(df_src[, tr, drop = FALSE])
      y_src_real[is.na(y_src_real)] <- 0
      y_src_real[y_src_real <= 0] <- 1e-6

      if (!input$opt_coda) {
        sc_src <- safe_scale_matrix(y_src_real)
        y_src_s <- sc_src$y
        y_val_s <- scale(y_val_raw, center = sc_src$center, scale = sc_src$scale)
        y_ctr_val <- sc_src$center
        y_scl_val <- sc_src$scale
        V_val <- matrix(0, nrow = 0, ncol = 0)
        N_tr_raw_val <- ncol(y_val_s)
      } else {
        ys_prop <- y_src_real / rowSums(y_src_real)
        yv_prop <- y_val_raw / rowSums(y_val_raw)
        N_tr_raw_val <- ncol(ys_prop)
        V_val <- compositions::ilrBase(D = N_tr_raw_val)
        y_src_s <- as.matrix(log(ys_prop) %*% V_val)
        y_val_s <- as.matrix(log(yv_prop) %*% V_val)
        y_ctr_val <- rep(0, ncol(y_val_s))
        y_scl_val <- rep(1, ncol(y_val_s))
      }

      ssa_val <- rep(0, N_val)
      ssa_src <- rep(0, nrow(df_src))
      ssa_grp_val <- rep(0, n_train_groups)

      conn_weights_val <- vapply(
        train_groups,
        function(g) {
          x <- input[[paste0("w_", g)]]
          if (is.null(x) || is.na(x) || x <= 0) 1.0 else as.numeric(x)
        },
        numeric(1)
      )

      size_active_val <- !is.null(input$size_col) && input$size_col != "None" && input$size_col %in% colnames(df_src)
      if (size_active_val) {
        ssa_src_raw <- df_src[[input$size_col]]
        ssa_src_raw[is.na(ssa_src_raw)] <- mean(ssa_src_raw, na.rm = TRUE)
        mean_s <- mean(ssa_src_raw, na.rm = TRUE)
        sd_s   <- sd(ssa_src_raw, na.rm = TRUE)
        if (!is.na(sd_s) && sd_s >= 1e-8) {
          ssa_src <- (ssa_src_raw - mean_s) / sd_s
          src_labels <- as.character(df_src[[gr_col]])
          ssa_grp_val <- sapply(train_groups, function(g) {
            x <- ssa_src[src_labels == g]
            if (length(x) == 0 || all(is.na(x))) 0 else mean(x, na.rm = TRUE)
          })
        }
      }

      use_levels_val <- !is.null(input$mix_cov) && input$mix_cov != "None" && input$pool_mode != "no_pool"
      group_id <- make_group_id(as.character(df_src[[gr_col]]), train_groups)

      use_b_src_val <- isTRUE(input$opt_b && input$opt_b_src && !input$opt_coda)
      ssa_grp_for_stan <- if (use_b_src_val) ssa_grp_val else rep(0, n_train_groups)

      stan_val_data <- list(
        N_mix = N_val,
        N_sources = nrow(y_src_s),
        N_tracers = ncol(y_src_s),
        N_groups = n_train_groups,
        y_mix = y_val_s,
        y_source = y_src_s,
        group_id = group_id,
        conn_weights = conn_weights_val,
        ssa_mix = as.array(ssa_val),
        ssa_source = as.array(ssa_src),
        ssa_group = as.array(ssa_grp_for_stan),
        y_center = as.array(as.numeric(y_ctr_val)),
        y_scale = as.array(as.numeric(y_scl_val)),
        holdout_id = 0L
      )

      if (use_ghost_val) {
        gh_anchor_val <- calc_ghost_anchor(y_src_s, group_id, n_train_groups, y_val_s)
        stan_val_data$ghost_anchor <- as.array(as.numeric(gh_anchor_val))
        stan_val_data$alpha_ghost  <- if (!is.null(input$ghost_alpha)) as.numeric(input$ghost_alpha) else 0.12
      }

      if (use_levels_val) {
        stan_val_data$N_mix_levels <- 1
        stan_val_data$mix_level_id <- as.array(rep(1, N_val))
      }

      if (input$opt_coda) {
        stan_val_data$N_tracers_raw <- N_tr_raw_val
        stan_val_data$V <- V_val
        stan_val_data$src_totals <- as.array(
          sapply(train_groups, function(g) {
            idx <- as.character(df_src[[gr_col]]) == g
            mean(rowSums(y_src_real[idx, , drop = FALSE], na.rm = TRUE))
          })
        )
      }

      h_mode_val <- if (!is.null(input$hotspot_mode)) input$hotspot_mode else "mean"
      if (h_mode_val == "k_scale") {
        stan_val_data$k_hotspot <- if (!is.null(input$k_hotspot) && !is.na(input$k_hotspot)) max(0.1, as.numeric(input$k_hotspot)) else 3.0
      }
      if (h_mode_val == "nested") {
        stan_val_data$sigma_hotspot <- if (!is.null(input$sigma_hotspot) && !is.na(input$sigma_hotspot)) max(0.1, as.numeric(input$sigma_hotspot)) else 1.5
      }

      model_id <- paste0(
        "VAL_",
        input$pool_mode,
        as.integer(input$opt_b),
        as.integer(use_b_src_val),
        input$cov_mode,
        as.integer(input$opt_coda),
        as.integer(use_levels_val),
        as.integer(input$opt_bias),
        as.integer(use_ghost_val),
        "_", h_mode_val
      )

      if (!is.null(values$compiled_list[[model_id]])) {
        mod <- values$compiled_list[[model_id]]
      } else {
        showNotification("Compiling model for validation (first time)...", type = "message")
        code <- generate_stan_code(
          pool_mode = input$pool_mode,
          b = input$opt_b,
          b_src = use_b_src_val,
          cov_mode = input$cov_mode,
          coda = input$opt_coda,
          has_groups = use_levels_val,
          bias = input$opt_bias,
          ghost = use_ghost_val,
          hotspot_mode = h_mode_val
        )
        values$compiled_list[[model_id]] <- stan_model(model_code = code)
        mod <- values$compiled_list[[model_id]]
      }

      fit_val <- tryCatch({
        if (input$calc_method == "vb") {
          vb(mod, data = stan_val_data, output_samples = 1000, seed = 123, init = "random", tol_rel_obj = 0.005)
        } else {
          sampling(mod, data = stan_val_data, iter = input$iter, chains = input$chains, init = 0.5)
        }
      }, error = function(e) {
        showNotification(paste("Validation error:", e$message), type = "error")
        return(NULL)
      })

      if (is.null(fit_val)) return(NULL)
      val_res <- rstan::extract(fit_val)

      if ("P_ind" %in% names(val_res)) {
        raw <- ensure_3d_samples(val_res$P_ind, n_out_groups)
        p_samples_array <- raw
        pred_means <- apply(raw, c(2, 3), median)
        pred_q025  <- apply(raw, c(2, 3), quantile, probs = 0.025)
        pred_q975  <- apply(raw, c(2, 3), quantile, probs = 0.975)
        pred_q25   <- apply(raw, c(2, 3), quantile, probs = 0.25)
        pred_q75   <- apply(raw, c(2, 3), quantile, probs = 0.75)
      } else if ("P_level" %in% names(val_res)) {
        raw <- val_res$P_level
        p_samples <- if (length(dim(raw)) == 3) matrix(raw[, 1, , drop = FALSE], nrow = dim(raw)[1], ncol = n_out_groups) else matrix(raw, ncol = n_out_groups)
        p_samples_array <- array(0, dim = c(nrow(p_samples), N_val, n_out_groups))
        for (v in seq_len(N_val)) p_samples_array[, v, ] <- p_samples
        pred_means <- matrix(rep(apply(p_samples, 2, median), N_val), nrow = N_val, byrow = TRUE)
        pred_q025  <- matrix(rep(apply(p_samples, 2, quantile, probs = 0.025), N_val), nrow = N_val, byrow = TRUE)
        pred_q975  <- matrix(rep(apply(p_samples, 2, quantile, probs = 0.975), N_val), nrow = N_val, byrow = TRUE)
        pred_q25   <- matrix(rep(apply(p_samples, 2, quantile, probs = 0.25), N_val), nrow = N_val, byrow = TRUE)
        pred_q75   <- matrix(rep(apply(p_samples, 2, quantile, probs = 0.75), N_val), nrow = N_val, byrow = TRUE)
      } else {
        raw <- matrix(val_res$P, ncol = n_out_groups)
        p_samples_array <- array(0, dim = c(nrow(raw), N_val, n_out_groups))
        for (v in seq_len(N_val)) p_samples_array[, v, ] <- raw
        pred_means <- matrix(rep(apply(raw, 2, median), N_val), nrow = N_val, byrow = TRUE)
        pred_q025  <- matrix(rep(apply(raw, 2, quantile, probs = 0.025), N_val), nrow = N_val, byrow = TRUE)
        pred_q975  <- matrix(rep(apply(raw, 2, quantile, probs = 0.975), N_val), nrow = N_val, byrow = TRUE)
        pred_q25   <- matrix(rep(apply(raw, 2, quantile, probs = 0.25), N_val), nrow = N_val, byrow = TRUE)
        pred_q75   <- matrix(rep(apply(raw, 2, quantile, probs = 0.75), N_val), nrow = N_val, byrow = TRUE)
      }

      pred_means <- matrix(pred_means, nrow = N_val, ncol = n_out_groups)
      colnames(pred_means) <- out_group_names

      values$val_data <- list(
        true = true_props_eval,
        pred = pred_means,
        q025 = matrix(pred_q025, nrow = N_val, ncol = n_out_groups),
        q975 = matrix(pred_q975, nrow = N_val, ncol = n_out_groups),
        q25  = matrix(pred_q25,  nrow = N_val, ncol = n_out_groups),
        q75  = matrix(pred_q75,  nrow = N_val, ncol = n_out_groups),
        samples = p_samples_array
      )

      pred_df_export <- as.data.frame(pred_means)
      colnames(pred_df_export) <- paste0("PRED_", out_group_names)
      values$val_export_full <- cbind(values$val_export, pred_df_export)

      showNotification("Validation complete! Metrics updated.", type = "message")
    })
  })

  output$download_val_data <- downloadHandler(
    filename = function() paste0("Validation_Results_", Sys.Date(), ".csv"),
    content = function(file) {
      req(values$val_export)
      df_out <- if (!is.null(values$val_export_full)) values$val_export_full else values$val_export
      write.csv(df_out, file, row.names = FALSE)
    }
  )

  output$val_plot <- renderPlot({
    req(values$val_data)
    true_df <- reshape2::melt(values$val_data$true)
    pred_df <- reshape2::melt(values$val_data$pred)
    plot_df <- data.frame(Source = true_df$Var2, True = true_df$value, Predicted = pred_df$value)

    ggplot(plot_df, aes(x = True, y = Predicted, color = Source)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey", linewidth = 1) +
      geom_point(size = 3, alpha = 0.7) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5, linetype = "dotted") +
      coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
      theme_bw(base_size = 18) +
      labs(title = "Model Accuracy (Virtual Mixtures)", x = "True Proportion (Known)", y = "Predicted Proportion (Model)")
  })

  output$val_metrics <- renderTable({
    req(values$val_data)
    t <- as.matrix(values$val_data$true)
    p <- as.matrix(values$val_data$pred)
    if (ncol(t) != ncol(p)) return(NULL)

    q025 <- as.matrix(values$val_data$q025); q975 <- as.matrix(values$val_data$q975)
    q25  <- as.matrix(values$val_data$q25);  q75  <- as.matrix(values$val_data$q75)

    mae  <- colMeans(abs(t - p), na.rm = TRUE)
    rmse <- sqrt(colMeans((t - p)^2, na.rm = TRUE))
    me   <- colMeans(p - t, na.rm = TRUE)

    naive_pred <- 1 / ncol(t)
    denom_nse <- colSums((t - naive_pred)^2, na.rm = TRUE)
    nse <- ifelse(denom_nse == 0, NA_real_, 1 - (colSums((t - p)^2, na.rm = TRUE) / denom_nse))

    w95 <- colMeans(q975 - q025, na.rm = TRUE)
    w50 <- colMeans(q75 - q25, na.rm = TRUE)
    eps <- 0.015
    p95_adj <- colMeans((t >= (q025 - eps)) & (t <= (q975 + eps)), na.rm = TRUE)
    p50_adj <- colMeans((t >= (q25 - eps)) & (t <= (q75 + eps)), na.rm = TRUE)

    true_dom <- t > 0.50; pred_dom <- p > 0.50
    hits <- colSums(true_dom & pred_dom, na.rm = TRUE)
    misses <- colSums(true_dom & !pred_dom, na.rm = TRUE)
    false_alarms <- colSums(!true_dom & pred_dom, na.rm = TRUE)

    denom_csi <- hits + misses + false_alarms
    csi <- ifelse(denom_csi == 0, NA_real_, hits / denom_csi)
    denom_hr <- hits + misses
    hr  <- ifelse(denom_hr == 0, NA_real_, hits / denom_hr)

    r2 <- sapply(seq_len(ncol(t)), function(i) {
      cc <- suppressWarnings(cor(t[, i], p[, i], use = "complete.obs"))
      if (is.na(cc)) NA_real_ else cc^2
    })

    crps_vals <- rep(NA_real_, ncol(t))
    if (!is.null(values$val_data$samples)) {
      crps_vals <- sapply(seq_len(ncol(t)), function(i) {
        tryCatch({
          src_samples <- values$val_data$samples[, , i]
          if (is.null(dim(src_samples))) src_samples <- matrix(src_samples, ncol = 1)
          mean(scoringRules::crps_sample(y = t[, i], dat = t(src_samples)), na.rm = TRUE)
        }, error = function(e) NA_real_)
      })
    }

    df_metrics <- data.frame(
      Source = colnames(t), MAE = mae, RMSE = rmse, ME = me, R2 = r2, CRPS = crps_vals, NSE = nse,
      `W50` = w50, `W95` = w95, `P50 (adj)` = p50_adj, `P95 (adj)` = p95_adj,
      `CSI (>50%)` = csi, `HR (>50%)` = hr, check.names = FALSE
    )

    mean_row <- data.frame(
      Source = "AVERAGE", MAE = mean(mae, na.rm = TRUE), RMSE = mean(rmse, na.rm = TRUE),
      ME = mean(abs(me), na.rm = TRUE), R2 = mean(r2, na.rm = TRUE), CRPS = mean(crps_vals, na.rm = TRUE),
      NSE = mean(nse, na.rm = TRUE), `W50` = mean(w50, na.rm = TRUE), `W95` = mean(w95, na.rm = TRUE),
      `P50 (adj)` = mean(p50_adj, na.rm = TRUE), `P95 (adj)` = mean(p95_adj, na.rm = TRUE),
      `CSI (>50%)` = mean(csi, na.rm = TRUE), `HR (>50%)` = mean(hr, na.rm = TRUE), check.names = FALSE
    )

    rbind(df_metrics, mean_row)
  }, digits = 3, na = "-")


  #==========================================================================
  # LOO-CV MODEL COMPARISON & DIAGNOSTICS
  #==========================================================================

  output$ui_multi_compare_selector <- renderUI({
    req(length(values$saved_models) > 0)
    choices <- names(values$saved_models)
    selectizeInput(
      "selected_models_loo", "Select Models to Compare:",
      choices = choices, selected = choices, multiple = TRUE,
      options = list(plugins = list("remove_button"))
    )
  })

  observeEvent(input$run_loo_compare, {
    req(length(input$selected_models_loo) > 0)
    withProgress(message = "Calculating LOO-CV for selected models...", value = 0, {
      n_models <- length(input$selected_models_loo)
      for (i in seq_along(input$selected_models_loo)) {
        m_name <- input$selected_models_loo[i]
        incProgress(1 / n_models, detail = paste("Evaluating", m_name))
        tryCatch({
          if (is.null(values$loo_objects[[m_name]])) {
            fit <- values$saved_models[[m_name]]
            is_vb_fit <- grepl("^\\[VB\\]", m_name)
            log_lik <- loo::extract_log_lik(fit, parameter_name = "log_lik", merge_chains = is_vb_fit)
            values$loo_objects[[m_name]] <- loo::loo(log_lik, cores = 1)
          }
        }, error = function(e) {
          showNotification(paste("Error calculating LOO for", m_name, ":", e$message), type = "error")
        })
      }
    })
  })

  output$loo_results_text <- renderPrint({
    loo_list <- values$loo_objects[names(values$loo_objects) %in% input$selected_models_loo]
    if (length(loo_list) > 0) {
      cat("=== LOO-CV Leaderboard (Approximate PSIS-LOO) ===\n\n")
      if (length(loo_list) == 1) {
        print(loo_list[[1]])
      } else {
        n_obs <- nrow(loo_list[[1]]$pointwise)
        if (n_obs == 1) {
          cat("WARNING: The dataset contains only one mixture (N_mix = 1).\n")
          scores <- sort(sapply(loo_list, function(x) x$estimates["elpd_loo", "Estimate"]), decreasing = TRUE)
          for (m in names(scores)) cat(sprintf("%s : %.2f\n", m, scores[m]))
        } else {
          comp <- loo::loo_compare(loo_list)
          cat("Negative values in 'elpd_diff' indicate worse predictive performance.\n\n")
          print(comp)
        }
      }
    } else {
      cat("Select models and click 'Run LOO-CV' to see the PSIS-LOO leaderboard.\n")
    }

    if (!is.null(values$exact_loo_results) && length(values$exact_loo_results) > 0) {
      cat("\n=======================================================\n")
      cat("=== EXACT LOO-CV LEADERBOARD                        ===\n")
      cat("=======================================================\n")
      for (m_name in names(values$exact_loo_results)) {
        ex_res <- values$exact_loo_results[[m_name]]
        cat(sprintf("Model: %s\n  elpd_loo:    %8.2f\n  se_elpd_loo: %8.2f\n-------------------------------------------------------\n",
                    m_name, ex_res$elpd_loo, ex_res$se_elpd_loo))
      }
    }
    cat("\nImportant: For no_pool, Exact LOO estimates prior-predictive performance\n")
    cat("for the left-out independent sample, because each mixture has its own P_ind.\n")
  })

  output$ui_pareto_selector <- renderUI({
    req(length(values$loo_objects) > 0)
    selectInput("pareto_model_select", "Select Model to inspect:", choices = names(values$loo_objects))
  })

  output$pareto_plot <- renderPlot({
    req(input$pareto_model_select, values$loo_objects[[input$pareto_model_select]])
    loo_obj <- values$loo_objects[[input$pareto_model_select]]
    k_values <- loo_obj$diagnostics$pareto_k
    if (is.null(k_values) || length(k_values) == 0) return(NULL)

    plot_df <- data.frame(Observation = seq_along(k_values), k = as.numeric(k_values))
    plot_df$k[is.na(plot_df$k)] <- 0
    plot_df$is_high <- plot_df$k > 0.7
    y_min <- min(0, plot_df$k, na.rm = TRUE)
    y_max <- max(c(1, plot_df$k), na.rm = TRUE) + 0.1

    ggplot(plot_df, aes(x = Observation, y = k)) +
      geom_hline(yintercept = 0.7, linetype = "dashed", color = "#e74c3c", linewidth = 1) +
      geom_hline(yintercept = 1.0, linetype = "solid", color = "#c0392b", linewidth = 0.5, alpha = 0.5) +
      geom_point(aes(color = is_high, size = is_high), alpha = 0.8) +
      scale_color_manual(values = c("FALSE" = "#2c3e50", "TRUE" = "#e74c3c")) +
      scale_size_manual(values = c("FALSE" = 3, "TRUE" = 5)) +
      geom_text(data = subset(plot_df, k > 0.7), aes(label = Observation), vjust = -1.5, color = "#c0392b", fontface = "bold", size = 5) +
      theme_bw(base_size = 16) +
      theme(legend.position = "none", panel.grid.minor = element_blank(), axis.line = element_line(linewidth = 1)) +
      labs(
        title = paste("Pareto k Diagnostics:", input$pareto_model_select),
        subtitle = "Red points (k > 0.7) indicate highly influential mixtures",
        x = "Mixture Index (Data point)", y = "Pareto shape k"
      ) +
      coord_cartesian(ylim = c(y_min, y_max))
  })

  observeEvent(input$run_exact_loo, {
    req(values$stan_data, values$compiled_model, values$last_model_name)
    n_mixtures <- values$stan_data$N_mix
    exact_elpd_i <- numeric(n_mixtures)

    withProgress(message = "Running Exact LOO-CV (slow)...", value = 0, {
      for (i in seq_len(n_mixtures)) {
        incProgress(0, detail = paste("Refitting target", i, "of", n_mixtures))
        cv_data <- values$stan_data
        cv_data$holdout_id <- i

        exact_elpd_i[i] <- tryCatch({
          fit_cv <- rstan::sampling(
            values$compiled_model, data = cv_data,
            iter = input$iter, chains = input$chains, cores = input$chains,
            init = 0.5,
            control = list(adapt_delta = 0.98, max_treedepth = 14),
            refresh = 0, show_messages = FALSE
          )
          all_log_liks <- rstan::extract(fit_cv, pars = "log_lik")$log_lik
          log_mean_exp(all_log_liks[, i])
        }, error = function(e) {
          warning(paste("Exact LOO failed for mixture", i, ":", e$message))
          NA_real_
        })
        incProgress(1 / n_mixtures)
      }
    })

    if (all(is.na(exact_elpd_i))) {
      showNotification("Exact LOO-CV failed for all folds.", type = "error")
      return()
    }

    total_exact_elpd <- sum(exact_elpd_i, na.rm = TRUE)
    se_exact_elpd <- if (n_mixtures > 1 && sum(!is.na(exact_elpd_i)) > 1) sqrt(n_mixtures * var(exact_elpd_i, na.rm = TRUE)) else NA_real_

    values$exact_loo_results[[values$last_model_name]] <- list(
      elpd_loo = total_exact_elpd,
      se_elpd_loo = se_exact_elpd,
      pointwise = exact_elpd_i
    )
    showNotification(sprintf("Exact LOO-CV complete. ELPD: %.2f", total_exact_elpd), type = "message")
  })

}

#==============================================================================
# 5. RUN APP
#==============================================================================

shinyApp(ui, server)