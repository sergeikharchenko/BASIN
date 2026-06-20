# ==============================================================================
# BASIN: Bayesian Sediment Source Appointment Tool
# ==============================================================================

# 1. Installation and loading of libraries
suppressPackageStartupMessages({
  suppressWarnings({
    if (!require("pacman", quietly = TRUE)) install.packages("pacman", quiet = TRUE)
    options(repos = c(CRAN = "https://cloud.r-project.org"))
    
    # Load packages silently (added 'zip' for cross-platform archive creation)
    pacman::p_load(shiny, rstan, tidyverse, caret, reshape2, dplyr, ggplot2, DT, shinythemes, data.table, MASS, zip, compositions, scoringRules)
  })
})

# Performance settings for Stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ==============================================================================
# 1. STAN CODE GENERATOR
# ==============================================================================
generate_stan_code <- function(h = FALSE, b = FALSE, cov_mode = "mix", coda = FALSE, has_groups = FALSE, bias = FALSE) {
  
  cov_mix <- cov_mode %in% c("mix", "full")
  cov_src <- cov_mode == "full"
  
  # --- 1. ДОБАВЛЯЕМ holdout_id В DATA ---
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
    
    int<lower=0, upper=N_mix> holdout_id;
\n"
  if (has_groups) data_blk <- paste0(data_blk, "\n    int<lower=1> N_mix_levels; \n    int<lower=1, upper=N_mix_levels> mix_level_id[N_mix];")
  
  if (coda) data_blk <- paste0(data_blk, "\n    int<lower=1> N_tracers_raw;\n    matrix[N_tracers_raw, N_tracers] V;")
  data_blk <- paste0(data_blk, "\n  }\n")
  
  param_blk <- "parameters {\n"
  if (h) {
    param_blk <- paste0(param_blk, "  simplex[N_groups] P_global;\n  simplex[N_groups] P_ind[N_mix];\n  real<lower=0> phi_ind;\n")
    if (has_groups) param_blk <- paste0(param_blk, "  simplex[N_groups] P_level[N_mix_levels];\n  real<lower=0> phi_level;\n")
  } else {
    if (has_groups) param_blk <- paste0(param_blk, "  simplex[N_groups] P_level[N_mix_levels];\n")
    else            param_blk <- paste0(param_blk, "  simplex[N_groups] P;\n")
  }
  
  param_blk <- paste0(param_blk, "  vector[N_tracers] mu_source[N_groups];\n")
  
  if (cov_src) {
    param_blk <- paste0(param_blk, "  cholesky_factor_corr[N_tracers] L_Omega_src[N_groups];\n  vector<lower=0>[N_tracers] L_sigma_src[N_groups];\n")
  } else {
    param_blk <- paste0(param_blk, "  vector<lower=0>[N_tracers] sigma_source[N_groups];\n")
  }
  
  if (b) {
    if (coda) param_blk <- paste0(param_blk, "  vector[N_tracers_raw] beta_sorting;\n")
    else      param_blk <- paste0(param_blk, "  vector[N_tracers] beta_sorting;\n")
  }
  
  if (cov_mix) param_blk <- paste0(param_blk, "  cholesky_factor_corr[N_tracers] L_Omega; vector<lower=0>[N_tracers] L_sigma;\n")
  else         param_blk <- paste0(param_blk, "  vector<lower=0>[N_tracers] sigma_mix;\n")
  
  if (bias)    param_blk <- paste0(param_blk, "  vector[N_tracers] delta_bias;\n")
  param_blk <- paste0(param_blk, "}\n")
  
  model_blk <- "model {\n"
  if (h) {
    model_blk <- paste0(model_blk, "  phi_ind ~ gamma(2, 0.1);\n  P_global ~ dirichlet(conn_weights);\n")
    if (has_groups) model_blk <- paste0(model_blk, "  phi_level ~ gamma(2, 0.1);\n  for(L in 1:N_mix_levels) P_level[L] ~ dirichlet(P_global * phi_level);\n")
  } else {
    if (has_groups) model_blk <- paste0(model_blk, "  for(L in 1:N_mix_levels) P_level[L] ~ dirichlet(conn_weights);\n")
    else            model_blk <- paste0(model_blk, "  P ~ dirichlet(conn_weights);\n")
  }
  
  if (cov_src) {
    model_blk <- paste0(model_blk, "  for (g in 1:N_groups) {\n    L_Omega_src[g] ~ lkj_corr_cholesky(2);\n    L_sigma_src[g] ~ cauchy(0, 1);\n    mu_source[g] ~ normal(0, 10);\n  }\n")
    model_blk <- paste0(model_blk, "  for (i in 1:N_sources) y_source[i]' ~ multi_normal_cholesky(mu_source[group_id[i]], diag_pre_multiply(L_sigma_src[group_id[i]], L_Omega_src[group_id[i]]));\n")
  } else {
    model_blk <- paste0(model_blk, "  for (g in 1:N_groups) mu_source[g] ~ normal(0, 10);\n")
    model_blk <- paste0(model_blk, "  for (i in 1:N_sources) y_source[i] ~ normal(mu_source[group_id[i]], sigma_source[group_id[i]]);\n")
  }
  
  if (b)       model_blk <- paste0(model_blk, "  beta_sorting ~ normal(0, 0.2);\n")
  if (cov_mix) model_blk <- paste0(model_blk, "  L_Omega ~ lkj_corr_cholesky(4); L_sigma ~ cauchy(0, 1);\n")
  if (bias)    model_blk <- paste0(model_blk, "  delta_bias ~ student_t(3, 0, 0.1);\n")
  
  model_blk <- paste0(model_blk, "\n  for (i in 1:N_mix) {\n")
  current_p_g <- if (h) { if (has_groups) "P_level[mix_level_id[i]][g]" else "P_ind[i][g]" } else { if (has_groups) "P_level[mix_level_id[i]][g]" else "P[g]" }
  if (h) {
    if (has_groups) model_blk <- paste0(model_blk, "    P_ind[i] ~ dirichlet(P_level[mix_level_id[i]] * phi_ind);\n")
    else            model_blk <- paste0(model_blk, "    P_ind[i] ~ dirichlet(P_global * phi_ind);\n")
    current_p_g <- "P_ind[i][g]"
  }
  
  model_blk <- paste0(model_blk, "    vector[N_tracers] mu_f;\n    vector[N_tracers] y_f = y_mix[i]';\n")
  
  if (coda) {
    model_blk <- paste0(model_blk, "    vector[N_tracers_raw] mu_pred_prop;\n    vector[N_tracers_raw] src_prop_matrix[N_groups];\n")
    model_blk <- paste0(model_blk, "    for (g in 1:N_groups) src_prop_matrix[g] = softmax(V * mu_source[g]);\n")
    model_blk <- paste0(model_blk, "    for (t in 1:N_tracers_raw) {\n      real val = 0;\n      for (g in 1:N_groups) {\n")
    if (b) model_blk <- paste0(model_blk, sprintf("        val += %s * src_prop_matrix[g, t] * exp(beta_sorting[t] * (ssa_mix[i] - mean(ssa_source)));\n", current_p_g))
    else   model_blk <- paste0(model_blk, sprintf("        val += %s * src_prop_matrix[g, t];\n", current_p_g))
    model_blk <- paste0(model_blk, "      }\n      mu_pred_prop[t] = val + 1e-9;\n    }\n")
    if (bias) model_blk <- paste0(model_blk, "    mu_f = (log(mu_pred_prop)' * V)' + delta_bias;\n")
    else      model_blk <- paste0(model_blk, "    mu_f = (log(mu_pred_prop)' * V)';\n")
  } else {
    model_blk <- paste0(model_blk, "    for (t in 1:N_tracers) {\n      real val = 0;\n")
    if (b) model_blk <- paste0(model_blk, sprintf("      for (g in 1:N_groups) val += %s * mu_source[g, t] * exp(beta_sorting[t] * (ssa_mix[i] - mean(ssa_source)));\n", current_p_g))
    else  model_blk <- paste0(model_blk, sprintf("      for (g in 1:N_groups) val += %s * mu_source[g, t];\n", current_p_g))
    if (bias) model_blk <- paste0(model_blk, "      mu_f[t] = val + delta_bias[t];\n    }\n")
    else      model_blk <- paste0(model_blk, "      mu_f[t] = val;\n    }\n")
  }
  
  # --- 2. ДОБАВЛЯЕМ МАСКИРОВАНИЕ ПРИ РАСЧЕТЕ ПРАВДОПОДОБИЯ ---
  if (cov_mix) model_blk <- paste0(model_blk, "    if (i != holdout_id) y_f ~ multi_student_t_cholesky(4, mu_f, diag_pre_multiply(L_sigma, L_Omega));\n")
  else         model_blk <- paste0(model_blk, "    if (i != holdout_id) y_f ~ student_t(4, mu_f, sigma_mix);\n")
  
  model_blk <- paste0(model_blk, "  }\n}\n")
  
  # --- 3. GENERATED QUANTITIES БЕЗ ИЗМЕНЕНИЙ (Считает log_lik для всех i) ---
  gen_blk <- "generated quantities {\n  vector[N_mix] log_lik;\n"
  gen_blk <- paste0(gen_blk, "  {\n")
  
  loop_body <- "\n    for (i in 1:N_mix) {\n"
  if (h) {
    loop_body <- paste0(loop_body, "      vector[N_tracers] mu_f;\n      vector[N_tracers] y_f = y_mix[i]';\n")
  } else {
    loop_body <- paste0(loop_body, "      vector[N_tracers] mu_f;\n      vector[N_tracers] y_f = y_mix[i]';\n")
  }
  
  if (coda) {
    loop_body <- paste0(loop_body, "      vector[N_tracers_raw] mu_pred_prop;\n      vector[N_tracers_raw] src_prop_matrix[N_groups];\n")
    loop_body <- paste0(loop_body, "      for (g in 1:N_groups) src_prop_matrix[g] = softmax(V * mu_source[g]);\n")
    loop_body <- paste0(loop_body, "      for (t in 1:N_tracers_raw) {\n        real val = 0;\n        for (g in 1:N_groups) {\n")
    if (b) loop_body <- paste0(loop_body, sprintf("          val += %s * src_prop_matrix[g, t] * exp(beta_sorting[t] * (ssa_mix[i] - mean(ssa_source)));\n", current_p_g))
    else   loop_body <- paste0(loop_body, sprintf("          val += %s * src_prop_matrix[g, t];\n", current_p_g))
    loop_body <- paste0(loop_body, "        }\n        mu_pred_prop[t] = val + 1e-9;\n      }\n")
    if (bias) loop_body <- paste0(loop_body, "      mu_f = (log(mu_pred_prop)' * V)' + delta_bias;\n")
    else      loop_body <- paste0(loop_body, "      mu_f = (log(mu_pred_prop)' * V)';\n")
  } else {
    loop_body <- paste0(loop_body, "      for (t in 1:N_tracers) {\n        real val = 0;\n")
    if (b) loop_body <- paste0(loop_body, sprintf("        for (g in 1:N_groups) val += %s * mu_source[g, t] * exp(beta_sorting[t] * (ssa_mix[i] - mean(ssa_source)));\n", current_p_g))
    else  loop_body <- paste0(loop_body, sprintf("        for (g in 1:N_groups) val += %s * mu_source[g, t];\n", current_p_g))
    if (bias) loop_body <- paste0(loop_body, "        mu_f[t] = val + delta_bias[t];\n      }\n")
    else      loop_body <- paste0(loop_body, "        mu_f[t] = val;\n      }\n")
  }
  
  if (cov_mix) loop_body <- paste0(loop_body, "      log_lik[i] = multi_student_t_cholesky_lpdf(y_f | 4, mu_f, diag_pre_multiply(L_sigma, L_Omega));\n")
  else         loop_body <- paste0(loop_body, "      log_lik[i] = student_t_lpdf(y_f | 4, mu_f, sigma_mix);\n")
  
  loop_body <- paste0(loop_body, "    }\n")
  gen_blk <- paste0(gen_blk, loop_body, "  }\n}\n")
  
  return(paste0(data_blk, param_blk, model_blk, gen_blk))
}

# ==============================================================================
# 2. USER INTERFACE (UI)
# ==============================================================================
ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel(
    title = HTML("<b>BASIN: Ba</b>ye<b>si</b>a<b>n</b> Sediment Source Apportiontment Tool"),
    windowTitle = "BASIN: Bayesian Sediment Source Apportiontment Tool"
  ),
  
  sidebarLayout(
    sidebarPanel(
      # Input section
      fileInput("file", "1. Upload data (CSV)", accept = ".csv"),
      uiOutput("ui_group_col"),
      uiOutput("ui_mix_label"),
      uiOutput("ui_mix_cov"),
      uiOutput("ui_size_col"),
      uiOutput("ui_tracer_select"),
      hr(),
      
      radioButtons("filter_mode", "Filtering mode (Particle Size / SSA):",
                   choices = c("Raw data only" = "raw", 
                               "Corrected only" = "corr", 
                               "Both (Hard filter)" = "both", 
                               "Any (Soft filter)" = "any"),
                   selected = "any", 
                   inline = TRUE),
      
      checkboxInput("opt_buddy", "1. Enable Geochemical Filter (Structure Preservation)", FALSE),
      checkboxInput("opt_hull", "2. PCA Convex Hull Penalty (Drop-out analysis)", FALSE),
      hr(),
      actionButton("run_filter", "2. Tracers filtering", class = "btn-warning", style="width:100%"),
      br(), br(),
      
      # Model Configuration
      h4("3. Model constructor:"),
      fluidRow(
        column(6, checkboxInput("opt_h", "Hierarchical", TRUE)),
        column(6, checkboxInput("opt_b", "Beta correction", FALSE))
      ),
      fluidRow(
        column(6, checkboxInput("opt_coda", "CoDA (CLR Transform)", FALSE)),
        column(6, checkboxInput("opt_bias", "Robust Bayesian Bias Absorption", FALSE))
      ),
      radioButtons("cov_mode", "Error Structure (Covariance):",
                   choices = c("Independent (No Covariance)" = "none",
                               "Mixture Covariance (Recommended)" = "mix",
                               "Full Bayes (Mixtures + Sources)" = "full"),
                   selected = "mix"),
      hr(),
      
      # MCMC Settings
      selectInput("calc_method", "4. Method:", 
                  choices = c("Markov Chains (MCMC)" = "mcmc", "Turbo (Fast VB) - experimental" = "vb")),
      helpText("MCMC is more accurate, VB is faster."),
      fluidRow(
        column(6, numericInput("iter", "Iterations:", 2000)),
        column(6, numericInput("chains", "Chains:", 4))
      ),
      hr(),
      
      textOutput("count_status"),
      br(),
      # Execution button
      actionButton("run_model", "5. Compile and run the model", class = "btn-primary", style="width:100%"),
      hr(),
      
      # Geographic/Connectivity Weights
      uiOutput("conn_inputs")
    ),
    
    mainPanel(
      tabsetPanel(
        # Tab 1: Raw Data inspection
        tabPanel("Data table", DTOutput("raw_table")),
        
        # Tab 2: Tracer selection & PCA
        tabPanel("Tracer selection", 
                 fluidRow(
                   column(5, 
                          h4("1. Final choice"),
                          helpText("Check or Uncheck variables before start"),
                          uiOutput("ui_final_tracers_chk"),
                          hr(),
                          h4("2. Filtering report"),
                          verbatimTextOutput("filter_log", placeholder = TRUE),
                          tags$style(type="text/css", "#filter_log { max-height: 300px; overflow-y: scroll; white-space: pre-wrap; }")
                   ),
                   column(7, 
                          h4("3. Multi-proxy Discriminant Visualization"),
                          tabsetPanel(
                            tabPanel("PCA (Variance Explorer)", plotOutput("pca_plot", height = "500px")),
                            tabPanel("LDA (Source Separation)", 
                                     plotOutput("lda_plot", height = "500px"),
                                     helpText("LDA maximizes the distance between source groups. Ellipses represent 95% confidence intervals.")
                            )
                          )
                   )
                 )
        ),
        
        # Tab 3: Model Results
        tabPanel("Results", 
                 plotOutput("density_plot", height = "500px"), 
                 hr(),
                 conditionalPanel(condition = "input.opt_b == true", plotOutput("beta_plot", height = "400px")),
                 hr(),
                 tableOutput("summary_table"),
                 hr(),
                 fluidRow(
                   column(6, downloadButton("download_ind_res", "Download Individual Results (CSV)", class = "btn-info", style="width:100%")),
                   column(6, downloadButton("download_chains", "Download Raw MCMC Chains (CSV)", class = "btn-warning", style="width:100%"))
                 ),
                 hr(),
                 h4("Individual Sample Posterior Distributions"),
                 fluidRow(
                   column(4, uiOutput("ui_show_sample_id")),
                   column(8, plotOutput("ind_density_plot", height = "350px"))
                 ),
                 br(),
                 downloadButton("download_all_plots", "Download All Density Plots (ZIP)", class = "btn-primary", style="width:100%")
        ),
        
        tabPanel("Reconstruction Analysis", 
                 uiOutput("reconstruction_ui")
        ),
        
        # Tab 4: Validation
        tabPanel("Validation (virtual mixtures)",
                 br(),
                 br(),
                 fluidRow(
                   column(3, 
                          radioButtons("val_gen_mode", "VM Generation Logic:",
                                       choices = c("1. Full Stoch. (Rand Props + MVN Tracers)" = "stoch_rand", 
                                                   "2. Semi-Stoch. (Simplex Props + MVN Tracers)" = "stoch_simp",
                                                   "3. Semi-Det. (Rand Props + Mean Tracers)" = "det_rand",
                                                   "4. Full Det. (Simplex Props + Mean Tracers)" = "det_simp"),
                                       selected = "stoch_rand")
                   ),
                   column(3, 
                          conditionalPanel(
                            condition = "input.val_gen_mode == 'stoch_rand' || input.val_gen_mode == 'det_rand'",
                            numericInput("val_n", "Number of VMs:", 20, min=5, max=500)
                          ),
                          conditionalPanel(
                            condition = "input.val_gen_mode == 'det_simp' || input.val_gen_mode == 'stoch_simp'",
                            numericInput("val_step", "Simplex Step (e.g. 0.05):", 0.05, min=0.01, max=0.50, step=0.01)
                          )
                   ),
                   column(3, 
                          selectInput("val_noise", "Data noise level:", 
                                      choices = c("Not (0%)"=0, "Low (5%)"=0.05, "Medium (15%)"=0.15, "High (25%)"=0.25), selected=0.05)
                   ),
                   column(3, 
                          actionButton("run_val", "1. Generate and estimate", class = "btn-success", style="width:100%"),
                          br(), br(),
                          downloadButton("download_val_data", "2. Download results (CSV)", style="width:100%")
                   )
                 ),
                 hr(),
                 h4("The model accuracy (True vs Predicted):"),
                 plotOutput("val_plot", height = "600px"),
                 h4("Errors & Performance Metrics:"),
                 tableOutput("val_metrics")
        ),
        
        # Tab 4.5: Model Comparison (LOO-CV)
        tabPanel("Model Comparison", 
                 br(),
                 h4("Bayesian Model Comparison & Diagnostics"),
                 helpText("Select multiple models to rank their predictive performance. The model with elpd_diff = 0 is the best."),
                 hr(),
                 fluidRow(
                   column(4, 
                          uiOutput("ui_multi_compare_selector"),
                          hr(),
                          actionButton("run_loo_compare", "Run LOO-CV", class = "btn-success", style="width:100%"),
                          hr(),
                          actionButton("run_exact_loo", "Run Exact LOO-CV (SLOW!)", class = "btn-warning", style="width:100%")
                   ),
                   column(8, 
                          verbatimTextOutput("loo_results_text"),
                          tags$style(type="text/css", "#loo_results_text { font-weight: bold; background-color: #f8f9fa; }")
                   )
                 ),
                 hr(),
                 h4("Pareto k Diagnostics (Influential Observations)"),
                 helpText("Points above the red dashed line (k > 0.7) strongly influence the model and might be outliers or poorly fitted mixtures."),
                 fluidRow(
                   column(4, uiOutput("ui_pareto_selector")),
                   column(8, plotOutput("pareto_plot", height = "400px"))
                 )
        ),
        
        # Tab 5: Help
        tabPanel("User Guide & Math", includeHTML("basin_help.html")) 
      )
    )
  ),
  
  hr(),
  tags$div(
    style = "text-align: center; color: #7f8c8d; padding-bottom: 20px; font-size: 0.9em;",
    tags$p(HTML("Developed by <b>Sergey Kharchenko</b>. Supported by <b>RSF</b>, grant No. 25-17-00143."))
  )
)

# ==============================================================================
# 3. SERVER LOGIC
# ==============================================================================
server <- function(input, output, session) {
  
  values <- reactiveValues(df = NULL, raw_df = NULL, groups = NULL, tracers = NULL, 
                           fit = NULL, stan_out = "", compiled_list = list(), level_names = NULL,
                           saved_models = list())
  
  observeEvent(list(input$opt_h, input$opt_b, input$cov_mode, input$opt_coda, input$calc_method, input$opt_buddy, input$opt_hull, input$opt_bias), {
    values$fit <- NULL
    values$stan_out <- "Settings have been changed. Re-run the model."
  })
  
  # 1. LOAD AND SAFE RAW DATA
  observeEvent(input$file, {
    values$raw_df <- data.table::fread(input$file$datapath, data.table=F, check.names=F)
    values$df <- values$raw_df
  })
  
  # 2. AUTO-SORTING
  observeEvent(list(input$group_col, input$mix_label), {
    req(values$raw_df, input$group_col, input$mix_label)
    df <- values$raw_df
    if(input$group_col %in% colnames(df)) {
      if(input$mix_label %in% df[[input$group_col]]) {
        # Divide sources/mixtures
        df_src <- df[df[[input$group_col]] != input$mix_label, ]
        df_mix <- df[df[[input$group_col]] == input$mix_label, ]
        
        # Sorting by abc
        df_src <- df_src[order(df_src[[input$group_col]]), ]
        values$df <- rbind(df_src, df_mix)
      } else {
        values$df <- df
      }
    }
  }, ignoreInit = TRUE)
  
  # Dynamic selectors rendering
  output$ui_group_col <- renderUI({ req(values$raw_df); selectInput("group_col", "Sources column:", choices = colnames(values$raw_df), selected = colnames(values$raw_df)[2]) })
  output$ui_mix_label <- renderUI({ req(input$group_col, values$raw_df); selectInput("mix_label", "Mixture label:", choices = unique(values$raw_df[[input$group_col]]), selected = tail(unique(values$raw_df[[input$group_col]]), 1)) })
  output$ui_mix_cov <- renderUI({ req(values$df); cols <- colnames(values$df); selectInput("mix_cov", "Mixture grouping covariate (e.g. Season/Trap):", choices = c("None", cols), selected = "None") })
  output$ui_size_col <- renderUI({ req(values$df); num_cols <- colnames(values$df)[sapply(values$df, is.numeric)]; selectInput("size_col", "SSA or SSA-proxy column (optional):", choices = c("None", num_cols), selected = "None") })
  
  output$ui_tracer_select <- renderUI({
    req(values$df)
    all_nums <- colnames(values$df)[sapply(values$df, is.numeric)]
    exclude <- c("ID", "id", "Id")
    if (!is.null(input$size_col) && input$size_col != "None") exclude <- c(exclude, input$size_col)
    if (!is.null(input$mix_cov) && input$mix_cov != "None") exclude <- c(exclude, input$mix_cov)
    cand <- setdiff(all_nums, exclude)
    selectizeInput("tracers_to_use", "Choose tracers:", choices = cand, selected = cand, multiple = TRUE)
  })
  
  output$ui_show_sample_id <- renderUI({
    req(values$fit, values$df, input$mix_label, input$group_col)
    df_mix <- values$df[values$df[[input$group_col]] == input$mix_label, ]
    mix_ids <- if("ID" %in% colnames(df_mix)) as.character(df_mix$ID) else paste0("Mix_", 1:nrow(df_mix))
    selectInput("show_sample_id", "Select Sample:", choices = mix_ids)
  })
  
  output$raw_table <- renderDT({ req(values$df); datatable(values$df, options = list(pageLength = 10, scrollX = TRUE)) })
  output$count_status <- renderText({ paste("Selected for analysis:", length(input$final_selected_tracers), "tracers") })
  
  output$conn_inputs <- renderUI({
    req(values$df, input$group_col, input$mix_label)
    grs <- unique(values$df[[input$group_col]]); grs <- grs[grs != input$mix_label]
    values$groups <- grs
    tagList(h4("Geographical weights:"), p(tags$small("e.g. Area * Erosion * IC")),
            lapply(grs, function(g) numericInput(paste0("w_", g), paste("Source:", g), value = 1.0, step = 0.1)))
  })
  
  # Filtering Logic
  values$suggested_tracers <- NULL
  values$filter_msg <- "Click Tracers filtering before running the model"
  
  # ADVANCED SSA-based filter
  observeEvent(input$run_filter, {
    req(values$df, input$tracers_to_use)
    
    withProgress(message = 'Analyzing tracers (Geochemical AI)...', {
      df <- values$df
      gr_col <- input$group_col
      m_lab <- input$mix_label
      ssa_col <- input$size_col
      
      candidates <- input$tracers_to_use
      src <- df[df[[gr_col]] != m_lab, ]
      mix <- df[df[[gr_col]] == m_lab, ]
      
      # --- 0. First check ---
      if(nrow(src) == 0 || nrow(mix) == 0) {
        showNotification("Error: Source or Mixture group is empty.", type="error")
        return()
      }
      
      # Computation SSA factor
      ssa_factor <- 1.0
      if (!is.null(ssa_col) && ssa_col != "None" && ssa_col %in% colnames(df)) {
        mean_ssa_mix <- mean(mix[[ssa_col]], na.rm = TRUE)
        mean_ssa_src <- mean(src[[ssa_col]], na.rm = TRUE)
        if (!is.na(mean_ssa_mix) && !is.na(mean_ssa_src) && mean_ssa_src > 0) {
          ssa_factor <- mean_ssa_mix / mean_ssa_src
        }
      }
      
      log_text <- c(paste("--- SSA-AWARE FILTERING (Mode:", input$filter_mode, ") ---"))
      log_text <- c(log_text, paste("[INFO] SSA enrichment factor:", round(ssa_factor, 2)))
      
      # --- 1. RANGE TEST (1D) ---
      pass_range <- c()
      for(t in candidates) {
        if(!is.numeric(df[[t]])) next
        
        mv_raw <- mean(mix[[t]], na.rm = TRUE)
        mv_norm <- mv_raw / ssa_factor
        s_min <- min(src[[t]], na.rm = TRUE) * 0.95
        s_max <- max(src[[t]], na.rm = TRUE) * 1.05
        
        check_raw  <- (mv_raw >= s_min && mv_raw <= s_max)
        check_norm <- (mv_norm >= s_min && mv_norm <= s_max)
        
        res_check <- switch(input$filter_mode,
                            "raw"  = check_raw,
                            "corr" = check_norm,
                            "both" = (check_raw && check_norm),
                            "any"  = (check_raw || check_norm))
        
        if (is.na(res_check)) res_check <- FALSE
        
        if(res_check) {
          pass_range <- c(pass_range, t)
          # SALVAGED logic
          if (!check_raw && check_norm) {
            log_text <- c(log_text, paste0("[RANGE] ", t, ": SALVAGED by SSA correction"))
          }
        } else {
          log_text <- c(log_text, paste0("[RANGE] ", t, ": REMOVED (Outside source polygon)"))
        }
      }
      
      # --- 2. BUDDY FILTER (Structural integrity) ---
      pass_buddy <- pass_range 
      if (input$opt_buddy && length(pass_range) > 2) {
        log_text <- c(log_text, "--- GEOCHEMICAL ASSOCIATIONS FILTER ---")
        src_cor <- suppressWarnings(cor(src[, pass_range], use = "pairwise.complete.obs", method = "spearman"))
        src_cor[is.na(src_cor)] <- 0
        diag(src_cor) <- 0
        bad_buddies <- c()
        
        for (t in pass_range) {
          buddies <- pass_range[abs(src_cor[t, pass_range]) > 0.85]
          if (length(buddies) < 2) next
          broken_links <- 0
          
          for (b in buddies) {
            coeffs <- tryCatch({ coef(lm(src[[t]] ~ src[[b]])) }, error = function(e) NULL)
            if (is.null(coeffs)) next
            
            err_raw  <- abs(mean(mix[[t]]) - (coeffs[1] + coeffs[2] * mean(mix[[b]]))) / (abs(coeffs[1] + coeffs[2] * mean(mix[[b]])) + 1e-9)
            err_norm <- abs((mean(mix[[t]])/ssa_factor) - (coeffs[1] + coeffs[2] * (mean(mix[[b]])/ssa_factor))) / (abs(coeffs[1] + coeffs[2] * (mean(mix[[b]])/ssa_factor)) + 1e-9)
            
            active_error <- switch(input$filter_mode,
                                   "raw"  = err_raw,
                                   "corr" = err_norm,
                                   "both" = max(err_raw, err_norm),
                                   "any"  = min(err_raw, err_norm))
            
            if (is.na(active_error) || active_error > 0.10) broken_links <- broken_links + 1
          }
          if (broken_links / length(buddies) >= 0.5) {
            bad_buddies <- c(bad_buddies, t)
            log_text <- c(log_text, paste0("[BUDDY] ", t, ": REMOVED (Broken ", broken_links, "/", length(buddies), " links)"))
          }
        }
        pass_buddy <- setdiff(pass_range, bad_buddies)
      }
      
      # --- 3. PCA CONVEX HULL (n-D Outlier Detection) ---
      pass_hull <- pass_buddy 
      if (input$opt_hull && length(pass_buddy) > 4) {
        log_text <- c(log_text, "--- PCA CONVEX HULL PENALTY ---")
        mix_norm <- mix
        for(tr_name in pass_buddy) mix_norm[[tr_name]] <- mix[[tr_name]] / ssa_factor
        
        calc_pca_dist <- function(tr_list, target_mix_df) {
          pca <- tryCatch(prcomp(src[, tr_list, drop=F], scale. = TRUE), error = function(e) NULL)
          if (is.null(pca)) return(NA)
          s_pcs <- predict(pca, src[, tr_list, drop=F])
          m_pcs <- predict(pca, target_mix_df[, tr_list, drop=F])
          d <- sqrt(sum((colMeans(m_pcs[,1:2, drop=F]) - colMeans(s_pcs[,1:2, drop=F]))^2))
          r <- mean(sqrt(rowSums(sweep(s_pcs[,1:2, drop=F], 2, colMeans(s_pcs[,1:2, drop=F]), "-")^2)))
          return(d / (r + 1e-9))
        }
        
        get_active_dist <- function(tr_list) {
          d_raw  <- calc_pca_dist(tr_list, mix)
          d_norm <- calc_pca_dist(tr_list, mix_norm)
          switch(input$filter_mode,
                 "raw"  = d_raw,
                 "corr" = d_norm,
                 "both" = max(d_raw, d_norm, na.rm=T),
                 "any"  = min(d_raw, d_norm, na.rm=T))
        }
        
        base_dist <- get_active_dist(pass_buddy)
        bad_hull <- c()
        if (!is.na(base_dist) && base_dist > 1.2) {
          for (t in pass_buddy) {
            new_dist <- get_active_dist(setdiff(pass_buddy, t))
            if (!is.na(new_dist) && new_dist < (base_dist * 0.95)) {
              bad_hull <- c(bad_hull, t)
              log_text <- c(log_text, paste0("[HULL] ", t, ": REMOVED (High PCA leverage)"))
            }
          }
        }
        pass_hull <- setdiff(pass_buddy, bad_hull)
      }
      
      # --- 4. KW + CORR ---
      pass_kw <- c()
      for(t in pass_hull) {
        pval <- tryCatch({ kruskal.test(as.formula(paste0("`", t, "` ~ `", gr_col, "`")), data = src)$p.value }, error = function(e) NA)
        if(!is.na(pval) && pval < 0.05) {
          pass_kw <- c(pass_kw, t)
        } else {
          log_text <- c(log_text, paste0("[KW-TEST] ", t, ": REMOVED (No discrimination)"))
        }
      }
      
      final_pass <- pass_kw
      if(length(final_pass) > 1) {
        cor_mat <- cor(src[, final_pass], use = "pairwise.complete.obs")
        cor_mat[is.na(cor_mat)] <- 0
        high_c_idx <- caret::findCorrelation(cor_mat, cutoff = 0.95)
        if(length(high_c_idx) > 0) {
          for(rt in final_pass[high_c_idx]) log_text <- c(log_text, paste0("[CORR] ", rt, ": REMOVED (Redundant)"))
          final_pass <- final_pass[-high_c_idx]
        }
      }
      
      log_text <- c(log_text, "------------------------", paste("PASSED:", length(final_pass), "tracers."))
      values$suggested_tracers <- final_pass
      values$filter_msg <- paste(log_text, collapse="\n")
      showNotification(paste("Complete. Suggested:", length(final_pass), "tracers."), type="message")
    })
  })
  
  output$filter_log <- renderText({ values$filter_msg })
  output$ui_final_tracers_chk <- renderUI({
    req(input$tracers_to_use); choices <- input$tracers_to_use
    selected <- if(!is.null(values$suggested_tracers)) values$suggested_tracers else choices
    checkboxGroupInput("final_selected_tracers", label = NULL, choices = choices, selected = selected, inline = TRUE)
  })
  outputOptions(output, "ui_final_tracers_chk", suspendWhenHidden = FALSE)
  
  # Plots
  output$pca_plot <- renderPlot({
    req(input$final_selected_tracers, input$group_col, input$mix_label)
    tr <- input$final_selected_tracers; df <- values$df; m_lab <- input$mix_label; gr_col <- input$group_col
    src <- df[df[[gr_col]] != m_lab, ]; mix <- df[df[[gr_col]] == m_lab, ]
    
    # PCA
    pca <- prcomp(src[, tr], scale. = TRUE)
    var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
    
    plot_df <- rbind(
      data.frame(predict(pca, src[, tr]), Group = src[[gr_col]]), 
      data.frame(predict(pca, mix[, tr]), Group = m_lab)
    )
    hull <- plot_df %>% filter(Group != m_lab) %>% group_by(Group) %>% slice(chull(PC1, PC2))
    
    # --- Loadings ---
    loadings <- as.data.frame(pca$rotation[, 1:2])
    loadings$Tracer <- rownames(loadings)
    
    scale_factor <- min(max(abs(plot_df$PC1)), max(abs(plot_df$PC2))) / max(abs(pca$rotation[,1:2])) * 0.8
    loadings$PC1 <- loadings$PC1 * scale_factor
    loadings$PC2 <- loadings$PC2 * scale_factor
    
    ggplot() + 
      geom_polygon(data=hull, aes(x=PC1, y=PC2, fill=Group), alpha=0.2) +
      
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
      
      geom_point(data=filter(plot_df, Group!=m_lab), aes(x=PC1, y=PC2, color=Group), size=5, alpha=0.8) +
      
      geom_segment(data=loadings, aes(x=0, y=0, xend=PC1, yend=PC2), 
                   arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.7, alpha = 0.6) +
      geom_text(data=loadings, aes(x=PC1*1.15, y=PC2*1.15, label=Tracer), 
                color = "black", fontface = "bold", size = 5) +
      
      geom_point(data=filter(plot_df, Group==m_lab), aes(x=PC1, y=PC2), 
                 shape=23, fill="red", size=8, stroke=2) +
      
      scale_color_viridis_d(option = "turbo") + 
      scale_fill_viridis_d(option = "turbo") +
      
      theme_bw(base_size = 18) + 
      theme(axis.line = element_line(linewidth=1.2), 
            axis.text = element_text(face="bold"),
            panel.grid.minor = element_blank()) +
      labs(title=paste("PCA Biplot (Mixing Polygon):", m_lab), 
           x=paste0("PC1 (", var_exp[1], "%)"), 
           y=paste0("PC2 (", var_exp[2], "%)"))
  })
  
  output$lda_plot <- renderPlot({
    req(input$final_selected_tracers, input$group_col, input$mix_label)
    tr <- input$final_selected_tracers; df <- values$df; m_lab <- input$mix_label; gr_col <- input$group_col
    src_clean <- df[df[[gr_col]] != m_lab, ] %>% drop_na(all_of(tr))
    if(length(unique(src_clean[[gr_col]])) < 2 || nrow(src_clean) <= length(tr)) return(NULL)
    
    lda_mod <- tryCatch({ MASS::lda(src_clean[, tr], grouping = src_clean[[gr_col]]) }, error = function(e) return(NULL))
    if(is.null(lda_mod)) return(NULL)
    
    lda_values <- predict(lda_mod)
    lda_plot_df <- data.frame(lda_values$x, Group = src_clean[[gr_col]])
    prop_lda <- lda_mod$svd^2 / sum(lda_mod$svd^2)
    
    # --- Scaling ---
    is_1d <- ncol(lda_values$x) == 1
    
    if(is_1d) {
      # Для 2 источников добавляем Jitter
      lda_plot_df$LD2 <- rnorm(nrow(lda_plot_df), 0, 0.05)
      colnames(lda_plot_df)[1] <- "LD1"
      x_label <- "LD1 (100%)"; y_label <- "Jitter"
      
      # Векторы лежат на оси X
      ld_loadings <- data.frame(LD1 = lda_mod$scaling[, 1], LD2 = 0, Tracer = rownames(lda_mod$scaling))
    } else {
      # Для 3+ источников полноценные 2 оси
      x_label <- paste0("LD1 (", round(prop_lda[1]*100, 1), "%)")
      y_label <- paste0("LD2 (", round(prop_lda[2]*100, 1), "%)")
      ld_loadings <- data.frame(LD1 = lda_mod$scaling[, 1], LD2 = lda_mod$scaling[, 2], Tracer = rownames(lda_mod$scaling))
    }
    
    lda_plot_df$SampleID <- if("ID" %in% colnames(src_clean)) src_clean$ID else seq_len(nrow(lda_plot_df))
    
    # Динамическое масштабирование стрелок
    max_val_points <- max(abs(lda_plot_df$LD1), if(!is_1d) abs(lda_plot_df$LD2) else 0)
    max_val_loadings <- max(abs(ld_loadings$LD1), abs(ld_loadings$LD2))
    scale_factor <- (max_val_points / (max_val_loadings + 1e-9)) * 0.8
    
    ld_loadings$LD1 <- ld_loadings$LD1 * scale_factor
    ld_loadings$LD2 <- ld_loadings$LD2 * scale_factor
    
    ggplot() + 
      # Эллипсы (Confidence intervals)
      stat_ellipse(data=lda_plot_df, aes(x=LD1, y=LD2, color=Group, fill=Group), 
                   geom = "polygon", alpha = 0.2, show.legend = FALSE) + 
      
      # Оси центра
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
      
      # Точки источников
      geom_point(data=lda_plot_df, aes(x=LD1, y=LD2, color=Group, fill=Group), size = 4, alpha = 0.8) + 
      geom_text(data=lda_plot_df, aes(x=LD1, y=LD2, label = SampleID, color=Group), 
                vjust = -1.2, size = 4, show.legend = FALSE) +
      
      # Векторы трассеров
      geom_segment(data=ld_loadings, aes(x=0, y=0, xend=LD1, yend=LD2), 
                   arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.7, alpha = 0.6) +
      geom_text(data=ld_loadings, aes(x=LD1*1.15, y=LD2*1.15, label=Tracer), 
                color = "black", fontface = "bold", size = 5) +
      
      theme_bw(base_size = 18) + 
      scale_color_viridis_d(option = "turbo") + 
      scale_fill_viridis_d(option = "turbo") +
      labs(title = "Linear Discriminant Analysis", x = x_label, y = y_label) + 
      theme(legend.position = "bottom",
            axis.line = element_line(linewidth=1.2), 
            panel.grid.minor = element_blank())
  })
  
  # --- MAIN MODEL EXECUTION ---
  observeEvent(input$run_model, {
    req(input$final_selected_tracers)
    if(length(input$final_selected_tracers) < 2) { showNotification("Error: Select at least 2 tracers!", type="error"); return() }
    if(input$opt_b && (is.null(input$size_col) || input$size_col == "None")) showNotification("Warning: Beta correction enabled but no SSA column selected.", type="warning")

    
    withProgress(message = 'Preparing and calculating...', value = 0, {
      
      tr <- input$final_selected_tracers; m_lab <- input$mix_label; gr_col <- input$group_col
      y_src_raw <- as.matrix(values$df[values$df[[gr_col]] != m_lab, tr])
      y_mix_raw <- as.matrix(values$df[values$df[[gr_col]] == m_lab, tr])
      y_src_raw[y_src_raw <= 0] <- 1e-6; y_mix_raw[y_mix_raw <= 0] <- 1e-6
      
      if (!input$opt_coda) {
        y_src <- scale(y_src_raw)
        y_mix <- scale(y_mix_raw, center = attr(y_src, "scaled:center"), scale = attr(y_src, "scaled:scale"))
        V_matrix <- matrix(0, nrow=0, ncol=0) # Заглушка
        N_tr_raw <- ncol(y_mix)
      } else { 
        ys_safe <- y_src_raw; ys_safe[ys_safe <= 0] <- 1e-6
        ym_safe <- y_mix_raw; ym_safe[ym_safe <= 0] <- 1e-6
        
        ys_prop <- ys_safe / rowSums(ys_safe)
        ym_prop <- ym_safe / rowSums(ym_safe)
        
        # ILR Transform
        N_tr_raw <- ncol(ys_prop)
        V_matrix <- compositions::ilrBase(D = N_tr_raw)
        
        y_src <- as.matrix(log(ys_prop) %*% V_matrix)
        y_mix <- as.matrix(log(ym_prop) %*% V_matrix)
      }
      y_src <- as.matrix(y_src); y_mix <- as.matrix(y_mix)
      y_src <- as.matrix(y_src); y_mix <- as.matrix(y_mix)
      
      # 3. SSA (Avoiding NaN / Inf)
      ssa_m <- rep(0, nrow(y_mix)); ssa_s <- rep(0, nrow(y_src))
      if (!is.null(input$size_col) && input$size_col != "None" && input$size_col %in% colnames(values$df)) {
        ssa_m_raw <- values$df[values$df[[gr_col]] == m_lab, input$size_col]
        ssa_s_raw <- values$df[values$df[[gr_col]] != m_lab, input$size_col]
        ssa_m_raw[is.na(ssa_m_raw)] <- mean(ssa_m_raw, na.rm=TRUE); ssa_s_raw[is.na(ssa_s_raw)] <- mean(ssa_s_raw, na.rm=TRUE)
        
        # Standardize SSA (Z-score)
        if(!all(is.na(ssa_s_raw))) {
          mean_s <- mean(ssa_s_raw); sd_s <- sd(ssa_s_raw) + 1e-9
          ssa_m <- (ssa_m_raw - mean_s) / sd_s
          ssa_s <- (ssa_s_raw - mean_s) / sd_s
        }
      }
      
      df_mix_rows <- values$df[values$df[[gr_col]] == m_lab, ]
      has_groups <- FALSE; N_mix_levels <- NULL; mix_level_id <- NULL
      if (!is.null(input$mix_cov) && input$mix_cov != "None" && input$mix_cov %in% colnames(df_mix_rows)) {
        has_groups <- TRUE; cov_factor <- as.factor(df_mix_rows[[input$mix_cov]])
        mix_level_id <- as.numeric(cov_factor); N_mix_levels <- length(levels(cov_factor)); values$level_names <- levels(cov_factor) 
      } else values$level_names <- c("Overall")
      
      model_id <- paste0(as.integer(input$opt_h), as.integer(input$opt_b), 
                         input$cov_mode, as.integer(input$opt_coda), as.integer(has_groups), as.integer(input$opt_bias))
      
      if(is.null(values$compiled_list[[model_id]])) {
        showNotification("Compiling Stan model (this may take time)...", type="message")
        code <- generate_stan_code(input$opt_h, input$opt_b, input$cov_mode, input$opt_coda, has_groups, bias = input$opt_bias)
        values$compiled_list[[model_id]] <- stan_model(model_code = code)
      }
      compiled_mod <- values$compiled_list[[model_id]]
      
      values$compiled_model <- compiled_mod
      
      # Data List
      s_data <- list(
        N_mix = nrow(y_mix), N_sources = nrow(y_src), 
        N_tracers = ncol(y_mix),      # Для ILR это будет D-1
        N_tracers_raw = N_tr_raw,     # Исходное количество D
        N_groups = length(values$groups),
        y_mix = y_mix, y_source = y_src, 
        V = V_matrix,                 # Передаем базис
        group_id = as.numeric(as.factor(values$df[[gr_col]][values$df[[gr_col]] != m_lab])),
        ssa_mix = as.array(ssa_m), ssa_source = as.array(ssa_s), 
        conn_weights = sapply(values$groups, function(g) input[[paste0("w_", g)]])
      )
      if(has_groups) { s_data$N_mix_levels <- N_mix_levels; s_data$mix_level_id <- as.array(mix_level_id) }
      
      s_data$holdout_id <- 0
      values$stan_data <- s_data
      
      logs <- capture.output({
        tryCatch({
          if(input$calc_method == "vb") {
            values$fit <- vb(compiled_mod, data = s_data, init = "random", output_samples = 2000, seed = 123, tol_rel_obj=0.005)
          } else {
            values$fit <- sampling(compiled_mod, data = s_data, iter = input$iter, chains = input$chains, 
                                   init = 0.5, control = list(adapt_delta = 0.99, max_treedepth = 15))
          }
          suppressWarnings(gc())
        }, error = function(e) { showNotification(paste("Stan Error:", e$message), type="error"); values$stan_out <- paste("ERROR:", e$message) })
      }, type = "output")
      
      # --- MODEL AUTOSAVING (for LOO-CV) ---
      if (!is.null(values$fit) && isS4(values$fit)) {
        
        opts <- c()
        
        if(input$opt_h) opts <- c(opts, "Hier") else opts <- c(opts, "Flat")
        
        cov_label <- switch(input$cov_mode,
                            "none" = "Cov:No",
                            "mix"  = "Cov:Mix",
                            "full" = "Cov:Full")
        opts <- c(opts, cov_label)
        
        if(input$opt_coda) opts <- c(opts, "ILR") else opts <- c(opts, "Scale")
        
        if(input$opt_b) opts <- c(opts, "+Beta")
        if(input$opt_bias) opts <- c(opts, "+Bias")
        
        config_str <- paste(opts, collapse="-")
        base_name <- paste0("[", toupper(input$calc_method), "] ", 
                            length(tr), "Tr | ", config_str)
        
        model_name <- base_name
        counter <- 1
        while(model_name %in% names(values$saved_models)) {
          model_name <- paste0(base_name, " (", counter, ")")
          counter <- counter + 1
        }
        
        values$saved_models[[model_name]] <- values$fit
        
        values$last_model_name <- model_name
        
        # =========================================================================
        # Project Folder and Auto-Export
        # =========================================================================
        safe_folder_name <- gsub("[^A-Za-z0-9_-]", "_", model_name)
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        run_folder <- file.path("BASIN_Workspaces", paste0(timestamp, "_", safe_folder_name))
        dir.create(run_folder, recursive = TRUE, showWarnings = FALSE)
        values$current_run_folder <- run_folder 
        
        tryCatch({
          write.csv(as.data.frame(values$fit), file.path(run_folder, "01_MCMC_Raw_Chains.csv"), row.names = FALSE)
          saveRDS(values$fit, file.path(run_folder, "02_Stanfit_Object.rds"))
          
          p_name <- if(input$opt_h) "P_global" else if (exists("has_groups") && has_groups) "P_level" else "P"
          samples <- rstan::extract(values$fit, pars=p_name)[[1]]
          if (length(dim(samples)) == 3) samples <- apply(samples, c(1,3), mean) 
          if (is.null(dim(samples)) || length(dim(samples)) < 2) samples <- matrix(samples, ncol = length(values$groups))
          colnames(samples) <- values$groups
          p_plot <- ggplot(reshape2::melt(samples), aes(x=value, fill=Var2)) + 
            geom_density(aes(y = after_stat(scaled)), alpha=0.6, linewidth=1.2) + 
            scale_fill_viridis_d(option="mako") + theme_classic() + labs(title="Global Posterior")
          ggsave(file.path(run_folder, "03_Global_Posterior_Density.png"), plot = p_plot, width = 8, height = 5)
          
          if (input$opt_b && "beta_sorting" %in% names(rstan::extract(values$fit))) {
            betas <- rstan::extract(values$fit)$beta_sorting
            if (is.null(dim(betas)) || length(dim(betas)) < 2) betas <- matrix(betas, ncol = 1)
            colnames(betas) <- input$final_selected_tracers
            b_plot <- ggplot(reshape2::melt(betas), aes(x=Var2, y=value, fill=Var2)) + 
              geom_boxplot(alpha=0.7, show.legend=FALSE) + geom_hline(yintercept=0, color="red", linetype="dashed") + theme_minimal()
            ggsave(file.path(run_folder, "04_Beta_Coefficients.png"), plot = b_plot, width = 8, height = 5)
          }
          
          meta <- data.frame(Model = model_name, Tracers = paste(input$final_selected_tracers, collapse=", "), Date = Sys.time())
          write.csv(meta, file.path(run_folder, "00_Run_Metadata.csv"), row.names = FALSE)
        }, error = function(e) { showNotification(paste("Auto-export failed:", e$message), type = "warning") })
        # =========================================================================
      }
      # ----------------------------------------------------
      
      values$stan_out <- paste(logs, collapse = "\n")
      
      tryCatch({
        output_dir <- "BASIN_Saved_Results"
        if(!dir.exists(output_dir)) dir.create(output_dir)
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        mix_ids <- if("ID" %in% colnames(df_mix_rows)) as.character(df_mix_rows$ID) else paste0("Mix_", 1:nrow(df_mix_rows))
        N_mix <- nrow(df_mix_rows); groups <- values$groups
        auto_df <- data.frame(Sample_ID = mix_ids, check.names = FALSE)
        
        if(input$opt_h) {
          raw_samples <- rstan::extract(values$fit, pars="P_ind")$P_ind
          if (length(dim(raw_samples)) == 2) dim(raw_samples) <- c(dim(raw_samples)[1], 1, dim(raw_samples)[2])
          for(i in seq_along(groups)) {
            src_data <- raw_samples[,,i, drop=FALSE]; auto_df[[paste0(groups[i], " (Mean)")]] <- as.numeric(colMeans(src_data)); auto_df[[paste0(groups[i], " (Median)")]] <- as.numeric(apply(src_data, 2, median))
          }
        } else {
          if (has_groups) {
            raw_samples <- rstan::extract(values$fit, pars="P_level")$P_level
            if (length(dim(raw_samples)) == 2) dim(raw_samples) <- c(dim(raw_samples)[1], 1, dim(raw_samples)[2])
            for(i in seq_along(groups)) {
              level_means <- colMeans(raw_samples[,,i, drop=FALSE]); level_meds <- apply(raw_samples[,,i, drop=FALSE], 2, median)
              auto_df[[paste0(groups[i], " (Mean)")]] <- level_means[mix_level_id]; auto_df[[paste0(groups[i], " (Median)")]] <- level_meds[mix_level_id]
            }
          } else {
            raw_samples <- rstan::extract(values$fit, pars="P")$P
            if (is.vector(raw_samples)) raw_samples <- matrix(raw_samples, ncol=length(groups))
            for(i in seq_along(groups)) {
              auto_df[[paste0(groups[i], " (Mean)")]] <- rep(mean(raw_samples[,i]), N_mix); auto_df[[paste0(groups[i], " (Median)")]] <- rep(median(raw_samples[,i]), N_mix)
            }
          }
        }
        write.csv(auto_df, file.path(output_dir, paste0("Results_", timestamp, ".csv")), row.names = FALSE)
      }, error = function(e) { showNotification(paste("Auto-save failed:", e$message), type = "warning") })
    })
  })
  
  # --- PLOTS & TABLES ---
  output$density_plot <- renderPlot({
    req(values$fit); if (!isS4(values$fit)) return(NULL)
    has_groups <- (!is.null(input$mix_cov) && input$mix_cov != "None")
    p_name <- if(input$opt_h) "P_global" else if (has_groups) "P_level" else "P"
    
    samples_list <- tryCatch(rstan::extract(values$fit), error = function(e) NULL)
    if(is.null(samples_list) || !(p_name %in% names(samples_list))) return(NULL)
    
    samples <- samples_list[[p_name]]; if (length(dim(samples)) == 3) samples <- apply(samples, c(1,3), mean) 
    if (is.null(dim(samples)) || length(dim(samples)) < 2) samples <- matrix(samples, ncol = length(values$groups))
    colnames(samples) <- values$groups; df_p <- reshape2::melt(samples)
    
    ggplot(df_p, aes(x=value, fill=Var2)) + geom_density(aes(y = after_stat(scaled)), alpha=0.6, linewidth=1.2) + 
      scale_x_continuous(limits=c(0, 1), breaks = seq(0, 1, 0.2)) + scale_y_continuous(limits=c(0, 1.05), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
      scale_fill_viridis_d(option="mako") + theme_classic(base_size=20) + theme(axis.line=element_line(linewidth=1.5)) + labs(x="Proportion", y="Normalized Density (Max=1)", fill="Source")
  })
  
  output$ind_density_plot <- renderPlot({
    req(values$fit, input$show_sample_id)
    df_mix <- values$df[values$df[[input$group_col]] == input$mix_label, ]
    mix_ids <- if("ID" %in% colnames(df_mix)) as.character(df_mix$ID) else paste0("Mix_", 1:nrow(df_mix))
    sample_idx <- which(mix_ids == input$show_sample_id); if (length(sample_idx) == 0) return(NULL)
    
    all_samples <- rstan::extract(values$fit); ind_data <- NULL
    if("P_ind" %in% names(all_samples)) ind_data <- all_samples$P_ind[, sample_idx, ]
    else if ("P_level" %in% names(all_samples)) ind_data <- all_samples$P_level[, as.numeric(as.factor(df_mix[[input$mix_cov]]))[sample_idx], ]
    else ind_data <- all_samples$P
    
    if(is.null(ind_data)) return(NULL)
    if(is.null(dim(ind_data))) ind_data <- matrix(ind_data, ncol = length(values$groups))
    colnames(ind_data) <- values$groups; df_p <- reshape2::melt(ind_data); df_p$Var2 <- factor(df_p$Var2, levels = values$groups)
    
    ggplot(df_p, aes(x=value, fill=Var2)) + geom_density(aes(y = after_stat(scaled)), alpha=0.6, linewidth=0.8) + 
      scale_x_continuous(limits=c(0, 1), breaks = seq(0, 1, 0.2)) + scale_y_continuous(limits=c(0, 1.05), expand = c(0, 0)) +
      scale_fill_viridis_d(option="turbo") + theme_classic(base_size=16) + labs(title=paste("Posterior Density for Sample:", input$show_sample_id), x="Proportion", y="Normalized Density (Max=1)", fill="Source")
  })
  
  output$download_all_plots <- downloadHandler(
    filename = function() { paste0("All_Sample_Plots_", Sys.Date(), ".zip") },
    content = function(file) {
      req(values$fit); temp_dir <- tempdir(); plot_folder <- file.path(temp_dir, "plots"); if(!dir.exists(plot_folder)) dir.create(plot_folder)
      df_mix <- values$df[values$df[[input$group_col]] == input$mix_label, ]
      mix_ids <- if("ID" %in% colnames(df_mix)) as.character(df_mix$ID) else paste0("Mix_", 1:nrow(df_mix))
      all_samples <- rstan::extract(values$fit)
      has_groups <- (!is.null(input$mix_cov) && input$mix_cov != "None"); cov_factor <- if(has_groups) as.numeric(as.factor(df_mix[[input$mix_cov]])) else NULL
      
      withProgress(message = 'Generating plots...', value = 0, {
        for (i in 1:length(mix_ids)) {
          ind_data <- NULL
          if("P_ind" %in% names(all_samples)) ind_data <- all_samples$P_ind[, i, ]
          else if ("P_level" %in% names(all_samples)) ind_data <- all_samples$P_level[, cov_factor[i], ]
          else ind_data <- all_samples$P
          
          if(is.null(dim(ind_data))) ind_data <- matrix(ind_data, ncol = length(values$groups))
          colnames(ind_data) <- values$groups; df_p <- reshape2::melt(ind_data); df_p$Var2 <- factor(df_p$Var2, levels = values$groups)
          p <- ggplot(df_p, aes(x=value, fill=Var2)) + geom_density(aes(y = after_stat(scaled)), alpha=0.6) + scale_x_continuous(limits=c(0, 1)) + scale_fill_viridis_d(option="turbo") + theme_classic() + labs(title=paste("Sample:", mix_ids[i]), x="Proportion", y="Density")
          ggsave(file.path(plot_folder, paste0(mix_ids[i], ".png")), plot = p, width = 6, height = 4, dpi = 150)
          incProgress(1/length(mix_ids))
        }
      })
      zip::zip(zipfile = file, files = list.files(plot_folder, full.names = FALSE), root = plot_folder)
    }
  )
  
  output$beta_plot <- renderPlot({
    req(values$fit, input$opt_b); samples_all <- rstan::extract(values$fit); if (!"beta_sorting" %in% names(samples_all)) return(NULL)
    betas <- samples_all$beta_sorting; if (is.null(dim(betas)) || length(dim(betas)) < 2) betas <- matrix(betas, ncol = 1)
    colnames(betas) <- input$final_selected_tracers; df_b <- reshape2::melt(betas)
    ggplot(df_b, aes(x=Var2, y=value, fill=Var2)) + geom_boxplot(alpha=0.7, show.legend=F) + geom_hline(yintercept=0, color="red", linetype="dashed") + theme_minimal(base_size=18) + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) + labs(title="Beta Coefficients (Particle Size Correction)", x="Tracer", y="Beta Value")
  })
  
  output$summary_table <- renderTable({
    req(values$fit)
    has_groups <- (!is.null(input$mix_cov) && input$mix_cov != "None")
    p_param <- if(input$opt_h) "P_global" else if(has_groups) "P_level" else "P"
    res_p <- as.data.frame(summary(values$fit, pars=p_param)$summary)
    samples_p <- rstan::extract(values$fit, pars=p_param)[[1]]
    if (length(dim(samples_p)) == 3) samples_p <- apply(samples_p, c(1,3), mean)
    if(is.null(dim(samples_p))) samples_p <- matrix(samples_p, ncol = length(values$groups))
    
    q25 <- apply(samples_p, 2, quantile, probs = 0.25); q75 <- apply(samples_p, 2, quantile, probs = 0.75)
    medians_p <- apply(samples_p, 2, median); median_p_corr <- medians_p / sum(medians_p)
    modes_p <- apply(samples_p, 2, function(x) { d <- density(x); d$x[which.max(d$y)] }); modes_p_corr <- modes_p / sum(modes_p)
    
    df_final <- data.frame(
      Source = values$groups, mean = res_p$mean[1:length(values$groups)], sd = res_p$sd[1:length(values$groups)],
      `2.5%` = res_p$`2.5%`[1:length(values$groups)], `25% (Q1)` = q25, median = medians_p, `median corr` = median_p_corr,
      mode = modes_p, `mode corr` = modes_p_corr, `75% (Q3)` = q75, `97.5%` = res_p$`97.5%`[1:length(values$groups)],
      `W50 (IQR)` = q75 - q25, check.names = FALSE
    )
    if("n_eff" %in% colnames(res_p)) { df_final$n_eff <- res_p$n_eff[1:length(values$groups)]; df_final$Rhat <- res_p$Rhat[1:length(values$groups)] }
    
    if(input$opt_h) {
      res_phi <- as.data.frame(summary(values$fit, pars="phi_ind")$summary); samples_phi <- rstan::extract(values$fit, pars="phi_ind")[[1]]
      q25phi <- quantile(samples_phi, 0.25); q75phi <- quantile(samples_phi, 0.75)
      df_phi <- data.frame(Source = "Phi (Precision)", mean = res_phi$mean, sd = res_phi$sd, `2.5%` = res_phi$`2.5%`, `25% (Q1)` = q25phi, median = median(samples_phi), `median corr` = NA, mode = density(samples_phi)$x[which.max(density(samples_phi)$y)], `mode corr` = NA, `75% (Q3)` = q75phi, `97.5%` = res_phi$`97.5%`, `W50 (IQR)` = q75phi - q25phi, check.names = FALSE)
      if("n_eff" %in% colnames(res_phi)) { df_phi$n_eff <- res_phi$n_eff; df_phi$Rhat <- res_phi$Rhat }
      df_final <- rbind(df_final, df_phi)
    }
    return(df_final)
  }, digits = 3, na = "-")
  
  output$download_ind_res <- downloadHandler(
    filename = function() { paste0("Individual_Results_Detailed_", Sys.Date(), ".csv") },
    content = function(file) {
      req(values$fit, values$df); if (!isS4(values$fit)) return(NULL)
      df_mix_rows <- values$df[values$df[[input$group_col]] == input$mix_label, ]
      mix_ids <- if ("ID" %in% colnames(df_mix_rows)) as.character(df_mix_rows$ID) else paste0("Mix_", 1:nrow(df_mix_rows))
      N_mix <- nrow(df_mix_rows); groups <- values$groups; final_df <- data.frame(Sample_ID = mix_ids, check.names = FALSE)
      
      tryCatch({
        has_groups <- (!is.null(input$mix_cov) && input$mix_cov != "None")
        if (input$opt_h) {
          raw_samples <- rstan::extract(values$fit, pars = "P_ind")$P_ind
          if (length(dim(raw_samples)) == 2) dim(raw_samples) <- c(dim(raw_samples)[1], 1, dim(raw_samples)[2])
          for (i in seq_along(groups)) {
            src_matrix <- matrix(raw_samples[, , i, drop = FALSE], nrow = dim(raw_samples)[1], ncol = N_mix)
            q25 <- as.numeric(apply(src_matrix, 2, quantile, probs = 0.25)); q75 <- as.numeric(apply(src_matrix, 2, quantile, probs = 0.75))
            final_df[[paste0(groups[i], " (Mean)")]] <- as.numeric(colMeans(src_matrix)); final_df[[paste0(groups[i], " (Median)")]] <- as.numeric(apply(src_matrix, 2, median))
            final_df[[paste0(groups[i], " (2.5%)")]] <- as.numeric(apply(src_matrix, 2, quantile, probs = 0.025)); final_df[[paste0(groups[i], " (97.5%)")]] <- as.numeric(apply(src_matrix, 2, quantile, probs = 0.975))
            final_df[[paste0(groups[i], " (25%)")]] <- q25; final_df[[paste0(groups[i], " (75%)")]] <- q75; final_df[[paste0(groups[i], " (W50)")]] <- q75 - q25
          }
        } else if (has_groups) {
          raw_samples <- rstan::extract(values$fit, pars = "P_level")$P_level
          if (length(dim(raw_samples)) == 2) dim(raw_samples) <- c(dim(raw_samples)[1], 1, dim(raw_samples)[2])
          cov_factor <- as.numeric(as.factor(df_mix_rows[[input$mix_cov]]))
          for (i in seq_along(groups)) {
            src_data <- raw_samples[,,i, drop=FALSE]
            final_df[[paste0(groups[i], " (Mean)")]] <- colMeans(src_data)[cov_factor]; final_df[[paste0(groups[i], " (Median)")]] <- apply(src_data, 2, median)[cov_factor]
          }
          final_df$Note <- "Level-Specific Model"
        } else {
          raw_samples <- rstan::extract(values$fit, pars = "P")$P
          if (is.vector(raw_samples)) raw_samples <- matrix(raw_samples, ncol = length(groups))
          for (i in seq_along(groups)) {
            src_vec <- raw_samples[, i]; q25 <- quantile(src_vec, 0.25); q75 <- quantile(src_vec, 0.75)
            final_df[[paste0(groups[i], " (Mean)")]] <- rep(mean(src_vec), N_mix); final_df[[paste0(groups[i], " (Median)")]] <- rep(median(src_vec), N_mix)
            final_df[[paste0(groups[i], " (2.5%)")]] <- rep(quantile(src_vec,0.025), N_mix); final_df[[paste0(groups[i], " (97.5%)")]] <- rep(quantile(src_vec,0.975), N_mix)
            final_df[[paste0(groups[i], " (25%)")]] <- rep(q25, N_mix); final_df[[paste0(groups[i], " (75%)")]] <- rep(q75, N_mix); final_df[[paste0(groups[i], " (W50)")]] <- rep(q75-q25, N_mix)
          }
          final_df$Note <- "Static Model"
        }
        write.csv(final_df, file, row.names = FALSE)
      }, error = function(e) { showNotification(paste("Download Error:", e$message), type = "error") })
    }
  )
  
  output$download_chains <- downloadHandler(
    filename = function() { paste0("MCMC_Raw_Chains_", Sys.Date(), ".csv") },
    content = function(file) { req(values$fit); write.csv(as.data.frame(values$fit), file, row.names = FALSE) }
  )
  
  # --- RECONSTRUCTION UI ---
  output$reconstruction_ui <- renderUI({
    if (is.null(values$fit)) return(tags$div(style="padding:50px; text-align:center; color:gray;", h3("No model results available.")))
    
    dynamic_height <- max(450, ceiling(length(input$final_selected_tracers) / 3) * 280)
    has_beta <- input$opt_b && "beta_sorting" %in% names(rstan::extract(values$fit))
    
    tagList(
      br(), h4("Geochemical Mass-Balance Reconstruction"), hr(),
      
      fluidRow(
        column(8, 
               if(has_beta) {
                 radioButtons("recon_mode", "Reconstruction Mode:",
                              choices = c("1. Conservative (Raw / No Beta correction)" = "raw", 
                                          "2. Corrected (Beta size-correction applied)" = "corr"),
                              selected = "raw", inline = TRUE)
               } else {
                 p(tags$b("Reconstruction Mode: "), "Conservative (Beta correction was disabled in the model)")
               }
        ),
        column(4, align = "right", 
               downloadButton("download_recon_csv", "Download Full Details", class = "btn-info"))
      ),
      br(),
      
      # Graphs and tables
      fluidRow(
        column(8, plotOutput("recon_plot", height = paste0(dynamic_height, "px"))), 
        column(4, h4("Error per Tracer (%)"), tableOutput("recon_metrics_table"))
      ),
      hr(), h4("Sample-specific reconstruction performance"), 
      plotOutput("recon_error_boxplot", height = paste0(dynamic_height, "px"))
    )
  })
  
  recon_calc <- reactive({
    req(values$fit, values$df, input$final_selected_tracers)
    tr <- input$final_selected_tracers; gr_col <- input$group_col; m_lab <- input$mix_label; groups <- values$groups
    has_groups <- (!is.null(input$mix_cov) && input$mix_cov != "None")
    p_param <- if(input$opt_h) "P_ind" else if (has_groups) "P_level" else "P"
    
    samples <- rstan::extract(values$fit)
    p_samples <- samples[[p_param]]
    if(input$opt_h && length(dim(p_samples)) == 2) dim(p_samples) <- c(dim(p_samples)[1], 1, dim(p_samples)[2])
    
    df_mix_rows <- values$df[values$df[[gr_col]] == m_lab, ]
    N_mix <- nrow(df_mix_rows)
    
    if (input$opt_h) medians <- apply(p_samples, c(2,3), median)
    else if (has_groups) medians <- apply(p_samples, c(2,3), median)[as.numeric(as.factor(df_mix_rows[[input$mix_cov]])), ]
    else medians <- matrix(rep(apply(p_samples, 2, median), N_mix), ncol=length(groups), byrow=T)
    
    medians_corr <- medians / rowSums(medians); df_src <- values$df[values$df[[gr_col]] != m_lab, ]
    src_means_mat <- matrix(0, nrow = length(groups), ncol = length(tr)); src_range <- numeric(length(tr))
    
    ## COMMENTED - Error vs. Full Source Range (Min-Max)
    #for(i in seq_along(tr)) {
    #  t_vals <- df_src[[tr[i]]]; src_range[i] <- max(t_vals, na.rm=T) - min(t_vals, na.rm=T)
    #  for(g in seq_along(groups)) src_means_mat[g, i] <- mean(df_src[df_src[[gr_col]] == groups[g], tr[i]], na.rm=T)
    #}
    
    ## Error vs. Source Means with Variation based on Latorre et al. (2021)
    for(i in seq_along(tr)) {
      # 1. Mean and SD for the sources
      src_means <- tapply(df_src[[tr[i]]], df_src[[gr_col]], mean, na.rm=TRUE)
      src_sds <- tapply(df_src[[tr[i]]], df_src[[gr_col]], sd, na.rm=TRUE)
      src_sds[is.na(src_sds)] <- 0 # Защита от NA, если в каком-то источнике ровно 1 проба
      
      # 2. Latorre et al. (2021) for normalization factor
      src_range[i] <- max(src_means + src_sds, na.rm=TRUE) - min(src_means - src_sds, na.rm=TRUE)
      
      for(g in seq_along(groups)) {
        src_means_mat[g, i] <- mean(df_src[df_src[[gr_col]] == groups[g], tr[i]], na.rm=T)
      }
    }
    
    obs_conc <- as.matrix(df_mix_rows[, tr])
    pred_raw <- matrix(0, nrow = N_mix, ncol = length(tr))
    pred_corr <- matrix(0, nrow = N_mix, ncol = length(tr))
    
    has_beta <- input$opt_b && "beta_sorting" %in% names(samples)
    
    if (has_beta) {
      beta_med <- apply(samples$beta_sorting, 2, median)
      ssa_m_raw <- df_mix_rows[[input$size_col]]; ssa_s_raw <- df_src[[input$size_col]]
      ssa_m_raw[is.na(ssa_m_raw)] <- mean(ssa_m_raw, na.rm=TRUE); ssa_s_raw[is.na(ssa_s_raw)] <- mean(ssa_s_raw, na.rm=TRUE)
      mean_s <- mean(ssa_s_raw); sd_s <- sd(ssa_s_raw) + 1e-9
      ssa_m_z <- (ssa_m_raw - mean_s) / sd_s; ssa_s_z <- (ssa_s_raw - mean_s) / sd_s
      mean_ssa_src_z <- mean(ssa_s_z) 
    }
    
    for(i in 1:N_mix) {
      for(t_idx in 1:length(tr)) {
        base_val <- sum(medians_corr[i, ] * src_means_mat[, t_idx])
        pred_raw[i, t_idx] <- base_val
        if (has_beta) pred_corr[i, t_idx] <- base_val * exp(beta_med[t_idx] * (ssa_m_z[i] - mean_ssa_src_z))
        else          pred_corr[i, t_idx] <- base_val
      }
    }
    
    res_list <- list()
    for(t_idx in seq_along(tr)) {
      o <- obs_conc[, t_idx]; p_r <- pred_raw[, t_idx]; p_c <- pred_corr[, t_idx]
      
      res_list[[t_idx]] <- data.frame(
        Sample = row.names(obs_conc), Tracer = tr[t_idx], Observed = o, 
        Pred_Raw = p_r, Pred_Corr = p_c, 
        Bias_Raw_pct = ((p_r - o) / (abs(o) + 1e-9)) * 100,
        Bias_Corr_pct = ((p_c - o) / (abs(o) + 1e-9)) * 100,
        sMAPE_Raw = (abs(p_r - o) / ((abs(p_r) + abs(o) + 1e-9) / 2)) * 100,
        sMAPE_Corr = (abs(p_c - o) / ((abs(p_c) + abs(o) + 1e-9) / 2)) * 100,
        NRMSE_src_Raw = (abs(p_r - o) / (src_range[t_idx] + 1e-9)) * 100,
        NRMSE_src_Corr = (abs(p_c - o) / (src_range[t_idx] + 1e-9)) * 100,
        Beta_Shift_pct = if(has_beta) ((p_c - p_r) / (abs(p_r) + 1e-9)) * 100 else 0
      )
    }
    df_final <- do.call(rbind, res_list)
    
    clamp_cols <- c("Bias_Raw_pct", "Bias_Corr_pct", "Beta_Shift_pct")
    for(c in clamp_cols) { df_final[[c]][df_final[[c]] > 500] <- 500; df_final[[c]][df_final[[c]] < -500] <- -500 }
    return(df_final)
  })
  
  # --- DYNAMIC DATA CHOISE FOR THE GRAPHS ---
  active_recon_data <- reactive({
    df <- recon_calc()
    mode <- input$recon_mode
    if(is.null(mode)) mode <- "corr" 
    
    res <- data.frame(Sample = df$Sample, Tracer = df$Tracer, Observed = df$Observed)
    
    if (mode == "corr") {
      res$Predicted <- df$Pred_Corr; res$Bias_pct <- df$Bias_Corr_pct
      res$sMAPE <- df$sMAPE_Corr; res$NRMSE_src <- df$NRMSE_src_Corr
    } else {
      res$Predicted <- df$Pred_Raw; res$Bias_pct <- df$Bias_Raw_pct
      res$sMAPE <- df$sMAPE_Raw; res$NRMSE_src <- df$NRMSE_src_Raw
    }
    return(res)
  })
  
  # --- GRAPHS and TABLE ---
  output$recon_error_boxplot <- renderPlot({ 
    req(active_recon_data())
    ggplot(active_recon_data(), aes(x = Tracer, y = Bias_pct, fill = Tracer)) + 
      geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1) + 
      geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.fill = "white") + 
      theme_minimal(base_size = 16) + scale_fill_viridis_d(option = "mako", guide = "none") + 
      labs(title = "Reconstruction Bias per Tracer", y = "Bias (%)", x = "") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x = element_blank()) + 
      coord_cartesian(ylim = c(-150, 150)) 
  })
  
  output$recon_metrics_table <- renderTable({ 
    active_recon_data() %>% group_by(Tracer) %>% 
      summarise(`Mean Bias (%)` = mean(Bias_pct), `sMAPE (%)` = mean(sMAPE), `Error vs Src Range (%)` = mean(NRMSE_src)) 
  }, digits = 1)
  
  output$recon_plot <- renderPlot({ 
    ggplot(active_recon_data(), aes(x = Observed, y = Predicted, color = Tracer)) + 
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") + 
      geom_point(size = 3, alpha = 0.7) + facet_wrap(~Tracer, scales = "free") + 
      theme_bw(base_size = 14) + 
      labs(title = "Geochemical Reconstruction", x = "Observed Concentration", y = "Calculated Concentration") 
  })
  
  output$download_recon_csv <- downloadHandler(
    filename = function() { paste0("Reconstruction_Details_", Sys.Date(), ".csv") }, 
    content = function(file) { write.csv(recon_calc(), file, row.names = FALSE) }
  )
  
  # --- VALIDATION BLOCK ---
  observeEvent(input$run_val, {
    req(values$df, input$final_selected_tracers)
    withProgress(message = 'Validation (Virtual Mixtures)...', value = 0, {
      
      tr <- input$final_selected_tracers; gr_col <- input$group_col; m_lab <- input$mix_label
      df_src <- values$df[values$df[[gr_col]] != m_lab, ]
      
      groups <- unique(df_src[[gr_col]]); n_groups <- length(groups)
      
      src_params <- list()
      for(g in groups) {
        sub_data <- as.matrix(df_src[df_src[[gr_col]] == g, tr]); sub_data[is.na(sub_data)] <- 0
        cov_mat <- tryCatch({ cm <- cov(sub_data); diag(cm) <- diag(cm) + 1e-6; cm }, error = function(e) diag(apply(sub_data, 2, var)) + 1e-6)
        src_params[[g]] <- list(mu = colMeans(sub_data), sigma = cov_mat)
      }
      
      is_simplex <- input$val_gen_mode %in% c("det_simp", "stoch_simp")
      is_stoch_tracers <- input$val_gen_mode %in% c("stoch_rand", "stoch_simp")
      
      # А. Apportionments generator
      if (is_simplex) {
        step <- input$val_step; if (is.null(step) || step <= 0) step <- 0.05
        seq_vals <- seq(0, 1, by = step); grid <- expand.grid(rep(list(seq_vals), n_groups - 1)); grid$last <- 1 - rowSums(grid)
        grid <- grid[grid$last >= -1e-9 & grid$last <= 1 + 1e-9, ]; grid[grid < 0] <- 0; grid <- t(apply(grid, 1, function(x) x / sum(x)))
        N_val <- nrow(grid); true_props <- as.matrix(grid); colnames(true_props) <- groups
        if (N_val > 2000) { true_props <- true_props[sample(1:N_val, 2000), ]; N_val <- 2000 }
      } else {
        N_val <- input$val_n; true_props <- matrix(0, nrow = N_val, ncol = n_groups); colnames(true_props) <- groups
        for(i in 1:N_val) { p_rand <- runif(n_groups); true_props[i, ] <- p_rand / sum(p_rand) }
      }
      
      # Б. Tracers generator
      y_val_raw <- matrix(0, nrow = N_val, ncol = length(tr)); colnames(y_val_raw) <- tr
      for(i in 1:N_val) {
        mix_vec <- numeric(length(tr)); p_current <- true_props[i, ]
        for(k in 1:n_groups) {
          if(is_stoch_tracers) {
            sim <- tryCatch(MASS::mvrnorm(1, src_params[[groups[k]]]$mu, src_params[[groups[k]]]$sigma), 
                            error = function(e) rnorm(length(tr), src_params[[groups[k]]]$mu, sqrt(diag(src_params[[groups[k]]]$sigma))))
          } else {
            sim <- src_params[[groups[k]]]$mu
          }
          mix_vec <- mix_vec + (sim * p_current[k])
        }
        noise <- as.numeric(input$val_noise)
        if(noise > 0) { mix_vec <- mix_vec * runif(length(mix_vec), 1 - noise, 1 + noise); mix_vec <- mix_vec + rnorm(length(mix_vec), 0, 0.001) }
        mix_vec[mix_vec <= 1e-6] <- 1e-6; y_val_raw[i, ] <- mix_vec
      }
      
      for(j in 1:ncol(y_val_raw)) if(sd(y_val_raw[,j]) < 1e-9) y_val_raw[,j] <- y_val_raw[,j] + rnorm(N_val, 0, 0.01)
      
      export_df <- data.frame(ID = paste0("Virtual_Mix_", 1:N_val)); export_df <- cbind(export_df, as.data.frame(y_val_raw))
      props_df <- as.data.frame(true_props); colnames(props_df) <- paste0("TRUE_", colnames(props_df))
      values$val_export <- cbind(export_df, props_df)
      
      y_src_real <- as.matrix(df_src[, tr]); y_src_real[y_src_real <= 0] <- 1e-6
      
      if(!input$opt_coda) {
        y_src_s <- scale(y_src_real)
        y_val_s <- scale(y_val_raw, center=attr(y_src_s, "scaled:center"), scale=attr(y_src_s, "scaled:scale"))
        V_val <- matrix(0, nrow=0, ncol=0)
        N_tr_raw_val <- ncol(y_val_s)
      } else { 
        ys_safe <- y_src_real
        yv_safe <- y_val_raw; yv_safe[yv_safe <= 0] <- 1e-6
        
        ys_prop <- ys_safe / rowSums(ys_safe)
        yv_prop <- yv_safe / rowSums(yv_safe)
        
        N_tr_raw_val <- ncol(ys_prop)
        V_val <- compositions::ilrBase(D = N_tr_raw_val)
        
        y_src_s <- as.matrix(log(ys_prop) %*% V_val)
        y_val_s <- as.matrix(log(yv_prop) %*% V_val) 
      }
      
      ssa_val <- rep(0, N_val); ssa_src <- rep(0, nrow(df_src))
      if(!is.null(input$size_col) && input$size_col != "None" && input$size_col %in% colnames(df_src)) {
        ssa_src_raw <- df_src[[input$size_col]]
        ssa_src_raw[is.na(ssa_src_raw)] <- mean(ssa_src_raw, na.rm = TRUE)
        if(!all(is.na(ssa_src_raw))) ssa_src <- (ssa_src_raw - mean(ssa_src_raw)) / (sd(ssa_src_raw) + 1e-9)
      }
      
      has_groups <- (!is.null(input$mix_cov) && input$mix_cov != "None")
      stan_val_data <- list(
        N_mix = N_val, N_sources = nrow(y_src_s), 
        N_tracers = ncol(y_src_s), N_tracers_raw = N_tr_raw_val,
        N_groups = n_groups,
        y_mix = y_val_s, y_source = y_src_s, 
        V = V_val,
        group_id = as.numeric(as.factor(df_src[[gr_col]])),
        ssa_mix = as.array(ssa_val), ssa_source = as.array(ssa_src), conn_weights = rep(1, n_groups)
      )
      if(has_groups) { stan_val_data$N_mix_levels <- 1; stan_val_data$mix_level_id <- as.array(rep(1, N_val)) }
      
      stan_val_data$holdout_id <- 0
      
      has_groups_val <- (!is.null(input$mix_cov) && input$mix_cov != "None")
      model_id <- paste0(as.integer(input$opt_h), as.integer(input$opt_b), 
                         input$cov_mode, as.integer(input$opt_coda), as.integer(has_groups_val), as.integer(input$opt_bias))
      
      if(!is.null(values$compiled_list[[model_id]])) {
        mod <- values$compiled_list[[model_id]]
      } else {
        showNotification("Compiling model for validation (first time)...", type="message")
        code <- generate_stan_code(input$opt_h, input$opt_b, input$cov_mode, input$opt_coda, has_groups_val, bias = input$opt_bias)
        values$compiled_list[[model_id]] <- stan_model(model_code = code)
        mod <- values$compiled_list[[model_id]]
      }
      
      fit_val <- tryCatch({
        if (input$calc_method == "vb") vb(mod, data = stan_val_data, output_samples = 1000, seed=123, init="random", tol_rel_obj=0.005)
        else sampling(mod, data = stan_val_data, iter = input$iter, chains = input$chains, init = 0.5)
      }, error = function(e) { showNotification(paste("Validation error:", e$message), type="error"); return(NULL) })
      if(is.null(fit_val)) return(NULL)
      
      val_res <- rstan::extract(fit_val)
      pred_means <- NULL; pred_q025 <- NULL; pred_q975 <- NULL; pred_q25 <- NULL; pred_q75 <- NULL
      p_samples_array <- NULL # Переменная для хранения 3D массива цепей
      
      if ("P_ind" %in% names(val_res)) {
        raw <- val_res$P_ind
        if (length(dim(raw)) == 2) { if (N_val == 1) dim(raw) <- c(dim(raw)[1], 1, dim(raw)[2]) else dim(raw) <- c(dim(raw)[1], N_val, n_groups) }
        p_samples_array <- raw # Формат: [iter, N_val, n_groups]
        pred_means <- apply(raw, c(2,3), median); pred_q025 <- apply(raw, c(2,3), quantile, probs = 0.025); pred_q975 <- apply(raw, c(2,3), quantile, probs = 0.975); pred_q25 <- apply(raw, c(2,3), quantile, probs = 0.25); pred_q75 <- apply(raw, c(2,3), quantile, probs = 0.75)
      } else if ("P_level" %in% names(val_res)) {
        raw <- val_res$P_level; p_samples <- if (length(dim(raw)) == 3) raw[, 1, ] else raw; if (is.vector(p_samples)) p_samples <- matrix(p_samples, ncol = n_groups)
        p_samples_array <- array(0, dim = c(nrow(p_samples), N_val, n_groups))
        for(v in 1:N_val) p_samples_array[, v, ] <- p_samples # Транслируем на N_val
        p_mean <- apply(p_samples, 2, median); pred_means <- matrix(rep(p_mean, N_val), nrow = N_val, byrow = TRUE); pred_q025 <- matrix(rep(apply(p_samples, 2, quantile, probs=0.025), N_val), nrow=N_val, byrow=TRUE); pred_q975 <- matrix(rep(apply(p_samples, 2, quantile, probs=0.975), N_val), nrow=N_val, byrow=TRUE); pred_q25 <- matrix(rep(apply(p_samples, 2, quantile, probs=0.25), N_val), nrow=N_val, byrow=TRUE); pred_q75 <- matrix(rep(apply(p_samples, 2, quantile, probs=0.75), N_val), nrow=N_val, byrow=TRUE)
      } else {
        raw <- val_res$P; if (is.vector(raw)) raw <- matrix(raw, ncol = n_groups)
        p_samples_array <- array(0, dim = c(nrow(raw), N_val, n_groups))
        for(v in 1:N_val) p_samples_array[, v, ] <- raw # Транслируем на N_val
        p_mean <- apply(raw, 2, median); pred_means <- matrix(rep(p_mean, N_val), nrow = N_val, byrow = TRUE); pred_q025 <- matrix(rep(apply(raw, 2, quantile, probs=0.025), N_val), nrow=N_val, byrow=TRUE); pred_q975 <- matrix(rep(apply(raw, 2, quantile, probs=0.975), N_val), nrow=N_val, byrow=TRUE); pred_q25 <- matrix(rep(apply(raw, 2, quantile, probs=0.25), N_val), nrow=N_val, byrow=TRUE); pred_q75 <- matrix(rep(apply(raw, 2, quantile, probs=0.75), N_val), nrow=N_val, byrow=TRUE)
      }
      
      if (is.vector(pred_means)) { pred_means <- t(as.matrix(pred_means)); pred_q025 <- t(as.matrix(pred_q025)); pred_q975 <- t(as.matrix(pred_q975)); pred_q25 <- t(as.matrix(pred_q25)); pred_q75 <- t(as.matrix(pred_q75)) }
      
      colnames(pred_means) <- groups
      values$val_data <- list(true = true_props, pred = pred_means, q025 = pred_q025, q975 = pred_q975, q25 = pred_q25, q75 = pred_q75, samples = p_samples_array)
      
      pred_df_export <- as.data.frame(pred_means); colnames(pred_df_export) <- paste0("PRED_", groups)
      values$val_export_full <- cbind(values$val_export, pred_df_export)
      showNotification("Validation complete! Metrics updated.", type="message")
    })
  })
  
  output$download_val_data <- downloadHandler(filename = function() { paste0("Validation_Results_", Sys.Date(), ".csv") }, content = function(file) { write.csv(if(!is.null(values$val_export_full)) values$val_export_full else values$val_export, file, row.names = FALSE) })
  
  output$val_plot <- renderPlot({
    req(values$val_data)
    true_df <- reshape2::melt(values$val_data$true); pred_df <- reshape2::melt(values$val_data$pred)
    plot_df <- data.frame(Source = true_df$Var2, True = true_df$value, Predicted = pred_df$value)
    ggplot(plot_df, aes(x = True, y = Predicted, color = Source)) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey", linewidth = 1) +
      geom_point(size = 3, alpha = 0.7) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5, linetype = "dotted") + coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) + theme_bw(base_size = 18) + labs(title = "Model Accuracy (Virtual Mixtures)", x = "True Proportion (Known)", y = "Predicted Proportion (Model)")
  })
  
  output$val_metrics <- renderTable({
    req(values$val_data)
    t <- as.matrix(values$val_data$true)
    p <- as.matrix(values$val_data$pred)
    q025 <- as.matrix(values$val_data$q025)
    q975 <- as.matrix(values$val_data$q975)
    q25 <- as.matrix(values$val_data$q25)
    q75 <- as.matrix(values$val_data$q75)
    
    mae <- colMeans(abs(t - p))
    rmse <- sqrt(colMeans((t - p)^2))
    me <- colMeans(p - t); naive_pred <- 1 / ncol(t)
    nse <- 1 - (colSums((t - p)^2) / colSums((t - naive_pred)^2))
    w95 <- colMeans(q975 - q025)
    w50 <- colMeans(q75 - q25)
    eps <- 0.015
    p95_adj <- colMeans((t >= (q025 - eps)) & (t <= (q975 + eps)))
    p50_adj <- colMeans((t >= (q25 - eps))  & (t <= (q75 + eps)))
    true_dom <- t > 0.50
    pred_dom <- p > 0.50
    hits <- colSums(true_dom & pred_dom)
    misses <- colSums(true_dom & !pred_dom)
    false_alarms <- colSums(!true_dom & pred_dom)
    denom_csi <- hits + misses + false_alarms
    csi <- ifelse(denom_csi == 0, NA, hits / denom_csi)
    denom_hr <- hits + misses
    hr <- ifelse(denom_hr == 0, NA, hits / denom_hr)
    
    r2 <- sapply(1:ncol(t), function(i) cor(t[, i], p[, i])^2)
    
    # --- Exact CRPS via raw MCMC ---
    crps_vals <- sapply(1:ncol(t), function(i) {
      src_samples <- values$val_data$samples[, , i] 
      if(is.null(dim(src_samples))) { src_samples <- matrix(src_samples, ncol = 1) }
      crps_vec <- scoringRules::crps_sample(y = t[, i], dat = t(src_samples))
      
      mean(crps_vec)
    })
    
    df_metrics <- data.frame(Source = colnames(t), MAE = mae, RMSE = rmse, ME = me, R2 = r2, CRPS = crps_vals, NSE = nse, `W50` = w50, `W95` = w95, `P50 (adj)` = p50_adj, `P95 (adj)` = p95_adj, `CSI (>50%)` = csi, `HR (>50%)` = hr, check.names = FALSE)
    mean_row <- data.frame(Source = "AVERAGE", MAE = mean(mae), RMSE = mean(rmse), ME = mean(abs(me)), R2 = mean(r2), CRPS = mean(crps_vals), NSE = mean(nse), `W50` = mean(w50), `W95` = mean(w95), `P50 (adj)` = mean(p50_adj), `P95 (adj)` = mean(p95_adj), `CSI (>50%)` = mean(csi, na.rm=TRUE), `HR (>50%)` = mean(hr, na.rm=TRUE), check.names = FALSE)
    rbind(df_metrics, mean_row)
  }, digits = 3, na = "-")
  
  # ==============================================================================
  # LOO-CV MODEL COMPARISON & DIAGNOSTICS LOGIC
  # ==============================================================================
  
  # Container for the loo-objects
  values$loo_objects <- list()
  
  # selector fir comparison
  output$ui_multi_compare_selector <- renderUI({
    req(length(values$saved_models) > 0)
    choices <- names(values$saved_models)
    selectizeInput("selected_models_loo", "Select Models to Compare:", 
                   choices = choices, 
                   selected = choices, # По умолчанию выбираем все
                   multiple = TRUE,
                   options = list(plugins = list('remove_button')))
  })
  
  # Run comparison
  observeEvent(input$run_loo_compare, {
    req(length(input$selected_models_loo) > 0)
    withProgress(message = 'Calculating LOO-CV for selected models...', value = 0, {
      n_models <- length(input$selected_models_loo)
      tryCatch({
        for (i in seq_along(input$selected_models_loo)) {
          m_name <- input$selected_models_loo[i]
          incProgress(1/n_models, detail = paste("Evaluating", m_name))
          if (is.null(values$loo_objects[[m_name]])) {
            fit <- values$saved_models[[m_name]]
            log_lik <- loo::extract_log_lik(fit, parameter_name = "log_lik", merge_chains = FALSE)
            values$loo_objects[[m_name]] <- loo::loo(log_lik, cores = 1)
          }
        }
      }, error = function(e) { showNotification(paste("Error calculating LOO:", e$message), type="error") })
    })
  })
  
  # --- DYNAMIC RESULTS OUTPUT ---
  output$loo_results_text <- renderPrint({
    # 1. PSIS-LOO output
    loo_list <- values$loo_objects[names(values$loo_objects) %in% input$selected_models_loo]
    
    if (length(loo_list) > 0) {
      cat("=== LOO-CV Leaderboard (Approximate PSIS-LOO) ===\n\n")
      if (length(loo_list) == 1) {
        print(loo_list[[1]])
      } else {
        n_obs <- nrow(loo_list[[1]]$pointwise)
        if (n_obs == 1) {
          cat("⚠️ WARNING: The dataset contains only one mixture (N_mix = 1).\n")
          scores <- sapply(loo_list, function(x) x$estimates["elpd_loo", "Estimate"])
          scores_sorted <- sort(scores, decreasing = TRUE)
          for (m in names(scores_sorted)) { cat(sprintf("%s : %.2f\n", m, scores_sorted[m])) }
        } else {
          comp <- loo::loo_compare(loo_list)
          cat("Negative values in 'elpd_diff' indicate worse predictive performance.\n\n")
          print(comp)
        }
      }
    } else {
      cat("Select models and click 'Calculate & Compare' to see PSIS-LOO leaderboard.\n")
    }
    
    # 2. EXACT LOO-CV Output
    if (length(values$exact_loo_results) > 0) {
      cat("\n=======================================================\n")
      cat("=== EXACT LOO-CV LEADERBOARD                        ===\n")
      cat("=======================================================\n")
      for (m_name in names(values$exact_loo_results)) {
        ex_res <- values$exact_loo_results[[m_name]]
        cat(sprintf("Model: %s\n", m_name))
        cat(sprintf("  elpd_loo:    %8.2f\n", ex_res$elpd_loo))
        cat(sprintf("  se_elpd_loo: %8.2f\n", ex_res$se_elpd_loo))
        cat("-------------------------------------------------------\n")
      }
      cat("Note: This is the mathematically exact elpd without approximations.\n")
    }
  })
  
  # --- PARETO k DIAGNOSTICS ---
  output$ui_pareto_selector <- renderUI({
    req(length(values$loo_objects) > 0)
    selectInput("pareto_model_select", "Select Model to inspect:", 
                choices = names(values$loo_objects))
  })
  
  output$pareto_plot <- renderPlot({
    req(input$pareto_model_select, values$loo_objects[[input$pareto_model_select]])
    
    loo_obj <- values$loo_objects[[input$pareto_model_select]]
    
    k_values <- loo_obj$diagnostics$pareto_k
    
    plot_df <- data.frame(
      Observation = seq_along(k_values),
      k = k_values,
      is_high = k_values > 0.7
    )
    
    ggplot(plot_df, aes(x = Observation, y = k)) +
      geom_hline(yintercept = 0.7, linetype = "dashed", color = "#e74c3c", linewidth = 1) +
      geom_hline(yintercept = 1.0, linetype = "solid", color = "#c0392b", linewidth = 0.5, alpha = 0.5) +
      
      geom_point(aes(color = is_high, size = is_high), alpha = 0.8) +
      scale_color_manual(values = c("FALSE" = "#2c3e50", "TRUE" = "#e74c3c")) +
      scale_size_manual(values = c("FALSE" = 3, "TRUE" = 5)) +
      
      geom_text(data = subset(plot_df, k > 0.7), 
                aes(label = Observation), 
                vjust = -1.5, color = "#c0392b", fontface = "bold", size = 5) +
      
      # Appearance
      theme_bw(base_size = 16) +
      theme(legend.position = "none",
            panel.grid.minor = element_blank(),
            axis.line = element_line(linewidth = 1)) +
      labs(title = paste("Pareto k Diagnostics:", input$pareto_model_select),
           subtitle = "Red points (k > 0.7) indicate highly influential mixtures",
           x = "Mixture Index (Data point)",
           y = "Pareto shape k") +
      coord_cartesian(ylim = c(min(0, min(k_values)), max(k_values) + 0.1)) # Динамический масштаб Y
  })
  
  # =========================================================================
  # EXACT LOO-CV MODULE
  # =========================================================================
  observeEvent(input$run_exact_loo, {
    req(values$stan_data, values$compiled_model)
    n_mixtures <- values$stan_data$N_mix
    exact_elpd_i <- numeric(n_mixtures)
    
    log_mean_exp <- function(x) {
      max_x <- max(x)
      max_x + log(mean(exp(x - max_x)))
    }
    
    withProgress(message = 'Running Exact LOO-CV...', value = 0, {
      for (i in 1:n_mixtures) {
        incProgress(0, detail = paste("Refitting target", i, "of", n_mixtures))
        
        cv_data <- values$stan_data
        
        cv_data$holdout_id <- i 
        
        tryCatch({
          fit_cv <- rstan::sampling(
            values$compiled_model, data = cv_data, 
            iter = input$iter, chains = input$chains, cores = input$chains,
            refresh = 0, show_messages = FALSE
          )
          
          all_log_liks <- rstan::extract(fit_cv, pars = "log_lik")$log_lik
          holdout_log_liks <- all_log_liks[, i]
          
          exact_elpd_i[i] <- log_mean_exp(holdout_log_liks)
        }, error = function(e) {
          warning(paste("Exact LOO failed for mixture", i, ":", e$message))
          exact_elpd_i[i] <- NA 
        })
        incProgress(1/n_mixtures)
      }
    })
    
    total_exact_elpd <- sum(exact_elpd_i, na.rm = TRUE)
    se_exact_elpd <- sqrt(n_mixtures * var(exact_elpd_i, na.rm = TRUE))
    
    values$exact_loo_results[[values$last_model_name]] <- list(
      elpd_loo = total_exact_elpd,
      se_elpd_loo = se_exact_elpd,
      pointwise = exact_elpd_i
    )
    showNotification(sprintf("Exact LOO-CV complete. ELPD: %.2f", total_exact_elpd), type = "message")
  })
  # =========================================================================
}

shinyApp(ui, server)