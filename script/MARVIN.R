##############################################
# Stomatal Conductance Model Fitting Script
# Multi-species version (looped by species × data.type)
##############################################

require(likelihood)

# === MODEL DEFINITIONS ===
Exponential <- function(A, B, psi) {
  A * exp(-B * psi)
}
define_parsE <- function(input_df) {
  list(
    parsE = list(A = max(input_df$gs), B = 1, sd = 2),
    par_loE = list(A = 0, B = 0.1, sd = 0.0005),
    par_highE = list(A = max(input_df$gs) * 2, B = 10, sd = 50)
  )
}

Logistic <- function(A, B, Xo, psi) {
  A / (1 + ((psi / Xo)^B))
}
define_parsL <- function(input_df) {
  list(
    parsL = list(A = 100, B = 2, Xo = 1, sd = 2),
    par_loL = list(A = 0, B = 0.1, Xo = 0, sd = 0.0005),
    par_highL = list(A = max(input_df$gs) * 2, B = 4, Xo = 5, sd = 50)
  )
}

Sigmoidal <- function(A, B, Xo, psi) {
  A / (1 + exp(-((psi - Xo) / B)))
}

define_parsS <- function(input_df) {
  # --- General case ---
  parsS <- list(
    A = max(input_df$gs),
    B = -0.5,
    Xo = mean(range(input_df$psi, na.rm = TRUE)),
    sd = 2
  )
  par_loS <- list(
    A = 0,
    B = -1.25,
    Xo = min(input_df$psi, na.rm = TRUE) - 0.1,
    sd = 0.0005
  )
  par_highS <- list(
    A = max(input_df$gs) * 1.5,
    B = 0,
    Xo = max(input_df$psi, na.rm = TRUE) + 0.1,
    sd = 50
  )

  # --- Special case for Trae (commented out, ready to activate) ---
  if (tolower(input_df$species[1]) == "trae") {
    parsS <- list(
      A = max(input_df$gs) * 0.9,
      B = -0.4,
      Xo = mean(range(input_df$psi, na.rm = TRUE)),
      sd = 2
    )
    par_loS <- list(
      A = 0.7 * max(input_df$gs),
      B = -0.9,
      Xo = min(input_df$psi, na.rm = TRUE) - 0.1,
      sd = 0.0005
    )
    par_highS <- list(
      A = max(input_df$gs) * 1.1,
      B = -0.1,
      Xo = max(input_df$psi, na.rm = TRUE) + 0.1,
      sd = 50
    )
  }

  #Return correctly named list
  list(parsS = parsS, par_loS = par_loS, par_highS = par_highS)
}

Exponential2 <- function(A, B, C, psi) {
  C + A * exp(-B * psi)
}
define_parsE2 <- function(input_df) {
  list(
    parsE2 = list(A = max(input_df$gs), B = 0.5, C = 0.1, sd = 2),
    par_loE2 = list(A = 0, B = 0.011, C = 0.001, sd = 0.0005),
    par_highE2 = list(A = max(input_df$gs) * 2, B = 20, C = 20, sd = 50)
  )
}

Linear <- function(A, B, psi) {
  B * psi + A
}

# === FITTING FUNCTIONS ===
do_the_thing_nonlinear <- function(
  input_df,
  model_type,
  pars1,
  par_lo1,
  par_hi1,
  model_name,
  output_plot_dir
) {
  # --- Special handling for Sigmoidal ---
  if (tolower(model_name) == "sigmoidal") {
    sig_params <- define_parsS(input_df) # uses the current definition without Trae special case
    pars1 <- sig_params$parsS
    par_lo1 <- sig_params$par_loS
    par_hi1 <- sig_params$par_highS
  }

  # --- Old special case for Trae (commented out) ---
  # if (tolower(model_name) == "sigmoidal" && tolower(input_df$species[1]) == "trae") {
  #   sig_params <- define_parsS(input_df)  # would use the Trae-specific bounds
  #   pars1   <- sig_params$parsS
  #   par_lo1 <- sig_params$par_loS
  #   par_hi1 <- sig_params$par_highS
  # }

  var <- list(psi = "psi", x = "gs", mean = "predicted", log = TRUE)
  res <- anneal(
    model = model_type,
    par = pars1,
    source_data = input_df,
    var = var,
    par_lo = par_lo1,
    par_hi = par_hi1,
    dep_var = "gs",
    pdf = dnorm,
    max_iter = 50000,
    show_display = FALSE,
    temp_red = 0.05,
    hessian = T
  )

  sterror <- res$std_errs
  A <- as.numeric(res$best_pars[1])
  B <- as.numeric(res$best_pars[2])
  C <- as.numeric(ifelse(length(res$best_pars) >= 3, res$best_pars[3], NA))
  D <- as.numeric(ifelse(length(res$best_pars) >= 4, res$best_pars[4], NA))

  # vcov_mat <- tryCatch(
  #   solve(res$hessian),
  #   error = function(e) NA
  # )

  vcov_mat <- res$var_covar_mat

  # ---- New calculations ----
  model_lower <- tolower(model_name)
  if (model_lower == "logistic") {
    Kmax <- A
  } else if (model_lower == "exponential") {
    Kmax <- A
  } else if (model_lower == "sigmoidal") {
    Kmax <- A / (1 + exp(C / B))
  } else if (model_lower == "exponential2") {
    Kmax <- C + A
  } else if (model_lower == "linear") {
    Kmax <- A
  } else {
    Kmax <- NA
  }

  if (model_lower == "logistic") {
    Kmax_0.1 <- A / (1 + (0.1 / C)^B)
  } else if (model_lower == "exponential") {
    Kmax_0.1 <- A * exp(-(B * 0.1))
  } else if (model_lower == "sigmoidal") {
    Kmax_0.1 <- A / (1 + exp(-(-0.1 - C) / B))
  } else if (model_lower == "exponential2") {
    Kmax_0.1 <- C + A * exp(-B * 0.1)
  } else if (model_lower == "linear") {
    Kmax_0.1 <- A + B * 0.1
  } else {
    Kmax_0.1 <- NA
  }

  if (model_lower == "logistic") {
    P50_0.1kmax <- C * (A / (Kmax_0.1 * 0.5) - 1)^(1 / B)
  } else if (model_lower == "exponential") {
    P50_0.1kmax <- (-1 / B) * log((Kmax_0.1 * 0.5) / A)
  } else if (model_lower == "sigmoidal") {
    P50_0.1kmax <- -B * log(A / (Kmax_0.1 * 0.5) - 1) + C
  } else if (model_lower == "exponential2") {
    P50_0.1kmax <- (-1 / B) * log(((Kmax_0.1 * 0.5) - C) / A)
  } else if (model_lower == "linear") {
    P50_0.1kmax <- (Kmax_0.1 / 2 - A) / B
  } else {
    P50_0.1kmax <- NA
  }

  if (model_lower == "linear") {
    psi_kleaf20 <- ((0.8 - 1) * (-A) / B) * -1
    psi_kleaf50 <- ((0.5 - 1) * (-A) / B) * -1
    psi_kleaf80 <- ((0.2 - 1) * (-A) / B) * -1
    psi_kleaf95 <- ((0.05 - 1) * (-A) / B) * -1
    #RWC_PLRC50 <- (50-A)/B
    RWC_fvfm50 <- (50 - A) / B
  } else if (model_lower == "logistic") {
    psi_kleaf20 <- C * (1 / 0.8 - 1)^(1 / B)
    psi_kleaf50 <- C * (1 / 0.5 - 1)^(1 / B)
    psi_kleaf80 <- C * (1 / 0.2 - 1)^(1 / B)
    psi_kleaf95 <- C * (1 / 0.05 - 1)^(1 / B)
    # RWC_PLRC50 <- ((A/50)-1)^(1/B)*C
    RWC_fvfm50 <- ((A / 50) - 1)^(1 / B) * C
  } else if (model_lower == "sigmoidal") {
    psi_kleaf20 <- -B * log(A / (0.8 * A) - 1) + C
    psi_kleaf50 <- -B * log(A / (0.5 * A) - 1) + C
    psi_kleaf80 <- -B * log(A / (0.2 * A) - 1) + C
    psi_kleaf95 <- -B * log(A / (0.05 * A) - 1) + C
    # RWC_PLRC50 <- -B*log((A/50)-1)+C
    RWC_fvfm50 <- -B * log((A / 50) - 1) + C
  } else if (model_lower == "exponential") {
    psi_kleaf20 <- log(0.8) / (-B)
    psi_kleaf50 <- log(0.5) / (-B)
    psi_kleaf80 <- log(0.2) / (-B)
    psi_kleaf95 <- log(0.05) / (-B)
    #RWC_PLRC50 <- (-1/B)* log((50)/A)
    RWC_fvfm50 <- (-1 / B) * log((50) / A)
  } else if (model_lower == "exponential2") {
    psi_kleaf20 <- log(((C * 0.8) + (0.8 * A) - C) / A) / (-B)
    psi_kleaf50 <- log(((C * 0.5) + (0.5 * A) - C) / A) / (-B)
    psi_kleaf80 <- log(((C * 0.2) + (0.2 * A) - C) / A) / (-B)
    psi_kleaf95 <- log(((C * 0.05) + (0.05 * A) - C) / A) / (-B)
    # RWC_PLRC50 <- (-1/B)* log((50-C)/A)
    RWC_fvfm50 <- (-1 / B) * log((50 - C) / A)
  } else {
    psi_kleaf20 <- psi_kleaf50 <- psi_kleaf80 <- psi_kleaf95 <- NA
  }
  # ---------------------------

  parvec <- c(
    Species = as.character(input_df$species[1]),
    data.type = as.character(input_df$data.type[1]),
    model = model_name,
    A = A,
    B = B,
    C = C,
    D = D,
    sterrorA = sterror[1],
    sterrorB = sterror[2],
    sterrorC = ifelse(length(sterror) >= 3, sterror[3], NA),
    sterrorD = ifelse(length(sterror) >= 4, sterror[4], NA),
    loglikeli = res$max_likeli,
    rsq = res$R2,
    slope = sum(res$source_data$predicted * res$source_data$gs) /
      sum(res$source_data$predicted^2),
    AIC = res$aic,
    AICcorr = res$aic_corr,
    N = length(res$source_data$gs),
    Kmax = Kmax,
    Kmax_at_0.1MPa = Kmax_0.1,
    P50_0.1kmax = P50_0.1kmax,
    psi_kleaf20 = psi_kleaf20,
    psi_kleaf50 = psi_kleaf50,
    psi_kleaf80 = psi_kleaf80,
    psi_kleaf95 = psi_kleaf95,
    # RWC_PLRC50 = RWC_PLRC50,
    RWC_fvfm50 = RWC_fvfm50
  )

  png(
    file.path(
      output_plot_dir,
      paste0(
        input_df$species[1],
        "_",
        input_df$data.type[1],
        "_",
        model_name,
        ".png"
      )
    ),
    width = 800,
    height = 600
  )
  plot(
    res$source_data$psi,
    res$source_data$gs,
    xlab = "Water Potential (-MPa)",
    ylab = "Stomatal Conductance (mmol m-2 s-1)",
    main = paste(input_df$species[1], input_df$data.type[1], model_name)
  )
  for_plotting <- cbind(res$source_data$psi, res$source_data$predicted)
  for_plotting <- for_plotting[order(for_plotting[, 1]), ]
  lines(for_plotting[, 1], for_plotting[, 2], col = "blue", lwd = 2)
  dev.off()

  return(list(
    Species = as.character(input_df$species[1]),
    data.type = as.character(input_df$data.type[1]),
    model = model_name,
    A = A,
    B = B,
    C= C,
    D=D,
    # C = ifelse(length(res$best_pars) >= 3, res$best_pars[3], NA),
    # D = ifelse(length(res$best_pars) >= 4, res$best_pars[4], NA),
    sterrorA = sterror[[1]],
    sterrorB = sterror[[2]],
    sterrorC = sterror[[3]],
    sterrorD = sterror[["sd"]],
    #sterrorC = ifelse(length(sterror) >= 3, sterror[3], NA),
    #sterrorD = ifelse(length(sterror) >= 4, sterror[4], NA),
    loglikeli = res$max_likeli,
    rsq = res$R2,
    slope = sum(res$source_data$predicted * res$source_data$gs) /
      sum(res$source_data$predicted^2),
    AIC = res$aic,
    AICcorr = res$aic_corr,
    N = length(res$source_data$gs),
    Kmax = Kmax,
    Kmax_at_0.1MPa = Kmax_0.1,
    P50_0.1kmax = P50_0.1kmax,
    psi_kleaf20 = psi_kleaf20,
    psi_kleaf50 = psi_kleaf50,
    psi_kleaf80 = psi_kleaf80,
    psi_kleaf95 = psi_kleaf95,
    # RWC_PLRC50 = RWC_PLRC50
    # RWC_fvfm50 = RWC_fvfm50
    #best_pars = res$best_pars,
    # std_errs = res$std_errs,
    # hessian = res$hessian,
    vcov = vcov_mat,
    # logLik = res$max_likeli,
    # AIC = res$aic,
    # AICc = res$aic_corr,
    # R2 = res$R2,
    # fitted = res$source_data$predicted,
    # data = res$source_data,
    model = model_name,
    species = input_df$species[1],
    data.type = input_df$data.type[1]
  ))
}


do_the_thing_linear <- function(
  input_df,
  model_type,
  model_name,
  output_plot_dir
) {
  var <- list(psi = "psi", x = "gs", mean = "predicted", log = TRUE)
  lm_fit <- lm(input_df$gs ~ input_df$psi)
  pars <- list(A = coef(lm_fit)[1], B = coef(lm_fit)[2], sd = 1)
  par_lo <- list(
    A = coef(lm_fit)[1] * 0.05,
    B = coef(lm_fit)[2] * 2,
    sd = 0.005
  )
  par_hi <- list(A = coef(lm_fit)[1] * 4, B = coef(lm_fit)[2] * 0.1, sd = 20)

  res <- anneal(
    model = model_type,
    par = pars,
    source_data = input_df,
    var = var,
    par_lo = par_lo,
    par_hi = par_hi,
    dep_var = "gs",
    pdf = dnorm,
    max_iter = 50000, # probably no need to specify this 50k is the default
    show_display = FALSE,
    temp_red = 0.05,
    hessian = T
  )

  sterror <- res$std_errs
  A <- as.numeric(res$best_pars[1])
  B <- as.numeric(res$best_pars[2])

  # vcov_mat <- tryCatch(
  #   solve(res$hessian),
  #   error = function(e) NA
  # )
  vcov_mat <- res$var_covar_mat

  # ---- New calculations ----
  Kmax_0.1 <- A + B * 0.1
  Kmax <- A
  P50_0.1kmax <- (Kmax_0.1 / 2 - A) / B
  psi_kleaf20 <- ((0.8 - 1) * (-A) / B) * -1
  psi_kleaf50 <- ((0.5 - 1) * (-A) / B) * -1
  psi_kleaf80 <- ((0.2 - 1) * (-A) / B) * -1
  psi_kleaf95 <- ((0.05 - 1) * (-A) / B) * -1
  # RWC_PLRC50 <- (50-A)/B
  RWC_fvfm50 <- (50 - A) / B
  # ---------------------------

  parvec <- c(
    Species = as.character(input_df$species[1]),
    data.type = as.character(input_df$data.type[1]),
    model = model_name,
    A = A,
    B = B,
    C = ifelse(length(res$best_pars) >= 3, res$best_pars[3], NA),
    D = ifelse(length(res$best_pars) >= 4, res$best_pars[4], NA),
    sterrorA = sterror[1],
    sterrorB = sterror[2],
    sterrorC = ifelse(length(sterror) >= 3, sterror[3], NA),
    sterrorD = ifelse(length(sterror) >= 4, sterror[4], NA),
    loglikeli = res$max_likeli,
    rsq = res$R2,
    slope = sum(res$source_data$predicted * res$source_data$gs) /
      sum(res$source_data$predicted^2),
    AIC = res$aic,
    AICcorr = res$aic_corr,
    N = length(res$source_data$gs),
    Kmax = Kmax,
    Kmax_at_0.1MPa = Kmax_0.1,
    P50_0.1kmax = P50_0.1kmax,
    psi_kleaf20 = psi_kleaf20,
    psi_kleaf50 = psi_kleaf50,
    psi_kleaf80 = psi_kleaf80,
    psi_kleaf95 = psi_kleaf95,
    # RWC_PLRC50 = RWC_PLRC50
    RWC_fvfm50 = RWC_fvfm50
  )

  png(
    file.path(
      output_plot_dir,
      paste0(
        input_df$species[1],
        "_",
        input_df$data.type[1],
        "_",
        model_name,
        ".png"
      )
    ),
    width = 800,
    height = 600
  )
  plot(
    res$source_data$psi,
    res$source_data$gs,
    xlab = "Water Potential (-MPa)",
    ylab = "Stomatal Conductance (mmol m-2 s-1)",
    main = paste(input_df$species[1], input_df$data.type[1], model_name)
  )
  for_plotting <- cbind(res$source_data$psi, res$source_data$predicted)
  for_plotting <- for_plotting[order(for_plotting[, 1]), ]
  lines(for_plotting[, 1], for_plotting[, 2], col = "blue", lwd = 2)
  dev.off()

  return(list(
    Species = as.character(input_df$species[1]),
    data.type = as.character(input_df$data.type[1]),
    model = model_name,
    A = A,
    B = B,
    C = C,
    D = D,
    #C = ifelse(length(res$best_pars) >= 3, res$best_pars[3], NA),
    #D = ifelse(length(res$best_pars) >= 4, res$best_pars[4], NA),
    sterrorA = sterror[[1]],
    sterrorB = sterror[[2]],
    sterrorC = sterror[[3]],
    sterrorD = sterror[["sd"]],
    #sterrorC = ifelse(length(sterror) >= 3, sterror[3], NA),
    #sterrorD = ifelse(length(sterror) >= 4, sterror[4], NA),
    loglikeli = res$max_likeli,
    rsq = res$R2,
    slope = sum(res$source_data$predicted * res$source_data$gs) /
      sum(res$source_data$predicted^2),
    AIC = res$aic,
    AICcorr = res$aic_corr,
    N = length(res$source_data$gs),
    Kmax = Kmax,
    Kmax_at_0.1MPa = Kmax_0.1,
    P50_0.1kmax = P50_0.1kmax,
    psi_kleaf20 = psi_kleaf20,
    psi_kleaf50 = psi_kleaf50,
    psi_kleaf80 = psi_kleaf80,
    psi_kleaf95 = psi_kleaf95,
    # RWC_PLRC50 = RWC_PLRC50
    # RWC_fvfm50 = RWC_fvfm50
    #best_pars = res$best_pars,
    # std_errs = res$std_errs,
    # hessian = res$hessian,
    vcov = vcov_mat,
    # logLik = res$max_likeli,
    # AIC = res$aic,
    # AICc = res$aic_corr,
    # R2 = res$R2,
    # fitted = res$source_data$predicted,
    # data = res$source_data,
    model = model_name,
    species = input_df$species[1],
    data.type = input_df$data.type[1]
  ))
}


# === MAIN EXECUTION ===
# MGB : load temporary data to test the fits
library(here)
library(dplyr)
#load(here("scof2012.rda")) # load scof2012 data


input_file <- "../results/KplantKroots_results/Kroot/Kroot_curve.csv"
output_dir <- "../results/KplantKroots_results/Kroot/Kroot_anneal_modified/"
plot_dir <- file.path(output_dir, "plots")
dir.create(output_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

df <- read.csv(input_file)
# df <- scof2012 |>
#   mutate(data.type = 1) |>
#   rename(gs = kl)

# Containers for results
all_results_df <- data.frame() # For Excel/table summaries
all_results_r <- list() # For full R objects (fitted values, Hessian, vcov, etc.)

species_list <- unique(df$species)
datatype_list <- unique(df$data.type)

for (sp in species_list) {
  for (dt in datatype_list) {
    subset_df <- df[df$species == sp & df$data.type == dt, ]

    if (nrow(subset_df) < 3) {
      next
    } # Skip if too few points

    # --- Exponential ---
    pE <- define_parsE(subset_df)
    fitE <- do_the_thing_nonlinear(
      subset_df,
      Exponential,
      pE$parsE,
      pE$par_loE,
      pE$par_highE,
      "Exponential",
      plot_dir
    )
    all_results_r[[paste(sp, dt, "Exponential", sep = "_")]] <- fitE
    all_results_df <- rbind(all_results_df, fitE$summary)

    # --- Exponential2 ---
    pE2 <- define_parsE2(subset_df)
    fitE2 <- do_the_thing_nonlinear(
      subset_df,
      Exponential2,
      pE2$parsE2,
      pE2$par_loE2,
      pE2$par_highE2,
      "Exponential2",
      plot_dir
    )
    all_results_r[[paste(sp, dt, "Exponential2", sep = "_")]] <- fitE2
    all_results_df <- rbind(all_results_df, fitE2$summary)

    # --- Logistic ---
    pL <- define_parsL(subset_df)
    fitL <- do_the_thing_nonlinear(
      subset_df,
      Logistic,
      pL$parsL,
      pL$par_loL,
      pL$par_highL,
      "Logistic",
      plot_dir
    )
    all_results_r[[paste(sp, dt, "Logistic", sep = "_")]] <- fitL
    all_results_df <- rbind(all_results_df, fitL$summary)

    # --- Sigmoidal ---
    pS <- define_parsS(subset_df)
    fitS <- do_the_thing_nonlinear(
      subset_df,
      Sigmoidal,
      pS$parsS,
      pS$par_loS,
      pS$par_highS,
      "Sigmoidal",
      plot_dir
    )
    all_results_r[[paste(sp, dt, "Sigmoidal", sep = "_")]] <- fitS
    all_results_df <- rbind(all_results_df, fitS$summary)

    # --- Linear ---
    fitLnr <- do_the_thing_linear(subset_df, Linear, "Linear", plot_dir)
    all_results_r[[paste(sp, dt, "Linear", sep = "_")]] <- fitLnr
    all_results_df <- rbind(all_results_df, fitLnr$summary)
  }
}

# === CLEAN COLUMN NAMES (if needed) ===
names(all_results_df) <- gsub("A\\.A", "A", names(all_results_df))
names(all_results_df) <- gsub("B\\.B", "B", names(all_results_df))
names(all_results_df) <- sub("\\.A$", "A", names(all_results_df))
names(all_results_df) <- sub("\\.B$", "B", names(all_results_df))

# MGB: check that the variance covariance matrix exists for each model type/species

## all should equal FALSE
vcov_exists <- sapply(seq_along(all_results_r), \(res) {
  is.null(all_results_r[[res]]$vcov)
})

# all should be TRUE
vcov_is_matrix <- sapply(seq_along(all_results_r), \(res) {
  is.matrix(all_results_r[[res]]$vcov)
})

# these model types should create a 4x4 matrix(Lin, Exp1 = 3x3)
# i.e.  3 parameters and error
mods_with_three_params <- c("Exponential2", "Sigmoidal", "Logistic")

# all should be true
vcov_matrix_dims_correct <- sapply(seq_along(all_results_r), \(res) {
  vcov_mat_dim <- dim(all_results_r[[res]]$vcov)
  mod_type <- all_results_r[[res]]$model

  if (mod_type %in% mods_with_three_params) {
    return(all(vcov_mat_dim == c(4, 4)))
  } else {
    return(all(vcov_mat_dim == c(3, 3)))
  }
})

# === SAVE RESULTS ===
write.csv(
  all_results_df,
  file.path(output_dir, "model_fitting_results_Kroot_anneal_modified.csv"),
  row.names = FALSE
)
saveRDS(all_results_r, file.path(output_dir, "model_fitting_results_full_Kroot.RDS"))
