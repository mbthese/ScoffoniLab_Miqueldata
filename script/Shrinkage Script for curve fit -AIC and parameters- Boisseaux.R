##############################################
# Stomatal Conductance Model Fitting Script
# Multi-species version (looped by species × data.type)
##############################################

require(likelihood)

# === MODEL DEFINITIONS ===
Exponential <- function (A,B,psi){ A * (exp(-B * psi)) }
define_parsE <- function(input_df){
  list(
    parsE     = list(A = max(input_df$gs), B = 0.1, sd = 2),
    par_loE   = list(A = 0, B = 0, sd = 0.0005),
    par_highE = list(A = max(input_df$gs)*2, B = 10, sd = 50)
  )
}

Logistic <- function (A,B,Xo,psi){ A / (1 + ((psi / Xo)^B)) }
define_parsL <- function(input_df){
  list(
    parsL = list(A = 0.8, B = -1, Xo = 60, sd = 0.05),
    par_loL = list(A = 0, B = -5, Xo = 10, sd = 0.0005),
    par_highL = list(A =100, B = 5, Xo = 80, sd = 0.5)
  )
}

Sigmoidal <- function(A, B, Xo, psi) {A / (1 + exp(-((psi - Xo) / B)))}
define_parsS <- function(input_df) {
  list(
    parsS     = list(A = max(input_df$gs), B = 0, Xo = 30, sd = 2),
    par_loS   = list(A = 0, B = 0, Xo = 10, sd = 0.0005),
    par_highS = list(A = max(input_df$gs) * 1.5, B = -30, Xo = 90, sd = 50))
}

Exponential2 <- function (A,B,C,psi){ C + A * (exp(-B * psi)) }
define_parsE2 <- function(input_df){
  list(
    parsE2     = list(A = max(input_df$gs), B = 0.5, C = 0.1, sd = 2),
    par_loE2   = list(A = 0, B = 0.011, C = 0.001, sd = 0.0005),
    par_highE2 = list(A = max(input_df$gs)*2, B = 20, C = 20, sd = 50)
  )
}

Linear <- function (A,B,psi){ B * psi + A }

# === FITTING FUNCTIONS ===
do_the_thing_nonlinear <- function(input_df, model_type, pars1, par_lo1, par_hi1, model_name, output_plot_dir) {
  # --- Special handling for Sigmoidal ---
  if (tolower(model_name) == "sigmoidal") {
    sig_params <- define_parsS(input_df)  # uses the current definition without Trae special case
    pars1   <- sig_params$parsS
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
  
  var <- list(psi="psi", x="gs", mean="predicted", log=TRUE)
  res <- anneal(model = model_type, par = pars1,
                source_data = input_df, var = var,
                par_lo = par_lo1, par_hi = par_hi1,
                dep_var = "gs", pdf = dnorm,
                max_iter = 50000, show_display = FALSE,
                temp_red = 0.05)
  
  sterror <- res$std_errs
  A <- as.numeric(res$best_pars[1])
  B <- as.numeric(res$best_pars[2])
  C <- as.numeric(ifelse(length(res$best_pars) >= 3, res$best_pars[3], NA))
  D <- as.numeric(ifelse(length(res$best_pars) >= 4, res$best_pars[4], NA))
  
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
  } else if (model_lower == "logistic") {
    psi_kleaf20 <- C * (1 / 0.8 - 1)^(1 / B)
    psi_kleaf50 <- C * (1 / 0.5 - 1)^(1 / B)
    psi_kleaf80 <- C * (1 / 0.2 - 1)^(1 / B)
    psi_kleaf95 <- C * (1 / 0.05 - 1)^(1 / B)
  } else if (model_lower == "sigmoidal") {
    psi_kleaf20 <- -B * log(A / (0.8 * A) - 1) + C
    psi_kleaf50 <- -B * log(A / (0.5 * A) - 1) + C
    psi_kleaf80 <- -B * log(A / (0.2 * A) - 1) + C
    psi_kleaf95 <- -B * log(A / (0.05 * A) - 1) + C
  } else if (model_lower == "exponential") {
    psi_kleaf20 <- log(0.8) / (-B)
    psi_kleaf50 <- log(0.5) / (-B)
    psi_kleaf80 <- log(0.2) / (-B)
    psi_kleaf95 <- log(0.05) / (-B)
  } else if (model_lower == "exponential2") {
    psi_kleaf20 <- log(((C * 0.8) + (0.8 * A) - C) / A) / (-B)
    psi_kleaf50 <- log(((C * 0.5) + (0.5 * A) - C) / A) / (-B)
    psi_kleaf80 <- log(((C * 0.2) + (0.2 * A) - C) / A) / (-B)
    psi_kleaf95 <- log(((C * 0.05) + (0.05 * A) - C) / A) / (-B)
  } else {
    psi_kleaf20 <- psi_kleaf50 <- psi_kleaf80 <- psi_kleaf95 <- NA
  }
  # ---------------------------
  
  parvec <- c(
    Species     = as.character(input_df$species[1]),
    data.type   = as.character(input_df$data.type[1]),
    model       = model_name,
    A           = A,
    B           = B,
    C           = C,
    D           = D,
    sterrorA    = sterror[1],
    sterrorB    = sterror[2],
    sterrorC    = ifelse(length(sterror) >= 3, sterror[3], NA),
    sterrorD    = ifelse(length(sterror) >= 4, sterror[4], NA),
    loglikeli   = res$max_likeli,
    rsq         = res$R2,
    slope       = sum(res$source_data$predicted * res$source_data$gs) / sum(res$source_data$predicted^2),
    AIC         = res$aic,
    AICcorr     = res$aic_corr,
    N           = length(res$source_data$gs),
    Kmax        = Kmax,
    Kmax_at_0.1MPa = Kmax_0.1,
    P50_0.1kmax = P50_0.1kmax,
    psi_kleaf20 = psi_kleaf20,
    psi_kleaf50 = psi_kleaf50,
    psi_kleaf80 = psi_kleaf80,
    psi_kleaf95 = psi_kleaf95
  )
  
  png(file.path(output_plot_dir, paste0(input_df$species[1], "_", input_df$data.type[1], "_", model_name, ".png")),
      width = 800, height = 600)
  plot(res$source_data$psi, res$source_data$gs,
       xlab = "Water Potential (-MPa)", ylab = "Stomatal Conductance (mmol m-2 s-1)",
       main = paste(input_df$species[1], input_df$data.type[1], model_name))
  for_plotting <- cbind(res$source_data$psi, res$source_data$predicted)
  for_plotting <- for_plotting[order(for_plotting[,1]),]
  lines(for_plotting[,1], for_plotting[,2], col="blue", lwd=2)
  dev.off()
  
  return(parvec)
}


do_the_thing_linear <- function(input_df, model_type, model_name, output_plot_dir) {
  var <- list(psi="psi", x="gs", mean="predicted", log=TRUE)
  lm_fit <- lm(input_df$gs ~ input_df$psi)
  pars <- list(A = coef(lm_fit)[1], B = coef(lm_fit)[2], sd = 1)
  par_lo <- list(A = coef(lm_fit)[1]*0.05, B = coef(lm_fit)[2]*2, sd = 0.005)
  par_hi <- list(A = coef(lm_fit)[1]*4, B = coef(lm_fit)[2]*0.1, sd = 20)
  
  res <- anneal(model = model_type, par = pars,
                source_data = input_df, var = var,
                par_lo = par_lo, par_hi = par_hi,
                dep_var = "gs", pdf = dnorm,
                max_iter = 50000, show_display = FALSE,
                temp_red = 0.05)
  
  sterror <- res$std_errs
  A <- as.numeric(res$best_pars[1])
  B <- as.numeric(res$best_pars[2])
  
  # ---- New calculations ----
  Kmax_0.1 <- A + B * 0.1
  Kmax <- A
  P50_0.1kmax <- (Kmax_0.1 / 2 - A) / B
  psi_kleaf20 <- ((0.8 - 1) * (-A) / B) * -1
  psi_kleaf50 <- ((0.5 - 1) * (-A) / B) * -1
  psi_kleaf80 <- ((0.2 - 1) * (-A) / B) * -1
  psi_kleaf95 <- ((0.05 - 1) * (-A) / B) * -1
  # ---------------------------
  
  parvec <- c(
    Species     = as.character(input_df$species[1]),
    data.type   = as.character(input_df$data.type[1]),
    model       = model_name,
    A           = A,
    B           = B,
    C           = ifelse(length(res$best_pars) >= 3, res$best_pars[3], NA),
    D           = ifelse(length(res$best_pars) >= 4, res$best_pars[4], NA),
    sterrorA    = sterror[1],
    sterrorB    = sterror[2],
    sterrorC    = ifelse(length(sterror) >= 3, sterror[3], NA),
    sterrorD    = ifelse(length(sterror) >= 4, sterror[4], NA),
    loglikeli   = res$max_likeli,
    rsq         = res$R2,
    slope       = sum(res$source_data$predicted * res$source_data$gs) / sum(res$source_data$predicted^2),
    AIC         = res$aic,
    AICcorr     = res$aic_corr,
    N           = length(res$source_data$gs),
    Kmax        = Kmax,
    Kmax_at_0.1MPa = Kmax_0.1,
    P50_0.1kmax = P50_0.1kmax,
    psi_kleaf20 = psi_kleaf20,
    psi_kleaf50 = psi_kleaf50,
    psi_kleaf80 = psi_kleaf80,
    psi_kleaf95 = psi_kleaf95
  )
  
  png(file.path(output_plot_dir, paste0(input_df$species[1], "_", input_df$data.type[1], "_", model_name, ".png")),
      width = 800, height = 600)
  plot(res$source_data$psi, res$source_data$gs,
       xlab = "Water Potential (-MPa)", ylab = "Stomatal Conductance (mmol m-2 s-1)",
       main = paste(input_df$species[1], input_df$data.type[1], model_name))
  for_plotting <- cbind(res$source_data$psi, res$source_data$predicted)
  for_plotting <- for_plotting[order(for_plotting[,1]),]
  lines(for_plotting[,1], for_plotting[,2], col="blue", lwd=2)
  dev.off()
  
  return(parvec)
}


# === MAIN EXECUTION ===
input_file <- "../results/shrinkage/shrinkage_curve_area.csv"
output_dir <- "../results/shrinkage/area/"
plot_dir <- file.path(output_dir, "plots")
dir.create(output_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

df <- read.csv(input_file)

all_results <- data.frame()

species_list <- unique(df$species)
datatype_list <- unique(df$data.type)

for (sp in species_list) {
  for (dt in datatype_list) {
    subset_df <- df[df$species == sp & df$data.type == dt, ]
    if (nrow(subset_df) < 3) next
    
   pE <- define_parsE(subset_df)
   all_results <- rbind(all_results, do_the_thing_nonlinear(subset_df, Exponential, pE$parsE, pE$par_loE, pE$par_highE, "Exponential", plot_dir))

  pE2 <- define_parsE2(subset_df)
  all_results <- rbind(all_results, do_the_thing_nonlinear(subset_df, Exponential2, pE2$parsE2, pE2$par_loE2, pE2$par_highE2, "Exponential2", plot_dir))

    pL <- define_parsL(subset_df)
    all_results <- rbind(all_results, do_the_thing_nonlinear(subset_df, Logistic, pL$parsL, pL$par_loL, pL$par_highL, "Logistic", plot_dir))

   pS <- define_parsS(subset_df)
   all_results <- rbind(all_results, do_the_thing_nonlinear(subset_df, Sigmoidal, pS$parsS, pS$par_loS, pS$par_highS, "Sigmoidal", plot_dir))
    
    all_results <- rbind(all_results, do_the_thing_linear(subset_df, Linear, "Linear", plot_dir))
  }
}

# === FIX COLUMN NAMES (A, B) ===
names(all_results) <- gsub("A\\.A", "A", names(all_results))
names(all_results) <- gsub("B\\.B", "B", names(all_results))
names(all_results) <- sub("\\.A$", "A", names(all_results))
names(all_results) <- sub("\\.B$", "B", names(all_results))

# === SAVE CLEANED RESULTS ===
write.csv(all_results, file.path(output_dir, "model_fitting_results.csv"), row.names = FALSE)
