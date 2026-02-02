##############################################
# RC Model Fitting Script
# Multi-species version (looped by species × data.type)
##############################################

# Notes:
# sd =   to model the variability of the data, assumed standard deviation of residuals between your observed gs values and the model predictions.
# It assumes the data has some error modeled as a normal distribution
# So sd gives the optimizer an initial guess of how much variability there is in the data.

require(likelihood)

# === MODEL DEFINITIONS ===

Logistic <- function (A,B,Xo,psi){ A / (1 + ((psi / Xo)^B)) } #if plotting with psi
define_parsL <- function(input_df){
  list(
    parsL     = list(A = 0.5, B = -1, Xo = 1, sd = 2),
    par_loL   = list(A = 0, B = -10, Xo = 0, sd = 0.0005),
    par_highL = list(A = max(input_df$gs)*2, B = 5, Xo = 20, sd = 50)
  )
}

Exponential <- function (A,B,psi){ A * (1- exp(-B * psi)) }
define_parsE <- function(input_df){
  list(
    parsE     = list(A = 1, B = 0.1, sd = 2),
    par_loE   = list(A = 1, B = 0.01, sd = 0.0005),
    par_highE = list(A = 1, B = 10, sd = 50)
  )
}

Sigmoidal <- function(A, B, Xo, psi) {
  A / (1 + exp(-((psi - Xo) / B)))  
}

define_parsS <- function(input_df) {
  list(
  parsS     = list(A = 1, B = 0.5, Xo = 20, sd = 2),
  par_loS   = list(A = 1, B = 0, Xo = 10, sd = 0.0005),
  par_highS = list(A = 1, B = 4, Xo = 90, sd = 50))
}

Exponential2 <- function (A,B,C,psi){ C + A * (1-exp(-B * psi)) }
define_parsE2 <- function(input_df){
  list(
    parsE2     = list(A = 1, B = 0.1, C = 0.1, sd = 2),
    par_loE2   = list(A = 1, B = 0.01, C = 0, sd = 0.0005),
    par_highE2 = list(A = 1, B = 10, C = 10, sd = 50)
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
  
  
  # --- Dehydration dataset ---
  max_RC <- data.frame(
    Species = c("Avsa","Chga","Hovu","Pegl","Trae","Zema"),
    Kmax = c(1.0036853, 0.9984288,1.0138406, 0.9967606, 0.9653397, 0.9919813)
  )
  
  # --- Recovery dataset ---
  # mean_FvFm_recovery <- data.frame(
  #   Species = c("Avsa","Chga","Hovu","Pegl","Trae","Zema"),
  #   Kmax = c(0.7483636, 0.7428, 0.7353333, 0.74, 0.62, 0.7058571)
  # )
  
  # --- Select dataset by uncommenting ---
  mean_RC <- max_RC
  # mean_FvFm <- mean_FvFm_recovery
  
  # --- Get Kmax dynamically for the current species ---
  current_species <- unique(input_df$species)
  if (length(current_species) != 1) stop("Input dataframe should contain only one species per call.")
  
  Kmax <- mean_RC$Kmax[mean_RC$Species == current_species]
  
  if (length(Kmax) != 1) stop(paste("Kmax not found for species:", current_species))
  
  
  if (model_lower == "linear") {
    psi_kleaf20 <- (0.8 - A) / B
    psi_kleaf50 <- (0.5  - A) / B
    psi_kleaf80 <- (0.2 - A) / B
    psi_kleaf95 <- (0.05 - A) / B
  } else if (model_lower == "weibull") {
    psi_kleaf20 <- B * (-log(1 - 0.8 * Kmax / A))^(1 / C)
    psi_kleaf50 <- B * (-log(1 - 0.5 * Kmax / A))^(1 / C)
    psi_kleaf80 <- B * (-log(1 - 0.2 * Kmax / A))^(1 / C)
    psi_kleaf95 <- B * (-log(1 - 0.05 * Kmax / A))^(1 / C)
    
  } else if (model_lower == "logistic") {
    psi_kleaf20 <- C * ((A / 0.80 - 1)^(1 / B))
    psi_kleaf50 <- C * ((A/ 0.50 - 1)^(1 / B))  # if only C * ((A/ 0.50 - 1)^(1 / B)) then i don't solve for y=  A/2 , with A=1 because A for the curve is 1.3, so i need to factor it in
    psi_kleaf80 <- C * ((A / 0.20 - 1)^(1 / B))
    psi_kleaf95 <- C * ((A / 0.05 - 1)^(1 / B))
  } else if (model_lower == "sigmoidal") {
    psi_kleaf20 <- C - B * log(0.80 / (1 - 0.80))
    psi_kleaf50 <- C - B * log(0.50 / (1 - 0.50))  # = C
    psi_kleaf80 <- C - B * log(0.20 / (1 - 0.20))
    psi_kleaf95 <- C - B * log(0.05 / (1 - 0.05))
  } else if (model_lower == "exponential") {
    psi_kleaf20 <- -(1 / B) * log(1 - 0.20*A)  
    psi_kleaf50 <- -(1 / B) * log(1 - 0.50*A)  
    psi_kleaf80 <- -(1 / B) * log(1 - 0.80*A)  
    psi_kleaf95 <- -(1 / B) * log(1 - 0.95*A) 
  } else if (model_lower == "exponential2") {
    psi_kleaf20 <- -(1/B) * log(((1 - ((0.2-C)/A))))
    psi_kleaf50 <- -(1/B) * log(((1 - ((0.5-C)/A))))
    psi_kleaf80 <--(1/B) * log(((1 - ((0.8-C)/A))))
    psi_kleaf95 <- -(1/B) * log(((1 - ((0.95-C)/A))))
    
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
    psi_kleaf20 = psi_kleaf20,
    psi_kleaf50 = psi_kleaf50,
    psi_kleaf80 = psi_kleaf80,
    psi_kleaf95 = psi_kleaf95
  )
  
  png(file.path(output_plot_dir, paste0(input_df$species[1], "_", input_df$data.type[1], "_", model_name, ".png")),
      width = 800, height = 600)
  plot(res$source_data$psi, res$source_data$gs,
       xlab = "RWC (%)", ylab = "RC",
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
  
  max_RC <- data.frame(
    Species = c("Avsa","Chga","Hovu","Pegl","Trae","Zema"),
    Kmax = c(1.0036853, 0.9984288,1.0138406, 0.9967606, 0.9653397, 0.9919813)
  )
  
  # --- Recovery dataset ---
  # mean_FvFm_recovery <- data.frame(
  #   Species = c("Avsa","Chga","Hovu","Pegl","Trae","Zema"),
  #   Kmax = c(0.7483636, 0.7428, 0.7353333, 0.74, 0.62, 0.7058571)
  # )
  
  # --- Select dataset by uncommenting ---
  mean_RC <- max_RC
  # mean_FvFm <- mean_FvFm_recovery
  
  # --- Get Kmax dynamically for the current species ---
  current_species <- unique(input_df$species)
  if (length(current_species) != 1) stop("Input dataframe should contain only one species per call.")
  
  Kmax <- mean_RC$Kmax[mean_RC$Species == current_species]
  
  if (length(Kmax) != 1) stop(paste("Kmax not found for species:", current_species))
  
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
    psi_kleaf20 = psi_kleaf20,
    psi_kleaf50 = psi_kleaf50,
    psi_kleaf80 = psi_kleaf80,
    psi_kleaf95 = psi_kleaf95
  )
  
  png(file.path(output_plot_dir, paste0(input_df$species[1], "_", input_df$data.type[1], "_", model_name, ".png")),
      width = 800, height = 600)
  plot(res$source_data$psi, res$source_data$gs,
       xlab = "RWC (%)", ylab = "RC",
       main = paste(input_df$species[1], input_df$data.type[1], model_name))
  for_plotting <- cbind(res$source_data$psi, res$source_data$predicted)
  for_plotting <- for_plotting[order(for_plotting[,1]),]
  lines(for_plotting[,1], for_plotting[,2], col="blue", lwd=2)
  dev.off()
  
  return(parvec)
}


# === MAIN EXECUTION ===
input_file <- "../results/RC/RWC/RC_curve.csv"
output_dir <- "../results/RC/RWC/"
plot_dir <- file.path(output_dir, "plots_RC_RWC")
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

     # pW <- define_parsW(subset_df)
     # all_results <- rbind(all_results, do_the_thing_nonlinear(subset_df, Weibull, pW$parsW, pW$par_loW, pW$par_highW, "Weibull", plot_dir))

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

