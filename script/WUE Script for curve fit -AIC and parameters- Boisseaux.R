##############################################
# Stomatal Conductance Model Fitting Script
# Multi-species version (looped by species × data.type)
##############################################

require(likelihood)

# === MODEL DEFINITIONS ===

Gaussian <- function(a, mu, sigma, psi) {
  a * exp(-(psi - mu)^2 / (2 * sigma^2))
}

define_parsG <- function(input_df){
  list(
    parsG     = list(a = max(input_df$gs), mu = mean(input_df$psi, na.rm = TRUE), sigma = 1, sd = 2),
    par_loG   = list(a = 0, mu = min(input_df$psi, na.rm = TRUE), sigma = 0.01, sd = 0.0005),
    par_highG = list(a = max(input_df$gs)*2, mu = max(input_df$psi, na.rm = TRUE), sigma = 10, sd = 50)
  )
}


# === FITTING FUNCTIONS ===
do_the_thing_nonlinear <- function(input_df, model_type, pars1, par_lo1, par_hi1, model_name, output_plot_dir) {

  var <- list(psi="psi", x="gs", mean="predicted", log=TRUE)
  res <- anneal(model = model_type, par = pars1,
                source_data = input_df, var = var,
                par_lo = par_lo1, par_hi = par_hi1,
                dep_var = "gs", pdf = dnorm,
                max_iter = 50000, show_display = FALSE,
                temp_red = 0.05)
  
  sterror <- res$std_errs
  a     <- as.numeric(res$best_pars[1])
  mu    <- as.numeric(res$best_pars[2])
  sigma <- as.numeric(res$best_pars[3])

  #output results
  parvec <- c(
    Species   = as.character(input_df$species[1]),
    data.type = as.character(input_df$data.type[1]),
    model     = "Gaussian",
    a         = a,
    mu        = mu,
    sigma     = sigma,
    sterror_a = sterror[1],
    sterror_mu = sterror[2],
    sterror_sigma = sterror[3],
    loglikeli = res$max_likeli,
    rsq       = res$R2,
    AIC       = res$aic,
    AICcorr   = res$aic_corr,
    N         = length(res$source_data$gs)
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
input_file <- "../results/WUE/WUE_curve.csv"
output_dir <- "../results/WUE/"
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
    
    pG <- define_parsG(subset_df)
    all_results <- rbind(all_results,do_the_thing_nonlinear(subset_df, Gaussian, pG$parsG, pG$par_loG, pG$par_highG, "Gaussian", plot_dir))
  }
}

# === FIX COLUMN NAMES (A, B) ===
names(all_results) <- gsub("A\\.A", "A", names(all_results))
names(all_results) <- gsub("B\\.B", "B", names(all_results))
names(all_results) <- sub("\\.A$", "A", names(all_results))
names(all_results) <- sub("\\.B$", "B", names(all_results))

# === SAVE CLEANED RESULTS ===
write.csv(all_results, file.path(output_dir, "model_fitting_results.csv"), row.names = FALSE)
