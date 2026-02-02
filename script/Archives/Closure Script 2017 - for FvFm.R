source("script/R function MB 2017 - for FvFm.R")
dd = read.csv("../results/fvfm/FVFM_recovery_curve.csv",header=TRUE)

dd$to_split = paste(dd$species, dd$data.type)
aggregate(dd, by =list(dd$to_split), FUN=length) -> speciesN
speciesN[which(speciesN$to_split > 5),1] -> species_list
#Trying to fit curves with less than 6 points returns -inf for the AICcor. Are those the binned data? Would it be possible to make the bins smaller so there are more points?

# first, just to make sure that things work as expected, let's print out
# the kleaf column for each species
for (ii in 1:length(species_list)){
  subset(dd, dd$to_split == species_list[ii], select =c(1:5)) -> data_by_sp
  print(data_by_sp)
}
#I couldn't get the list format to save as a usable .csv file, so I changed everything to dataframes

# use these kleaf's as inputs to the likelihood function


#Log-Logistic
pdf(file = '../results/fvfm/FVFM_recovery_LogLogistic_fits.pdf', width = 8.5, height = 5, onefile = T)
modelfitting_results[FALSE, ] -> modelfitting_results  # Empty results data frame

for (ii in 1:length(species_list)) {
  subset(dd, dd$to_split == species_list[ii], select = c(1:5)) -> data_by_sp
  define_parsLL(data_by_sp) -> par_estimates
  parsLL = par_estimates[[1]]
  par_loLL = par_estimates[[2]]
  par_highLL = par_estimates[[3]]
  LogLogistic_fits = do_the_thing_nonlinear(data_by_sp, LogLogistic, parsLL, par_loLL, par_highLL)
  rbind(modelfitting_results, as.data.frame(LogLogistic_fits)) -> modelfitting_results
  print(ii)
}
dev.off()
modelfitting_results$Kmax <- modelfitting_results$A.A
modelfitting_results$psi_kleaf20 <- modelfitting_results$C.Xo * (1 / 0.8 - 1)^(1 / modelfitting_results$B.B)
modelfitting_results$psi_kleaf50 <- modelfitting_results$C.Xo * (1 / 0.5 - 1)^(1 / modelfitting_results$B.B)
modelfitting_results$psi_kleaf80 <- modelfitting_results$C.Xo * (1 / 0.2 - 1)^(1 / modelfitting_results$B.B)
modelfitting_results$psi_kleaf95 <- modelfitting_results$C.Xo * (1 / 0.05 - 1)^(1 / modelfitting_results$B.B)


write.csv(modelfitting_results, file = "../results/fvfm/FVFM_recovery_LogLogistic_fits.csv")
cat("Log-Logistic has ", nrow(modelfitting_results), " fitted curves.\n")

# # Weibull
# pdf(file = '../results/FVFM_recovery_Weibull_fits.pdf', width = 8.5, height = 5, onefile = T)
# modelfitting_results[FALSE, ] -> modelfitting_results  # Empty results data frame
# 
# for (ii in 1:length(species_list)) {
#   subset(dd, dd$to_split == species_list[ii], select = c(1:5)) -> data_by_sp
#   define_parsW(data_by_sp) -> par_estimates
#   parsW = par_estimates[[1]]
#   par_loW = par_estimates[[2]]
#   par_highW = par_estimates[[3]]
#   Weibull_fits = do_the_thing_nonlinear(data_by_sp, Weibull, parsW, par_loW, par_highW)
#   rbind(modelfitting_results, as.data.frame(Weibull_fits)) -> modelfitting_results
#   print(ii)
# }
# dev.off()
# modelfitting_results$Kmax <- modelfitting_results$A.A
# modelfitting_results$psi_kleaf20 <- modelfitting_results$B.B * (-log(0.8))^(1 / modelfitting_results$C.C)
# modelfitting_results$psi_kleaf50 <- modelfitting_results$B.B * (-log(0.5))^(1 / modelfitting_results$C.C)
# modelfitting_results$psi_kleaf80 <- modelfitting_results$B.B * (-log(0.2))^(1 / modelfitting_results$C.C)
# modelfitting_results$psi_kleaf95 <- modelfitting_results$B.B * (-log(0.05))^(1 / modelfitting_results$C.C)
# 
# 
# write.csv(modelfitting_results, file = "../results/FVFM_recovery_Weibull_fits.csv")
# cat("Weibull has ", nrow(modelfitting_results), " fitted curves.\n")


