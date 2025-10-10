#1

data<-read.csv("../results/GEx_A_curve_avsa.csv", stringsAsFactors = FALSE, strip.white = TRUE, na.strings = "")
library(likelihood)

########
# Assuming your data frame is called "data" and the column with the species codes is called "species"

# Get unique species codes
unique_species <- unique(data$species)

# Subset the data frame to include only the columns "psi" and "gs" and the unique species codes
subset_data <- subset(data, select = c("species", "psi", "gs"))

# View the first few rows of the subsetted data frame
head(subset_data)
#subset_data_CRMULL <- subset(subset_data, select = )


# Assuming your data frame is called "data" and the column with the species codes is called "Species"

# Split the data frame by unique species codes
split_data <- split(subset_data, subset_data$species)

# Create separate data frames for each species code and assign them to individual variables
list2env(split_data, envir = .GlobalEnv)

# Now, each species data frame is stored as a separate variable with the species code as its name.
# For example, if you have a species code "ABC", you can access its data frame using the variable name "ABC".

#data

#original: SAKI_test <- read.csv(".csv", header=TRUE)

model <-function (A,B,Xo,psi){A/(1+((psi/Xo)^B))}
likelisum<-function(pars,psi,observed)
{
  A<-pars[[1]]
  B<-pars[[2]]
  Xo<-pars[[3]]
  sd<-pars[[4]]
  predicted <-model(A,B,Xo,psi)
  likelisum<-sum(dnorm(observed,predicted,sd,log=T))
}
pars <- list(A=100,B=2,Xo=1,sd=2)
res<-optim(pars,likelisum,method="SANN",control=list(fnscale=-1,maxit=5000),hessian=T,psi=Avsa$psi,observed=Avsa$gs)

pars <-list(A=res$par[1], B=res$par[2], Xo=res$par[3], sd=res$par[4])
res

res<-optim(pars,likelisum,method="Nelder-Mead",control=list(fnscale=-1,maxit=5000),hessian=T,psi=Avsa$psi,observed=Avsa$gs)
expected<-model(res$par[1], res$par[2], res$par[3], Avsa$psi)
res
slope <-sum(expected*Avsa$gs)/sum(Avsa$gs^2)
slope
rsq <- 1-(sum((Avsa$gs-expected)^2)/sum((Avsa$gs-mean(Avsa$gs))^2))
resid<-(expected-Avsa$gs)
hist(resid)
sterror<-sqrt(diag(solve(-1*res$hessian)))
cat("parameters",res$par,"\n loglikeli",res$value,"\n rsq",rsq,"\n slope",slope,"\n sterror",sterror)

##PLOT  

# Create a new JPEG device and set the width and height in pixels

jpeg("AVSA.jpg", width = 800, height = 600)

# Create the plot

plot3 <- read.csv("../results/FVFM_dehydration_curve.csv", header=TRUE)

plot(plot3$psi, plot3$FvFm, type = "p", col = "magenta",xlab = "RWC", ylab = "Fv/Fm", cex.axis = 2.0, cex.lab= 2.0)

lines(plot3$psi, plot3$AVSA, type = "l", col = "magenta", lwd = 3)

points(plot3$psi,plot3$FvFm, type = "p", col = "forestgreen")

lines(plot3$psi, plot3$AVSA, type = "l", col = "forestgreen", lwd = 3)

plot3


# Close the JPEG device and save the file

dev.off()

