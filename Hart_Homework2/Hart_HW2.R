# Amanda Hart
# Homework 2 FISH 559 

setwd("/Users/ahart2/Research/TMB_Workshop_Fall2018/ClassWork/Hart_Homework2")
library(gtools)

##### Read in data #####
data <- read.table("Home2.dat", skip=2)
colnames(data) <- c("Year", "Sex_M", "Sex_F") # numbers of males and females

########## PART 1 #################
# Conventional estimate
###################################
sex_ratio1 <- sum(data[,"Sex_M"])/sum(data[,"Sex_M"] + data[,"Sex_F"])
sex_ratio1

########## PART 2 #################
# Fit model in TMB
###################################
require(TMB)

# Read in data here
data <- data[-which(data$Sex_M==0 | data$Sex_F==0),] # remove rows with 0s
fishCount <- data[,"Sex_M"] + data[,"Sex_F"]
Years <- data$Year
Males <- data$Sex_M
Females <- data$Sex_F

# Create list of data objects, names of list items must match DATA objects in Cpp code
ModelData <- list(Years = Years, Males = Males, Females = Females, fishCount = fishCount)

# Create list of parameters and provide initial values (may include parameter vector, e.g. param_vec = rep(0,5))
ModelParameters <- list(dummy=0, logitp_ratio=logit(0.5)) # must be a list even if only 1 item in list


# Compile Cpp code
compile("Hart_HW2.cpp") # file must be in working directory or provide full file path name
dyn.load(dynlib("Hart_HW2"))

# Use map function to specify which parameters to estimate, those that are not estimated are fixed at initial values and must have a factor(NA) in map list
ModelMap <- list(dummy = factor(NA)) # rep(factor(NA),5) for a parameter vector of length 5
# ModelMap <- list(p_ratio = factor(NA)) # rep(factor(NA),5) for a parameter vector of length 5

# Construct objective function to optimize based on data, parameters, and Cpp code
Model <- MakeADFun(data = ModelData, parameters = ModelParameters, DLL="Hart_HW2",silent=T,map = ModelMap) # silent=T silences a bunch of extra print statements

# Fit model to data using structure provided by MakeADFun() function call
fit <- nlminb(Model$par, Model$fn, Model$gr, control=list(rel.tol=1e-12,eval.max=100000,iter.max=1000)) # notice the inclusion of the lower and upper bound vectors
# fit <- nlminb(Model$par, Model$fn, Model$gr, control=list(eval.max=100000,iter.max=1000)) # no bounds on parameters

##### Fitted model results #####
# Best parameter estimates
best <- inv.logit(Model$env$last.par.best)
print(best) 

# Report parameter estimates & std error 
rep <- sdreport(Model)  
print(summary(rep))


########## PART 3 #################
# Fit to data with varying p_ratio
###################################
# Given 25 yr data with annual sample 100 fish and beta distributed p (p~beta(2,1)) each year
# fit to 1000 data sets and estimate p_ratio
# compare to 95% CI from PART 2

Nsim <- 1000
Nfish <- 100
Nyear <- 25

require(TMB)
# Compile Cpp code
compile("Hart_HW2.cpp") # file must be in working directory or provide full file path name
dyn.load(dynlib("Hart_HW2"))

RatioResults <- matrix(NA, ncol=2, nrow=Nsim)
colnames(RatioResults) <- c("p_ratio", "std_error")

# Run data simulations and fit model
for(isim in 1:Nsim){
  ##### Simulate data #####
  Males <- rep(NA,25)
  Females <- rep(NA,25)
  Years <- seq(1, 25, by=1)
  for(iyr in 1:Nyear){
    pbar <- rbeta(1,2,1) # Pick pbar from beta distribution
    Males[iyr] <- rbinom(1, size=Nfish, prob=pbar) # Use pbar to calculate number of males
    Females[iyr] <- 100 - Males[iyr]
  }
  
  # Create data, parameter & map lists 
  ModelData <- list(Years = Years, Males = Males, Females = Females, fishCount = Nfish)
  ModelParameters <- list(dummy=0, logitp_ratio=logit(0.5))
  ModelMap <- list(dummy = factor(NA)) # rep(factor(NA),5) for a parameter vector of length 5
  
  # Construct objective function to optimize based on data, parameters, and Cpp code
  Model <- MakeADFun(data = ModelData, parameters = ModelParameters, DLL="Hart_HW2",silent=T,map = ModelMap) # silent=T silences a bunch of extra print statements
  
  # Fit model to data using structure provided by MakeADFun() function call
  fit <- nlminb(Model$par, Model$fn, Model$gr, control=list(rel.tol=1e-12,eval.max=100000,iter.max=1000)) # notice the inclusion of the lower and upper bound vectors
  
  ##### Fitted model results #####
  # Best parameter estimates
  best <- inv.logit(Model$env$last.par.best)
  
  # Report parameter estimates & std error 
  rep <- sdreport(Model)  
  
  # store p_ratio, std. error
  RatioResults[isim,] <- summary(rep)[2,]
}

hist(x=RatioResults[,"p_ratio"], 
     xlab="Sex ratio (p)",
     main="Histogram of sex ratio",
     ylim=c(0,200),
     xlim=c(LowerCI-0.1,UpperCI+0.1))
abline(v=UpperCI, col="red", lwd=2)
abline(v=LowerCI, col="red", lwd=2)

# !!! now just need to figure out CI (make likelihood profile)






