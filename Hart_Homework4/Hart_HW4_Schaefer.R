# Amanda Hart
# Homework 4 FISH 559

# bzero ==1165.16

setwd("/Users/ahart2/Research/TMB_Homeworks_UW/Hart_Homework4")

require(TMB)

# Read in data 
data <- read.csv("HWK4.csv", header = TRUE)
Catch <- data[1:30,2:3]
colnames(Catch) <- c("Catch", "Year")
Catch_vals <- Catch$Catch
Catch_vals <- as.numeric(as.character(Catch_vals))
Catch_yr <- Catch$Year
Catch_yr <- as.numeric(as.character(Catch_yr))
Survey <- data[32:37,2:4]
colnames(Survey) <- c("Abundance", "SE", "Year")
Survey_vals <- Survey$Abundance
Survey_vals <- as.numeric(as.character(Survey_vals))
Survey_yr <- Survey$Year
Survey_yr <- as.numeric(as.character(Survey_yr))
Survey_SE <- Survey$SE
Survey_SE <- as.numeric(as.character(Survey_SE))


# Create list of data objects, names of list items must match DATA objects in Cpp code
ModelData <- list(survey_vals = Survey_vals, survey_yr = Survey_yr, survey_SE = Survey_SE, catch_obs = Catch_vals, catch_yr = Catch_yr)

# Create list of parameters and provide initial values (may include parameter vector, e.g. param_vec = rep(0,5))
ModelParameters <- list(dummy=0, Bzero=1000, logr_growth = -1.609438, logF_y = c(rep(-2.0, length(Catch[,"Year"])-1))) # must be a list even if only 1 item in list
  # log_r is sensitive to starting value near the bounds

# Compile Cpp code
compile("Hart_HW4_Schaefer.cpp") # file must be in working directory or provide full file path name
dyn.load(dynlib("Hart_HW4_Schaefer"))


# Use map function to specify which parameters to estimate, those that are not estimated are fixed at initial values and must have a factor(NA) in map list
ModelMap <- list(dummy = factor(NA)) # rep(factor(NA),5) for a parameter vector of length 5

# Construct objective function to optimize based on data, parameters, and Cpp code
Model <- MakeADFun(data = ModelData, parameters = ModelParameters, DLL="Hart_HW4_Schaefer",silent=T,map = ModelMap) # silent=T silences a bunch of extra print statements


# ###################################### Run TestSchaefer.cpp
# compile("TestSchaefer.cpp") # file must be in working directory or provide full file path name
# dyn.load(dynlib("TestSchaefer"))
# 
# 
# # Use map function to specify which parameters to estimate, those that are not estimated are fixed at initial values and must have a factor(NA) in map list
# ModelMap <- list(dummy = factor(NA)) # rep(factor(NA),5) for a parameter vector of length 5
# 
# # Construct objective function to optimize based on data, parameters, and Cpp code
# Model <- MakeADFun(data = ModelData, parameters = ModelParameters, DLL="TestSchaefer",silent=T,map = ModelMap) # silent=T silences a bunch of extra print statements
# ###################################





# Set bounds on different parameters, length of this vector must equal number of estimate parameters
  # dummy=0, Bzero=220, logr_growth bounded from (-10, 0) = exp(4.539993e-05,1),    F_y bounded from logF_y (-20 to -0..01) = exp(2.061154e-09,0.9900498 ) 
lowbnd <- c(500,-10,rep(-20,length(Catch[,"Year"])-1)) # rep( 0.1, 5) for a parameter vector of length 5 with lower bound 0.1, syntax for upper bound is the same
uppbnd <- c(1500,0,rep(-0.01,length(Catch[,"Year"])-1))

# Fit model to data using structure provided by MakeADFun() function call
# eval.max = max number evaluations of objective function
# iter.max = max number of iterations allowed
# rel.tol = relative tolerance 
fit <- nlminb(Model$par, Model$fn, Model$gr, control=list(rel.tol=1e-12,eval.max=100000,iter.max=1000), lower=lowbnd,upper=uppbnd) # notice the inclusion of the lower and upper bound vectors
# fit <- nlminb(Model$par, Model$fn, Model$gr, control=list(eval.max=100000,iter.max=1000)) # no bounds on parameters

##### Fitted model results #####
# Best parameter estimates
best <- Model$env$last.par.best
print(best) 

# Report parameter estimates & std error
rep <- sdreport(Model)  
print(summary(rep))

# Print objective function
print(Model$report()$obj_fun)

# # print objective (likelihood)
# fit$objective
# 
# # Check for Hessian
# VarCo <- solve(Model$he())
# print(sqrt(diag(VarCo)))

# Get reported info & predicted data
# Predicted <- Model$report()$PredictedVariableNameHere 
obj_fun <- Model$report()$obj_fun
obj_fun
rNLL <- Model$report()$rNLL
rNLL
BzeroNLL <- Model$report()$BzeroNLL
BzeroNLL
biomass <- Model$report()$biomass
biomass
catch_pred <- Model$report()$catch_pred
catch_pred


