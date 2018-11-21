# Amanda Hart
# Homework 4 FISH 559

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

# Specify model (non-estimated) parameter values
rho_val <- 0.995
M_val <- 0.2
w_val <- 0.5


# Create list of data objects, names of list items must match DATA objects in Cpp code
ModelData <- list(survey_vals = Survey_vals, survey_yr = Survey_yr, survey_SE = Survey_SE, catch_obs = Catch_vals, catch_yr = Catch_yr, rho = rho_val, w_dat = w_val, M_dat = M_val)

# Create list of parameters and provide initial values (may include parameter vector, e.g. param_vec = rep(0,5))
# ModelParameters <- list(dummy=0, Bzero=220, h_steep = 0.5, F_y = c(rep(0.5, length(Catch[,"Year"])-1))) # this may need to be of length 29 not 30, (30th is F for final projection year), if shorten to 29, don't store final recruitment & biomass in .cpp file
ModelParameters <- list(dummy=0, Bzero=1000, h_steep = 0.5, logF_y = c(rep(-2.0, length(Catch[,"Year"])-1))) 

# Compile Cpp code
compile("Hart_HW4_Deriso.cpp") # file must be in working directory or provide full file path name
dyn.load(dynlib("Hart_HW4_Deriso"))


# Use map function to specify which parameters to estimate, those that are not estimated are fixed at initial values and must have a factor(NA) in map list
ModelMap <- list(dummy = factor(NA)) # rep(factor(NA),5) for a parameter vector of length 5

# Construct objective function to optimize based on data, parameters, and Cpp code
Model <- MakeADFun(data = ModelData, parameters = ModelParameters, DLL="Hart_HW4_Deriso",silent=T,map = ModelMap) # silent=T silences a bunch of extra print statements

# Set bounds on different parameters, length of this vector must equal number of estimate parameters
  # dummy=0, Bzero=220, h_steep = 0.5, F_y = c(rep(0.5, length(Catch[,"Year"])))
lowbnd <- c(500,0.2,rep(-20,length(Catch[,"Year"])-1)) # rep( 0.1, 5) for a parameter vector of length 5 with lower bound 0.1, syntax for upper bound is the same
uppbnd <- c(1500,1,rep(-0.01,length(Catch[,"Year"])-1)) # no bounds for mapped params!!!!!!!!!!

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

# print objective (likelihood)
fit$objective

# Get reported info & predicted data
# Predicted <- Model$report()$PredictedVariableNameHere 
obj_fun <- Model$report()$obj_fun
obj_fun
biomass_pred <- Model$report()$biomass
biomass_pred
catch_pred <- Model$report()$catch_pred
catch_pred
recruitment <- Model$report()$recruitment
recruitment
# catch_obs <- Model$report()$catch_obs
# catch_obs
# Catch_vals


################################ MCMC ###############################
library(tmbstan)

# MCMC <- tmbstan(obj = Model,
#                 iter = 10000, # 1000 iterations
#                 init = "last.par.best", # ParamList, # "par" uses defaults from model object, alternatively use "last.par.best"
#                 chains = 1, # length(ParamList),
#                 warmup = 1000,
#                 thin = 10,
#                 lower = lowbnd,
#                 upper = uppbnd)
# 
# # randomly generate our own lists, pick 1 F at random for entire vector
# 
# 
# # library(adnuts)
# # 
# # # start MCMC chain at multiple locations (in this case, multiple values of h_steep)
# h_list <- c(best["h_steep"], seq(0.1, 0.9, by = 0.1)) # best estimated parameter and a series of other points
# numParamList <- 2
# ParamList <- vector("list", length(numParamList))
# for(i in 1:length(numParamList)){
#   # pick h from between bounds (uniform distribution)
#   h_val <- runif(1, 0.2, 1)
#   # pick Bzero from uniform distribution
#   Bzero_val <- runif(1, 500, 1500)
#   # pick F from uniform distribution
#   F_temp <- runif(1, -20, -0.01)
#   F_vals <- rep(F_temp, length(Catch[,"Year"])-1)
# 
# 
#   # Set different starting h values
#   tempParamList <- list(dummy=0, Bzero=Bzero_val, h_steep = h_val, logF_y = F_vals) # F_y may need to be updated based on the above
#     # what do you do with parameters that don't have priors?
# 
#   # fit model with stan
#   # Add to list that is passed to MCMC function
#   ParamList[[i]] <- tempParamList # if you don't add the double brackets, bad things happen
#   # Param list needs to be a list of lists conatining initial parameter vectors (sample initial parameter values from prior)
# }

numParamList <- 1

make_init <- function(chainID = 1){
  list(dummy=0, Bzero=runif(1, 500, 1500), h_steep = runif(1, 0.2, 1), logF_y = rep(runif(1, -20, -0.01), length(Catch[,"Year"])-1)) 
}

ParamList <- lapply(1:numParamList, function(id) make_init(chainID = id)) # function(id) make_init(chainID = id) is the same as function(id){make_init(chainID = id)} and is the function applied to list items in the list of length numParamList



MCMC_Schaefer <- tmbstan(obj = Model,
                         iter = 20000, # 1000 iterations
                         init = ParamList, # "par" uses defaults from model object, alternatively use "last.par.best"
                         chains = length(ParamList),
                         warmup = 2000,
                         thin = 10,
                         lower = lowbnd,
                         upper = uppbnd)

print(MCMC_Schaefer)
summary(MCMC_Schaefer)
plot(MCMC_Schaefer)
# pairs(MCMC_Schaefer)
 
# # To resolve this error:
# Error in .Call("is_Null_NS", ns) : 
#   "is_Null_NS" not resolved from current namespace (rstan)
# # I ran this line of code
# devtools::install_github("stan-dev/rstan/rstan/rstan@bugfix/issue%23535")





 
 
 
 # # pick h from prior distribution
 # temp_logith <- logit((h_steep - 0.2)/0.8);
 # logith_sample <- dnorm(temp_logith, 0.51, 4);
 # h_sample <- inv.logit(logith_sample)*0.8 + 0.2
# 
# Seeds <- seq(1, length(ParamList))
# MCMC <- sample_tmb(Model, iter = 1000, init = ParamList,  chains = length(ParamList), seeds = Seeds, 
#                    warmup = 100,
#                    thin = 0.1) # initial param vector wrong length
# 
# MCMC <- sample_tmb(Model, iter = 1000, init = ModelParameters,  chains = 4, seeds = Seeds, 
#                    warmup = 100,
#                    thin = 0.1) # initial param vector wrong length
# 
# MCMC <- sample_tmb(Model, iter = 1000, init = ModelParameters,  chains = 3, seeds = Seeds, 
#                    warmup = 100,
#                    thin = 0.1) # length init != number of chains
# 
# MCMC <- sample_tmb(Model, 
#                    iter = 10000, 
#                    init = NULL,  # I really should have something here, default uses model values
#                    chains = 3, 
#                    seeds = Seeds, # check that this is correct
#                    warmup = 1000, # too few and it causes issues
#                    thin = 1)  # needs to be a number >= 1
# 
# # warmup = burn-in
# # thin = thinning rate to apply
# # algorithm default = NUTS
# 
# # do I understand  chains vs. init? I think these things should be equal in length



