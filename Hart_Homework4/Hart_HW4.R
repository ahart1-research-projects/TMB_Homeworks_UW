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
ModelData <- list(survey_vals = Survey_vals, survey_yr = Survey_yr, survey_SE = Survey_SE, 
                  catch_obs = Catch_vals, catch_yr = Catch_yr, rho = rho_val, w_dat = w_val, 
                  M_dat = M_val, Nproj = 0, catch_proj = 1000)

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

numParamList <- 1

make_init <- function(chainID = 1){
  list(dummy=0, Bzero=runif(1, 500, 1500), h_steep = runif(1, 0.2, 1), logF_y = rep(runif(1, -20, -0.01), length(Catch[,"Year"])-1)) 
}

ParamList <- lapply(1:numParamList, function(id) make_init(chainID = id)) # function(id) make_init(chainID = id) is the same as function(id){make_init(chainID = id)} and is the function applied to list items in the list of length numParamList



MCMC_Deriso <- tmbstan(obj = Model,
                         iter = 20000, # 2000 iterations
                         init = ParamList, # "par" uses defaults from model object, alternatively use "last.par.best"
                         chains = length(ParamList),
                         warmup = 2000,
                         thin = 10,
                         lower = lowbnd,
                         upper = uppbnd)

print(MCMC_Deriso)
  # summary(MCMC_Deriso)
  # plot(MCMC_Deriso)
  # traceplot(MCMC_Deriso, pars=names(Model$par))
  # stan_plot(MCMC_Deriso, pars=names(Model$par)) # hard to look at with scale difference
stan_trace(MCMC_Deriso, pars=names(Model$par)) # breaks it down by parameter =)
# stan_hist(MCMC_Deriso, pars=names(Model$par)) # similar to below but with bar plots
stan_dens(MCMC_Deriso, pars=names(Model$par)) # This is the best thing ever (as far as plots go, the distributions are wrong)
stan_scat(MCMC_Deriso, pars=c("Bzero", "h_steep")) # Don't know why you would do this but you can plot 2 against eachother
stan_diag(MCMC_Deriso) # requires more than 1 chain, what is this?
# stan_rhat(MCMC_Deriso, pars=names(Model$par)) # ick
stan_ess(MCMC_Deriso, pars=names(Model$par)) # plot effective sample size
stan_mcse(MCMC_Deriso, pars=names(Model$par))
# stan_ac(MCMC_Deriso, pars=names(Model$par)) # ick

# pairs(MCMC_Deriso)
 
# # To resolve this error:
# Error in .Call("is_Null_NS", ns) : 
#   "is_Null_NS" not resolved from current namespace (rstan)
# # I ran this line of code
# devtools::install_github("stan-dev/rstan/rstan/rstan@bugfix/issue%23535")

MCMC_Bzero <- extract(MCMC_Deriso)$Bzero
MCMC_h_steep <- extract(MCMC_Deriso)$h_steep
MCMC_logFy <- extract(MCMC_Deriso)$logF_y




######################################################################################################
########################### Part C with projections ##################################################
######################################################################################################
# Compile Cpp code
compile("Hart_HW4_Deriso.cpp") # file must be in working directory or provide full file path name
dyn.load(dynlib("Hart_HW4_Deriso"))

CheckTargetB <- function(catch_proj = 1000, Nproj = 1){
  # Note that only catch_proj is passed in as a parameter, all other data required to fit model is assumed to be global (i.e. you run the data lines at the start of this script)
  # projection over 10 years (Nproj = 10)
  
  # Return: calculated difference between biomass in last year and biomass target (here target = 0.5*Bzero)
  
  # Create list of data objects, names of list items must match DATA objects in Cpp code
  ModelData <- list(survey_vals = Survey_vals, survey_yr = Survey_yr, 
                    survey_SE = Survey_SE, catch_obs = Catch_vals, catch_yr = Catch_yr, 
                    rho = rho_val, w_dat = w_val, M_dat = M_val, Nproj = Nproj, catch_proj = catch_proj)
  
  CheckMeetRef <- NULL
  for(imcmc in 1:length(MCMC_Bzero)){
    # # Create list of parameters and provide initial values (may include parameter vector, e.g. param_vec = rep(0,5))
    # ModelParameters <- list(dummy=0, Bzero=MCMC_Bzero[imcmc], h_steep = MCMC_h_steep[imcmc], logF_y = MCMC_logFy[imcmc,]) 
    # 
    # 
    # # Use map function to specify which parameters to estimate, those that are not estimated are fixed at initial values and must have a factor(NA) in map list
    # ModelMap <- list(dummy = factor(NA)) # rep(factor(NA),5) for a parameter vector of length 5
    # 
    # Construct objective function to optimize based on data, parameters, and Cpp code
    Model <- MakeADFun(data = ModelData, parameters = ModelParameters, DLL="Hart_HW4_Deriso",silent=T,map = ModelMap) # silent=T silences a bunch of extra print statements
    
    Model$fn(c(MCMC_Bzero[imcmc], MCMC_h_steep[imcmc], MCMC_logFy[imcmc,]))
    
    # # Set bounds on different parameters, length of this vector must equal number of estimate parameters
    # # dummy=0, Bzero=220, h_steep = 0.5, F_y = c(rep(0.5, length(Catch[,"Year"])))
    # lowbnd <- c(500, 0.2, rep(-20,length(Catch[,"Year"])-1)) 
    # uppbnd <- c(1500, 1, rep(-0.01,length(Catch[,"Year"])-1)) # no bounds for mapped params!!!!!!!!!!
    # 
    # # Fit model to data using structure provided by MakeADFun() function call
    # fit <- nlminb(Model$par, Model$fn, Model$gr, control=list(rel.tol=1e-12,eval.max=100000,iter.max=1000), lower=lowbnd,upper=uppbnd) # notice the inclusion of the lower and upper bound vectors
    # 
    # ##### Fitted model results #####
    # # Best parameter estimates
    # best <- Model$env$last.par.best
    # # Report parameter estimates & std error
    # rep <- sdreport(Model)  
    
    # Get reported info & predicted data
    obj_fun <- Model$report()$obj_fun
    biomass_pred <- Model$report()$biomass
    catch_pred <- Model$report()$catch_pred
    recruitment <- Model$report()$recruitment
    
    ##### Check if biomass in final projection year >= 0.5*Bzero (half of K)
    CheckMeetRef[imcmc] <- biomass_pred[length(biomass_pred)] >= 0.5*best["Bzero"] 
  }
  
  PropMeetRef <- length(which(CheckMeetRef == TRUE))/length(CheckMeetRef)
  PropDiff <- 0.6 - PropMeetRef # we want biomass to be >= 0.5*Bzero (CheckMeetRef = TRUE) 60% of the time
  
  return(PropDiff)
}

# Find F for which biomass will recover to 0.5*Bzero during 10 year projection (set by Nproj)
  # If this returns an error that "f() values at end points not of opposite sign" then there is no catch within the interval for which the stock recovers within Nproj number of years
Result_projCatch_10yr <- uniroot(CheckTargetB, interval= c(0, 2000), Nproj = 10)




#####################################################################################################################################
###### Part D use Result_projCatch_10yr to do MCMC projections 
#####################################################################################################################################
# !!!!!!!!!! this next bit shoud take the MCMC param values and then do & save the biomass projections for plotting
Nproj <- 10
catch_proj <= Result_projCatch_10yr$root

# set up storage for biomass vectors
BiomassResults <- matrix(NA,nrow=length(MCMC_Bzero),ncol=length(Catch_yr)+Nproj) # Matrix with rows for MCMC calls, and columns = to number of years

ModelData <- list(survey_vals = Survey_vals, survey_yr = Survey_yr, survey_SE = Survey_SE, 
                  catch_obs = Catch_vals, catch_yr = Catch_yr, rho = rho_val, w_dat = w_val, 
                  M_dat = M_val, Nproj = 0, catch_proj = 1000)

for(imcmc in 1:length(MCMC_Bzero)){
  # Construct objective function to optimize based on data, parameters, and Cpp code
  Model <- MakeADFun(data = ModelData, parameters = ModelParameters, DLL="Hart_HW4_Deriso",silent=T,map = ModelMap) # silent=T silences a bunch of extra print statements
  
  Model$fn(c(MCMC_Bzero[imcmc], MCMC_h_steep[imcmc], MCMC_logFy[imcmc,]))
 
  # Get reported info & predicted data
  obj_fun <- Model$report()$obj_fun
  biomass_pred <- Model$report()$biomass
  catch_pred <- Model$report()$catch_pred

  # Save biomass vector
  BiomassResults[imcmc,] <- biomass_pred
}


# Calculate quantiles for BiomassResults & plot
StandardizedBiomass <- matrix(NA, ncol = ncol(BiomassResults), nrow = nrow(BiomassResults))
for(irow in 1:nrow(BiomassResults)){
  StandardizedBiomass[irow,] <- BiomassResults[irow,]/MCMC_Bzero[irow]
}

quant <- matrix(NA,nrow=5,ncol = length(Catch_yr)+Nproj)
for(iyr in 1:ncol(StandardizedBiomass)){
  quant[,iyr] <- quantile(StandardizedBiomass[,iyr], probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
}
Years <- 1:(length(Catch_yr)+Nproj)

# Plot of biomass posterior
ymax <- max(quant)
plot(Years,quant[3,],xlab="Year",ylab="Biomass",ylim=c(0,ymax),type='l', main="Deriso Biomass Timeseries")
xx <- c(Years,rev(Years))
yy <- c(quant[1,],rev(quant[5,]))
polygon(xx,yy,col="gray10")
xx <- c(Years,rev(Years))
yy <- c(quant[2,],rev(quant[4,]))
polygon(xx,yy,col="gray90")
lines(Years,quant[3,],lwd=3,lty=3)

##### Compare distributions of B_y/Bzero over time #######################
# first proj year
hist(StandardizedBiomass[,31], main="Histogram of B_y/K in first projection year (2002)",
     xlab = "B_y/K")

# last proj year
hist(StandardizedBiomass[,ncol(StandardizedBiomass)], main="Histogram of B_y/K in last projection year (2012)",
     xlab = "B_y/K")






