# Amanda Hart
# Homework 4 FISH 559

####################################################################################################################################
###### Fit without projection 
####################################################################################################################################
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
ModelData <- list(survey_vals = Survey_vals, survey_yr = Survey_yr, survey_SE = Survey_SE, catch_obs = Catch_vals, catch_yr = Catch_yr, Nproj = 0, catch_proj = 100)

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

# Set bounds on different parameters, length of this vector must equal number of estimate parameters
  # dummy=0, Bzero=220, logr_growth bounded from (-10, 0) = exp(4.539993e-05,1),    F_y bounded from logF_y (-20 to -0..01) = exp(2.061154e-09,0.9900498 ) 
lowbnd <- c(500,-10,rep(-20,length(Catch[,"Year"])-1)) # rep( 0.1, 5) for a parameter vector of length 5 with lower bound 0.1, syntax for upper bound is the same
uppbnd <- c(1500,0,rep(-0.01,length(Catch[,"Year"])-1))

# Fit model to data using structure provided by MakeADFun() function call

fit <- nlminb(Model$par, Model$fn, Model$gr, control=list(rel.tol=1e-12,eval.max=100000,iter.max=1000), lower=lowbnd,upper=uppbnd) # notice the inclusion of the lower and upper bound vectors

##### Fitted model results #####
# Best parameter estimates
best <- Model$env$last.par.best
print(best) 

# Report parameter estimates & std error
rep <- sdreport(Model)  
print(summary(rep))

# Get reported info & predicted data
# Predicted <- Model$report()$PredictedVariableNameHere 
obj_fun <- Model$report()$obj_fun
obj_fun
rNLL <- Model$report()$rNLL
rNLL
BzeroNLL <- Model$report()$BzeroNLL
BzeroNLL
biomass_pred <- Model$report()$biomass
biomass_pred
catch_pred <- Model$report()$catch_pred
catch_pred

####################################################################################################################################
###### Do MCMC 
####################################################################################################################################
library(tmbstan)

numParamList <- 3

make_init <- function(chainID = 1){
  list(dummy=0, Bzero=runif(1, 500, 1500), logr_growth = runif(1, -10, 0), logF_y = rep(runif(1, -20, -0.01), length(Catch[,"Year"])-1)) 
}

ParamList <- lapply(1:numParamList, function(id) make_init(chainID = id)) # function(id) make_init(chainID = id) is the same as function(id){make_init(chainID = id)} and is the function applied to list items in the list of length numParamList



MCMC_Schaefer <- tmbstan(obj = Model,
                       iter = 20000, # 2000 iterations
                       init = ParamList, # "par" uses defaults from model object, alternatively use "last.par.best"
                       chains = length(ParamList),
                       warmup = 2000,
                       thin = 100,
                       lower = lowbnd,
                       upper = uppbnd)

##### Store MCMC results for projections
MCMC_Bzero <- extract(MCMC_Schaefer)$Bzero
MCMC_logr_growth <- extract(MCMC_Schaefer)$logr_growth
MCMC_logFy <- extract(MCMC_Schaefer)$logF_y

##### MCMC plots
print(MCMC_Schaefer)
# summary(MCMC_Deriso)
# plot(MCMC_Deriso)
# traceplot(MCMC_Deriso, pars=names(Model$par))
# stan_plot(MCMC_Deriso, pars=names(Model$par)) # hard to look at with scale difference
stan_trace(MCMC_Schaefer, pars=names(Model$par)) # breaks it down by parameter =)
# stan_hist(MCMC_Deriso, pars=names(Model$par)) # similar to below but with bar plots
stan_dens(MCMC_Schaefer, pars=names(Model$par)) # This is the best thing ever (as far as plots go, the distributions are wrong)
stan_scat(MCMC_Schaefer, pars=c("Bzero", "r_growth")) # Don't know why you would do this but you can plot 2 against eachother
stan_diag(MCMC_Schaefer) # requires more than 1 chain, what is this?
# stan_rhat(MCMC_Deriso, pars=names(Model$par)) # ick
stan_ess(MCMC_Deriso, pars=names(Model$par)) # plot effective sample size
stan_mcse(MCMC_Deriso, pars=names(Model$par))
# stan_ac(MCMC_Deriso, pars=names(Model$par)) # ick

# # To resolve this error:
# Error in .Call("is_Null_NS", ns) : 
#   "is_Null_NS" not resolved from current namespace (rstan)
# # I ran this line of code
# devtools::install_github("stan-dev/rstan/rstan/rstan@bugfix/issue%23535")





#####################################################################################################################################
###### Part C with projections 
#####################################################################################################################################
# Compile Cpp code
compile("Hart_HW4_Schaefer.cpp") 
dyn.load(dynlib("Hart_HW4_Schaefer"))

CheckTargetB <- function(catch_proj = 100, Nproj = 1){
  # Note that only catch_proj is passed in as a parameter, all other data required to fit model is assumed to be global (i.e. you run the data lines at the start of this script)
  # projection over 10 years (Nproj = 10)
  
  # Return: calculated difference between biomass in last year and biomass target (here target = 0.5*Bzero)
  
  # Create list of data objects, names of list items must match DATA objects in Cpp code
  ModelData <- list(survey_vals = Survey_vals, survey_yr = Survey_yr, survey_SE = Survey_SE, 
                    catch_obs = Catch_vals, catch_yr = Catch_yr, Nproj = Nproj, catch_proj = catch_proj)
  
  CheckMeetRef <- NULL
  for(imcmc in 1:length(MCMC_Bzero)){
    # # Create list of parameters and provide initial values (may include parameter vector, e.g. param_vec = rep(0,5))
    # ModelParameters <- list(dummy=0, Bzero=MCMC_Bzero[imcmc], logr_growth = MCMC_logr_growth[imcmc], logF_y = MCMC_logFy[imcmc,]) 
    # # log_r is sensitive to starting value near the bounds
    # 
    # # Use map function to specify which parameters to estimate, those that are not estimated are fixed at initial values and must have a factor(NA) in map list
    # ModelMap <- list(dummy = factor(NA)) # rep(factor(NA),5) for a parameter vector of length 5
    # 
    # Construct objective function to optimize based on data, parameters, and Cpp code
    Model <- MakeADFun(data = ModelData, parameters = ModelParameters, DLL="Hart_HW4_Schaefer",silent=T,map = ModelMap) # silent=T silences a bunch of extra print statements
    
    Model$fn(c(MCMC_Bzero[imcmc],  MCMC_logr_growth[imcmc], MCMC_logFy[imcmc,]))
    # # Set bounds on different parameters, length of this vector must equal number of estimate parameters
    # # dummy=0, Bzero=220, logr_growth bounded from (-10, 0) = exp(4.539993e-05,1),    F_y bounded from logF_y (-20 to -0..01) = exp(2.061154e-09,0.9900498 ) 
    # lowbnd <- c(500,-10,rep(-20,length(Catch[,"Year"])-1)) # rep( 0.1, 5) for a parameter vector of length 5 with lower bound 0.1, syntax for upper bound is the same
    # uppbnd <- c(1500,0,rep(-0.01,length(Catch[,"Year"])-1))
    # 
    # # Fit model to data using structure provided by MakeADFun() function call
    # fit <- nlminb(Model$par, Model$fn, Model$gr, control=list(rel.tol=1e-12,eval.max=100000,iter.max=1000), lower=lowbnd,upper=uppbnd) # notice the inclusion of the lower and upper bound vectors
    # 
    # ##### Fitted model results #####
    # # Best parameter estimates
    # best <- Model$env$last.par.best
    # # Report parameter estimates & std error
    # rep <- sdreport(Model)  
    # 
    # Get reported info & predicted data
    obj_fun <- Model$report()$obj_fun
    biomass_pred <- Model$report()$biomass
    catch_pred <- Model$report()$catch_pred
    
    ##### Check if biomass in final projection year >= 0.5*Bzero (half of K)
    CheckMeetRef[imcmc] <- biomass_pred[length(biomass_pred)] >= 0.5*best["Bzero"] 
  }
  
  PropMeetRef <- length(which(CheckMeetRef == TRUE))/length(CheckMeetRef)
  PropDiff <- 0.6 - PropMeetRef # we want biomass to be >= 0.5*Bzero (CheckMeetRef = TRUE) 60% of the time
  
  return(PropDiff)
}

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
                  catch_obs = Catch_vals, catch_yr = Catch_yr, Nproj = Nproj, catch_proj = catch_proj)

for(imcmc in 1:length(MCMC_Bzero)){
  # Construct objective function to optimize based on data, parameters, and Cpp code
  Model <- MakeADFun(data = ModelData, parameters = ModelParameters, DLL="Hart_HW4_Schaefer",silent=T,map = ModelMap) # silent=T silences a bunch of extra print statements
  
  Model$fn(c(MCMC_Bzero[imcmc],  MCMC_logr_growth[imcmc], MCMC_logFy[imcmc,]))

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
plot(Years,quant[3,],xlab="Year",ylab="Biomass",ylim=c(0,ymax),type='l', main="Schaefer Biomass Timeseries")
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
