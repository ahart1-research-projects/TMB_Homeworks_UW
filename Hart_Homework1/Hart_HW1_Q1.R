# Amanda Hart
# Homework 1 FISH 559 Question 1

##### Read in data #####
setwd("/Users/ahart2/Research/TMB_Workshop_Fall2018/ClassWork/Hart_Homework1")

data <- read.table("HOME1A.TXT")
colnames(data) <- c("StreamID", "Ignore", "StreamData") 


##### Fit model with LME #####
library(nlme)
StreamData <- data[,"StreamData"]
StreamID <- data[,"StreamID"]
check_data <- list(StreamData, StreamID)
check <- lme(StreamData ~ 1, data = check_data, random = ~ 1 | StreamID, method = "ML")
# random = ~1|stream
summary(check)
coef(check) # prints random effect values

LowerCI_LME <- intervals(check)[1]$fixed[1]
UpperCI_LME <- intervals(check)[1]$fixed[3]

# Correct parameter values
# beta 78.62872
# between sigma_b 14.34159
# within sigma_e 6.083881
# log likelihood -66.65819


CalcLikelihood <- function(beta=NULL, sigma_b=NULL, sigma_e=NULL){
  # Args:
    # beta = value for mean stream density
    # sigma_b = value for between stream variation in avg density
    # sigma_e = value for residual variation in stream density
  
  # Defined here
    # data = data table containing a column labeld "StreamID", and "StreamData"
      data <- read.table("HOME1A.TXT")
      colnames(data) <- c("StreamID", "Ignore", "StreamData")
    # nobs = number of observations for each stream (assumed to be equal for all streams as is the case with this dataset)
      nobs <- 3
    # integrate_Lower = value for lower bound of integration
      integrate_Lower <- -5
    # integrate_Upper = value for upper bound of integration
      integrate_Upper <- 5
    # Ncalls = value for number of function calls, default = 1
      Ncalls <- 500 #????? I don't know where this number is supposed to come from so I picked a number which will make h small
      
  # Returns: negative log likelihood
  
  # read data inside functioin, args=param
  
  # Info from data
  Nstreams <- length(unique(data[,"StreamID"]))
  StreamIDList <- c(unique(data[,"StreamID"]))
  
  # Random effect (bi) steps
  b_step <- seq(from=integrate_Lower, to=integrate_Upper, by=(integrate_Upper-integrate_Lower)/(Ncalls-1))
  
  # Calculations for simpson's rule
  simpson_h <- (integrate_Upper - integrate_Lower)/(Ncalls - 1) # (bmax-bmin)/(ncalls-1)
  SimpsonInterval <- c(integrate_Lower, integrate_Lower+simpson_h, integrate_Upper) # contains values of random effect bi to evaluate (between -5 and 5)
  
  # Start likelihood calculations
  Likelihood <- rep(NA,Nstreams)
  for(istream in 1:Nstreams){ # loop over streams
    SimpsonsLike <- rep(NA,length(b_step)) # rep(NA,length(SimpsonInterval))
    IstreamRows <- which(data[,"StreamID"]==StreamIDList[istream]) # which rows contain observations from this stream
    
    for(ibval in 1:length(b_step)){ # loop over each value of bi
      prob_b_given_sigmab <- 1
      for(iobs in 1:nobs){ # loop over observations in istream (could change this so nobs changes with stream if number of obs at each stream changes)
        temp_prob <- exp(-1*((data[IstreamRows[iobs],"StreamData"] - beta - b_step[ibval]*sigma_b)^2)/(2*(sigma_e^2)))/((sigma_e)*((2*pi)^0.5))
        # print(temp_prob) #!!! these are all super small numbers
        prob_b_given_sigmab <- prob_b_given_sigmab*temp_prob # multiply all the temp_probs together for a given stream(istream)
        # print(prob_b_given_sigmab) #!!! rounds to zero
      }
      
      # Calculate function value at all values of bi and store in SimpsonsLike vector
      SimpsonsLike[ibval] = (exp(-1*(b_step[ibval]^2)/2)/((2*pi)^0.5))*prob_b_given_sigmab
    }
    
    # !!! this is still super wrong, it gives tiny values that tend to round to zero
    # Do multiplication for simpson's rule
    SimpsonsLike[1] <- SimpsonsLike[1] # first and last multiplied by 1
    SimpsonsLike[length(SimpsonsLike)] <- SimpsonsLike[length(SimpsonsLike)]
    SimpsonsLike[seq(2, length(SimpsonsLike)-1, by=2)] <- 4*SimpsonsLike[seq(2, length(SimpsonsLike)-1, by=2)] # multiply evens by 4 except the last number
    SimpsonsLike[seq(3, length(SimpsonsLike)-1, by=2)] <- 2*SimpsonsLike[seq(3, length(SimpsonsLike)-1, by=2)] # multiply odds by 2 except the last number
    SimpsonsLike <- SimpsonsLike/3
    
    #print(SimpsonsLike)
    # print(sum(SimpsonsLike))
    
    # Calculate likelihood component for single stream (integrate over all possible values of bi)
    Like_1_stream <- simpson_h*(sum(SimpsonsLike)) 
    
    # Calculate likelihood across all streams
    Likelihood[istream] <- Like_1_stream # store each stream likelihood in a Likelihood vector
    print(Likelihood)
  }
  # Calculate total likelihood
  TotLikelihood <- sum(log(Likelihood)) # sum of log value = product of non-log
  print(TotLikelihood)
  # Calculate negative log likelihood
  NLL <- -1*(TotLikelihood)
  print(NLL)
}



##### Calculate parameters #####
library(stats4)
param_calc <- mle(CalcLikelihood, start=list(beta=78, sigma_b=15, sigma_e=6)) # then vary beta but fix other things
summary(param_calc) # correct parameter estimates to 2 decimal points
BestNLL <- (summary(param_calc)@m2logL)/2 
BestBeta <- coef(summary(param_calc))[1]

##### Calculate likelihood with given parameters #####
CalcLikelihood(beta = 78.626037, sigma_b = 14.347004, sigma_e = 6.084041)


##### Create likelihood profile #####
# Fix beta at different values and estimate negative log likelihood
beta_list <- seq(60,100,by=1)
for(isim in 1:length(beta_list)){
  param_calc <- mle(CalcLikelihood, start=list(sigma_b=15, sigma_e=6), fixed=list(beta=beta_list[isim])) 
  LogLikelihood[isim] <- (summary(param_calc)@m2logL)/2 
}

# Standardize negative log likelihood
LogLikelihood <- LogLikelihood-BestNLL 

NLLPlotx <- beta_list
NLLPloty <- LogLikelihood

# 95% CI bounds for likelihood 
LowerCI <- BestNLL - (1.92) 
UpperCI <- BestNLL + 1.92
lowerCI <- lm(NLLPlotx[6:7]~NLLPloty[6:7])
beta_lowerCI <- coef(lowerCI)[2]*1.92+coef(lowerCI)[1] #= beta*coef(lowerCI)[2]
upperCI <- lm(NLLPlotx[32:33]~NLLPloty[32:33])
beta_upperCI <- coef(upperCI)[2]*1.92+coef(upperCI)[1]

# Plot Loglikelihood and CIs
plot(x=NLLPlotx, y=NLLPloty, type="l",
     main="Negative Log Likelihood Profile",
     xlab="Beta",
     ylab="Negative log likelihood")
#lines(NLLPlotx[which(NLLPlotx >= beta_lowerCI & NLLPlotx <= beta_upperCI)], NLLPloty[which(NLLPlotx >= beta_lowerCI & NLLPlotx <= beta_upperCI)], col="red", lwd=5)
abline(v=LowerCI_LME, col="blue", lwd=2)
abline(v=UpperCI_LME, col="blue", lwd=2)
abline(v=beta_lowerCI, col="red", lwd=2)
abline(v=beta_upperCI, col="red", lwd=2)
legend("topright", legend=c("LME 95% CI", "Likelihood 95% CI"), fill=c("blue", "red"))

LowerCI_LME
beta_lowerCI
UpperCI_LME
beta_upperCI




# Questions I have ??????
# should I really standardize log likelihood for profile so -NLL - -BestNLL?
# how do you find 95% CI, what am I comparing these to?
# how do I construct a likelihood profile? Look back at TMB workshop example
# fix github so my TMB templates are online
# fix this so best NLL-1.92 not best beta - 1.92
   # to plot need to figure out what beta value corresponds with likelihood














##### This is the original version of the function but it makes it difficult to run through mle() so above it is changed to only have parameters as arguments
CalcLikelihood <- function(data=NULL, nobs=NULL, integrate_Lower=NULL, integrate_Upper=NULL, Ncalls=1, beta=NULL, sigma_b=NULL, sigma_e=NULL){
  # Args:
  # data = data table containing a column labeld "StreamID", and "StreamData"
  # nobs = number of observations for each stream (assumed to be equal for all streams as is the case with this dataset)
  # integrate_Lower = value for lower bound of integration
  # integrate_Upper = value for upper bound of integration
  # Ncalls = value for number of function calls, default = 1
  # beta = value for mean stream density
  # sigma_b = value for between stream variation in avg density
  # sigma_e = value for residual variation in stream density
  
  # Returns: likelihood
  
  # Info from data
  Nstreams <- length(unique(data[,"StreamID"]))
  StreamIDList <- c(unique(data[,"StreamID"]))
  
  # Calculations for simpson's rule
  simpson_h <- (integrate_Upper - integrate_Lower)/(Ncalls - 1) # (bmax-bmin)/(ncalls-1)
  SimpsonInterval <- c(integrate_Lower, integrate_Lower+simpson_h, integrate_Upper)
  
  # Start likelihood calculations
  Likelihood <- NULL
  for(istream in 1:Nstreams){ # loop over streams
    SimpsonsLike <- rep(NA,length(SimpsonInterval))
    IstreamRows <- which(data[,"StreamID"]==StreamIDList[istream]) # which rows contain observations from this stream
    
    for(isimpson in 1:length(SimpsonInterval)){ # calculate values at three intervals for simpson's rule integration
      prob_b_given_sigmab <- 1
      
      for(iobs in 1:nobs){ # loop over observations in istream (could change this so nobs changes with stream if number of obs at each stream changes)
        temp_prob <- exp(-1*(data[IstreamRows[iobs],"StreamData"] - beta - SimpsonInterval[isimpson]*sigma_b)^2/(2*sigma_e^2))/(sigma_e*((2*pi)^0.5))
        prob_b_given_sigmab <- prob_b_given_sigmab*temp_prob # multiply all the temp_probs together for a given stream(istream)
      }
      # Calculate Simpson's interval value
      SimpsonsLike[isimpson] = (exp(-1*(SimpsonInterval[isimpson]^2)/2)/((2*pi)^0.5))*prob_b_given_sigmab
      
      # Calculate likelihood component for single stream
      Like_1_stream <- simpson_h*(1/3*SimpsonsLike[1] + 4/3*SimpsonsLike[2] + 1/3*SimpsonsLike[3])
    }
    # Calculate likelihood across all streams
    Likelihood <- c(Likelihood,Like_1_stream) # store each stream likelihood in a Likelihood vector
  }
}


