# Amanda Hart
# Homework 2 FISH 559 

setwd("/Users/ahart2/Research/TMB_Workshop_Fall2018/ClassWork/Hart_Homework2")

##### Read in data #####
data <- read.table("Home2.dat", skip=2)
colnames(data) <- c("Year", "Sex_M", "Sex_F") # numbers of males and females

########## PART 1 #################
# Show conventional estimator = MLE
###################################

##### Conventional estimate
sex_ratio1_1 <- sum(data[,"Sex_M"])/sum(data[,"Sex_M"] + data[,"Sex_F"])
sex_ratio1_1

##### Maximum likelihood estimate
LikeFunction <- function(log_p_est = 0.5){
  # Args:
    # log_p_est = log value for estimated sex ratio p (restricts this to + values for optimization purposes)
  
  # Returns: Negative log likelihood value
  
  p_est <- exp(log_p_est)
  # Read in data
  data <- read.table("Home2.dat", skip=2)
  colnames(data) <- c("Year", "Sex_M", "Sex_F") # numbers of males and females

  M_est <- rep(NA, length(data[,"Year"]))
  # Likelihood <- rep(NA, length(data[,"Year"]))
  for(yr in 1:length(data[,"Year"])){
    fishCount <- (data[yr,"Sex_M"] + data[yr,"Sex_F"])

    #Likelihood[yr] <- dbinom(data[yr,"Sex_M"], size = fishCount, prob = p_est) # data[yr,"Sex_M] successful trials given fishCount independent trials and p_est probability
    M_est[yr] <- rbinom(1, size=fishCount, prob=p_est) # returns 1 value from data of size fishCount, and probability of sucess = p_est
    # rbinom() will return a warning (NAs produced) if the solver provides a negative value for p_est (no negative probabilities possible)
  }
  print(M_est)
  # Calculate sum of squares
  Likelihood <- (data[,"Sex_M"] - M_est)^2
  print(Likelihood)
  
  # Calculate total & negative log likelihoods
  Likelihood <- sum(Likelihood)
  NLL <- -1*log(Likelihood)
  print(p_est)
  print(NLL)
  return(NLL)
}

LikeFunction(log_p_est=log(0.5)) # Test function works 

sex_ratio1_2 <- mle(LikeFunction, start=list(log_p_est=log(0.5)))
summary(sex_ratio1_2) 


#pdf of binom then find max likelihood


##### TO DO #####
# Check syntax in Part 1 for dbinom (I think this gives me the likelihood of picking data[,"Sex_M"] number of male fish given a total sample of M+F fish and a probability of getting a male = p_est)
    # also figure out why it gives me an error ( the size argument is an integer and it doesn't cause problems if there are no observations)
    # rbinom gives NAs
    # something is super wrong with rbinom() doing the likelihood with dbinom() didn't produce the same answer, maybe size argument is wrong
