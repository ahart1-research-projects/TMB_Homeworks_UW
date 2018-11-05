# Amanda Hart
# Homework 3 FISH 559 

# This script contains population dynamics functions

# Calculate Z_mortality (total mortality at age) #!!!!!! need to calculate by sex
CalcZMortality <- function(FishingMortality = 0,
                           M_par = NULL,
                           SelectivityAtAge = NULL){
  # Args
    # FishingMortality: Number between 0 and 1, corresponds to fishing mortality
    # M_par: Number between 0 and 1, represents natural mortality, no default
    # SelectivityAtAge: Matrix of proportions corresponding to selectivity at age and sex (columnames = "Selectivity_F" and "Selectivity_M"), no default
  
  # Returns: Matrix containing columns for female and male total mortality (Z) at age
  
  ZtempList <- NULL
  Z_mortality <- matrix(NA, length(Ages),2)
  colnames(Z_mortality) <- c("Z_M", "Z_F")
  sexes <- c("F", "M")
  
  for(isex in sexes){
    for(iage in 1:length(Ages)){
      Z_mortality[iage,paste("Z",isex, sep="_")] <- M_par + SelectivityAtAge[iage,paste("Selectivity",isex, sep="_")]*FishingMortality
    }
  }
  
  return(Z_mortality)
}


# Calculate NAA
CalcNAA <- function(Z_mortality = NULL,
                    Ages = NULL){
  # Args:
    # Z_mortality: Matrix containing columns for female and male total mortality (Z) at age, no default
    # Ages: Vector of ages
  
  # Return: Matrix of abundance at age by sex
  
  NAA <- matrix(NA, length(Ages), 2)
  colnames(NAA) <- c("NAA_M", "NAA_F")
  sexes <- c("F", "M")
  
  for(isex in sexes){
    for(iage in 1:length(Ages)){
      # NAA for age 0 (first age class)
      if(iage == 1){
        NAA[iage,paste("NAA",isex, sep="_")] <- 0.5
      } else if(iage == max(Ages)){ # NAA for last age
        NAA[iage,paste("NAA",isex, sep="_")] <- NAA[max(Ages)-1,paste("NAA",isex, sep="_")]*exp(-Z_mortality[iage-1,paste("Z",isex, sep="_")])
      } else{
        # NAA for intermediate ages
        NAA[iage,paste("NAA",isex, sep="_")] <- NAA[iage-1,paste("NAA",isex, sep="_")]*exp(-Z_mortality[iage-1,paste("Z",isex, sep="_")])
      }
    }
  }
  return(NAA)
}

  
# Calculate YPR !!!!!!! need to make sure isex matches for all 3
CalcYPR <- function(FishingMortality = 0,
                    WeightAtAge = NULL,
                    SelectivityAtAge = NULL,
                    NAA = NULL,
                    Z_mortality = NULL,
                    Ages = NULL){
  # Args
    # FishingMortality: Number corresponding to fishing mortalities, default = 0
    # WeightAtAge: Vector of numbers corresponding to weight at age, no default
    # SelectivityAtAge: Matrix of proportions corresponding to selectivity at age and sex (columnames = "Sel_F" and "Sel_M"), no default
    # NAA: Matrix of abundance at age by sex given FishingMortality, no default
    # Z_mortality: Matrix containing columns for female and male total mortality (Z) at age, no default
    # Ages: Vector of ages
  
  # Return: Number corresponding to Yield Per Recruit 
  
  sexes <- c("F", "M")
  # Calculate Yield Per Recruit (YPR)
  tempYPRList <- NULL
  for(isex in sexes){
    for(iage in 1:length(Ages)){
      tempYPR <- WeightAtAge[iage,paste("Weight",isex, sep="_")]*SelectivityAtAge[iage,paste("Selectivity",isex, sep="_")]*FishingMortality*NAA[iage,paste("NAA",isex, sep="_")]*(1-exp(-Z_mortality[iage,paste("Z",isex, sep="_")]))/Z_mortality[iage,paste("Z",isex, sep="_")]
      tempYPRList <- c(tempYPRList, tempYPR)
    }
  }
  #print(tempYPRList)
  YPR <- sum(tempYPRList)
  
  return(YPR)
}


# Calculate SBPR
CalcSBPR <- function(NAA = NULL,
                     FecundityAtAge = NULL,
                     Ages = NULL){
  # Args:
    # NAA: Matrix of abundance at age by sex given FishingMortality, no default
    # FecundityAtAge: Vector of numbers corresponding to fecundity at age, no default
    # Ages: Vector of ages
    
  # Return: Number for SBPR 
  
  SBPRList <- NULL
  for(iage in 1:length(Ages)){
    SBPRtemp <- NAA[iage,"NAA_F"]*FecundityAtAge[iage]
    SBPRList <- c(SBPRList, SBPRtemp)
  }
  SBPR <- sum(SBPRList)
  
  return(SBPR)
}


# Calculate Recruitment based on settings
CalcRecruit <- function(FishingMortality = 0,
                        StockRecruitForm = "BevertonHolt",
                        SBPR = NULL,
                        SBPRzero = NULL,
                        Steepness = NULL,
                        Rzero = NULL){
  # Args
    # FishingMortality: Number corresponding to fishing mortalities, default = 0
    # StockRecruitForm: String specifying form of stock recruit relationship, default = "BevertonHolt"
      # = "BevertonHolt"
      # = "Ricker"
      # = "PellaTomlinson
    # SBPR: Number corresponding to spawner biomass per recruit, no default
    # SBPRzero: Spawnerbiomass per recruit when fishing mortality = 0, default = 0
    # Steepness: Number for steepness, no default
    # Rzero: Number, recruitment at zero fishing mortality, no default
  
  # Return: Number, recruitment calculated based on the StockRecruitForm
  
  if(StockRecruitForm == "BevertonHolt"){
    Recruit <- Rzero*(0.8*SBPR*Steepness - SBPRzero*(0.2 - 0.2*Steepness))/((Steepness-0.2)*SBPR)
  } else if(StockRecruitForm == "Ricker"){
    Recruit <- 0.8*SBPRzero*Rzero*(log(SBPR) - log(SBPRzero) + log(Steepness/0.2)/0.8)/(log(Steepness/0.2)*SBPR)
  } else if(StockRecruitForm == "PellaTomlinson"){
    Recruit <- (Rzero*SBPRzero/SBPR)*((1 + (0.2 - 0.2^(gamma+1))/(Steepness-0.2) - SBPRzero*(0.2 - 0.2^(gamma+1))/(SBPR*(Steepness - 0.2)))^(1/gamma))
  }
  
  return(Recruit)
}


# Calculate Yield
CalcYield <- function(YPR = NULL,
                      Recruit = NULL){
  # Args: 
    # Recruit: Number corresponding to recruitment, no default
    # YPR: Number corresponding to Yield Per Recruit, no default
  
  # Return: Number corresponding to Yield
  
  Yield <- YPR*Recruits
  
  return(Yield)
}


# Calculate SSB
CalcSSB <- function(SBPR = NULL,
                    Recruit = NULL){
  # Args:
    # SBPR: Number corresponding to spawner biomass per recruit, no default
    # Recruit: Number corresponding to recruitment, no default
    
  # Returns: Number for spawning stock biomass
  
  SSB <- SBPR*Recruit
  
  return(SSB)
}

