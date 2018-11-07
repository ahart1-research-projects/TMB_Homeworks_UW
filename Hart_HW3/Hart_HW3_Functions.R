# Amanda Hart
# Homework 3 FISH 559 

# This script contains population dynamics functions

# Calculate Z_mortality (total mortality at age) #!!!!!! need to calculate by sex
CalcZMortality <- function(FishingMortality = 0,
                           M_par = NULL,
                           SelectivityAtAge = NULL,
                           Ages = NULL){
  # Args
    # FishingMortality: Number between 0 and 1, corresponds to fishing mortality
    # M_par: Number between 0 and 1, represents natural mortality, no default
    # SelectivityAtAge: Matrix of proportions corresponding to selectivity at age and sex (columnames = "Selectivity_F" and "Selectivity_M"), no default
    # Ages: Vector of fish ages
  
  # Returns: Matrix containing columns for female and male total mortality (Z) at age
  
  Z_mortality <- matrix(NA, length(Ages),2)
  colnames(Z_mortality) <- c("Z_M", "Z_F")
  sexes <- c("F", "M")
  
  for(isex in sexes){
    for(iage in 1:length(Ages)){
      # print(SelectivityAtAge[iage,paste("Selectivity",isex, sep="_")])
      Z_mortality[iage,paste("Z",isex, sep="_")] <- as.numeric(M_par + SelectivityAtAge[iage,paste("Selectivity",isex, sep="_")]*FishingMortality)
    }
  }
  
  return(Z_mortality)
}


# Calculate NAA
CalcNAA <- function(FishingMortality = 0,
                    Z_mortality = NULL,
                    Ages = NULL){
  # Args:
    # FishingMortality: Number between 0 and 1, corresponds to fishing mortality
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
      } else if(iage == length(Ages)){ # NAA for last age
        NAA[iage,paste("NAA",isex, sep="_")] <- NAA[iage-1,paste("NAA",isex, sep="_")]*exp(-Z_mortality[iage-1,paste("Z",isex, sep="_")])/(1-exp(-Z_mortality[iage,paste("Z",isex, sep="_")]))
      } else if(1 < iage & iage < length(Ages)){
        # NAA for intermediate ages
        NAA[iage,paste("NAA",isex, sep="_")] <- NAA[iage-1,paste("NAA",isex, sep="_")]*exp(-Z_mortality[iage-1,paste("Z",isex, sep="_")])
      }
    }
  }
  print(paste("NAA calculated at Fishing Mortality =", FishingMortality, sep=" "))
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
      #print(tempYPRList)
    }
  }
  #print(tempYPRList)
  YPR <- sum(tempYPRList)
  
  return(YPR)
}


# Calculate SBPR
CalcSBPR <- function(FishingMortality = 0,
                     NAA = NULL,
                     FecundityAtAge = NULL,
                     Ages = NULL){
  # Args:
    # FishingMortality: Number between 0 and 1, corresponds to fishing mortality
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
  
  print(paste("SBPR calculated at Fishing Mortality =", FishingMortality, sep=" "))
  return(SBPR)
}


# Calculate Recruitment based on settings
CalcRecruit <- function(FishingMortality = 0,
                        StockRecruitForm = "BevertonHolt",
                        SBPR = NULL,
                        SBPRzero = NULL,
                        Steepness = NULL,
                        Rzero = NULL,
                        gamma_par = 1){
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
    # gamma_par: Number, parameter for Pella-Tomlinson stock recruit form, default = 1
  
  # Return: Number, recruitment calculated based on the StockRecruitForm
  
  if(StockRecruitForm == "BevertonHolt"){
    # Recruit <- Rzero*(0.8*SBPR*Steepness - SBPRzero*(0.2 - 0.2*Steepness))/((Steepness-0.2)*SBPR)
    alpha <- SBPRzero*(0.2-0.2*Steepness)/(0.8*Steepness)
    beta <- (0.2 - Steepness)/(-0.8*Steepness*Rzero)
    Recruit <- (SBPR - alpha)/(beta*SBPR)
    print(paste("Doesn't require gamma_par parameter:", gamma_par, sep=" "))
  } else if(StockRecruitForm == "Ricker"){
    # Recruit <- 0.8*SBPRzero*Rzero*(log(SBPR) - log(SBPRzero) + log(Steepness/0.2)/0.8)/(log(Steepness/0.2)*SBPR)
    alpha <- 1/(SBPRzero*exp(-(log(Steepness/0.2)/0.8)))
    beta <- log(Steepness/0.2)/(0.8*SBPRzero*Rzero)
    Recruit <- (log(alpha) + log(SBPR*Recruit))/(beta*SBPR)
    print(paste("Doesn't require gamma_par parameter:", gamma_par, sep=" "))
  } else if(StockRecruitForm == "PellaTomlinson"){
    # Recruit <- (Rzero*SBPRzero/SBPR)*((1 + (0.2 - 0.2^(gamma_par+1))/(Steepness-0.2) - SBPRzero*(0.2 - 0.2^(gamma_par+1))/(SBPR*(Steepness - 0.2)))^(1/gamma_par))
    alpha <- 1/SBPRzero
    beta <- (h-0.2)/(0.2 - 0.2^(gamma_par + 1))
    Recruit <- (Rzero*SBPRzero/SBPR)*((1 + 1/beta - 1/(alpha*beta*SBPR))^(1/gamma_par))
  }
  
  return(Recruit)
}

# Calculate Rzero and SBPRzero
CalcSBPRzero <- function(M_par = NULL,
                         SelectivityAtAge = NULL,
                         WeightAtAge = NULL,
                         FecundityAtAge = NULL,
                         gamma_par = 1,
                         Ages = NULL,
                         Steepness = NULL,
                         Rzero = NULL,
                         StockRecruitForm = "BevertonHolt"){
  # Args: 
  # FishingMortality: Number corresponding to fishing mortalities, default = 0
  # M_par: Number between 0 and 1, represents natural mortality, no default
  # SelectivityAtAge: Matrix of proportions corresponding to selectivity at age and sex (columnames = "Selectivity_F" and "Selectivity_M"), no default
  # WeightAtAge: Vector of numbers corresponding to weight at age, no default
  # FecundityAtAge: Vector of numbers corresponding to fecundity at age, no default
  # gamma_par: Number, parameter for Pella-Tomlinson stock recruit relationship, default = 1
  # Ages: Vector of ages
  # Steepness: Number for steepness, no default
  # Rzero: Number, recruitment at zero fishing mortality, no default
  # StockRecruitForm: String specifying form of stock recruit relationship, default = "BevertonHolt"
  # = "BevertonHolt"
  # = "Ricker"
  # = "PellaTomlinson
  
  # Return: SBPRzero
  
  FishingMortality <- 0
  Z_mortality <- CalcZMortality(FishingMortality = 0, M_par = M_par, SelectivityAtAge = SelectivityAtAge, Ages = Ages)
  NAA <- CalcNAA(Z_mortality = Z_mortality, Ages = Ages)
  SBPRzero <- CalcSBPR(NAA = NAA, FecundityAtAge = FecundityAtAge, Ages = Ages)
  
  return(SBPRzero)
}


# Calculate Yield
CalcYield <- function(FishingMortality = 0,
                      M_par = NULL,
                      SelectivityAtAge = NULL,
                      WeightAtAge = NULL,
                      FecundityAtAge = NULL,
                      gamma_par = 1,
                      Ages = NULL,
                      Steepness = NULL,
                      Rzero = NULL,
                      StockRecruitForm = "BevertonHolt"){
  # Args: 
    # FishingMortality: Number corresponding to fishing mortalities, default = 0
    # M_par: Number between 0 and 1, represents natural mortality, no default
    # SelectivityAtAge: Matrix of proportions corresponding to selectivity at age and sex (columnames = "Selectivity_F" and "Selectivity_M"), no default
    # WeightAtAge: Vector of numbers corresponding to weight at age, no default
    # FecundityAtAge: Vector of numbers corresponding to fecundity at age, no default
    # gamma_par: Number, parameter for Pella-Tomlinson stock recruit relationship, default = 1
    # Ages: Vector of ages
    # Steepness: Number for steepness, no default
    # Rzero: Number, recruitment at zero fishing mortality, no default
    # StockRecruitForm: String specifying form of stock recruit relationship, default = "BevertonHolt"
      # = "BevertonHolt"
      # = "Ricker"
      # = "PellaTomlinson
    
  # Return: List containing calculated Yield with associated FishingMortality, Recruitment, SBPR, NAA, YPR, SBPRzero, and Z_mortality
  
  Ages <- seq(0,length(FecundityAtAge)-1)
  
  # First calculate SBPRzero
  SBPRzero <- CalcSBPRzero(M_par = M_par, SelectivityAtAge = SelectivityAtAge, WeightAtAge = WeightAtAge, FecundityAtAge = FecundityAtAge, 
                           gamma_par = gamma_par, Ages = Ages, Steepness = Steepness, Rzero = Rzero, StockRecruitForm = StockRecruitForm)
  print(SBPRzero)
  # Calculate total mortality at age for both sexes
  Z_mortality <- CalcZMortality(FishingMortality = FishingMortality, M_par = M_par, SelectivityAtAge = SelectivityAtAge, Ages = Ages)
  # Calculate abundance at age for both sexes
  NAA <- CalcNAA(FishingMortality = FishingMortality, Z_mortality = Z_mortality, Ages = Ages)
  # Calculate Yield Per Recruit
  YPR <- CalcYPR(FishingMortality = FishingMortality, WeightAtAge = WeightAtAge, SelectivityAtAge = SelectivityAtAge, NAA = NAA, Z_mortality = Z_mortality, Ages = Ages)
  # Calculate Spawner Biomas Per Recruit
  SBPR <- CalcSBPR(FishingMortality = FishingMortality, NAA = NAA, FecundityAtAge = FecundityAtAge, Ages = Ages)
  # Calculate recruitment
  # print(FishingMortality)
  # print(SBPR)
  # print(SBPRzero)
  # print(Rzero)
  # print(Steepness)
  Recruit <- CalcRecruit(FishingMortality = FishingMortality, StockRecruitForm = StockRecruitForm, SBPR = SBPR, SBPRzero = SBPRzero, Steepness = Steepness, Rzero = Rzero, gamma_par=gamma_par)
  # Calculate yield
  Yield <- YPR*Recruit
  
  Result <- NULL
  Result$Yield <- Yield
  Result$Recruit <- Recruit
  # print(Recruit)
  # print(SBPR)
  Result$FishingMortality <- FishingMortality
  Result$SBPR <- SBPR
  Result$NAA <- NAA
  Result$Z_mortality <- Z_mortality
  Result$YPR <- YPR
  Result$SBPRzero <- SBPRzero
  
  return(Result)
}


# Calculate SSB
CalcSSB <- function(FishingMortality = 0,
                    SBPR = NULL,
                    Recruit = NULL){
  # Args:
    # FishingMortality: Number between 0 and 1, corresponds to fishing mortality
    # SBPR: Number corresponding to spawner biomass per recruit, no default
    # Recruit: Number corresponding to recruitment, no default
    
  # Returns: Number for spawning stock biomass
  
  SSB <- SBPR*Recruit
  
  print(paste("SSB calculated at Fishing Mortality =", FishingMortality, sep=" "))
  return(SSB)
}






