# Amanda Hart
# Homework 3 FISH 559 


setwd("/Users/ahart2/Research/TMB_Homeworks_UW/Hart_HW3")

data <- read.table("HOME3.TXT", header=TRUE)
colnames(data) <- c("Age", "Fecundity", "Weight_F", "Weight_M", "Selectivity_F", "Selectivity_M")





# Put a function here
HW3Function <- function(FishingMortality = seq(0,1, by=0.01),
                        M_par = NULL,
                        FecundityAtAge = NULL,
                        WeightAtAgeMatrix = NULL,
                        SelectivityAtAgeMatrix = NULL,
                        Rzero = NULL,
                        Steepness = NULL,
                        StockRecruitForm = "BevertonHolt",
                        gamma_par = NULL,
                        DoPlot = TRUE){
  # Args
    # FishingMortality: Vector of fishing mortalities to assess
    # M_par: Number between 0 and 1, represents natural mortality, no default
    # FecundityAtAge: Vector of numbers corresponding to fecundity at age, no default
    # WeightAtAgeMatrix: Matrix of numbers corresponding to weight at age by sex (colnames = "Weight_F", "Weight_M"), no default
    # SelectivityAtAgeMatrix: Matrix of proportions corresponding to selectivity at age by sex (colnames = "Selectivity_F", "Selectivity_M"), no default
    # Rzero: Number, recruitment at zero fishing mortality, no default
    # Steepness: Number, no default
    # StockRecruitForm: String specifying form of stock recruit relationship, default = "BevertonHolt"
      # = "BevertonHolt"
      # = "Ricker"
      # = "PellaTomlinson
    # gamma_par: number, no default
    # DoPlot: Boolean, specifies if plots are produced, default = TRUE
  
  # Return:
  
  Ages <- seq(0,length(FecundityAtAge)-1)
  

  
  #################################################################################################################
  # Storage vector
  SBPR <- rep(NA, length(FishingMortality))
  Recruit <- rep(NA, length(FishingMortality))
  Yield <- rep(NA, length(FishingMortality))
  SSB <- rep(NA, length(FishingMortality))
  
  # Loop over fishing efforts
  for(ifish in 1:length(FishingMortality)){
    # # Calculate total mortality at age for both sexes
    # Z_mortality <- CalcZMortality(FishingMortality = FishingMortality[ifish], M_par = M_par, SelectivityAtAge = SelectivityAtAgeMatrix)
    # # Calculate abundance at age for both sexes
    # NAA <- CalcNAA(Z_mortality = Z_mortality, Ages = Ages)
    # # Calculate Yield Per Recruit
    # YPR[ifish] <- CalcYPR(FishingMortality = FishingMortality[ifish], WeightAtAge = WeightAtAgeMatrix, SelectivityAtAge = SelectivityAtAgeMatrix, NAA = NAA, Z_mortality = Z_mortality, Ages = Ages)
    # # Calculate Spawner Biomas Per Recruit
    # SBPR[ifish] <- CalcSBPR(NAA = NAA, FecundityAtAge = FecundityAtAge, Ages = Ages)
    # # Calculate recruitment
    # Recruit[ifish] <- CalcRecruit(FishingMortality = FishingMortality[ifish], StockRecruitForm = "BevertonHolt", SBPR = SBPR[ifish], SBPRzero = SBPRzero, Steepness = Steepness, Rzero = Rzero)
    
    # Calculate yield & store corresponding SBPR & Recruitment
    Result <- CalcYield(FishingMortality = FishingMortality[ifish], M_par = M_par, SelectivityAtAge = SelectivityAtAgeMatrix, WeightAtAge = WeightAtAgeMatrix, FecundityAtAge = FecundityAtAge, 
                        gamma_par = gamma_par, Ages = Ages, Steepness = Steepness, Rzero = Rzero, StockRecruitForm = StockRecruitForm)
    SBPR[ifish] <- Result$SBPR
    #print(Result$Z_mortality)
    Recruit[ifish] <- Result$Recruit
    # print(SBPR)
    # print(SBPRzero)
    # print(Rzero)
    # print(Recruit)
    Yield[ifish] <- Result$Yield
    NAA <- Result$NAA
    #print(NAA)
    #print(Result$YPR)
    # Calculate Spawning Stock Biomass
    print("here")
    print(Result$SBPR)
    print(Result$Recruit)
    SSB[ifish] <- CalcSSB(SBPR = Result$SBPR, Recruit = Result$Recruit, FishingMortality = FishingMortality[ifish])
    print(SSB)
  }
  
  # Calculate MSY and FMSY
  FMSY <- uniroot(CalcFMSY,interval=c(0,1), StockRecruitForm = StockRecruitForm, gamma_par = gamma_par, Steepness = Steepness, Rzero = Rzero) #, Steepness = Steepness)
  MSYResults <- CalcYield(FishingMortality = FMSY$root, M_par = M_par, SelectivityAtAge = SelectivityAtAgeMatrix, WeightAtAge = WeightAtAgeMatrix, FecundityAtAge = FecundityAtAge,
                          gamma_par = gamma_par, Ages = Ages, Steepness = Steepness, Rzero = Rzero, StockRecruitForm = StockRecruitForm)
  MSY <- MSYResults$Yield
  print(MSY)
  print(max(Result$Yield))

  # Calculate Fcrash
  # Fcrash <- uniroot(CalcFcrash, interval=c(0,1), StockRecruitForm = StockRecruitForm, Steepness = Steepness, Rzero = Rzero)
  # Fcrash <- Fcrash$root
  print("FMSY")
  print(FMSY$root)
  print("MSY")
  print(MSY)
  # print("Fcrash")
  # print(Fcrash)

  
  if(DoPlot == TRUE){
    par(mfrow=c(2,2)) # this makes a 4X4 plot
    
    # Subset data
    plotRecruit <- Recruit[which(Recruit >=0)]
    plotSSB <- SSB[which(Recruit >= 0)]
    plotYield <- Yield[which(Recruit >= 0)]
    plotFishingMortality <- FishingMortality[which(Recruit >= 0)]
    
    # Plot
    plot(x = plotSSB, y = plotRecruit, main = "Recruitment vs. SSB", xlab = "SSB", ylab = "Recruitment")
    xreplacement <- seq(0:max(SSB))
    yreplacement <- seq(0:max(SSB))
    lines(x=xreplacement, y=yreplacement)
    plot(x = plotSSB, y = plotYield, main = "Yield vs. SSB", xlab = "SSB", ylab = "Yield")
    abline(h = MSY, col = "red")
    plot(x = plotFishingMortality, y = plotYield, main = "Yield vs. Fishing Mortality", xlab = "Fishing Mortality", ylab = "Yield")
    abline(v = FMSY$root, col="red")
    # abline(v = Fcrash, col = "green")
    abline(h = MSY, col = "blue")
    legend("topright", legend=c("FMSY", "MSY", "Fcrash"), fill=c("red", "blue", "green"))
    plot(x = plotFishingMortality, y = plotSSB, main = "SSB vs. Fishing Mortality", xlab = "Fishing Mortality", ylab = "SSB")
    mtext(StockRecruitForm, outer = TRUE, cex = 1.5, side = 3, line = -1.5)
  }
}

# Test with single value of F
HW3Function(FishingMortality = 0.1, M_par = 0.15, FecundityAtAge = data$Fecundity, WeightAtAgeMatrix = WeightAtAgeMatrix,
            SelectivityAtAgeMatrix = SelectivityAtAgeMatrix, Rzero = 1, Steepness = 0.5, StockRecruitForm = "BevertonHolt", 
            gamma_par = 0.2, DoPlot = FALSE)

# Run across a range of F
HW3Function(FishingMortality = seq(0,1,by=0.001), M_par = 0.2, FecundityAtAge = data$Fecundity, 
            WeightAtAgeMatrix = WeightAtAgeMatrix, SelectivityAtAgeMatrix = SelectivityAtAgeMatrix, 
            Rzero = 1, Steepness = 0.5, StockRecruitForm = "Ricker", gamma_par = 0.2, 
            DoPlot = TRUE)


CalcFMSY <- function(F_par, StockRecruitForm = "BevertonHolt", gamma_par = 1, Steepness = NULL, Rzero = NULL){
  # Args:
  # F_par is fishing mortality to vary using uniroot
  # StockRecruitForm is the stock recruit relationship
  # gamma_par is parameter needed if StockRecruitForm = "PellaTomlinson", default = 1
  # Steepness of stock recruit relationship
  # Rzero = recruitment at 0 fishing
  
  Result1 <- CalcYield(FishingMortality = (F_par + 0.001), M_par = 0.15, SelectivityAtAge = SelectivityAtAgeMatrix, WeightAtAge = WeightAtAgeMatrix, FecundityAtAge = data$Fecundity, 
                       gamma_par = gamma_par, Ages = Ages, Steepness = Steepness, Rzero = Rzero, StockRecruitForm = StockRecruitForm)
  
  Result2 <- CalcYield(FishingMortality = (F_par - 0.001), M_par = 0.15, SelectivityAtAge = SelectivityAtAgeMatrix, WeightAtAge = WeightAtAgeMatrix, FecundityAtAge = data$Fecundity, 
                       gamma_par = gamma_par, Ages = Ages, Steepness = Steepness, Rzero = Rzero, StockRecruitForm = StockRecruitForm)
  print(Result1$Yield)
  print(Result2$Yield)
  
  Difference <- (Result2$Yield - Result1$Yield)/0.002
  return(Difference)
}

uniroot(CalcFMSY,interval=c(0,1), StockRecruitForm = "PellaTomlinson", gamma_par = 1, Steepness = 0.5, Rzero=1)

#!!!!!!!!!!!! this doesn't work for Ricker or BevertonHolt
CalcFcrash <- function(F_par, StockRecruitForm = "BevertonHolt", Steepness = NULL, Rzero = NULL, gamma_par = 1){
  # Args:
  # F_par is fishing mortality to vary using uniroot
  # StockRecruitForm is the stock recruit relationship
  # Steepness is slope of the stock-recruit relationship
  # Rzero is recruitment under no fishing
  # gamma_par: optional parameter passed when 
  
  if(StockRecruitForm == "BevertonHolt"){ # This works
    Result <- CalcYield(FishingMortality = (F_par), M_par = 0.15, SelectivityAtAge = SelectivityAtAgeMatrix, WeightAtAge = WeightAtAgeMatrix, FecundityAtAge = data$Fecundity, 
                        gamma_par = gamma_par, Ages = Ages, Steepness = Steepness, Rzero = Rzero, StockRecruitForm = "BevertonHolt")
    alpha <- Result$SBPRzero*(0.2-0.2*Steepness)/(0.8*Steepness)
    
    print(Result$SBPR)
    print(alpha)
    
    Difference <- Result$SBPR - alpha
  } else if (StockRecruitForm == "Ricker"){ # this is set up such that calculating at F=1 and F=0 are both negative differences so there is no zero between them, rethink this
    Result <- CalcYield(FishingMortality = (F_par), M_par = 0.15, SelectivityAtAge = SelectivityAtAgeMatrix, WeightAtAge = WeightAtAgeMatrix, FecundityAtAge = data$Fecundity, 
                        gamma_par = gamma_par, Ages = Ages, Steepness = Steepness, Rzero = Rzero, StockRecruitForm = "Ricker")
    alpha <- Rzero/(Result$SBPRzero*exp(-(log(Steepness/0.2)/0.8)))
    
    print(Result$SBPR)
    print(alpha)
    
    Difference <- Result$SBPR - alpha
  } else if (StockRecruitForm == "PellaTomlinson"){
    Result <- CalcYield(FishingMortality = (F_par), M_par = 0.15, SelectivityAtAge = SelectivityAtAgeMatrix, WeightAtAge = WeightAtAgeMatrix, FecundityAtAge = data$Fecundity, 
                        gamma_par = gamma_par, Ages = Ages, Steepness = Steepness, Rzero = Rzero, StockRecruitForm = "Ricker")
    alpha <- 1/Result$SBPRzero
    
    print(Result$SBPR)
    print(alpha)
    
    Difference <- Result$SBPR - alpha
  }
  
  return(Difference)
}

uniroot(CalcFcrash, interval=c(0,1), StockRecruitForm = "BevertonHolt", Steepness = 0.5, Rzero = 1, gamma_par = 1)












# max yield 0.02662691 as calculated


# Test CalcYield Equation
ResultTest <- CalcYield(FishingMortality = 0.5, M_par = M_par, SelectivityAtAge = SelectivityAtAgeMatrix, WeightAtAge = WeightAtAgeMatrix, FecundityAtAge = FecundityAtAge, 
                        Ages = Ages, Steepness = Steepness, Rzero = Rzero, StockRecruitForm = "BevertonHolt")


SelectivityAtAgeMatrix <- cbind(data$Selectivity_F, data$Selectivity_M)
colnames(SelectivityAtAgeMatrix) <- c("Selectivity_F", "Selectivity_M")

WeightAtAgeMatrix <- cbind(data$Weight_F, data$Weight_M)
colnames(WeightAtAgeMatrix) <- c("Weight_F", "Weight_M")


# Test
M_par <- 0.15
Steepness  <- 0.5
Rzero <- 1






# CalcMSY <- function(i){
#   # Args: FishingMortality = FishingMortality, Yield = Yield
#     # FishingMortality: Vector of fishing mortalities to assess
#     # Yield: Vector of yields across FishingMortalities 
#   
#   # Returns: MSY value where derivative of yield is zero
#   
#   h_interval <- FishingMortality[2] - FishingMortality[1]
#   #for(i in 1:length(FishingMortality)-1){
#     YieldDeriv <- (Yield[i+1] - Yield[i])/(2*h_interval)
#   #}
# }







# # Storage vectors
# Recruits <- rep(NA, length(FishingMortality))
# Yield <- rep(NA, length(FishingMortality))
# YPR <- rep(NA, length(FishingMortality)) # yield per recruit
# # Z_mortality <- rep(NA, length(FishingMortality)) # total mortality matrix, iage rows, ifish = columns
# # NAA # Abundance at age, iage = rows, ifish = columns
# 
# # Need to Calculate
# SBPR # 1 valueper F value (length = length(FishingMortality))
# SBPRzero # SBPR at zero fishing
# 
# # If statements to calculate R, SSB, Yield
# for(ifish in 1:length(FishingMortality)){
#   
# }


