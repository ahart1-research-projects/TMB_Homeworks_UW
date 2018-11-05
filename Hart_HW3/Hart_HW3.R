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
                        gamma = NULL){
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
    # gamma: number, no default
  
  # Return:
  
  
  
}

Ages <- seq(0,length(FecundityAtAge)-1)

SelectivityAtAgeMatrix <- cbind(data$Selectivity_F, data$Selectivity_M)
colnames(SelectivityAtAgeMatrix) <- c("Selectivity_F", "Selectivity_M")

WeightAtAgeMatrix <- cbind(data$Weight_F, data$Weight_M)
colnames(WeightAtAgeMatrix) <- c("Weight_F", "Weight_M")


# Test
M_par <- 0.2
Steepness  <- 0.3
Rzero <- 1000



# First calculate SBPRzero
Z_mortality <- CalcZMortality(FishingMortality = 0, M_par = M_par, SelectivityAtAge = SelectivityAtAgeMatrix)
NAA <- CalcNAA(Z_mortality = Z_mortality, Ages = Ages)
SBPRzero <- CalcSBPR(NAA = NAA, FecundityAtAge = data$Fecundity, Ages = Ages)

#################################################################################################################
# Storage vector
YPR <- rep(NA, length(FishingMortality))
SBPR <- rep(NA, length(FishingMortality))
Recruit <- rep(NA, length(FishingMortality))
Yield <- rep(NA, length(FishingMortality))
SSB <- rep(NA, length(FishingMortality))

# Loop over fishing efforts
for(ifish in 1:length(FishingMortality)){
  # Calculate total mortality at age for both sexes
  Z_mortality <- CalcZMortality(FishingMortality = FishingMortality[ifish], M_par = M_par, SelectivityAtAge = SelectivityAtAgeMatrix)
  # Calculate abundance at age for both sexes
  NAA <- CalcNAA(Z_mortality = Z_mortality, Ages = Ages)
  # Calculate Yield Per Recruit
  YPR[ifish] <- CalcYPR(FishingMortality = FishingMortality[ifish], WeightAtAge = WeightAtAgeMatrix, SelectivityAtAge = SelectivityAtAgeMatrix, NAA = NAA, Z_mortality = Z_mortality, Ages = Ages)
  # Calculate Spawner Biomas Per Recruit
  SBPR[ifish] <- CalcSBPR(NAA = NAA, FecundityAtAge = FecunditiyAtAge, Ages = Ages)
  # Calculate recruitment
  Recruit[ifish] <- CalcRecruit(FishingMortality = FishingMortality[ifish], StockRecruitForm = "BevertonHolt", SBPR = SBPR[ifish], SBPRzero = SBPRzero, Steepness = Steepness, Rzero = Rzero)
  # Calculate yield
  Yield[ifish] <- CalcYield(YPR = YPR[ifish], Recruit = Recruit[ifish])
  # Calculate Spawning Stock Biomass
  SSB[ifish] <- CalcSSB(SBPR = SBPR[ifish], Recruit = Recruit[ifish])
}







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












# par(mfrow=c(2,2)) # this makes a 4X4 plot

# !!!!!!!!!! SSB shouldn't be negative (probably because recruitment is negative)
# !!!!!!!!!! figure out how to plot replacement line
plot(x = SSB, y = Recruit, main = "Recruitment vs. SSB")
plot(x = SSB, y = Yield, main = "Yield vs. SSB")
plot(x = FishingMortality, y = Yield, main = "Yield vs. Fishing Mortality")
plot(x = FishingMortality, y = SSB, main = "SSB vs. Fishing Mortality")
