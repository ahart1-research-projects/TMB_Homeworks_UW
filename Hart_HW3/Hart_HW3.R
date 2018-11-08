# Amanda Hart
# Homework 3 FISH 559 


setwd("/Users/ahart2/Research/TMB_Homeworks_UW/Hart_HW3")

data <- read.table("HOME3.TXT", header=TRUE)
colnames(data) <- c("Age", "Fecundity", "Weight_F", "Weight_M", "Selectivity_F", "Selectivity_M")


SelectivityAtAgeMatrix <- cbind(data$Selectivity_F, data$Selectivity_M)
colnames(SelectivityAtAgeMatrix) <- c("Selectivity_F", "Selectivity_M")

WeightAtAgeMatrix <- cbind(data$Weight_F, data$Weight_M)
colnames(WeightAtAgeMatrix) <- c("Weight_F", "Weight_M")



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
  
  # Return: Matrix containing StockRecruitForm, MSY, FMSY, Fcrash
  
  FinalResults <- matrix(NA, length(FishingMortality),10)
  colnames(FinalResults) <- c("StockRecruitForm", "Steepness", "MSY", "FMSY", "Fcrash", "SSB0", "SSB_MSY", "FishingMortality", "SSB", "Recruitment")
  FinalResults[,"StockRecruitForm"] <- StockRecruitForm
  
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
    Recruit[ifish] <- Result$Recruit
    Yield[ifish] <- Result$Yield
    NAA <- Result$NAA
    SSB[ifish] <- Result$SSB
  }
  FinalResults[,"FishingMortality"] <- FishingMortality
  FinalResults[,"SSB"] <- SSB
  FinalResults[,"Recruitment"] <- Recruit
  FinalResults[,"SSB0"] <- Result$SBPRzero*Rzero
  FinalResults[,"Steepness"] <- Steepness
  
  # Calculate MSY and FMSY
  FMSY <- uniroot(CalcFMSY,interval=c(0,1), StockRecruitForm = StockRecruitForm, gamma_par = gamma_par, Steepness = Steepness, Rzero = Rzero) #, Steepness = Steepness)
  MSYResults <- CalcYield(FishingMortality = FMSY$root, M_par = M_par, SelectivityAtAge = SelectivityAtAgeMatrix, WeightAtAge = WeightAtAgeMatrix, FecundityAtAge = FecundityAtAge,
                          gamma_par = gamma_par, Ages = Ages, Steepness = Steepness, Rzero = Rzero, StockRecruitForm = StockRecruitForm)
  MSY <- MSYResults$Yield
  FinalResults[,"MSY"] <- MSY
  FinalResults[,"FMSY"] <- FMSY$root
  FinalResults[,"SSB_MSY"] <- MSYResults$SBPR*MSYResults$Recruit

  # Calculate Fcrash
  Fcrash <- uniroot(CalcFcrash, interval=c(0,1), StockRecruitForm = StockRecruitForm, Steepness = 0.5, Rzero = 1, gamma_par = 1)
  Fcrash <- Fcrash$root
  FinalResults[,"Fcrash"] <- Fcrash

  print(min(Recruit[which(Recruit >= 0)])) # This prints the minimum Recruitment greater than or equal to zero
  print(FishingMortality[which(Recruit == min(Recruit[which(Recruit >= 0)]))]) # prints fishing mortality corresponding to min recruitment
  print(Fcrash)
  
  if(DoPlot == TRUE){
    par(mfrow=c(2,2)) # this makes a 4X4 plot
    
    # Subset data
    plotRecruit <- Recruit[which(Recruit >=0)]
    plotSSB <- SSB[which(Recruit >= 0)]
    plotYield <- Yield[which(Recruit >= 0)]
    plotFishingMortality <- FishingMortality[which(Recruit >= 0)]
    
    # Plot
    plot(x = plotSSB, y = plotRecruit, main = "Recruitment vs. SSB", xlab = "SSB", ylab = "Recruitment")
    xreplacement <- seq(0,max(plotSSB))
    yreplacement <- seq(0,max(plotSSB))
    lines(x=xreplacement, y=yreplacement)
    plot(x = plotSSB, y = plotYield, main = "Yield vs. SSB", xlab = "SSB", ylab = "Yield")
    abline(h = MSY, col = "blue")
    plot(x = plotFishingMortality, y = plotYield, main = "Yield vs. Fishing Mortality", xlab = "Fishing Mortality", ylab = "Yield")
    abline(v = FMSY$root, col="red")
    abline(v = Fcrash, col = "green")
    abline(h = MSY, col = "blue")
    legend("topright", legend=c("FMSY", "MSY", "Fcrash"), fill=c("red", "blue", "green"))
    plot(x = plotFishingMortality, y = plotSSB, main = "SSB vs. Fishing Mortality", xlab = "Fishing Mortality", ylab = "SSB",
         xlim = c(0,max(plotFishingMortality)+0.1))
    abline(v = Fcrash, col = "green")
    legend("topright", legend=c("Fcrash"), fill=c("green"))
    mtext(StockRecruitForm, outer = TRUE, cex = 1.5, side = 3, line = -1.5)
  }
  
  return(FinalResults)
}

# Test with single value of F
HW3Function(FishingMortality = c(0,0.1,0.2), M_par = 0.15, FecundityAtAge = data$Fecundity, WeightAtAgeMatrix = WeightAtAgeMatrix,
            SelectivityAtAgeMatrix = SelectivityAtAgeMatrix, Rzero = 1, Steepness = 0.5, StockRecruitForm = "BevertonHolt", 
            gamma_par = 1, DoPlot = TRUE)

# Run across a range of F
Question3Results <- HW3Function(FishingMortality = c(seq(0,1,by=0.001)), M_par = 0.15, FecundityAtAge = data$Fecundity, 
                               WeightAtAgeMatrix = WeightAtAgeMatrix, SelectivityAtAgeMatrix = SelectivityAtAgeMatrix, 
                               Rzero = 1, Steepness = 0.5, StockRecruitForm = "BevertonHolt", gamma_par = 1, 
                               DoPlot = TRUE)



###########################################################################################
##### Question 4 
###########################################################################################
FormList <- c("BevertonHolt", "Ricker", "PellaTomlinson")
SteepnessList <- seq(0.25, 0.95, by = 0.05)

for(iform in FormList){
  FormResult <- NULL
  plotTable <- matrix(NA, length(SteepnessList), 2)
  colnames(plotTable) <- c("Steepness", "SSBratio")
  par(mfrow=c(1,4))
  for(isteep in 1:length(SteepnessList)){
    TempResult <- HW3Function(FishingMortality = c(seq(0,1,by=0.001)), M_par = 0.15, FecundityAtAge = data$Fecundity, 
                              WeightAtAgeMatrix = WeightAtAgeMatrix, SelectivityAtAgeMatrix = SelectivityAtAgeMatrix, 
                              Rzero = 1, Steepness = SteepnessList[isteep], StockRecruitForm = iform, gamma_par = 1, 
                              DoPlot = FALSE)
    TempResult <- as.data.frame(TempResult)
    SSBzeroRatio <- as.numeric(as.character(TempResult$SSB_MSY))/as.numeric(as.character((TempResult$SSB0)))
    TempResult <- cbind(TempResult, SSBzeroRatio)
    plotTable[isteep, "Steepness"] <- as.numeric(as.character(unique(TempResult$Steepness)))
    plotTable[isteep, "SSBratio"] <- as.numeric(as.character(unique(TempResult$SSBzeroRatio)))
    
    FormResult <- rbind(FormResult, TempResult)
  }
  # Plot Functional Form
  plot(0,0,xlim = c(0,1),ylim = c(0,1),type = "n")
  color <- rainbow(length(SteepnessList))
  for(isteep in 1:length(SteepnessList)){
    tempy <- as.numeric(as.character(FormResult[which(FormResult[,"Steepness"]==SteepnessList[isteep]),"Recruitment"]))
    tempy <- tempy[which(tempy >=0)]
    tempx <- as.numeric(as.character(FormResult[which(FormResult[,"Steepness"]==SteepnessList[isteep]),"SSB"]))
    tempx <- tempx[which(tempy >=0)]
    
    lines(x=tempx, y=tempy, col = color[isteep])
    lines(x = tempx, y = tempy, col = color[isteep])
  }
  # legend("bottomright",legend = SteepnessList, fill = color)
  
  # Plot SSB_MSY/SSB0 ratio information
  colnames(plotTable) <- c("Steepness", "S(Fmsy)/S(0)")
  plotTable[,"S(Fmsy)/S(0)"] <- round(plotTable[,"S(Fmsy)/S(0)"],4)
  
  library(gridExtra)
  library(grid)
  theme1 <- ttheme_minimal(core=list(bg_params = list(fill = c(color, rep("white", length(SteepnessList))), col=NA), fg_params=list(fontface=3)),
                           colhead=list(fg_params=list(col="black", fontface=4L)))
  plot(1,1,type="n", axes=FALSE, ann=FALSE)
  grid.table(plotTable, theme = theme1) # need to work on this grid since it isn't in the right place

}














