# Amanda Hart
# Homework 1 FISH 559 Question 2

setwd("/Users/ahart2/Research/TMB_Workshop_Fall2018/ClassWork/Hart_Homework1")

data <- read.table("HOME1B.TXT", header=TRUE)
colnames(data) <- c("SSB", "Recruitment", "Species", "SBPR")

# Fit linear random effects to 11 species data
# ln(alpha) = random effect
# beta = species-specific fixed effect
# noise around stock-recruit = log-normal
# variation around stock-recruit same for all species



# Ricker stock-recruit 
# Recruit_sp_yr = alpha_sp*SSB_sp_yr*exp(-1*beta_sp*SSB_sp_yr)/SBPR_sp_yr

# Log of Ricker stock-recruit
# log(Recruit_sp_yr) = log(alpha_sp) - log(SBPR_sp_yr) + log(SSB_sp_yr) - beta_sp*SSB_sp_yr

# SBPR0 = SBPR
# SSB_yr = SSB


##### Read in data #####
data <- read.table("HOME1B.TXT", header=TRUE)
colnames(data) <- c("SSB", "Recruitment", "SpeciesID", "SBPR")
data[,"SpeciesID"] <- as.factor(data[,"SpeciesID"])

SSB <- data[,"SSB"]
Recruitment <- data[,"Recruitment"]
SpeciesID <- data[,"SpeciesID"]
SBPR <- data[,"SBPR"]

RecData <- list(SSB, Recruitment, SpeciesID, SBPR)


##### Fit mixed effect model #####
library(nlme)
# recruitment depends on fixed effects for data provided in RecData with a random effect for species ID
result <- lme((log(Recruitment) + log(SBPR) - log(SSB))  ~ SSB:SpeciesID, random = ~ 1 | SpeciesID, data=data, method = "ML")
 # !!!! look at syntax for SpeciesID, use -1   (~ SSB:SpeciesID-1 this gives slope only, no intercept) # here I want the intercept
# SSB:SpeciesID-1 estimates beta by species
summary(result)

coef(result) # gives log alpha by species
resid(result) # give residuals for each row of data
result$logLik
result$data # data you gave model
beta_list <- result$coefficients$fixed # gives - beta by species as a fixed effect
alpha_list <- as.list(result$coefficients$random) # gives log alpha by species
alpha_list <- c(alpha_list[[1]])

##### Use estimated alphas and betas to predict recruitment #####
alpha_list <- exp(alpha_list)
beta_list <- -1*beta_list

data <- cbind(data,rep(NA,nrow(data)))
colnames(data)[ncol(data)] <- "R_pred"

for(irow in 1:nrow(data)){
  sp <- as.numeric(data[irow,"SpeciesID"])
  R_pred <- alpha_list[sp]*data[irow,"SSB"]*exp(-1*beta_list[sp]*data[irow,"SSB"])/data[irow,"SBPR"]
  data[irow,"R_pred"] <- R_pred
}

# Plot R & predicted R by species
for(isp in 1:length(unique(data[,"SpeciesID"]))){
  plot(y =data[which(data[,"SpeciesID"]==isp),"Recruitment"], x=data[which(data[,"SpeciesID"]==isp),"SSB"],  col="red", main=paste("Species",isp, sep=" "), ylab="Recruitment", xlab="SSB")
  lines(y=data[which(data[,"SpeciesID"]==isp),"R_pred"], x=data[which(data[,"SpeciesID"]==isp),"SSB"], col="black")
  legend("topright", legend=c("Observed R", "Predicted R"), fill=c("red", "black"))
}


##### Calculate recruitment deviations #####
RecDevs <- data$R_pred - data$Recruitment
names(RecDevs) <- data[,"SpeciesID"]

# Plot RecDevs by species
for(isp in 1:length(unique(data[,"SpeciesID"]))){
  plot(RecDevs[which(names(RecDevs)==isp)], main=paste("Species", isp, sep=" "), ylab="Rec Devs")
  abline(h=0)
}








##### TO DO #####
# plot residuals (of all data)
# look at documentation and see what residuals are for

