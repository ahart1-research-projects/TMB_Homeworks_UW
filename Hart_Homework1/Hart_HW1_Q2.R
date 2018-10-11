# Amanda Hart
# Homework 1 FISH 559 Question 2

setwd("/Users/ahart2/Research/TMB_Homeworks_UW/Hart_Homework1")

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
summary(result)

coef(result) # gives log alpha by species
Residuals <- resid(result) # give residuals for each row of data
data <- cbind(data,Residuals) # give residuals for each row of data
colnames(data)[ncol(data)] <- c("Residuals")
result$logLik
# result$data # data you gave model
beta_list <- result$coefficients$fixed # gives (- beta) by species as a fixed effect
logalpha_list <- as.list(result$coefficients$random) # gives log alpha by species
logalpha_list <- c(logalpha_list[[1]])

##### Calculate recruitment deviations #####
RecDevs <- data$Residuals
names(RecDevs) <- data[,"SpeciesID"]

# Plot RecDevs and QQ plot by species
for(isp in 1:length(unique(data[,"SpeciesID"]))){
  plot(RecDevs[which(names(RecDevs)==isp)], main=paste("Species", isp, sep=" "), ylab="Rec Devs")
  abline(h=0)
  
  # Check rec devs normally distributed
  qqnorm(RecDevs[which(names(RecDevs)==isp)], main=paste("Normal Q-Q Plot Species", isp, sep=" "))
  qqline(RecDevs[which(names(RecDevs)==isp)])
    # Rec devs appear normally distributed for most species. 
    # Species 5 has isn't quite normally distributed but this is expected due to small sample size.
}

# Summarize residual plots
plot(result,form=resid(.,type="p")~fitted(.)|SpeciesID,abline=0,pch=16)

# Boxplot residuals by species
boxplot(split(residuals(result),data$SpeciesID),ylab="Residual",xlab="Species",csi=0.2)

# qq plot for random effects
## plot fixed effects and species differ from eachother (not equal between species)
# plot f ixed effects

# what does augPred() do?
# is the below appropriate given what andre said yesterday about accounting for the random effects structrure?
# also plot random effects? see slide 21 how do I find random effect values? what does coef(result) give me

predict(model, newdata=range of values to predict at)

##### The below is not the right  way to do these projections ######
##### Use estimated alphas and betas to predict recruitment #####
alpha_list <- exp(logalpha_list)
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

















##### TO DO #####
# plot residuals (of all data)
# look at documentation and see what residuals are for

