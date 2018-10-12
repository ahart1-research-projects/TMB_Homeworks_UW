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
beta_list <- -1*(beta_list[-1])
logalpha_list <- as.list(result$coefficients$random) # gives log alpha by species
logalpha_list <- c(logalpha_list[[1]])
alpha_list <- exp(logalpha_list)

##### Calculate recruitment deviations #####
RecDevs <- data$Residuals
names(RecDevs) <- data[,"SpeciesID"]

# Plot RecDevs and QQ plot by species
for(isp in 1:length(unique(data[,"SpeciesID"]))){
  plot(RecDevs[which(names(RecDevs)==isp)], main=paste("Species", isp, sep=" "), ylab="Rec Devs")
  abline(h=0)
  
  # Check rec devs normally distributed
  qqnorm(RecDevs[which(names(RecDevs)==isp)], main=paste("Normal Q-Q Plot Rec Devs Species", isp, sep=" "))
  qqline(RecDevs[which(names(RecDevs)==isp)])
}

print("Rec devs generally appears normally distributed. Species 5 has a small sample size so the rec devs are not as obviously normally distributedand in some instances small and large values do not exactly follow a normal distribution (e.g. species 8).")

# Summarize residual plots
plot(result,form=resid(.,type="p")~fitted(.)|SpeciesID,abline=0,pch=16)

# Boxplot residuals by species
boxplot(split(residuals(result),data$SpeciesID),ylab="Residual",xlab="Species",csi=0.2, main = "Recruitment Deviations by Species")
print("Rec devs are fairly well centered around zero. However, variability changes between species (species 3-6 have lower variability), which violates our modeling assumption of equal variance across all species. Accounting for species-specific variability may improve model fit.")
qqnorm(residuals(result), main=paste("Normal Q-Q Plot All Rec Devs"))
qqline(residuals(result))

# QQ plot for random effects (log alpha)
qqnorm(logalpha_list, main="Normal Q-Q Plot log alpha")
qqline(logalpha_list)
print("We assumed logalpha was normally distributed. Estimated values for logalpha generally follow a normal distribution, but several values do not. Given a larger sample size I expect these points would appear closer to the normal distribution.")

# QQ plot for fixed effects (beta)
qqnorm(beta_list, main = "Normal Q-Q Plot beta")
qqline(beta_list)
print("Similarly to estimated values of log alpha, estimated beta values appear roughly normally distibuted but there are several points which deviate from this distribution.")


# Predicted Recruitment
newdata <- NULL
for(isp in 1:11){
  tempnewdata <- seq(min(data[which(data[,"SpeciesID"]==isp),"SSB"]), max(data[which(data[,"SpeciesID"]==isp),"SSB"]),0.1)
  newdataR <- cbind(tempnewdata, rep(isp, length(tempnewdata)), rep(unique(data[,"SBPR"])[isp], length(tempnewdata)))
  colnames(newdataR) <- c("SSB", "SpeciesID","SBPR")
  newdataR <- as.data.frame(newdataR)
  newdata <- rbind(newdata, newdataR)
}
newdata$SpeciesID <- as.factor(newdata$SpeciesID) # species must be a factor in order for predict to work

newdata <- cbind(newdata, predict(result,newdata=newdata))
colnames(newdata)[ncol(newdata)] <- "Predictions"

temp <- newdata[,"Predictions"] + log(newdata[,"SSB"]) - log(newdata[,"SBPR"]) # calculate predicted recruitment given model structure
temp <- exp(temp)
newdata <- cbind(newdata, temp)
colnames(newdata)[ncol(newdata)] <- "R_pred"

for(isp in 1:11){
  plot(newdata[which(newdata[,"SpeciesID"]==isp),"SSB"], newdata[which(newdata[,"SpeciesID"]==isp),"R_pred"],
       main=paste("Species",isp, sep=" "), ylab="Predicted Recruitment", xlab="SSB", type="l")
}


plot(newdata$SSB,newdata$R_pred, col=newdata$SpeciesID)
legend("bottomright", legend=levels(newdata$SpeciesID), fill=unique(newdata$SpeciesID))



