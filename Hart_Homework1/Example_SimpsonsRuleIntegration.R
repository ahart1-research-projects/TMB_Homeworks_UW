# Try coding Lecture 2 Simpson's rule integral correctly

Ncalls <- 3
h = (1-0)/(Ncalls-1)

b_step <- c(0,0+h,1)# seq(from=0, to=1, by=(1-0)/Ncalls) #!!! can't be a sequence, needs to be starting value + h... up until last call
print(b_step)

ResultVector <- rep(NA, length(b_step))
for(i in 1:length(b_step)){
  ResultVector[i] <- exp(-1*((b_step[i])^2))
}
print(ResultVector)

# ResultVector[1] <- ResultVector[1] # first and last multiplied by 1
# ResultVector[length(ResultVector)] <- ResultVector[length(ResultVector)]
# ResultVector[seq(2, length(ResultVector)-1, by=2)] <- 4*ResultVector[seq(2, length(ResultVector)-1, by=2)] # multiply evens by 4 except the last number
# #ResultVector[seq(3, length(ResultVector)-1, by=2)] <- 2*ResultVector[seq(3, length(ResultVector)-1, by=2)] # multiply odds by 2 except the last number
# ResultVector <- ResultVector/3
# 
# Likelihood <- h*(sum(ResultVector))
# print(Likelihood)

Likelihood_take2 <- h*(1/3*ResultVector[1] + 4/3*ResultVector[2] + 1/3*ResultVector[3])
print(Likelihood_take2)

