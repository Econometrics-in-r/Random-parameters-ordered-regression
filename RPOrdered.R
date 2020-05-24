rm(list = ls())
require(stats4)
require(maxLik)
require(randtoolbox)

# reading and storing data in a dataframe
dataset <- read.csv(file.choose(),header=T)

#Sample Size
N <- nrow(dataset) 

#Dependent variable (severity of crashes in this example); change the variable as required
DVar <- dataset$SEV

#Number of ordered levels
S <- NROW(unique(DVar))

#Matrix of indicators 
I=matrix(rep.int(0,N*S),ncol = S)
for (j in 1:S) { I[(DVar==j),j]=1 }

# #standardization; comment out if not used
# namesVar = c("AADT") # name of variables to be standardized
# meanVar = NULL
# sdVar = NULL
# for (namei in namesVar){
#   meanVar = c(meanVar,mean(dataset[,namei]))
#   sdVar = c(sdVar,sd(dataset[,namei]))            
#   names(meanVar)[length(meanVar)] = namei
#   names(sdVar)[length(sdVar)] = namei
#   dataset[,namei] =  (dataset[,namei] - meanVar[namei])/sdVar[namei]
# }

# Halton Draws 
preparedraws=function()
{
  d=1
  while(d<(length(normaldraws)+1))
  {
    draws1[,normaldraws[d]]<<- qnorm(draws1[,normaldraws[d]])
    d=d+1
  }
}

Ndraws=500      # set number of draws 
dimensions=2    # define number of random parameters in the model

# generate draws (using Halton)
draws1=as.matrix(halton(Ndraws*N,dimensions))

# assign names to individual sets of draws - need one entry per dimension
colnames(draws1)=c("HRbeta1","HRbeta2")
# define whether any draws should be transformed to Normals, which is also needed for e.g. lognormals (leave empty if not)
normaldraws=c("HRbeta1","HRbeta2")

# preparing draws for estimation - this may take a while
preparedraws()

# placeholder for probabilities
prob = matrix(NA,nrow = N*Ndraws, ncol = S)

# fixing parameters across grouped observations i.e. grouped random parameters
# Comment out if there is no panel
block = length(unique(dataset[,'ID']))
ngroup = length(unique(dataset[,'Group']))
for (i in 1:Ndraws){
  tempInd = ((i-1)*block*ngroup) + (1:block)
  for (ii in 2:ngroup){
    draws1[tempInd+(ii-1)*block,] = draws1[tempInd,]
  }
}

## data preparation
# separating the variables with fixed parameters 
dataF =  as.matrix(data.frame(log(dataset$Length),dataset$NumLane,dataset$curv))
# separating the variables with random parameters 
dataR = as.matrix(data.frame(log(dataset$AADT),dataset$LWIDTH))

dataR2=NULL
for(i in 1:Ndraws){
  dataR2=rbind(dataR2,dataR)
}

draws1 = draws1[,1:dimensions]

# Likelihood function; this has been written for only 3 ordered categories (S=3), adjust as needed
LL <- function(params){  
  Fbeta <- params[1:3] # Fixed parameters in the mean Function (excluding the constant)
  MRbeta <- params[4:5]  # Mean of Random parameters in the mean function
  SDRbeta <- params[6:7]  # Std of Random parameters in the mean function
  cutpoint1 = params[8] #first threshold of the ordered model
  cutpoint2 = params[9] #second threshold of the ordered model
  
  # vector of indipendent variables with fixed parameters
  offset = rep.int(dataF%*%as.matrix(Fbeta,ncol=1),Ndraws)
  # simulating random parameters from their means and standard deviation
  beta = t( t(draws1)*SDRbeta + MRbeta )
  # constructing the mean function
  mu <- exp(offset+rowSums(dataR2*beta))
  
  # cumulative probability functions for logistic distribution (ordered logit)
  prob[,1] <- plogis(cutpoint1 - mu)
  prob[,2] <- plogis(cutpoint2 - mu) - prob[,1]
  prob[,3] <- 1 - prob[,1] - prob[,2]
  
  # cumulative probability functions for normal distribution (ordered probit)
  # prob[,1] <- pnorm(cutpoint1 - mu)
  # prob[,2] <- pnorm(cutpoint2 - mu) - prob[,1]
  # prob[,3] <- 1 - prob[,1] - prob[,2]

  # simulated (log)probabilities of each ordered category
  PR1 <- log(rowMeans(matrix(prob[,1], ncol = Ndraws)))
  PR2 <- log(rowMeans(matrix(prob[,2], ncol = Ndraws)))
  PR3 <- log(rowMeans(matrix(prob[,3], ncol = Ndraws)))
  
  # simulated loglikelihood for fractional split model; this has been written for only 3 ordered categories, adjust as needed
  loglik <- sum(I[,1]*PR1+I[,2]*PR2+I[,3]*PR3)
  
  return(loglik)
}

# initial values for optimization
init <- c(2,-5.5,0.66,#fixed parameters
          0.17,-0.14,#mean of random parameters
          0.05,0.08,#standard deviation of random parameters
          0.1,1)#thresholds 

# optimization (maximization of likelihood function)
fit1 <- maxLik(LL,start=init,method="BFGS")

summary(fit1)






