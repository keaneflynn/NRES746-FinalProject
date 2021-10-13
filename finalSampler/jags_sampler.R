
###
### NRES 746 Final Project
###

### Morgan Byrne
### Keane Flynn
### James Golden

library("rjags")
library('jagsUI')
library('runjags')
library('coda')
library('MCMCvis')
testjags()




cat(

"

  model{
  
  ### Likelihood
  
  for(i in 1:nobs){
    logit(psi[i]) <- beta0 + beta1*bl[i] + siteff[sitevec[i]]    # Define deterministic function
    y[i] ~ dbern(psi[i])     # Define stochastic function
  }
  
  ### Random effect
  
  for(s in 1:nsites){
    siteff[s] ~ dnorm(0,siteprec)     # Hierarchical, site effects
  }
  
  ### Priors
  
  beta1 ~ dnorm(0,0.0001)     # Stochastic prior, low precision high variance, slope related to body length
  
  baseprob ~ dbeta(1,1)     # Stochastic prior
  
  beta0 <- log(baseprob/(1-baseprob))     # Deterministic function of baseprob, the intercept
  
  siteprec ~ dgamma(0.1,0.1)     # Hyperprior for site precision/variance
  
  }

    ", file="JAGSSamplerNRES746.txt"
)

### Data packages

data <- read.csv("ProjectData.csv")

data2 <- na.omit(data)

data2$site[data2$site == "18.1"] = 1
data2$site[data2$site == "18.2"] = 2
data2$site[data2$site == "18.3"] = 3
data2$site[data2$site == "18.4"] = 4
data2$site[data2$site == "18.5"] = 5
data2$site[data2$site == "18.6"] = 6

datalist <- list(
  y = data2$drift_forage,
  bl = scale(data2$length_mm)[,1],     # Fish body length scaled with mean=0 and sd=1
  nobs = nrow(data2),
  nsites = as.numeric(length(unique(data2$site))),
  sitevec = data2$site
)

### Initial values: beta1 and baseprob

inits <- function(){
  list(
    
    beta1 = rnorm(1,0,2),     # Because there is one slope value, we can use 1 because its a scalar
    # Can make a vector if we have multiple slopes/parameters
    
    baseprob = runif(1,0,1),     # Same as beta listed above
    
    siteprec = rgamma(1,0.1,0.1)     # Starting value for siteprec
    
  )
}

inits()

### Parameters to save

paramsave <- c("beta0","beta1","siteprec","baseprob")

### Run the model in JAGS

?jags

drift <- jags(data = datalist, parameters.to.save = paramsave, inits = inits,
              model.file = "JAGSSamplerNRES746.txt", n.chains = 3, n.adapt = 1000,
              n.iter = 30000, n.burnin = 35000, parallel = TRUE)

### Look at results of MCMC
summary(drift)
drift$summary
MCMC <- drift$samples

class(MCMC)

### Trace plots and density plots

plot(MCMC)
plot(MCMC[,"beta0"], ylab = "Likelihood", main = "Beta0")
plot(MCMC[,"baseprob"], ylab = "Likelihood", main = "Baseprob")
plot(MCMC[,"beta1"], ylab = "Likelihood", main = "Beta1")
plot(MCMC[,"deviance"], ylab = "Likelihood", main = "Deviance")

### Run a convergence diagnostic: Coda package

gelman.diag(MCMC) # If the upper confidence bound if below 1.1, good convergence!


### Other visualization: MCMCvis package

library(MCMCvis)
save.image("drift_forage.RData")

