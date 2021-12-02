###############################################################################
#### Steelhead Foraging Behavior Analysis
#### James Golden, Morgan Byrne, Keane Flynn
#### University of Nevada, Reno
##############################################################################

###
### Required Packages
### 
install.packages('jagsUI', 'runjags', 'coda', 'MCMCvis')
library('jagsUI', 'runjags', 'coda', 'MCMCvis')

###############################################################################
###
### Drift Foraging Model 
###
###############################################################################

cat(
  
  "

  model{
  
  ### Likelihood
  
  for(i in 1:nobs){
    logit(psi[i]) <- beta0 + beta1*bl[i] + beta2*volume[i] + beta3*dist[i] + beta4*nnd[i] #+ siteff[sitevec[i]]    # Define deterministic function
    y[i] ~ dbern(psi[i])     # Define stochastic function
    
    #~~~~
    # Bayes p
    #~~~~
    sim[i] ~ dbern(psi[i])
    
    SEEobs.all[i] <- (y[i] - psi[i])^2 # observed squared error
    SEEsim.all[i] <- (sim[i] - psi[i])^2 # simulated squared error
  }
  
  #~~~~
  # Bayes p (cont.)
  #~~~~
  
  SEEobs <- sum(SEEobs.all) # observed
  SEEsim <- sum(SEEsim.all) # simulated
  
  p.val <- step(SEEsim - SEEobs)
  
  
  ### Random effect
  
  # for(s in 1:nsites){
  #   siteff[s] ~ dnorm(0,siteprec)     # Hierarchical, site effects
  # }
  
   ### Priors
  
  beta1 ~ dnorm(0,1/2.25)     # Stochastic prior, low precision high variance, slope related to body length
  
  baseprob ~ dbeta(1,1)     # Stochastic prior
  
  beta0 <- log(baseprob/(1-baseprob))     # Deterministic function of baseprob, the intercept
  
  #siteprec ~ dgamma(0.0000001,0.0000001)     # Hyperprior for site precision/variance
  
  beta2 ~ dnorm(0,1/2.25)     #Low precision high variance for beta2, slope related to volume
  
  beta3 ~ dnorm(0,1/2.25)    #Slope related to distance traveled
  
  beta4 ~ dnorm(0, 1/2.25) #Slope Related to nearest neighbor

  
  }

    ", file="DriftForageJAGS_NRES746.txt"
)

### Data packages

data <- read.csv("~/Dropbox/746/Final Project/746 JAGS Analysis/Fish Analysis/ProjectData.csv")

data2 <- na.omit(data)

data2 <- data2[-which(data2$distance_cm_per_sec > 100),]
data2  = data2[-which(data2$volume > 500000),]
data2$site[data2$site == "18.1"] = 1
data2$site[data2$site == "18.2"] = 2
data2$site[data2$site == "18.3"] = 3
data2$site[data2$site == "18.4"] = 4
data2$site[data2$site == "18.5"] = 5
data2$site[data2$site == "18.6"] = 6

hist(data2$volume)

hist(dgamma(1000, 1/2.25, 1/2.25))

datalist <- list(
  y = data2$drift_forage,
  bl = scale(data2$length_mm)[,1],     # Fish body length scaled with mean=0 and sd=1
  volume = scale(data2$volume)[,1],
  dist = scale(data2$distance_cm_per_sec)[,1],
  nobs = nrow(data2),
  nsites = as.numeric(length(unique(data2$site))),
  sitevec = data2$site,
  nnd = scale(data2$NND_cm)[,1]
)

### Initial values: beta1 and baseprob

inits <- function(){
  list(
    
    beta1 = rnorm(1,0,2),     # Because there is one slope value, we can use 1 because its a scalar
    # Can make a vector if we have multiple slopes/parameters
    
    baseprob = runif(1,0,1),     # Same as beta listed above
    
    #siteprec = rgamma(1,0.1,0.1),     # Starting value for siteprec
    
    beta2 = rnorm(1,0,2), #starting value for beta2
    
    beta3 = rnorm(1,0,1), #starting value for beta3
    
    beta4 = rnorm(1,0,2) #starting value for beta4 
  )
}

inits()

### Parameters to save

paramsave <- c("beta1","beta2","beta3", "beta4", "siteprec", "siteff" ,'p.val',"SEEobs","SEEsim")

### Run the model in JAGS


set.seed(2021)
drift2 <- jags(data = datalist, parameters.to.save = paramsave, inits = inits,
              model.file = "DriftForageJAGS_NRES746.txt", n.chains = 3,
              n.iter = 100000, n.burnin = 50000, parallel = TRUE)

###############################################################################
###
### Benthic Foraging Model 
###
###############################################################################

cat(
  
  "

  model{
  
  ### Likelihood
  
  for(i in 1:nobs){
    logit(psi[i]) <- beta0 + beta1*bl[i] + beta2*volume[i] + beta3*dist[i] + beta4*nnd[i] # + siteff[sitevec[i]]     Define deterministic function
    y[i] ~ dbern(psi[i])     # Define stochastic function
    #~~~~
    # Bayes p
    #~~~~
    sim[i] ~ dbern(psi[i])
    
    SEEobs.all[i] <- (y[i] - psi[i])^2 # observed squared error
    SEEsim.all[i] <- (sim[i] - psi[i])^2 # simulated squared error
  }
  
   #~~~~
  # Bayes p (cont.)
  #~~~~
  
  SEEobs <- sum(SEEobs.all) # observed
  SEEsim <- sum(SEEsim.all) # simulated
  
  p.val <- step(SEEsim - SEEobs)
  
  ### Random effect
  
  # for(s in 1:nsites){
  #   siteff[s] ~ dnorm(0,siteprec)     # Hierarchical, site effects
  # }
  
    ### Priors
  
  beta1 ~ dnorm(0,1/2.25)     # Stochastic prior, low precision high variance, slope related to body length
  
  baseprob ~ dbeta(1,1)     # Stochastic prior
  
  beta0 <- log(baseprob/(1-baseprob))     # Deterministic function of baseprob, the intercept
  
  #siteprec ~ dgamma(0.1,0.1)     # Hyperprior for site precision/variance
  
  beta2 ~ dnorm(0,1/2.25)     #Low precision high variance for beta2, slope related to volume
  
  beta3 ~ dnorm(0,1/2.25)    #Slope related to distance traveled
  
  beta4 ~ dnorm(0, 1/2.25) #Slope Related to nearest neighbor
  
  }

    ", file="BenthicForageJAGS_NRES746.txt"
)

### Data packages

data <- read.csv("~/Dropbox/746/Final Project/746 JAGS Analysis/Fish Analysis/ProjectData.csv")

data2 <- na.omit(data)

data2 <- data2[-which(data2$distance_cm_per_sec > 100),]
data2  = data2[-which(data2$volume > 500000),]
data2$site[data2$site == "18.1"] = 1
data2$site[data2$site == "18.2"] = 2
data2$site[data2$site == "18.3"] = 3
data2$site[data2$site == "18.4"] = 4
data2$site[data2$site == "18.5"] = 5
data2$site[data2$site == "18.6"] = 6

datalist <- list(
  y = data2$benthic_forage,
  bl = scale(data2$length_mm)[,1],     # Fish body length scaled with mean=0 and sd=1
  volume = scale(data2$volume)[,1],
  dist = scale(data2$distance_cm_per_sec)[,1],
  nobs = nrow(data2),
  nsites = as.numeric(length(unique(data2$site))),
  sitevec = data2$site,
  nnd = scale(data2$NND_cm)[,1]
)
dist = scale(data2$distance_cm_per_sec)[,1]
nnd = scale(data2$NND_cm)[,1]
hist(dist)
hist(nnd)
summary(dist)
which.max(dist)
which(dist > 50)
data2
### Initial values: beta1 and baseprob

inits <- function(){
  list(
    
    beta1 = rnorm(1,0,2),     # Because there is one slope value, we can use 1 because its a scalar
    # Can make a vector if we have multiple slopes/parameters
    
    baseprob = runif(1,0,1),     # Same as beta listed above
    
    siteprec = rgamma(1,0.1,0.1),     # Starting value for siteprec
    
    beta2 = rnorm(1,0,2), #starting value for beta2
    
    beta3 = rnorm(1,0,1), #starting value for beta3
    
    beta4 = rnorm(1,0,2) #starting value for beta4 
  )
}

inits()

### Parameters to save

paramsave <- c("beta1","beta2","beta3", "beta4", "siteprec", "siteff" ,'p.val',"SEEobs","SEEsim")

### Run the model in JAGS


set.seed(2021)
benthicforage <- jags(data = datalist, parameters.to.save = paramsave, inits = inits,
               model.file = "BenthicForageJAGS_NRES746.txt", n.chains = 3, n.adapt = 1000,
               n.iter = 100000, n.burnin = 50000, parallel = TRUE)



###############################################################################
###
### Search Foraging Model 
###
###############################################################################


cat(
  
  "

  model{
  
  ### Likelihood
  
  for(i in 1:nobs){
    logit(psi[i]) <- beta0 + beta1*bl[i] + beta2*volume[i] + beta3*dist[i] + beta4*nnd[i] + siteff[sitevec[i]]    # Define deterministic function
    y[i] ~ dbern(psi[i])     # Define stochastic function
    #~~~~
    # Bayes p
    #~~~~
    sim[i] ~ dbern(psi[i])
    
    SEEobs.all[i] <- (y[i] - psi[i])^2 # observed squared error
    SEEsim.all[i] <- (sim[i] - psi[i])^2 # simulated squared error
  }
  
   #~~~~
  # Bayes p (cont.)
  #~~~~
  
  SEEobs <- sum(SEEobs.all) # observed
  SEEsim <- sum(SEEsim.all) # simulated
  
  p.val <- step(SEEsim - SEEobs)
  
  
  ### Random effect
  
  for(s in 1:nsites){
    siteff[s] ~ dnorm(0,siteprec)     # Hierarchical, site effects
  }
  
  ### Priors
  
  beta1 ~ dnorm(0,1/2.25)     # Stochastic prior, low precision high variance, slope related to body length
  
  baseprob ~ dbeta(1,1)     # Stochastic prior
  
  beta0 <- log(baseprob/(1-baseprob))     # Deterministic function of baseprob, the intercept
  
  siteprec ~ dgamma(0.1,0.1)     # Hyperprior for site precision/variance
  
  beta2 ~ dnorm(0,1/2.25)     #Low precision high variance for beta2, slope related to volume
  
  beta3 ~ dnorm(0,1/2.25)    #Slope related to distance traveled
  
  beta4 ~ dnorm(0, 1/2.25) #Slope Related to nearest neighbor
  
  }

    ", file="SearchForageJAGS_NRES746.txt"
)

### Data packages

data <- read.csv("~/Dropbox/746/Final Project/746 JAGS Analysis/Fish Analysis/ProjectData.csv")

data2 <- na.omit(data)

data2 <- data2[-which(data2$distance_cm_per_sec > 100),]
data2  = data2[-which(data2$volume > 500000),]
data2$site[data2$site == "18.1"] = 1
data2$site[data2$site == "18.2"] = 2
data2$site[data2$site == "18.3"] = 3
data2$site[data2$site == "18.4"] = 4
data2$site[data2$site == "18.5"] = 5
data2$site[data2$site == "18.6"] = 6

datalist <- list(
  y = data2$search_forage,
  bl = scale(data2$length_mm)[,1],     # Fish body length scaled with mean=0 and sd=1
  volume = scale(data2$volume)[,1],
  dist = scale(data2$distance_cm_per_sec)[,1],
  nobs = nrow(data2),
  nsites = as.numeric(length(unique(data2$site))),
  sitevec = data2$site,
  nnd = scale(data2$NND_cm)[,1]
)
hist(volume)

### Initial values: beta1 and baseprob

inits <- function(){
  list(
    
    beta1 = rnorm(1,0,2),     # Because there is one slope value, we can use 1 because its a scalar
    # Can make a vector if we have multiple slopes/parameters
    
    baseprob = runif(1,0,1),     # Same as beta listed above
    
    siteprec = rgamma(1,0.1,0.1),     # Starting value for siteprec
    
    beta2 = rnorm(1,0,2), #starting value for beta2
    
    beta3 = rnorm(1,0,1), #starting value for beta3
    
    beta4 = rnorm(1,0,2) #starting value for beta4 
  )
}

inits()

### Parameters to save

paramsave <- c("beta1","beta2","beta3", "beta4", "siteprec", "siteff" ,'p.val',"SEEobs","SEEsim")

### Run the model in JAGS


set.seed(2021)
searchforage <- jags(data = datalist, parameters.to.save = paramsave, inits = inits,
                     model.file = "SearchForageJAGS_NRES746.txt", n.chains = 3, n.adapt = 1000,
                     n.iter = 100000, n.burnin = 50000,parallel = TRUE)



###############################################################################
###
### Movement Foraging Model 
###
###############################################################################

cat(
  
  "

  model{
  
  ### Likelihood
  
  for(i in 1:nobs){
    logit(psi[i]) <- beta0 + beta1*bl[i] + beta2*volume[i] + beta3*dist[i] + beta4*nnd[i] + siteff[sitevec[i]]    # Define deterministic function
    y[i] ~ dbern(psi[i])     # Define stochastic function
     #~~~~
    # Bayes p
    #~~~~
    sim[i] ~ dbern(psi[i])
    
    SEEobs.all[i] <- (y[i] - psi[i])^2 # observed squared error
    SEEsim.all[i] <- (sim[i] - psi[i])^2 # simulated squared error
  }
  
   #~~~~
  # Bayes p (cont.)
  #~~~~
  
  SEEobs <- sum(SEEobs.all) # observed
  SEEsim <- sum(SEEsim.all) # simulated
  
  p.val <- step(SEEsim - SEEobs)
  
  ### Random effect
  
  for(s in 1:nsites){
    siteff[s] ~ dnorm(0,siteprec)     # Hierarchical, site effects
  }
  
    ### Priors
  
  beta1 ~ dnorm(0,1/2.25)     # Stochastic prior, low precision high variance, slope related to body length
  
  baseprob ~ dbeta(1,1)     # Stochastic prior
  
  beta0 <- log(baseprob/(1-baseprob))     # Deterministic function of baseprob, the intercept
  
  siteprec ~ dgamma(0.1,0.1)     # Hyperprior for site precision/variance
  
  beta2 ~ dnorm(0,1/2.25)     #Low precision high variance for beta2, slope related to volume
  
  beta3 ~ dnorm(0,1/2.25)    #Slope related to distance traveled
  
  beta4 ~ dnorm(0, 1/2.25) #Slope Related to nearest neighbor
  
  }

    ", file="MovementJAGS_NRES746.txt"
)

### Data packages

data <- read.csv("~/Dropbox/746/Final Project/746 JAGS Analysis/Fish Analysis/ProjectData.csv")

data2 <- na.omit(data)

data2 <- data2[-which(data2$distance_cm_per_sec > 100),]
data2  = data2[-which(data2$volume > 500000),]
data2$site[data2$site == "18.1"] = 1
data2$site[data2$site == "18.2"] = 2
data2$site[data2$site == "18.3"] = 3
data2$site[data2$site == "18.4"] = 4
data2$site[data2$site == "18.5"] = 5
data2$site[data2$site == "18.6"] = 6

datalist <- list(
  y = data2$movement ,
  bl = scale(data2$length_mm)[,1],     # Fish body length scaled with mean=0 and sd=1
  volume = scale(data2$volume)[,1],
  dist = scale(data2$distance_cm_per_sec)[,1],
  nobs = nrow(data2),
  nsites = as.numeric(length(unique(data2$site))),
  sitevec = data2$site,
  nnd = scale(data2$NND_cm)[,1]
)
dist = scale(data2$distance_cm_per_sec)[,1]
nnd = scale(data2$NND_cm)[,1]
hist(dist)
hist(nnd)
summary(dist)
which.max(dist)
which(dist > 50)
data2
### Initial values: beta1 and baseprob

inits <- function(){
  list(
    
    beta1 = rnorm(1,0,2),     # Because there is one slope value, we can use 1 because its a scalar
    # Can make a vector if we have multiple slopes/parameters
    
    baseprob = runif(1,0,1),     # Same as beta listed above
    
    siteprec = rgamma(1,0.1,0.1),     # Starting value for siteprec
    
    beta2 = rnorm(1,0,2), #starting value for beta2
    
    beta3 = rnorm(1,0,1), #starting value for beta3
    
    beta4 = rnorm(1,0,2) #starting value for beta4 
  )
}

inits()

### Parameters to save

paramsave <- c("beta1","beta2","beta3", "beta4", "siteprec", "siteff" ,'p.val',"SEEobs","SEEsim")

### Run the model in JAGS


set.seed(2021)
movement <- jags(data = datalist, parameters.to.save = paramsave, inits = inits,
                 model.file = "MovementJAGS_NRES746.txt", n.chains = 3, n.adapt = 1000,
                 n.iter = 100000, n.burnin = 50000, parallel = TRUE)



###############################################################################
###
### Surface Foraging Model 
###
###############################################################################

cat(
  
  "

  model{
  
  ### Likelihood
  
  for(i in 1:nobs){
    logit(psi[i]) <- beta0 + beta1*bl[i] + beta2*volume[i] + beta3*dist[i] + beta4*nnd[i] + siteff[sitevec[i]]    # Define deterministic function
    y[i] ~ dbern(psi[i])     # Define stochastic function
     #~~~~
    # Bayes p
    #~~~~
    sim[i] ~ dbern(psi[i])
    
    SEEobs.all[i] <- (y[i] - psi[i])^2 # observed squared error
    SEEsim.all[i] <- (sim[i] - psi[i])^2 # simulated squared error
  }
  
   #~~~~
  # Bayes p (cont.)
  #~~~~
  
  SEEobs <- sum(SEEobs.all) # observed
  SEEsim <- sum(SEEsim.all) # simulated
  
  p.val <- step(SEEsim - SEEobs)
  
  ### Random effect
  
  for(s in 1:nsites){
    siteff[s] ~ dnorm(0,siteprec)     # Hierarchical, site effects
  }
  
    ### Priors
  
  beta1 ~ dnorm(0,1/2.25)     # Stochastic prior, low precision high variance, slope related to body length
  
  baseprob ~ dbeta(1,1)     # Stochastic prior
  
  beta0 <- log(baseprob/(1-baseprob))     # Deterministic function of baseprob, the intercept
  
  siteprec ~ dgamma(0.1,0.1)     # Hyperprior for site precision/variance
  
  beta2 ~ dnorm(0,1/2.25)     #Low precision high variance for beta2, slope related to volume
  
  beta3 ~ dnorm(0,1/2.25)    #Slope related to distance traveled
  
  beta4 ~ dnorm(0, 1/2.25) #Slope Related to nearest neighbor
  
  }

    ", file="SurfaceStrikeJAGS_NRES746.txt"
)

### Data packages

data <- read.csv("~/Dropbox/746/Final Project/746 JAGS Analysis/Fish Analysis/ProjectData.csv")

data2 <- na.omit(data)

data2 <- data2[-which(data2$distance_cm_per_sec > 100),]
data2  = data2[-which(data2$volume > 500000),]
data2$site[data2$site == "18.1"] = 1
data2$site[data2$site == "18.2"] = 2
data2$site[data2$site == "18.3"] = 3
data2$site[data2$site == "18.4"] = 4
data2$site[data2$site == "18.5"] = 5
data2$site[data2$site == "18.6"] = 6

datalist <- list(
  y = data2$surface_strike,
  bl = scale(data2$length_mm)[,1],     # Fish body length scaled with mean=0 and sd=1
  volume = scale(data2$volume)[,1],
  dist = scale(data2$distance_cm_per_sec)[,1],
  nobs = nrow(data2),
  nsites = as.numeric(length(unique(data2$site))),
  sitevec = data2$site,
  nnd = scale(data2$NND_cm)[,1]
)
dist = scale(data2$distance_cm_per_sec)[,1]
nnd = scale(data2$NND_cm)[,1]
hist(dist)
hist(nnd)
summary(dist)
which.max(dist)
which(dist > 50)
data2
### Initial values: beta1 and baseprob

inits <- function(){
  list(
    
    beta1 = rnorm(1,0,2),     # Because there is one slope value, we can use 1 because its a scalar
    # Can make a vector if we have multiple slopes/parameters
    
    baseprob = runif(1,0,1),     # Same as beta listed above
    
    siteprec = rgamma(1,0.1,0.1),     # Starting value for siteprec
    
    beta2 = rnorm(1,0,2), #starting value for beta2
    
    beta3 = rnorm(1,0,1), #starting value for beta3
    
    beta4 = rnorm(1,0,2) #starting value for beta4 
  )
}

inits()

### Parameters to save

paramsave <- c("beta1","beta2","beta3", "beta4","siteff","p.val","SEEobs","SEEsim")

### Run the model in JAGS


set.seed(2021)
surfacestrike <- jags(data = datalist, parameters.to.save = paramsave, inits = inits,
                      model.file = "SurfaceStrikeJAGS_NRES746.txt", n.chains = 3, n.adapt = 1000,
                      n.iter = 100000, n.burnin = 50000, parallel = TRUE)


###############################################################################
###
### Attack Foraging Model 
###
###############################################################################

cat(
  
  "

  model{
  
  ### Likelihood
  
  for(i in 1:nobs){
    logit(psi[i]) <- beta0 + beta1*bl[i] + beta2*volume[i] + beta3*dist[i] + beta4*nnd[i] + siteff[sitevec[i]]    # Define deterministic function
    y[i] ~ dbern(psi[i])     # Define stochastic function
    #~~~~
    # Bayes p
    #~~~~
    sim[i] ~ dbern(psi[i])
    
    SEEobs.all[i] <- (y[i] - psi[i])^2 # observed squared error
    SEEsim.all[i] <- (sim[i] - psi[i])^2 # simulated squared error
  }
  
   #~~~~
  # Bayes p (cont.)
  #~~~~
  
  SEEobs <- sum(SEEobs.all) # observed
  SEEsim <- sum(SEEsim.all) # simulated
  
  p.val <- step(SEEsim - SEEobs)
  
  ### Random effect
  
  for(s in 1:nsites){
    siteff[s] ~ dnorm(0,siteprec)     # Hierarchical, site effects
  }
  
   ### Priors
  
  beta1 ~ dnorm(0,1/2.25)     # Stochastic prior, low precision high variance, slope related to body length
  
  baseprob ~ dbeta(1,1)     # Stochastic prior
  
  beta0 <- log(baseprob/(1-baseprob))     # Deterministic function of baseprob, the intercept
  
  siteprec ~ dgamma(0.1,0.1)     # Hyperprior for site precision/variance
  
  beta2 ~ dnorm(0,1/2.25)     #Low precision high variance for beta2, slope related to volume
  
  beta3 ~ dnorm(0,1/2.25)    #Slope related to distance traveled
  
  beta4 ~ dnorm(0, 1/2.25) #Slope Related to nearest neighbor
  
  }

    ", file="AttackJAGS_NRES746.txt"
)

### Data packages

data <- read.csv("~/Dropbox/746/Final Project/746 JAGS Analysis/Fish Analysis/ProjectData.csv")

data2 <- na.omit(data)

data2 <- data2[-which(data2$distance_cm_per_sec > 100),]
data2  = data2[-which(data2$volume > 500000),]
data2
data2$site[data2$site == "18.1"] = 1
data2$site[data2$site == "18.2"] = 2
data2$site[data2$site == "18.3"] = 3
data2$site[data2$site == "18.4"] = 4
data2$site[data2$site == "18.5"] = 5
data2$site[data2$site == "18.6"] = 6

datalist <- list(
  y = data2$surface_strike,
  bl = scale(data2$length_mm)[,1],     # Fish body length scaled with mean=0 and sd=1
  volume = scale(data2$volume)[,1],
  dist = scale(data2$distance_cm_per_sec)[,1],
  nobs = nrow(data2),
  nsites = as.numeric(length(unique(data2$site))),
  sitevec = data2$site,
  nnd = scale(data2$NND_cm)[,1]
)
dist = scale(data2$distance_cm_per_sec)[,1]
nnd = scale(data2$NND_cm)[,1]
hist(dist)
hist(nnd)
summary(dist)
which.max(dist)
which(dist > 50)
data2
data2[which(data2$volume > 500000),]
volume = scale(data2$volume)[,1]
hist(bl)
hist(volume)
hist(dist)
hist(nnd)
### Initial values: beta1 and baseprob

inits <- function(){
  list(
    
    beta1 = rnorm(1,0,2),     # Because there is one slope value, we can use 1 because its a scalar
    # Can make a vector if we have multiple slopes/parameters
    
    baseprob = runif(1,0,1),     # Same as beta listed above
    
    siteprec = rgamma(1,0.1,0.1),     # Starting value for siteprec
    
    beta2 = rnorm(1,0,2), #starting value for beta2
    
    beta3 = rnorm(1,0,1), #starting value for beta3
    
    beta4 = rnorm(1,0,2) #starting value for beta4 
  )
}

inits()

### Parameters to save

paramsave <- c("beta1","beta2","beta3", "beta4", "siteprec", "siteff" ,'p.val',"SEEobs","SEEsim")

### Run the model in JAGS


set.seed(2021)
attack <- jags(data = datalist, parameters.to.save = paramsave, inits = inits,
               model.file = "AttackJAGS_NRES746.txt", n.chains = 3, n.adapt = 1000,
               n.iter = 100000, n.burnin = 50000, parallel = TRUE)


##########################################################################################################
###
### Diagnostics 
###
##########################################################################################################

###
### MCMC Objects
### 

driftMCMC <- drift2$samples
benthicforageMCMC <- benthicforage$samples
searchforageMCMC <- searchforage$samples
movementMCMC <- movement$samples
surfacestrikeMCMC <- surfacestrike$samples
attackMCMC <- attack$samples

###
### Trace Plots/Density Plots
###

### Drift Forage 

# Betas

drift_forageMCMCtrace = MCMCtrace(driftMCMC,params = c("beta1","beta2","beta3","beta4"), type = 'trace',
                                  ISB = FALSE, exact = TRUE, post_zm = TRUE, pdf = FALSE,
                                  Rhat = TRUE, n.eff = TRUE, 
                                  main_tr = c("Beta1 - Body Length",
                                              "Beta2 - Volume"," Beta3 - Distance", 
                                              "Beta4 - Nearest Neighbor"), 
                                  main_den = c("Beta1 - Body Length",
                                               "Beta2 - Volume"," Beta3 - Distance"
                                               ,"Beta4 - Nearest Neighbor"))


drift_forageMCMCdens = MCMCtrace(driftMCMC,params = c("beta1","beta2","beta3","beta4"), type = 'density',
                                  ISB = FALSE, exact = TRUE, post_zm = TRUE, pdf = FALSE,ind = TRUE,
                                  Rhat = TRUE, n.eff = TRUE, 
                                  main_tr = c("Beta1 - Body Length",
                                              "Beta2 - Volume"," Beta3 - Distance", 
                                              "Beta4 - Nearest Neighbor"), 
                                  main_den = c("Beta1 - Body Length",
                                               "Beta2 - Volume"," Beta3 - Distance"
                                               ,"Beta4 - Nearest Neighbor"))

### Benthic Forage 

# Betas

benthic_forageMCMCtrace = MCMCtrace(benthicforageMCMC,params = c("beta1","beta2","beta3","beta4"),
                                    ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, 
                                    type = 'trace',
                                    Rhat = TRUE, n.eff = TRUE,
                                    main_tr = c("Beta1 - Body Length",
                                                "Beta2 - Volume"," Beta3 - Distance", 
                                                "Beta4 - Nearest Neighbor"), 
                                    main_den = c("Beta1 - Body Length",
                                                 "Beta2 - Volume"," Beta3 - Distance"
                                                 ,"Beta4 - Nearest Neighbor"))

benthic_forageMCMCdens = MCMCtrace(benthicforageMCMC,params = c("beta1","beta2","beta3","beta4"), type = 'density',
                                    ISB = FALSE, exact = TRUE, post_zm = TRUE, pdf = FALSE,ind = TRUE,
                                    Rhat = TRUE, n.eff = TRUE, 
                                    main_tr = c("Beta1 - Body Length",
                                                "Beta2 - Volume"," Beta3 - Distance", 
                                                "Beta4 - Nearest Neighbor"), 
                                    main_den = c("Beta1 - Body Length",
                                                 "Beta2 - Volume"," Beta3 - Distance"
                                                 ,"Beta4 - Nearest Neighbor"))

### Search Forage 

# Betas
search_forageMCMCtrace = MCMCtrace(searchforageMCMC,params = c("beta1","beta2","beta3","beta4"),
                                   ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, 
                                   type = 'trace',
                                   Rhat = TRUE, n.eff = TRUE,
                                   main_tr = c("Beta1 - Body Length",
                                               "Beta2 - Volume"," Beta3 - Distance", 
                                               "Beta4 - Nearest Neighbor"), 
                                   main_den = c("Beta1 - Body Length",
                                                "Beta2 - Volume"," Beta3 - Distance"
                                                ,"Beta4 - Nearest Neighbor"))

search_forageMCMCdens = MCMCtrace(searchforageMCMC,params = c("beta1","beta2","beta3","beta4"), type = 'density',
                                   ISB = FALSE, exact = TRUE, post_zm = TRUE, pdf = FALSE,ind = TRUE,
                                   Rhat = TRUE, n.eff = TRUE, 
                                   main_tr = c("Beta1 - Body Length",
                                               "Beta2 - Volume"," Beta3 - Distance", 
                                               "Beta4 - Nearest Neighbor"), 
                                   main_den = c("Beta1 - Body Length",
                                                "Beta2 - Volume"," Beta3 - Distance"
                                                ,"Beta4 - Nearest Neighbor"))

# Site Effect 

search_forageMCMCtracesiteff = MCMCtrace(searchforageMCMC,params = c("siteff[1]","siteff[2]","siteff[3]","siteff[4]","siteff[5]","siteff[6]"), type = 'trace',
                                  ISB = FALSE, exact = TRUE, post_zm = TRUE, pdf = FALSE,
                                  Rhat = TRUE, n.eff = TRUE)

search_forageMCMCdenssiteff = MCMCtrace(searchforageMCMC,params = c("siteff[1]","siteff[2]","siteff[3]","siteff[4]","siteff[5]","siteff[6]"), type = 'density',
                                  ISB = FALSE, exact = TRUE, post_zm = TRUE, pdf = FALSE,
                                  Rhat = TRUE, n.eff = TRUE, ind = TRUE)


### Movement Forage 

# Betas 

movement_forageMCMCtrace = MCMCtrace(movementMCMC,params = c("beta1","beta2","beta3","beta4"),
                                     ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, 
                                     type = 'trace',
                                     Rhat = TRUE, n.eff = TRUE,
                                     main_tr = c("Beta1 - Body Length",
                                                 "Beta2 - Volume"," Beta3 - Distance", 
                                                 "Beta4 - Nearest Neighbor"), 
                                     main_den = c("Beta1 - Body Length",
                                                  "Beta2 - Volume"," Beta3 - Distance"
                                                  ,"Beta4 - Nearest Neighbor"))


movement_forageMCMCdens = MCMCtrace(movementMCMC,params = c("beta1","beta2","beta3","beta4"), type = 'density',
                                     ISB = FALSE, exact = TRUE, post_zm = TRUE, pdf = FALSE,ind = TRUE,
                                     Rhat = TRUE, n.eff = TRUE, 
                                     main_tr = c("Beta1 - Body Length",
                                                 "Beta2 - Volume"," Beta3 - Distance", 
                                                 "Beta4 - Nearest Neighbor"), 
                                     main_den = c("Beta1 - Body Length",
                                                  "Beta2 - Volume"," Beta3 - Distance"
                                                  ,"Beta4 - Nearest Neighbor"))
# Site Effect

movement_forageMCMCtracesiteff = MCMCtrace(movementMCMC,params = c("siteff[1]","siteff[2]","siteff[3]","siteff[4]","siteff[5]","siteff[6]"), type = 'trace',
                                  ISB = FALSE, exact = TRUE, post_zm = TRUE, pdf = FALSE,
                                  Rhat = TRUE, n.eff = TRUE)

movement_forageMCMCdenssiteff = MCMCtrace(movementMCMC,params = c("siteff[1]","siteff[2]","siteff[3]","siteff[4]","siteff[5]","siteff[6]"), type = 'density',
                                  ISB = FALSE, exact = TRUE, post_zm = TRUE, pdf = FALSE,
                                  Rhat = TRUE, n.eff = TRUE, ind = TRUE)

### Surface Strike 

# Betas 

surfacestrikeMCMCtrace = MCMCtrace(surfacestrikeMCMC,params = c("beta1","beta2","beta3","beta4"), type = 'trace',
                                   ISB = FALSE, exact = TRUE, post_zm = TRUE, pdf = FALSE,
                                   Rhat = TRUE, n.eff = TRUE, 
                                   main_tr = c("Beta1 - Body Length",
                                               "Beta2 - Volume"," Beta3 - Distance", 
                                               "Beta4 - Nearest Neighbor"), 
                                   main_den = c("Beta1 - Body Length",
                                                "Beta2 - Volume"," Beta3 - Distance"
                                                ,"Beta4 - Nearest Neighbor"))

surfacestrikeMCMCdens = MCMCtrace(surfacestrikeMCMC,params = c("beta1","beta2","beta3","beta4"), type = 'density',
                                   ISB = FALSE, exact = TRUE, post_zm = TRUE, pdf = FALSE,ind = TRUE,
                                   Rhat = TRUE, n.eff = TRUE, 
                                   main_tr = c("Beta1 - Body Length",
                                               "Beta2 - Volume"," Beta3 - Distance", 
                                               "Beta4 - Nearest Neighbor"), 
                                   main_den = c("Beta1 - Body Length",
                                                "Beta2 - Volume"," Beta3 - Distance"
                                                ,"Beta4 - Nearest Neighbor"))
# Site Effect 


surfacestrikeMCMCtracesiteff = MCMCtrace(surfacestrikeMCMC,params = c("siteff[1]","siteff[2]","siteff[3]","siteff[4]","siteff[5]","siteff[6]"), type = 'trace',
                                  ISB = FALSE, exact = TRUE, post_zm = TRUE, pdf = FALSE,
                                  Rhat = TRUE, n.eff = TRUE)

surfacestrikeMCMCdenssiteff = MCMCtrace(surfacestrikeMCMC,params = c("siteff[1]","siteff[2]","siteff[3]","siteff[4]","siteff[5]","siteff[6]"), type = 'density',
                                  ISB = FALSE, exact = TRUE, post_zm = TRUE, pdf = FALSE,
                                  Rhat = TRUE, n.eff = TRUE, ind = TRUE)

### Attack Forage

# Betas

attackMCMCtrace = MCMCtrace(attackMCMC,params = c("beta1","beta2","beta3","beta4"), type = 'trace',
                            ISB = FALSE, exact = TRUE, post_zm = TRUE, pdf = FALSE,
                            Rhat = TRUE, n.eff = TRUE, 
                            main_tr = c("Beta1 - Body Length",
                                        "Beta2 - Volume"," Beta3 - Distance", 
                                        "Beta4 - Nearest Neighbor"), 
                            main_den = c("Beta1 - Body Length",
                                         "Beta2 - Volume"," Beta3 - Distance"
                                         ,"Beta4 - Nearest Neighbor"))

attackMCMCdens = MCMCtrace(attackMCMC,params = c("beta1","beta2","beta3","beta4"), type = 'density',
                            ISB = FALSE, exact = TRUE, post_zm = TRUE, pdf = FALSE,ind = TRUE,
                            Rhat = TRUE, n.eff = TRUE, 
                            main_tr = c("Beta1 - Body Length",
                                        "Beta2 - Volume"," Beta3 - Distance", 
                                        "Beta4 - Nearest Neighbor"), 
                            main_den = c("Beta1 - Body Length",
                                         "Beta2 - Volume"," Beta3 - Distance"
                                         ,"Beta4 - Nearest Neighbor"))

# Site Effect 

attackMCMCtracesiteff = MCMCtrace(attackMCMC,params = c("siteff[1]","siteff[2]","siteff[3]","siteff[4]","siteff[5]","siteff[6]"), type = 'trace',
                                  ISB = FALSE, exact = TRUE, post_zm = TRUE, pdf = FALSE,
                                  Rhat = TRUE, n.eff = TRUE)

attackMCMCdenssiteff = MCMCtrace(attackMCMC,params = c("siteff[1]","siteff[2]","siteff[3]","siteff[4]","siteff[5]","siteff[6]"), type = 'density',
                                  ISB = FALSE, exact = TRUE, post_zm = TRUE, pdf = FALSE,
                                  Rhat = TRUE, n.eff = TRUE, ind = TRUE)



### 
### Gelman-Rubin Diagnostics
###

### Drift Forage 

gelman.diag(driftMCMC)
gelman.plot(driftMCMC)

### Benthic Forage

gelman.diag(benthicforageMCMC)
gelman.plot(benthicforageMCMC)

### Search Forage

gelman.diag(searchforageMCMC)
gelman.plot(searchforageMCMC)

### Movement Forage 

gelman.diag(movementMCMC)
gelman.plot(movementMCMC)

### Surface Strike 

gelman.diag(surfacestrikeMCMC)
gelman.plot(surfacestrikeMCMC)

### Attack Forage

gelman.diag(attackMCMC)
gelman.plot(attackMCMC)


###
### Posterior Predictive Checks 
###
 
### Drift Forage 

plot(drift2$sims.list$SEEobs,drift2$sims.list$SEEsim,xlab="SSEobs",ylab="SSEsim", main = "Drift Forage")
abline(0,1,lwd = 2, col = "red")

### Benthic Forage 

plot(benthicforage$sims.list$SEEobs,benthicforage$sims.list$SEEsim,xlab="SSEobs",ylab="SSEsim", main = "Benthic Forage")
abline(0,1,lwd = 2, col = "red")

### Search Forage 

plot(searchforage$sims.list$SEEobs,searchforage$sims.list$SEEsim,xlab="SSEobs",ylab="SSEsim", main = "Search Forage")
abline(0,1,lwd = 2, col = "red")

### Movement Forage

plot(movement$sims.list$SEEobs,movement$sims.list$SEEsim,xlab="SSEobs",ylab="SSEsim", main = "Movement Forage")
abline(0,1,lwd = 2, col = "red")

### Surface Strike

plot(surfacestrike$sims.list$SEEobs,surfacestrike$sims.list$SEEsim,xlab="SSEobs",ylab="SSEsim", main = "Surface Forage")
abline(0,1,lwd = 2, col = "red")

### Attack 

plot(attack$sims.list$SEEobs,attack$sims.list$SEEsim,xlab="SSEobs",ylab="SSEsim", main = "Attack Forage")
abline(0,1,lwd = 2, col = "red")


##########################################################################################################
###
### Parameter Estimates 
###
##########################################################################################################

###
### Model Output Summaries
### 

### Drift Forage 

summary(drift2)
driftsum = drift2$summary
driftsum

### Benthic Forage 

summary(benthicforage)
benthicforagesum = benthicforage$summary
benthicforagesum

### Search Forage 

summary(searchforage)
searchforagesum = searchforage$summary
searchforagesum
### Movement Forage 

summary(movement)
movementsum = movement$summary
movementsum

### Surface Strike 

summary(surfacestrike)
surfacestrikesum = surfacestrike$summary
surfacestrikesum

### Attack

summary(attack)
attacksum = attack$summary
attacksum

##########################################################################################################
###
### Important Figures 
###
##########################################################################################################

### Bayesian Catepillar Plots

# Betas 

temp.layout <- layout(matrix(c(1,2,3,4,5,6),3,2,byrow=TRUE), respect = TRUE,widths=c(lcm(8),lcm(8)), heights=c(lcm(8),lcm(8),lcm(8)))
layout.show(temp.layout)

par(mar=c(0, 4, 2, 0) + 0.1)

drift.plot = MCMCplot(object = driftMCMC, ci = c(2.5, 97.5),
                      params = c('beta1','beta2','beta3','beta4'), main = "Drift Forage",
                      col = "black", labels = c("Body Length", "Volume", "Distance", "Nearest Neighbor"),
                      xlim = c(-1.5,0.5), guide_lines = TRUE)

par(mar=c(0, 1, 2, 0) + 0.1)
search.plot = MCMCplot(object = searchforageMCMC,ci = c(2.5, 97.5),
                       params = c('beta1','beta2','beta3','beta4'),
                       main = "Search Forage", col = "red",
                       labels = c("Body Length", "Volume", "Distance", "Nearest Neighbor"),
                       xlim = c(-1.5,1.0), guide_lines = TRUE)

par(mar=c(0, 1, 2, 0) + 0.1)
movement.plot = MCMCplot(object = movementMCMC, ci = c(2.5, 97.5),
                         params = c('beta1','beta2','beta3','beta4'), 
                         main = "Movement Forage", col = "blue",
                         labels = c("Body Length", "Volume", "Distance", "Nearest Neighbor"),
                         xlim = c(-.5,2.0), guide_lines = TRUE)

par(mar=c(2, 4, 0, 0) + 0.1)
benthicforage.plot = MCMCplot(object = benthicforageMCMC, ci = c(2.5, 97.5),
                              params = c('beta1','beta2','beta3','beta4'), 
                              main = "Benthic Forage", col = "purple",
                              labels = c("Body Length", "Volume", "Distance", "Nearest Neighbor"),
                              xlim = c(-4,2), guide_lines = TRUE)

par(mar=c(2, 1, 0, 0) + 0.1)
surfacestrike.plot = MCMCplot(object = surfacestrikeMCMC, ci = c(2.5, 97.5),
                              params = c('beta1','beta2','beta3','beta4'), 
                              main = "Surface Strike", col = "gold", 
                              labels = c("Body Length", "Volume", "Distance", "Nearest Neighbor"),
                              xlim = c(-5,5), guide_lines = TRUE)

par(mar=c(2, 1, 0, 0) + 0.1)
attack.plot = MCMCplot(object = attackMCMC, ci = c(2.5, 97.5),
                       params = c('beta1','beta2','beta3','beta4'), 
                       main = "Attack", col = "cyan", 
                       labels = c("Body Length", "Volume", "Distance", "Nearest Neighbor"), guide_lines = TRUE)

mtext("Confidence Around Mean Estimates for Each Model", side = 3, line = -7, outer = TRUE,
      cex = 1.2, adj = 0.6)
mtext("Parameters", side = 2, line = -11, outer = TRUE,
      cex = 1.2)

# Site Effect

temp.layout <- layout(matrix(c(1,2,3,4),2,2,byrow=TRUE), respect = TRUE,widths=c(lcm(8),lcm(8)), heights=c(lcm(8),lcm(8),lcm(8)))
layout.show(temp.layout)


par(mar=c(0, 1, 2, 0) + 0.1)
search.plot = MCMCplot(object = searchforageMCMC,ci = c(2.5, 97.5),
                       params = c('siteff'),
                       main = "Search Forage", col = "red",
                       labels = c("Site 1", "Site 2", "Site 3", "Site 4", "Site 5", "Site 6"),
                       xlim = c(-6,6), guide_lines = TRUE)

par(mar=c(0, 1, 2, 0) + 0.1)
movement.plot = MCMCplot(object = movementMCMC, ci = c(2.5, 97.5),
                         params = c('siteff'), 
                         main = "Movement Forage", col = "blue",
                         labels = c("Site 1", "Site 2", "Site 3", "Site 4", "Site 5", "Site 6"),
                         xlim = c(-6,6), guide_lines = TRUE)


par(mar=c(2, 1, 0, 0) + 0.1)
surfacestrike.plot = MCMCplot(object = surfacestrikeMCMC, ci = c(2.5, 97.5),
                              params = c('siteff'), 
                              main = "Surface Strike", col = "gold", 
                              labels = c("Site 1", "Site 2", "Site 3", "Site 4", "Site 5", "Site 6"),
                              xlim = c(-6,6), guide_lines = TRUE)

par(mar=c(2, 1, 0, 0) + 0.1)
attack.plot = MCMCplot(object = attackMCMC, ci = c(2.5, 97.5),
                       params = c('siteff'), 
                       main = "Attack", col = "cyan", 
                       labels = c("Site 1", "Site 2", "Site 3", "Site 4", "Site 5", "Site 6"),
                       xlim = c(-6,6), guide_lines = TRUE)

mtext("Confidence Around Site Effect", side = 3, line = -17, outer = TRUE,
      cex = 1.2, adj = 0.5)
mtext("Sites", side = 2, line = -13, outer = TRUE,
      cex = 1.2)

dev.off()

### P-Values 

p.vals = c(driftsum[5,1],
           benthicforagesum[5,1],
           searchforagesum[12,1],
           movementsum[12,1],
           surfacestrikesum[11,1],
           attacksum[12,1])
p.vals = round(p.vals,2)
p.vals = as.data.frame(as.numeric(p.vals))
models = c("Drift", "Benthic", "Search", "Movement","Surface", "Attack")
data.frame(models)
k = cbind(models,p.vals)
colnames(k) = c("Models", "P-Values")

library("kableExtra")

p.tab = k %>% 
  kbl(caption = "Bayesian P-Values for Each Model") %>% 
  kable_classic(full_width = FALSE, html_font = "Cambria")
p.tab

### Summary Tables 

Drift.Forage = driftsum[c("beta1","beta2","beta3","beta4","deviance"),c("mean","sd","2.5%","97.5%","Rhat")]
Benthic.Forage = benthicforagesum[c("beta1","beta2","beta3","beta4","deviance"),c("mean","sd","2.5%","97.5%","Rhat")]
Search.Forage = searchforagesum[c("beta1","beta2","beta3","beta4","deviance","siteff[1]","siteff[2]","siteff[3]","siteff[4]","siteff[5]","siteff[6]"),c("mean","sd","2.5%","97.5%","Rhat")]
Movement.Forage = movementsum[c("beta1","beta2","beta3","beta4","deviance","siteff[1]","siteff[2]","siteff[3]","siteff[4]","siteff[5]","siteff[6]"),c("mean","sd","2.5%","97.5%","Rhat")]
Surface.Strike = surfacestrikesum[c("beta1","beta2","beta3","beta4","deviance","siteff[1]","siteff[2]","siteff[3]","siteff[4]","siteff[5]","siteff[6]"),c("mean","sd","2.5%","97.5%","Rhat")]
Attack = attacksum[c("beta1","beta2","beta3","beta4","deviance","siteff[1]","siteff[2]","siteff[3]","siteff[4]","siteff[5]","siteff[6]"),c("mean","sd","2.5%","97.5%","Rhat")]

drift.t =  Drift.Forage %>% 
  kbl(caption = "Drift Forage Summary Table") %>% 
  kable_classic(full_width = FALSE, html_font = "Cambria")

benthic.t = Benthic.Forage %>% 
  kbl(caption = "Benthic Forage Summary Table") %>% 
  kable_classic(full_width = FALSE, html_font = "Cambria")

search.t = Search.Forage %>% 
  kbl(caption = "Search Forage Summary Table") %>% 
  kable_classic(full_width = FALSE, html_font = "Cambria")

movement.t = Movement.Forage %>% 
  kbl(caption = "Movement Forage Summary Table") %>% 
  kable_classic(full_width = FALSE, html_font = "Cambria")

surface.t = Surface.Strike %>% 
  kbl(caption = "Surface Forage Summary Table") %>% 
  kable_classic(full_width = FALSE, html_font = "Cambria")

Attack.t = Attack %>% 
  kbl(caption = "Attack Forage Summary Table") %>% 
  kable_classic(full_width = FALSE, html_font = "Cambria")

par(mfrow = c(6,1))
drift.t
benthic.t
search.t
movement.t
surface.t
Attack.t


















