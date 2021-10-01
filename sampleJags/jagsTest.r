library(dplyr)
df <- read.csv("Downloads/R-Program/NRES_746/NRES746-FinalProject/jagsTrial/finalDataFrame.csv") %>% na.omit()

dataTable <- list(
  y = df$drift_forage,
  bl = scale(df$length_mm)[,1],
  nobs = length(df)
)

inits <- function(){
  list(
    slope=rnorm(1,0,2),
    baseP=runif(1)
  )
}
initsTrial <- inits

parToSave <- c("baseP","slope")

output <- jags(data=dataTable,parameters.to.save=parToSave,inits=inits,
    model.file = "Downloads/R-Program/NRES_746/NRES746-FinalProject/jagsTrial/jagsFile.txt",
     n.iter=1000,n.cores=T,n.chains=2)

mcmc <- output$samples
plot(mcmc[,"baseP"])
gelman.diag(mcmc)