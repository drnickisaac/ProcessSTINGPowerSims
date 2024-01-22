#' defineModel_SS
#'
#' @details Defines the Nimble model for one species.
#' @param inclPhenology should the model account for seasonal variation?
#' @param inclPanTrap should the model include pan trap data?
#' @param incl2ndTransect should the model include data from the second transect walk?
#' @return a set of code
#' @import nimble
#' @export

defineModel_SS <- function(incl2ndTransect = TRUE,
                           inclPanTrap = TRUE,
                           inclPhenology = TRUE,
                           scalePheno = TRUE){

  modelcode <- nimbleCode({
    ######################### state model
    for(j in 1:nsite){
        for(t in 1:nyear){
          linPred[j,t] <- alpha.s + Trend * t #+ eta[j]
          log(lambda[j,t]) <- linPred[j,t]
          cloglog(psi[j,t]) <- linPred[j,t]
          z[j,t] ~ dbern(psi[j,t]) # True occupancy status
    }}

    ######################### state model priors
    #for(j in 1:nsite) {eta[j] ~ dnorm(0, sd=sd.eta)} # site-level random effect
    #sd.eta ~ T(dt(0, 1, 1), 0, 10) # constrained
    Trend ~ dnorm(0, tau=0.0001)
    alpha.s ~ dnorm(0, tau=0.0001)

    ######################### Obs model
    for(k in 1:nvisit) {
          ##### pan traps
          if(inclPanTrap){
            y1[k] ~ dbin(size = nT[k], prob = Py[k]) # Observed data
            Py[k] <- z[site[k], year[k]] * p1[k]
            logit(p1[k]) <- alpha.0 + alpha.1 * f_JD[JulDate[k]]
            #p1[k] <- alpha.p * pThin[k]
          } # Should I add site + year effects to detectability?

          ##### transects
          y2[k] ~ dpois(expectCount[k]) # Observed counts. Might need a NegBin here or Zero-inflated
          if(incl2ndTransect){
            y3[k] ~ dpois(expectCount[k]) # Observed counts. Might need a NegBin here or Zero-inflated
          }
          #expectCount[k] <- Multiplier * lambda[site[k], year[k]] * pThin[k]
          log(expectCount[k]) <- linPred[site[k], year[k]] * log(p2[k])
          logit(p2[k]) <- gamma.0 + gamma.1 * f_JD[JulDate[k]]

          #### shared: phenology
          #if(inclPhenology){
          #  logit(pThin[k]) <- phScale * f_JD[JulDate[k]] # assuming no site effect here
          #} else {
          #  pThin[k] <- 1
          #}
    }

    ######################### Obs model priors
    if(inclPanTrap){
      alpha.0 ~ dnorm(-2, tau = 0.0001) # logit detection probability per pan trap at peak phenology.
      alpha.1 ~ T(dt(0, 1, 1), 0, Inf) # scaling parameter for detection on pan trap
    }
    gamma.0 ~ dnorm(-2, tau = 0.0001) # intercept of detection probability GLM on transects
    gamma.1 ~ T(dt(0, 1, 1), 0, Inf) # slope of detection probability GLM on transects
    #Multiplier ~ T(dt(0, 1, 1), 0, Inf) # expected count per unit of lambda

    ######################### Seasonality shared effect
    if(inclPhenology){
         for (d in 1:365){
          f_JD[d] <- 1/((2*3.14159265359)^0.5 * beta2) * exp(-((d - (beta1))^2 / (2* beta2^2)))
          # could simplify this and evaluate only for dates in the dataset
         }
      if(scalePheno){
      # rescale the phenology curve so that it has a maximum value of 1
        phScale <- 1/max(f_JD[1:365])
        for (d in 1:365){
          f_JD[d] <- f_JD[d] * phScale
      }}
      beta1 ~ dunif(50, 300) # peak detectability/activity. Not constrained to fall within the field season (c(100, 250))
      beta2 ~ T(dt(0, 1, 1), 0, 500) # Half Cauchy. Stdev of phenology. At sd=500 the curve is entirely flat
      #phScale ~ T(dt(0, 1, 1), 0, Inf) # Half Cauchy
      #phScale <- 1/max(f_JD[1:365])
      } else {
        f_JD[1:365] <- 1/365
      }

    #########################  derived parameters
    for(t in 1:nyear){
      #psi.fs[t] <- mean(z[1:nsite],t)
      mu.lambda[t] <- mean(lambda[1:nsite,t])
    }
  })
  return(modelcode)
}

