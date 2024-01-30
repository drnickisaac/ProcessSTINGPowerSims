#' defineModel_SS
#'
#' @details Defines the Nimble model for one species.
#' @param inclPhenology should the model account for seasonal variation?
#' @param inclPanTrap should the model include pan trap data?
#' @param incl2ndTransect should the model include data from the second transect walk?
#' @param inclStateRE should there be a site-level random effect in the state model?
#' @return a set of code
#' @import nimble
#' @export

defineModel_SS <- function(incl2ndTransect = TRUE,
                           inclPanTrap = TRUE,
                           inclPhenology = TRUE,
                           inclStateRE = FALSE){

  modelcode <- nimbleCode({
    ######################### state model
    for(j in 1:nsite){
        for(t in 1:nyear){
          if(inclStateRE){
            linPred[j,t] <- lam.0 + Trend * (t-1) + eta[j]
          } else {
            linPred[j,t] <- lam.0 + Trend * (t-1)
          }
          log(lambda[j,t]) <- linPred[j,t]
          cloglog(psi[j,t]) <- linPred[j,t]
          z[j,t] ~ dbern(psi[j,t]) # True occupancy status
    }}

    ######################### state model priors
    if(inclStateRE){
      for(j in 1:nsite) {eta[j] ~ dnorm(0, sd=sd.eta)} # site-level random effect
      sd.eta ~ T(dt(0, 1, 1), 0, 10) # constrained
    }
    Trend ~ dnorm(0, tau=0.001)
    lam.0 ~ dnorm(0, tau=1/2.72)

    ######################### Obs model
    for(k in 1:nvisit) {
      ##### pan traps
      if(inclPanTrap){
        y1[k] ~ dbin(size = nT[k], prob = Py[k]) # Observed data
        Py[k] <- z[site[k], year[k]] * p1[k]
        if(inclPhenology){
          logit(p1[k]) <- alpha.0 + alpha.1 * (f_JD[JulDate[k]] - max(f_JD[1:365]))
        } else {
          logit(p1[k]) <- alpha.0
        }
      }

      ##### transects
      y2[k] ~ dpois(expectCount[k]) # Observed counts. Might need a NegBin here or Zero-inflated
      if(incl2ndTransect){
        y3[k] ~ dpois(expectCount[k]) # Observed counts. Might need a NegBin here or Zero-inflated
      }
      log(expectCount[k]) <- linPred[site[k], year[k]] * log(p2[k])
      if(inclPhenology){
        logit(p2[k]) <- gamma.0 + gamma.1 * (f_JD[JulDate[k]] - max(f_JD[1:365]))
      } else {
        logit(p2[k]) <- gamma.0
      }
    }

    ######################### Seasonality shared effect
    if(inclPhenology){
      for (d in 1:365){
        f_JD[d] <- 1/((2*3.14159265359)^0.5 * beta2) * exp(-((d - (beta1))^2 / (2* beta2^2)))
      }
    }
    ######################### Obs model priors
    if(inclPanTrap){
      #alpha.0 ~ dnorm(-2, tau = 0.0001) # logit detection probability per pan trap at peak phenology (or mean across year).
      alpha.0 ~ dnorm(-2, tau = 1/2.72) # logit detection probability per pan trap at peak phenology (or mean across year).
      if(inclPhenology){
        alpha.1 ~ T(dt(0, 1, 1), 0, Inf) # constrained to be positive
        #alpha.1 ~ dnorm(2, tau = 0.0001)
        } # scaling parameter for detection on pan trap
    }
    gamma.0 ~ dnorm(-2, tau = 1/2.72) # detection probability GLM on transects at peak phenology
    #gamma.0 ~ dnorm(-2, tau = 0.0001) # detection probability GLM on transects at peak phenology


    if(inclPhenology){
      #gamma.1 ~ dnorm(2, tau = 0.0001) # detection probability GLM on transects at peak phenology
      gamma.1 ~ T(dt(0, 1, 1), 0, Inf) # # constrained to be positive
      beta1 ~ dunif(100, 250) # peak detectability/activity. Not constrained to fall within the field season (c(100, 250))
      beta2 ~ T(dt(0, 1, 1), 10, 200) # Half Cauchy. Stdev of phenology. At sd=500 the curve is entirely flat
    }
    #########################  derived parameters
    for(t in 1:nyear){
      psi.fs[t] <- mean(z[1:nsite,t])
      mu.lambda[t] <- mean(lambda[1:nsite,t])
    }
  })
  return(modelcode)
}

