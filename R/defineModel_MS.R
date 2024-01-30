#' defineModel_MS
#'
#' @details Defines the Nimble model.
#' @param inclPhenology should the model account for seasonal variation?
#' @param inclPanTrap should the model include pan trap data?
#' @param incl2ndTransect should the model include data from the second transect walk?
#' @param inclStateRE should there be a site-level random effect in the state model?
#' @return a set of code
#' @import nimble
#' @export

defineModel_MS <- function(incl2ndTransect = TRUE,
                           inclPanTrap = TRUE,
                           inclPhenology = TRUE,
                           inclStateRE = FALSE){

  modelcode <- nimbleCode({
    ######################### state model
    for(i in 1:nsp){
      for(j in 1:nsite){
        for(t in 1:nyear){
          if(inclStateRE){
            linPred[i,j,t] <- lam.0[i] + spTr[i] * (t-1) + eta[j]
          } else {
            linPred[i,j,t] <- lam.0[i] + spTr[i] * (t-1)
          }
          log(lambda[i,j,t]) <- linPred[i,j,t]
          cloglog(psi[i,j,t]) <- linPred[i,j,t]
          z[i,j,t] ~ dbern(psi[i,j,t]) # True occupancy status
    }}}

    ######################### state model priors
    for(i in 1:nsp) {
      lam.0[i] ~ dnorm(0, tau=0.0001)
      spTr[i] ~ dnorm(Trend, tau=tau.trend)
    }
    tau.trend ~ T(dt(0, 1, 1), 0, Inf) # Half Cauchy

    if(inclStateRE){
      for(j in 1:nsite) {eta[j] ~ dnorm(0, sd=sd.eta)} # site-level random effect
      sd.eta ~ T(dt(0, 1, 1), 0, 10) # constrained
    }

    ######################### Obs model
    for(i in 1:nsp){
      for(k in 1:nvisit){
          ##### pan traps
          if(inclPanTrap){
            y1[i,k] ~ dbin(size = nT[k], prob = Py[i,k]) # Observed data
            Py[i,k] <- z[i, site[k], year[k]] * p1[i,k]
            if(inclPhenology){
              logit(p1[i,k]) <- alpha.0[i] + alpha.1[i]* (f_JD[i,JulDate[k]] - max(f_JD[i,1:365]))
            } else {
              logit(p1[i,k]) <- alpha.0[i]
            }
          } # Should I add site + year effects to detectability?

          ##### transects
          y2[i,k] ~ dpois(expectCount[i,k]) # Observed counts. Might need a NegBin here or Zero-inflated
          if(incl2ndTransect){
            y3[i,k] ~ dpois(expectCount[i,k]) # Observed counts. Might need a NegBin here or Zero-inflated
          }
          log(expectCount[i,k]) <- linPred[i, site[k], year[k]] * log(p2[i,k])
          if(inclPhenology){
            logit(p2[i,k]) <- gamma.0[i] + gamma.1[i]* (f_JD[i,JulDate[k]] - max(f_JD[i,1:365]))
          } else {
            logit(p2[i,k]) <- gamma.0[i]
          }
    }}

    ######################### Obs model priors
    for(i in 1:nsp){
      ######################### Obs model priors
      if(inclPanTrap){
        alpha.0[i] ~ dnorm(-2, tau = 0.0001) # logit detection probability per pan trap at peak phenology (or mean across year).
        #alpha.0[i] ~ dnorm(-2, tau = 1/2.72) # logit detection probability per pan trap at peak phenology (or mean across year).
        if(inclPhenology){
          #alpha.1[i] ~ T(dt(0, 1, 1), 0, Inf) # constrained to be positive
          alpha.1[i] ~ dnorm(2, tau = 0.0001)
        } # scaling parameter for detection on pan trap
      }
      #gamma.0[i] ~ dnorm(-2, tau = 1/2.72) # detection probability GLM on transects at peak phenology
      gamma.0[i] ~ dnorm(-2, tau = 0.0001) # detection probability GLM on transects at peak phenology
    }

    ######################### Seasonality shared effect
    if(inclPhenology){
        for(i in 1:nsp){
          for (d in 1:365){
            f_JD[i,d] <- 1/((2*3.14159265359)^0.5 * beta2[i]) * exp(-((d - (beta1[i]))^2 / (2* beta2[i]^2)))
            }
          gamma.1[i] ~ dnorm(2, tau = 0.0001)
          beta1[i] ~ dunif(100, 250) # peak detectability/activity. Constrained to fall within the field season
          beta2[i] ~ T(dt(0, 1, 1), 10, 200) # Half Cauchy. Stdev of phenology. At sd=500 the curve is entirely flat
      }
    }
    #########################  derived parameters
    for(i in 1:nsp){
      for(t in 1:nyear){
        psi.fs[i,t] <- mean(z[i,1:nsite,t])
        mu.lambda[i,t] <- mean(lambda[i,1:nsite, t])
    }}
  })
  return(modelcode)
}

