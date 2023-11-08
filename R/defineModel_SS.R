#' defineModel_SS
#'
#' @details Defines the Nimble model for one species.
#' @return a set of code
#' @import nimble
#' @export

defineModel_SS <- function(inclPhenology = TRUE,
                        inclPanTrap = TRUE){

  modelcode <- nimbleCode({
    ######################### state model
    for(j in 1:nsite){
        for(t in 1:nyear){
          linPred[j,t] <- alpha.s + Trend * t + eta[j]
          log(lambda[j,t]) <- linPred[j,t]
          cloglog(psi[j,t]) <- linPred[j,t]
          z[j,t] ~ dbern(psi[j,t]) # True occupancy status
    }}

    ######################### state model priors
    for(j in 1:nsite) {eta[j] ~ dnorm(0, sd=sd.eta)} # site-level random effect
    sd.eta ~ T(dt(0, 1, 1), 0, 10) # constrained
    Trend ~ dnorm(0, tau=0.0001)

    ######################### Obs model
    for(k in 1:nvisit) {
          ##### pan traps
          if(inclPanTrap){
            y1[k] ~ dbin(size = nT[k], prob = Py[k]) # Observed data
            Py[k] <- z[site[k], year[k]] * p1[k]
            p1[k] <- alpha.p * pThin[k]
          } # Should I add site + year effects to detectability?

          ##### transects
          y2[k] ~ dpois(lambdaThin[k]) # Observed counts. Might need a NegBin here or Zero-inflated
          y3[k] ~ dpois(lambdaThin[k]) # Observed counts. Might need a NegBin here or Zero-inflated
          lambdaThin[k] <- Multiplier * lambda[site[k], year[k]] * pThin[k]

          #### shared: phenology
          if(inclPhenology){
            logit(pThin[k]) <- phScale * f_JD[JulDate[k]] # assuming no site effect here
          } else {
            pThin[k] <- 1
          }
    }

    ######################### Obs model priors
    alpha.p ~ dunif(0, 1) # alpha.p is detection probability per pan trap at peak phenology. replace with beta dist
    Multiplier ~ T(dt(0, 1, 1), 0, Inf)

    ######################### Seasonality shared effect
    if(inclPhenology){
          beta1 ~ dunif(100, 250) # peak detectability/activity. Constrained to fall within the field season
          beta2 ~ T(dt(0, 1, 1), 0, 500) # Half Cauchy. Stdev of phenology. At sd=500 the curve is entirely flat
          phScale ~ T(dt(0, 1, 1), 0, Inf) # Half Cauchy

        for (d in 1:365){
          f_JD[d] <- 1/((2*3.14159265359)^0.5 * beta2) * exp(-((d - (beta1))^2 / (2* beta2^2)))
          # could simplify this and evaluate only for dates in the dataset
        }
      }

    #########################  derived parameters
    for(t in 1:nyear){
      #psi.fs[t] <- mean(z[1:nsite],t)
      mu.lambda[t] <- mean(lambda[1:nsite,t])
    }
  })
  return(modelcode)
}

