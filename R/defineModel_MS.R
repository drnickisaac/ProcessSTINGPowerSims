#' defineModel_MS
#'
#' @details Defines the Nimble model.
#' @return a set of code
#' @import nimble
#' @export

defineModel_MS <- function(inclPhenology = TRUE,
                        inclPanTrap = TRUE){

  modelcode <- nimbleCode({
    ######################### state model
    for(i in 1:nsp){
      for(j in 1:nsite){
        for(t in 1:nyear){
          linPred[i,j,t] <- alpha.s[i] + Trend * t  + spTr[i] * t + eta[j]
          log(lambda[i,j,t]) <- linPred[i,j,t]
          cloglog(psi[i,j,t]) <- linPred[i,j,t]
          z[i,j,t] ~ dbern(psi[i,j,t]) # True occupancy status
    }}}

    ######################### state model priors
    for(i in 1:nsp) {
      alpha.s[i] ~ dnorm(0, tau=0.0001)
      spTr[i] ~ dnorm(0, tau=tau.trend)
    }
    for(j in 1:nsite) {eta[j] ~ dnorm(0, sd=sd.eta)} # site-level random effect
    sd.eta ~ T(dt(0, 1, 1), 0, 10) # constrained
    tau.trend ~ T(dt(0, 1, 1), 0, Inf) # Half Cauchy
    Trend ~ dnorm(0, tau=0.0001)

    ######################### Obs model
    for(i in 1:nsp){
      for(k in 1:nvisit){
          ##### pan traps
          if(inclPanTrap){
            y1[i,k] ~ dbin(size = nT[k], prob = Py[i,k]) # Observed data
            Py[i,k] <- z[i, site[k], year[k]] * p1[i,k]
            p1[i,k] <- alpha.p[i] * pThin[i,k]
          } # Should I add site + year effects to detectability?

          ##### transects
          y2[i,k] ~ dpois(lambdaThin[i,k]) # Observed counts. Might need a NegBin here or Zero-inflated
          y3[i,k] ~ dpois(lambdaThin[i,k]) # Observed counts. Might need a NegBin here or Zero-inflated
          lambdaThin[i,k] <- Multiplier * lambda[i, site[k], year[k]] * pThin[i,k]

          #### shared: phenology
          if(inclPhenology){
            logit(pThin[i,k]) <- phScale[i] * f_JD[i,JulDate[k]] # assuming no site effect here
          } else {
            pThin[i,k] <- 1
          }
    }}

    ######################### Obs model priors
    for(i in 1:nsp){alpha.p[i] ~ dunif(0, 1)} # alpha.p is detection probability per pan trap at peak phenology. replace with beta dist
    Multiplier ~ T(dt(0, 1, 1), 0, Inf)

    ######################### Seasonality shared effect
    if(inclPhenology){
        for(i in 1:nsp){
          beta1[i] ~ dunif(100, 250) # peak detectability/activity. Constrained to fall within the field season
          beta2[i] ~ T(dt(0, 1, 1), 0, 500) # Half Cauchy. Stdev of phenology. At sd=500 the curve is entirely flat
          phScale[i] ~ T(dt(0, 1, 1), 0, Inf) # Half Cauchy

        for (d in 1:365){
          f_JD[i,d] <- 1/((2*3.14159265359)^0.5 * beta2[i]) * exp(-((d - (beta1[i]))^2 / (2* beta2[i]^2)))
          # could simplify this and evaluate only for dates in the dataset
        }
      }
    }
    #########################  derived parameters
    #for(i in 1:nsp){
    #  for(t in 1:nyr){
    #  psi.fs[i,t] <- mean(z[i,1:nsite,t])
    #  mu.lambda[i,t] <- mean(lambda[i,1:nsite, t])
    #}}
  })
  return(modelcode)
}

