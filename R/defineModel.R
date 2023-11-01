#' defineModel
#'
#' @details Defines the Nimble model.
#' @return a set of code
#' @import nimble
#' @export

defineModel <- function(inclPhenology = TRUE,
                        inclPanTrap = TRUE){

  modelcode <- nimbleCode({
    ######################### state model
    for(i in 1:nsp) {
      for(j in 1:nsite) {
        linPred[i,j] <- alpha.s[i] + eta[j]
        log(lambda[i,j]) <- linPred[i,j]
        cloglog(psi[i,j]) <- linPred[i,j]
        z[i,j] ~ dbern(psi[i,j]) # True occupancy status
        # need to add in a temporal component
    }}

    ######################### state model priors
    for(i in 1:nsp) {alpha.s[i] ~ dnorm(0, 0.0001)}
    for(j in 1:nsite) {eta[j] ~ dnorm(0, tau.eta)} # site-level random effect
    tau.eta ~ T(dt(0, 1, 1), 0, Inf) # Half Cauchy

    ######################### Obs model
    for(i in 1:nsp){
      for(k in 1:nvisit) {
          ##### pan traps
          if(inclPanTrap){
            y1[i,k] ~ dbin(size = nT[k], prob = Py[i,k]) # Observed data
            Py[i,k] <- z[i,site[k]] * p1[i,k]
            p1[i,k] <- alpha.p[i] * pThin[i,k]
          } # NEED TO ADD A SITE TERM IN HERE to the detectability element

          ##### transects
          y2[i,k] ~ dpois(lambdaThin[i,k]) # Observed counts. Might need a NegBin here or Zero-inflated
          y3[i,k] ~ dpois(lambdaThin[i,k]) # Observed counts. Might need a NegBin here or Zero-inflated
          lambdaThin[i,k] <- Multiplier * lambda[i, site[k]] * pThin[i,k]

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
    for(i in 1:nsp){
      psi.fs[i] <- mean(z[i,1:nsite])
      mu.lambda[i] <- mean(lambda[i,1:nsite])
    }
  })
  return(modelcode)
}

