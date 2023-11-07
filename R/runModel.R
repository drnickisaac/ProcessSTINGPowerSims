#' runModel
#'
#' @details Runs the model.
#' @param dataConstants dataframe produced by ProcessSimDatFile()
#' @param obsData dataframe produced by ProcessSimDatFile()
#' @param dataSumm dataframe produced by ProcessSimDatFile()
#' @param useNimble option to bypass the model fitting in Nimble (just for testing the code)
#' @param n.iter number of iterations for the Nimble model
#' @param inclPhenology should the model account for seasonal variation?
#' @param inclPanTrap should the model include pan trap data?
#' @param maxSp maximum number of species. If set at 1, the model is run sequentially for every species in the dataset
#' @param parallelize should the chains be run as separate processes on different cores?
#' @return a set of year effects
#' @export
#' @import nimble
#' @import pbmcapply
#' @import parallel

runModel <- function(dataConstants,
                     obsData,
                     dataSumm,
                     useNimble = TRUE,
                     n.iter = 1000,
                     inclPhenology = TRUE,
                     inclPanTrap = TRUE,
                     maxSp = 9999,
                     parallelize = FALSE){

  if(useNimble) {
      # truncate the dataset if there are too many species
      if(dim(obsData$y1)[1] > maxSp){
        obsData <- lapply(obsData, function(x) x[1:maxSp,])
        dataConstants$nsp <- maxSp
        dataSumm$occMatrix <- dataSumm$occMatrix[1:maxSp,,]
        dataSumm$naiveOcc <- dataSumm$naiveOcc[1:maxSp]
        dataSumm$reportingRate <- dataSumm$reportingRate[1:maxSp]
  }

  ###################################################################

  if(maxSp >1){ # Multispecies option

    # step 1 define the model code
      modelcode <- defineModel_MS(inclPhenology = inclPhenology, inclPanTrap = inclPanTrap) # this makes no difference!

    # step 2 create an operational from from NIMBLE/BUGS code
    model <- nimbleModel(code = modelcode,
                         constants = dataConstants,
                         data = obsData,
                         inits = list(z = dataSumm$occMatrix,
                                      alpha.s = cloglog(dataSumm$naiveOcc),
                                      alpha.p = dataSumm$reportingRate, # replace with reportingRate_1 when I can calculate it
                                      beta1 = rep(180, dataConstants$nsp),
                                      beta2 = rep(50, dataConstants$nsp),
                                      phScale = rep(1, dataConstants$nsp),
                                      Multiplier = 1,
                                      tau.eta = abs(rt(1, 1)),
                                      tau.trend = abs(rt(1, 1)),
                                      Trend = rnorm(n=1))
    )

    # step 3 build an MCMC object using buildMCMC(). we can add some customization here
    occMCMC <- buildMCMC(model,
                       monitors = c("Trend"),
                       thin = 3,
                       useConjugacy = FALSE) # useConjugacy controls whether conjugate samplers are assigned when possible

  # step 3 before compiling the MCMC object we need to compile the model first
  Cmodel <- compileNimble(model) # NJBI: I don't understand why this step is necessary

  # now the MCMC (project = NIMBLE model already associated with a project)
  CoccMCMC <- compileNimble(occMCMC, project = model)
  # instantaneous

  # and now we can use either $run or runMCMC() on the compiled model object.
  if(parallelize){
    av_cores <- parallel::detectCores() - 1
    runMCMC_samples <- pbmcapply::pbmclapply(1:3, function(i)
                              runMCMC(
                                mcmc = CoccMCMC,
                                nburnin = n.iter/2,
                                niter = n.iter,
                                nchains = 1, samplesAsCodaMCMC = T),
                              mc.cores = av_cores)

  } else {
    runMCMC_samples <- runMCMC(CoccMCMC,
                             nburnin = n.iter/2,
                             niter = n.iter,
                             nchains = 3, samplesAsCodaMCMC = T)
  }

  ############################################ end multispecies

    } else { # sequential single-species option
      # step 1 define the model code
      modelcode <- defineModel_SS(inclPhenology = inclPhenology, inclPanTrap = inclPanTrap)

      # step 2 create an operational from from NIMBLE/BUGS code
      model <- nimbleModel(code = modelcode,
                           constants = dataConstants#,
                           #data = obsData,
                           #inits = list(z = dataSumm$occMatrix,
                                        #alpha.s = cloglog(dataSumm$naiveOcc),
                                        #alpha.p = dataSumm$reportingRate, # replace with reportingRate_1 when I can calculate it
                                        #beta1 = rep(180, dataConstants$nsp),
                                        #beta2 = rep(50, dataConstants$nsp),
                                        #phScale = rep(1, dataConstants$nsp),
                                        #Multiplier = 1)
      )

      # step 3 build an MCMC object using buildMCMC(). we can add some customization here
      occMCMC <- buildMCMC(model,
                           monitors = c("Trend"),
                           thin = 3,
                           useConjugacy = FALSE) # useConjugacy controls whether conjugate samplers are assigned when possible

      # step 3 before compiling the MCMC object we need to compile the model first
      Cmodel <- compileNimble(model) # NJBI: I don't understand why this step is necessary

      # now the MCMC (project = NIMBLE model already associated with a project)
      CoccMCMC <- compileNimble(occMCMC, project = model)
      # instantaneous

      lapply(1:nsp, function(i){
        # might have to use nimbleMCMC rather than runMCMC for this.
        # and now we can use either $run or runMCMC() on the compiled model object.
        runMCMC_samples <- runMCMC(CoccMCMC,
                                   data = lapply(obsData, function(x) x[,i]), # this might throw an error
                                   inits = list(z = dataSumm$occMatrix[i],
                                     alpha.s = cloglog(dataSumm$naiveOcc)[i],
                                     alpha.p = dataSumm$reportingRate[i], # replace with reportingRate_1 when I can calculate it
                                     beta1 = rep(180, dataConstants$nsp),
                                     beta2 = rep(50, dataConstants$nsp),
                                     phScale = rep(1, dataConstants$nsp),
                                     Multiplier = 1,
                                     tau.eta = abs(rt(1, 1)),
                                     tau.trend = abs(rt(1, 1)),
                                     Trend = rnorm(n=1)),
                                   nburnin = n.iter/2, niter = n.iter, nchains = 3, samplesAsCodaMCMC = T)

      })
    }

  #####################################################################

  }
  else {
    # for simplicity, let's just report the annual fitted total count
    totalObs <- rowSums(sapply(obsData, colSums))
    mod <- with(dataConstants, glm(totalObs ~ factor(year) + site, family = "poisson"))

    yearEff <- c(1, exp(coef(mod)[grep("factor", names(coef(mod)))]))
    names(yearEff) <- paste0("Year", 1:dataConstants$nyear)
    return(yearEff)
  }
}

