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
#' @return a set of year effects
#' @import nimble
#' @export

# NEED  TO FIGURE OUT HOW TO DO THIS FOR A SINGLE SPECIES OPTION

runModel <- function(dataConstants,
                     obsData,
                     dataSumm,
                     useNimble = FALSE,
                     n.iter = 1000,
                     inclPhenology = TRUE,
                     inclPanTrap = TRUE,
                     maxSp = 9999){

  if(useNimble) {
    # warning("this is only being written on the fly!")

    # truncate the dataset if there are too many species
    if(dim(obsData$y1)[1] > maxSp){
      obsData <- lapply(obsData, function(x) x[1:maxSp,])
      dataConstants$nsp <- maxSp
      dataSumm$occMatrix <- dataSumm$occMatrix[1:maxSp,]
    }

    ###################################################################
    # step 1 define the model code
    modelcode <- defineModel_MS(inclPhenology = inclPhenology, inclPanTrap = inclPanTrap) # this makes no difference!

    # step 2 create an operational from from NIMBLE/JAGS/BUGS code
    model <- nimbleModel(code = modelcode,
                         constants = dataConstants,
                         data = obsData,
                         inits = list(z = dataSumm$occMatrix,
                                      alpha.s = cloglog(dataSumm$naiveOcc),
                                      alpha.p = dataSumm$reportingRate, # replace with reportingRate_1 when I can calculate it
                                      beta1 = rep(180, dataConstants$nsp),
                                      beta2 = rep(50, dataConstants$nsp),
                                      phScale = rep(1, dataConstants$nsp),
                                      Multiplier = 1)
    )

    # step 3 build an MCMC object using buildMCMC(). we can add some customization here
    occMCMC <- buildMCMC(model,
                       monitors = c("mu.lambda", "psi.fs",
                                    #'alpha.s', "beta.s",
                                    #'alpha.p', "phScale","Multiplier",
                                    #"beta1", "beta2"
                                    Trend),
                       thin = 3,
                       useConjugacy = FALSE) # useConjugacy controls whether conjugate samplers are assigned when possible

  # step 3 before compiling the MCMC object we need to compile the model first
  Cmodel <- compileNimble(model) # NJBI: I don't understand why this step is necessary

  # now the MCMC (project = NIMBLE model already associated with a project)
  CoccMCMC <- compileNimble(occMCMC, project = model)
  # instantaneous

  # and now we can use either $run or runMCMC() on the compiled model object.
  runMCMC_samples <- runMCMC(CoccMCMC, nburnin = n.iter/2, niter = n.iter, nchains = 3, samplesAsCodaMCMC = T)

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
