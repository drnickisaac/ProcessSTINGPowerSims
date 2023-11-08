#' runModel
#'
#' @details Runs the model.
#' @param dataConstants dataframe produced by ProcessSimDatFile()
#' @param obsData dataframe produced by ProcessSimDatFile()
#' @param dataSummStats dataframe produced by ProcessSimDatFile()
#' @param useNimble option to bypass the model fitting in Nimble (just for testing the code)
#' @param inclPhenology should the model account for seasonal variation?
#' @param inclPanTrap should the model include pan trap data?
#' @param multiSp should the model be run as a multispecies model, or many single-species models?
#' @param parallelize should the chains be run as separate processes on different cores?
#' @param n.iter number of iterations for the Nimble model. Default is 1000.
#' @param n.burn number of iterations for the burn-in. If `NULL` (the default), then it will be set to `n.iter/2`.
#' @param n.thin thinning for the MCMC chains. Defaults to 5
#' @param n.chain number of MCMC chains. Defaults to 3
#' @param maxSp maximum number of species to be modelled
#' @return a set of year effects
#' @export
#' @import nimble
#' @import pbmcapply
#' @import parallel
#' @import coda

runModel <- function(dataConstants,
                     obsData,
                     dataSummStats,
                     useNimble = TRUE,
                     inclPhenology = TRUE,
                     inclPanTrap = TRUE,
                     multiSp = TRUE,
                     parallelize = FALSE,
                     n.iter = 1000,
                     n.burn = NULL,
                     n.thin = 5,
                     n.chain = 3,
                     maxSp = 9999){

  ###################################################################

if(useNimble) {
  if(is.null(n.burn)) n.burn = n.iter/2

  ###################################################################
  # truncate the dataset if there are too many species
  if(dim(obsData$y1)[1] > maxSp){
    obsData <- lapply(obsData, function(x) x[1:maxSp,])
    dataConstants$nsp <- maxSp
    dataSummStats$occMatrix <- dataSummStats$occMatrix[1:maxSp,,]
    dataSummStats$naiveOcc <- dataSummStats$naiveOcc[1:maxSp]
    dataSummStats$reportingRate <- dataSummStats$reportingRate[1:maxSp]
  }

  ###################################################################

  if(multiSp == TRUE){ # Multispecies option

    # step 1 define the model code
    modelcode <- defineModel_MS(inclPhenology = inclPhenology, inclPanTrap = inclPanTrap) # this makes no difference!

    # step 2 create an operational from from NIMBLE/BUGS code
    model <- nimbleModel(code = modelcode,
                         constants = dataConstants,
                         data = obsData,
                         inits = list(z = dataSummStats$occMatrix,
                                      alpha.s = cloglog(dataSummStats$naiveOcc),
                                      alpha.p = dataSummStats$reportingRate, # replace with reportingRate_1 when I can calculate it
                                      beta1 = rep(180, dataConstants$nsp),
                                      beta2 = rep(50, dataConstants$nsp),
                                      phScale = rep(1, dataConstants$nsp),
                                      Multiplier = 1,
                                      tau.eta = 0.2,
                                      eta = rnorm(n=dataConstants$nsite, mean=0, sd=2),
                                      tau.trend = abs(rt(1, 1)),
                                      Trend = rnorm(n=1))
    )

    # step 3 build an MCMC object using buildMCMC(). we can add some customization here
    occMCMC <- buildMCMC(model,
                       monitors = c("Trend"),
                       thin = n.thin,
                       useConjugacy = FALSE) # useConjugacy controls whether conjugate samplers are assigned when possible

    # step 3 before compiling the MCMC object we need to compile the model first
    Cmodel <- compileNimble(model)

    # now the MCMC (project = NIMBLE model already associated with a project)
    CoccMCMC <- compileNimble(occMCMC, project = model)

    # and now we can use either $run or runMCMC() on the compiled model object.
    if(parallelize){
      av_cores <- parallel::detectCores() - 1
      runMCMC_samples <- pbmcapply::pbmclapply(1:n.chain, function(i)
                                runMCMC(
                                  mcmc = CoccMCMC,
                                  nburnin = n.burn,
                                  niter = n.iter,
                                  nchains = 1, samplesAsCodaMCMC = T),
                                mc.cores = av_cores)

    } else {
      runMCMC_samples <- runMCMC(CoccMCMC,
                               nburnin = n.burn,
                               niter = n.iter,
                               nchains = n.chain, samplesAsCodaMCMC = T)
    }
    yearEff <- runMCMC_samples

  ############################################ end multispecies

    } else { # sequential single-species option
      # step 1 define the model code
      modelcode <- defineModel_SS(inclPhenology = inclPhenology, inclPanTrap = inclPanTrap)

      # step 2 create an operational from from NIMBLE/BUGS code
      model <- nimbleModel(code = modelcode,
                           inits = list(beta1 = 180,
                                        beta2 = 50,
                                        phScale = 1,
                                        Multiplier = 1,
                                        tau.eta = 0.2,
                                        eta = rnorm(n=dataConstants$nsite, mean=0, sd=2),
                                        Trend = rnorm(n=1)),
                           constants = dataConstants[!names(dataConstants) %in% "nsp"])

      # step 3 build an MCMC object using buildMCMC(). we can add some customization here
      occMCMC <- buildMCMC(model,
                           monitors = c("mu.lambda",
                                        "Trend"),
                           useConjugacy = FALSE) # useConjugacy controls whether conjugate samplers are assigned when possible

      # step 3 before compiling the MCMC object we need to compile the model first
      Cmodel <- compileNimble(model)

      # now the MCMC (project = NIMBLE model already associated with a project)
      CoccMCMC <- compileNimble(occMCMC, project = model)

  #######

      single_species_model <- function(sp, spDat, dataSummStats, n.iter, Cmodel, CoccMCMC){

        # add the data
        Cmodel$setData(spDat)

        # finish initialization
        Cmodel$setInits(list(z = dataSummStats$occMatrix[sp],
                        alpha.s = cloglog(dataSummStats$naiveOcc)[sp],
                        alpha.p = dataSummStats$reportingRate[sp] # replace with reportingRate_1 when I can calculate it
                        ))

        # and now we can use $run on the compiled model object.
        samplesList <- list()
        for(i in 1:n.chain){
          CoccMCMC$run(niter = n.iter,
                       nburnin = n.burn,
                       chain = i,
                       thin = n.thin,
                       reset = TRUE,
                       time = TRUE)
          tmp <- as.matrix(CoccMCMC$mvSamples)
          samplesList[[i]] <- tmp
        }
        samplesList <- coda::as.mcmc.list(lapply(samplesList, as.mcmc))
      }

  ####### run the model for each species

    if(parallelize){
        av_cores <- parallel::detectCores() - 1
        yearEff <- pbmcapply::pbmclapply(1:maxSp, function(i){
          single_species_model(sp=i,
                               spDat=lapply(obsData, function(x) x[i,]),
                               dataSummStats=dataSummStats,
                               n.iter=n.iter,
                               Cmodel, CoccMCMC)
          },
          mc.cores = av_cores
        )
      } else {
        yearEff <- lapply(1:maxSp, function(i){
          single_species_model(sp=i,
                               spDat=lapply(obsData, function(x) x[i,]),
                               dataSummStats=dataSummStats,
                               n.iter=n.iter,
                               Cmodel, CoccMCMC)
          }
        )
      }
      names(yearEff) <- dimnames(dataSummStats$occMatrix)[[1]][1:maxSp]
    }

  #####################################################################

  }
  else {
    # for simplicity, let's just report the annual fitted total count
    totalObs <- rowSums(sapply(obsData, colSums))
    mod <- with(dataConstants, glm(totalObs ~ factor(year) + site, family = "poisson"))

    yearEff <- c(1, exp(coef(mod)[grep("factor", names(coef(mod)))]))
    names(yearEff) <- paste0("Year", 1:dataConstants$nyear)
  }

  return(yearEff)
}

