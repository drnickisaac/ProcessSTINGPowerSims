#' runModel
#'
#' @details Runs the model.
#' @param dataConstants dataframe produced by ProcessSimDatFile()
#' @param obsData dataframe produced by ProcessSimDatFile()
#' @param dataSumm$stats dataframe produced by ProcessSimDatFile()
#' @param useNimble option to bypass the model fitting in Nimble (just for testing the code)
#' @param inclPhenology should the model account for seasonal variation?
#' @param inclPanTrap should the model include pan trap data?
#' @param multiSp should the model be run as a multispecies model, or many single-species models?
#' @param parallelize should the chains be run as separate processes on different cores?
#' @param allPars if `TRUE` then all model parameters are monitored. If `FALSE`, just `mu.lambda` and `Trend`.
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
                     dataSumm,
                     useNimble = TRUE,
                     inclPhenology = TRUE,
                     inclPanTrap = TRUE,
                     multiSp = TRUE,
                     parallelize = FALSE,
                     allPars = FALSE,
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
      dataSumm$occMatrix <- dataSumm$occMatrix[1:maxSp,,]
      dataSumm$stats <- dataSumm$stats[1:maxSp,]
      dataConstants$nsp <- maxSp
      print(paste('Warning: only the first', maxSp, 'will be used in modelling: others will be ignored'))
    }

    ###################################################################

    if(multiSp == TRUE){ # Multispecies option

      # step 1 define the model code
      modelcode <- defineModel_MS(inclPhenology = inclPhenology, inclPanTrap = inclPanTrap) # this makes no difference!

      # step 2 create an operational from from NIMBLE/BUGS code
      model <- nimbleModel(code = modelcode,
                           constants = dataConstants,
                           data = obsData,
                           inits = list(z = dataSumm$occMatrix,
                                        alpha.s = cloglog(dataSumm$stats$naiveOcc),
                                        alpha.p = dataSumm$stats$reportingRate, # replace with reportingRate_1 when I can calculate it
                                        beta1 = rep(180, dataConstants$nsp),
                                        beta2 = rep(50, dataConstants$nsp),
                                        #phScale = rep(1, dataConstants$nsp),
                                        Multiplier = 1,
                                        sd.eta = 2,
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
                           constants = dataConstants[!names(dataConstants) %in% "nsp"],
                           data = lapply(obsData, function(x) x[1,]), # values for species 1
                           inits = list(z = dataSumm$occMatrix[1,,], # value for species 1
                                        alpha.s = cloglog(dataSumm$stats$naiveOcc)[1], # value for species 1
                                        alpha.0 = dataSumm$stats$reportingRate[1], # value for species 1
                                        alpha.1 = 1,
                                        gamma.0 = -2,
                                        gamma.1 = 1,
                                        beta1 = 180,
                                        beta2 = 50,
                                        #phScale = 1,
                                        #Multiplier = 1,
                                        #sd.eta = 2,
                                        #eta = rnorm(n=dataConstants$nsite, mean=0, sd=2),
                                        Trend = rnorm(n=1))
                           )

      # step 3 build an MCMC object using buildMCMC(). we can add some customization here
      params <- c("mu.lambda","Trend")
      if(allPars) params <- c(params,
                              'alpha.s',
                              'alpha.0','alpha.1',
                              'gamma.0','gamma.1',
                              #"sd.eta"
                              #"phScale","Multiplier","alpha.p"
                              "beta1", "beta2")

      occMCMC <- buildMCMC(model,
                           monitors = params,
                           useConjugacy = FALSE) # useConjugacy controls whether conjugate samplers are assigned when possible

      # step 3 before compiling the MCMC object we need to compile the model first
      Cmodel <- compileNimble(model)

      # now the MCMC (project = NIMBLE model already associated with a project)
      CoccMCMC <- compileNimble(occMCMC, project = model)

      ####################################################################################

      single_species_model <- function(sp, spDat, dataSumm,
                                       n.iter, n.burn, n.thin,
                                       Cmodel, CoccMCMC){

        # add the data for the species of interest
        Cmodel$setData(spDat)

        # finish initialization
        Cmodel$setInits(list(z = dataSumm$occMatrix[sp,,],
                             alpha.s = cloglog(dataSumm$stats$naiveOcc)[sp],
                             alpha.p = dataSumm$stats$reportingRate[sp] # replace with reportingRate_1 when I can calculate it
        ))

        # test whether the model is fully initialised
        if(is.na(Cmodel$calculate())) {stop("model not fully initialized")}
        Cmodel$initializeInfo()

        # and now we can use $run on the compiled model object.
        samplesList <- list()
        for(i in 1:n.chain){
          CoccMCMC$run(niter = n.iter,
                       nburnin = n.burn,
                       chain = i,
                       thin = n.thin,
                       reset = TRUE)
          samplesList[[i]] <- as.matrix(CoccMCMC$mvSamples)
        }
        samplesList <- coda::as.mcmc.list(lapply(samplesList, as.mcmc))
      }

      ####################################################################################

      ####### run the model for each species

      if(parallelize){
        av_cores <- parallel::detectCores() - 1
        yearEff <- pbmcapply::pbmclapply(1:maxSp, function(i){
          single_species_model(sp=i,
                               spDat=lapply(obsData, function(x) x[i,]),
                               dataSumm=dataSumm,
                               n.iter = n.iter,
                               n.burn = n.burn,
                               n.thin = n.thin,
                               Cmodel, CoccMCMC)
        },
        mc.cores = av_cores
        )
      } else {
        yearEff <- lapply(1:maxSp, function(i){
          single_species_model(sp=i,
                               spDat=lapply(obsData, function(x) x[i,]),
                               dataSumm = dataSumm,
                               n.iter = n.iter,
                               n.burn = n.burn,
                               n.thin = n.thin,
                               Cmodel, CoccMCMC)
        }
        )
      }
      names(yearEff) <- dimnames(dataSumm$occMatrix)[[1]][1:maxSp]
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

