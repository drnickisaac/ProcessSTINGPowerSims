#' ProcessSimDatFile
#'
#' @details
#' Processes one file produced by Nacho's simulations. Extracts metadata, formats the data for the model, then runs the model. The file must be in RDS format.
#'
#' @param filename the name of the datafile
#' @param inPath location where the `filename` can be found
#' @param outPath location where the results should be written. If `NULL` then results are returned rather than written.
#' @param useNimble option to bypass the model fitting in Nimble (just for testing the code)
#' @param inclPanTrap should the model include pan trap data?
#' @param incl2ndTransect should the model include data from the second transect walk?
#' @param inclPhenology should the model account for seasonal variation?
#' @param inclStateRE should there be a site-level random effect in the state model?
#' @param multiSp should the model be run for each species separately, or in a single model?
#' @param parallelize option to parallelize across MCMC chains
#' @param allPars if `TRUE` then all model parameters are monitored. If `FALSE`, just `mu.lambda` and `Trend`.
#' @param n.iter number of iterations for the Nimble model. Default is 1000.
#' @param n.burn number of iterations for the burn-in. If `NULL` (the default), then it will be set to `n.iter/2`.
#' @param n.thin thinning for the MCMC chains. Defaults to 5
#' @param n.chain number of MCMC chains. Defaults to 3
#' @param minSite the threshold minimum number of sites for a species to be considered for modelling
#' @param maxSite defines a limit on the number of sites in the database
#' @param maxSp defines the maximum number of species to model. Species with numbers greater than this are ignored
#' @return if `outpath` is NULL then a list comprising model output and metadata. Otherwise nothing
#' @export


ProcessSimDatFile <- function(filename,
                              inPath = ".",
                              outPath = NULL,
                              useNimble = TRUE,
                              inclPhenology = TRUE,
                              incl2ndTransect = TRUE,
                              inclPanTrap = TRUE,
                              inclStateRE = FALSE,
                              multiSp = FALSE,
                              parallelize = FALSE,
                              allPars = FALSE,
                              n.iter = 1000,
                              n.burn = NULL,
                              n.thin = 5,
                              n.chain = 3,
                              minSite = 1,
                              maxSite = 999,
                              maxSp = 9999){

#####################################################################


  #Read in the data; check the file exists then open it
  if(grepl("\\.zip", inPath)){
    con <- unz(inPath, filename, open = "rb")
    fileNamePath <- gzcon(con)
    indata <- readRDS(fileNamePath)
  } else {
    fileNamePath <- file.path(inPath, filename)
    if(file.exists(fileNamePath))
      indata <- readRDS(fileNamePath)
    else
      stop(paste0(fileNamePath, "not found"))
  }

  print(paste("Successfully read in", filename))

  # format the data (includes removing species found on few sites)
  formattedData <- formatData(indata, minSite = minSite, maxSite = maxSite,
                              inclPanTrap = inclPanTrap,
                              incl2ndTransect = incl2ndTransect)

  # if appropriate, limit the number of species
  numSp <- min(maxSp, formattedData$md$settings$sp_modelled)
  formattedData$md$settings$maxSp <- numSp

  # SUMMARISE the data
  dataSumm <- with(formattedData, summariseData(obsData, dataConstants))

  # run the model
  modelEff <- runModel(formattedData$dataConstants,
                       formattedData$obsData,
                       dataSumm = dataSumm,
                       useNimble = useNimble,
                       inclPanTrap = inclPanTrap,
                       incl2ndTransect = incl2ndTransect,
                       inclPhenology = inclPhenology,
                       inclStateRE = inclStateRE,
                       multiSp = multiSp,
                       parallelize = parallelize,
                       allPars = allPars,
                       n.iter = n.iter,
                       n.burn = n.burn,
                       n.thin = n.thin,
                       n.chain = n.chain,
                       maxSp = numSp)

  #clean up the name of the file. Remove the folder names and file suffix
  name <- gsub(filename, patt = "\\.rds", repl = "")
  name <- gsub(name, patt = "scenario_", repl = "")
  name <- strsplit(name, "/")[[1]]
  name <- name[length(name)]
  name <- paste0(name, "_Res_",numSp,"Sp_")
  if(useNimble) {
    name <- paste0(name, n.iter,"it.rds")
  } else {name <- paste0(name, "GLM.rds")}

  # finish up an complete the job
  output <- list(dataFileName = filename,
                 outFileName = name,
                 md = formattedData$md,
                 dataSummStats = dataSumm$stats,
                 modelEff = modelEff,
                 runSettings = list(data= c(incl2ndTransect = incl2ndTransect,
                                           inclPanTrap = inclPanTrap,
                                           minSite = minSite,
                                           numSp = numSp),
                                    model = c(
                                           useNimble = useNimble,
                                           inclPhenology = inclPhenology,
                                           inclStateRE = inclStateRE,
                                           multiSp = multiSp,
                                           parallelize = parallelize,
                                           n.iter = n.iter,
                                           n.burn = n.burn,
                                           n.thin = n.thin,
                                           n.chain = n.chain)),
                spDat = list(nVis = sapply(formattedData$obsData, function(x) rowSums(x>0)),
                             sumObs = sapply(formattedData$obsData, rowSums)),
                Version = packageVersion("ProcessSTINGPowerSims"),
                timestamp  = format(Sys.time(), "%y-%m-%d %H:%M:%S"))

  if(!is.null(outPath)){
    if(!dir.exists(outPath))
      dir.create(outPath)
    saveRDS(file = file.path(outPath, name), object = output)
    return(NULL)
  } else {
    return(output)
  }
}
