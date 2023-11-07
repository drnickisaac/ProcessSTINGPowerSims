#' ProcessSimDatFile
#'
#' @details
#' Processes one file produced by Nacho's simulations. Extracts metadata, formats the data for the model, then runs the model
#'
#' @param filename the name of the datafile
#' @param inPath location where the `filename` can be found
#' @param outPath location where the results should be written. If `NULL` then results are returned rather than written.
#' @param useNimble option to bypass the model fitting in Nimble (just for testing the code)
#' @param n.iter number of iterations for the Nimble model
#' @param inclPhenology should the model account for seasonal variation?
#' @param inclPanTrap should the model include pan trap data?
#' @param maxSp defines the maximum number of species to model. Species with numbers greater than this are ignored
#' @param multiSp should the model be run for each species separately, or in a single model?
#' @param parallelize option to parallelize across MCMC chains
#' @return if `outpath` is NULL then a list comprising model output and metadata. Otherwise nothing
#' @export


ProcessSimDatFile <- function(filename,
                              inPath = ".",
                              outPath = NULL,
                              useNimble = TRUE,
                              n.iter = 1000,
                              inclPhenology = TRUE,
                              inclPanTrap = TRUE,
                              maxSp = 9999,
                              multiSp = FALSE,
                              parallelize = FALSE){

  #Read in the data
  fileNamePath <- file.path(inPath, filename)
  if(file.exists(fileNamePath))
     indata <- readRDS(fileNamePath)
  else
    stop(paste0(fileNamePath, "not found"))

  # format the data
  formattedData <- formatData(indata)

  # if appropriate, limit the number of species
  formattedData$md$maxSp <- min(maxSp, formattedData$md$sp_obs)

  # SUMMARISE the data
  dataSumm <- with(formattedData, summariseData(obsData, dataConstants))

  # run the model
  modelEff <- runModel(formattedData$dataConstants,
                       formattedData$obsData,
                       dataSumm = dataSumm,
                       useNimble = useNimble,
                       maxSp = formattedData$md$maxSp,
                       multiSp = multiSp,
                       parallelize = parallelize)

  if(!multiSp){
  # for the single species option ...
  # collate the trends for each species into a multispecies average
  }

  # finish up an complete the job
  name <- gsub(filename, patt = "\\.rds", repl = "")
  output <- list(name = name,
                md = formattedData$md,
                dataSumm = dataSumm,
                modelEff = modelEff,
                timestamp  = format(Sys.time(), "%y%m%d%H%M%S"))

  if(!is.null(outPath)){
    if(!dir.exists(outPath))
      dir.create(outPath)
    saveRDS(file = file.path(outPath, paste0(name, "_Result.rds")), object = output)
  } else {
    return(output)
  }
}
