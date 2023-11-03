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
#' @return if `outpath` is NULL then a list comprising model output and metadata. Otherwise nothing
#' @export


ProcessSimDatFile <- function(filename, inPath = ".",
                              outPath = NULL,
                              useNimble = FALSE,
                              n.iter = 1000,
                              inclPhenology = TRUE,
                              inclPanTrap = TRUE,
                              maxSp = 9999){

  indata <- readRDS(file.path(inPath, filename))

  # get the metadata
  md <- formatMetadata(indata)
  if(maxSp < 9999) md$maxSp <- maxSp

  # format the data
  formattedData <- formatData(indata)

  # SUMMARISE the data
  dataSumm <- with(formattedData, summariseData(obsData, dataConstants))

  if(maxSp == 1){
    # for the single species option:
    #1) need to reduce `dataSumm` and `formattedData$obsData` to single species
    #2) loop runModel to run over each species in term
    #3) collate the trends for each species into a multispecies average
  }

  # run the model
  modelEff <- runModel(formattedData$dataConstants, formattedData$obsData,
                       dataSumm = dataSumm,
                       useNimble = useNimble, maxSp = maxSp)

  # finish up an complete the job
  output <- list(name = gsub(filename, patt = "\\.rds", repl = ""),
                md=md,
                dataSumm = dataSumm,
                modelEff = modelEff,
                timestamp  = format(Sys.time(), "%y%m%d%H%M%S"))

  if(!is.null(outPath)){
    saveRDS(file = file.path(outPath, filename), object = output)
  } else {
    return(output)
  }
}
