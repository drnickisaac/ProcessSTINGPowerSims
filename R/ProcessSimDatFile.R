#' ProcessSimDatFile
#'
#' @details
#' Processes one file produced by Nacho's simulations. Extracts metadata, formats the data for the model, then runs the model. The file must be in RDS format.
#'
#' @param filename the name of the datafile
#' @param inPath location where the `filename` can be found
#' @param outPath location where the results should be written. If `NULL` then results are returned rather than written.
#' @param useNimble option to bypass the model fitting in Nimble (just for testing the code)
#' @param inclPhenology should the model account for seasonal variation?
#' @param inclPanTrap should the model include pan trap data?
#' @param multiSp should the model be run for each species separately, or in a single model?
#' @param parallelize option to parallelize across MCMC chains
#' @param allPars if `TRUE` then all model parameters are monitored. If `FALSE`, just `mu.lambda` and `Trend`.
#' @param n.iter number of iterations for the Nimble model. Default is 1000.
#' @param n.burn number of iterations for the burn-in. If `NULL` (the default), then it will be set to `n.iter/2`.
#' @param n.thin thinning for the MCMC chains. Defaults to 5
#' @param n.chain number of MCMC chains. Defaults to 3
#' @param minSite the threshold minimum number of sites for a species to be considered for modelling
#' @param maxSp defines the maximum number of species to model. Species with numbers greater than this are ignored
#' @return if `outpath` is NULL then a list comprising model output and metadata. Otherwise nothing
#' @import reshape2
#' @export


ProcessSimDatFile <- function(filename,
                              inPath = ".",
                              outPath = NULL,
                              useNimble = TRUE,
                              inclPhenology = TRUE,
                              inclPanTrap = TRUE,
                              multiSp = FALSE,
                              parallelize = FALSE,
                              allPars = FALSE,
                              n.iter = 1000,
                              n.burn = NULL,
                              n.thin = 5,
                              n.chain = 3,
                              minSite = 1,
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

  # now restrict the data to species that occur on at least `minSite` sites
  if(minSite > 1){
    site_sp <- reshape2::acast(indata, siteID ~ species,
                               value.var = "obs",
                               fun = function(x) length(x) > 0, fill = 0)
    sp_n_Site <- colSums(site_sp)
    sp2incl <- names(sp_n_Site[sp_n_Site > minSite])
    indata <- subset(indata, species %in% sp2incl)
    nExcl <- length(sp_n_Site) - length(sp2incl)
    print(paste('Note:',nExcl,'species out of', length(sp_n_Site), 'have been excluded because they occur on fewer than', minSite, 'sites'))
    print(paste('We proceed to modelling with', length(sp2incl), 'species'))
  }

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
                       multiSp = multiSp,
                       parallelize = parallelize,
                       allPars = allPars,
                       n.iter = n.iter,
                       n.burn = n.burn,
                       n.thin = n.thin,
                       n.chain = n.chain,
                       maxSp = formattedData$md$maxSp)

  #clean up the name of the file. Remove the folder names and file suffix
  name <- gsub(filename, patt = "\\.rds", repl = "")
  name <- gsub(name, patt = "scenario_", repl = "")
  name <- strsplit(name, "/")[[1]]
  name <- name[length(name)]
  name <- paste0(name, "_Res_",maxSp,"Sp_",n.iter,"it.rds")

  # finish up an complete the job
  output <- list(dataFileName = filename,
                 outFileName = name,
                 md = formattedData$md,
                 dataSummStats = dataSumm$stats,
                 modelEff = modelEff,
                 runSettings = c(inclPhenology = inclPhenology,
                                   inclPanTrap = inclPanTrap,
                                   multiSp = multiSp,
                                   parallelize = parallelize,
                                   n.iter = n.iter,
                                   n.burn = n.burn,
                                   n.thin = n.thin,
                                   n.chain = n.chain),
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
