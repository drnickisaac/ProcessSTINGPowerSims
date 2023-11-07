#' ProcessSimFiles
#'
#' @details
#' Takes the name of a folder and processes them in a large batch
#'
#' @param inPath location of folder containing many simulation runs
#' @param outPath location of folder to store results
#' @param collateOutput should the output be stored for collation all together (`TRUE`) or written as separate files (`FALSE`)
#' @param useNimble option to bypass the model fitting in Nimble (just for testing the code)
#' @param n.iter number of iterations for the Nimble model
#' @param inclPhenology should the model account for seasonal variation?
#' @param inclPanTrap should the model include pan trap data?
#' @param multiSp should the model be run for species separately, or in a multispecies model?
#' @param maxSp defines the maximum number of species to model. Species with numbers greater than this are ignored
#' @param maxFiles defines the maximum number of files to work with. Additional files are ignored
#' @param parallelize option to parallelize across MCMC chains
#' @return Nothing
#' from reshape2 @import dcast
#' @export


ProcessSimFiles <- function(inPath = ".",
                            outPath = "output",
                            collateOutput = TRUE,
                            useNimble = TRUE,
                            n.iter = 1000,
                            inclPhenology = TRUE,
                            inclPanTrap = TRUE,
                            multiSp = FALSE,
                            maxSp = 9999,
                            maxFiles = 9999,
                            parallelize = FALSE){

  simfiles <- list.files(inPath)
  if(length(simfiles) > maxFiles) simfiles <- simfiles[1:maxFiles]

  if(!dir.exists(outPath)) dir.create(outPath)

  # parallelize this bit
  if(collateOutput){
    # we're going to take all the outputs together and format them into something pretty
    output <- lapply(simfiles, ProcessSimDatFile,
                     inPath = inPath, outPath = NULL,
                     useNimble = useNimble,
                     maxSp = maxSp,
                     parallelize = parallelize,
                     multiSp = multiSp)

    # extract the metadata
    names <- sapply(output, function(x) x$name)
    md <- cbind(names, t(sapply(output, function(x) x$md)))
    yrEff <- data.frame(cbind(names, t(sapply(output, function(x) x$modelEff))))

    timestamp <- format(Sys.time(), "%Y%m%d%H%M")
    write.csv(md, file = file.path(outPath, paste0("md",timestamp,".csv")))
    write.csv(yrEff, file = file.path(outPath, paste0("yrEffects",timestamp,".csv")))

  } else {
    sapply(simfiles, ProcessSimDatFile,
           inPath = inPath, outPath = outPath,
           maxSp = maxSp,
           parallelize = parallelize)
  }
}
