#' ProcessSimFolder
#'
#' @details
#' Takes the name of a folder and processes them in a large batch. The folder can be a zipped archive.
#' The files must be RDS format.
#'
#' @param inPath location of folder containing many simulation runs
#' @param outPath location of folder to store results
#' @param collateOutput should the output be stored for collation all together (`TRUE`) or written as separate files (`FALSE`)
#' @param useNimble option to bypass the model fitting in Nimble (just for testing the code)
#' @param inclPhenology should the model account for seasonal variation?
#' @param inclPanTrap should the model include pan trap data?
#' @param multiSp should the model be run for species separately, or in a multispecies model?
#' @param parallelize option to parallelize across MCMC chains
#' @param n.iter number of iterations for the Nimble model. Default is 1000.
#' @param n.burn number of iterations for the burn-in. If `NULL` (the default), then it will be set to `n.iter/2`.
#' @param n.thin thinning for the MCMC chains. Defaults to 5
#' @param n.chain number of MCMC chains. Defaults to 3
#' @param maxSp defines the maximum number of species to model. Species with numbers greater than this are ignored
#' @param maxFiles defines the maximum number of files to work with. Additional files are ignored
#' @return Nothing
#' @export


ProcessSimFolder <- function(inPath = ".",
                            outPath = "output",
                            collateOutput = FALSE,
                            useNimble = TRUE,
                            inclPhenology = TRUE,
                            inclPanTrap = TRUE,
                            multiSp = FALSE,
                            parallelize = FALSE,
                            n.iter = 1000,
                            n.burn = NULL,
                            n.thin = 5,
                            n.chain = 3,
                            maxSp = 9999,
                            maxFiles = 9999){

  ##############################################

  # first job: extract the names of the files
  # this will be different depending on whether inPath refers to a folder or zipped archive
  if(grepl(".zip", inPath))
    simfiles <- as.character(unzip(inPath, list = TRUE)$Name)
  else
    simfiles <- list.files(inPath)

  # in either case, restrict to the set that are rds files
  simfiles <- dir_files[grepl("\\.rds", dir_files)]

  # create an output directory
  if(!dir.exists(outPath)) dir.create(outPath)

  output <- lapply(simfiles[1:maxFiles],
                   ProcessSimDatFile,
                   zip = zip,
                   inPath = inPath,
                   outPath = ifelse(collateOutput, NULL, outPath),
                   useNimble = useNimble,
                   inclPhenology = inclPhenology,
                   inclPanTrap = inclPanTrap,
                   parallelize = parallelize,
                   n.iter = n.iter,
                   n.burn = n.burn,
                   n.thin = n.thin,
                   n.chain = n.chain,
                   multiSp = multiSp,
                   maxSp = maxSp)

  if(collateOutput){
    # we're going to take all the outputs together and format them into something pretty

     # extract the metadata
    names <- sapply(output, function(x) x$name)
    md <- cbind(names, t(sapply(output, function(x) x$md)))
    yrEff <- data.frame(cbind(names, t(sapply(output, function(x) x$modelEff))))

    timestamp <- format(Sys.time(), "%Y%m%d%H%M")
    write.csv(md, file = file.path(outPath, paste0("md",timestamp,".csv")))
    write.csv(yrEff, file = file.path(outPath, paste0("yrEffects",timestamp,".csv")))
  }
}
