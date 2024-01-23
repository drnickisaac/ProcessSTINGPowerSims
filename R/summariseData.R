#' summariseData
#'
#' @details Summary statistics on the input data
#' @param obsData dataframe produced by ProcessSimDatFile()
#' @param dataConstants dataframe produced by ProcessSimDatFile()
#' @return summary statistics
#' @import dplyr
#' @import reshape2
#' @export

summariseData <- function(obsData, dataConstants,
                          inclPanTrap = TRUE,
                          incl2ndTransect = TRUE){

  # calculate for each visit whether the species was observed across 1-3 data types (depending on what exists)
  obs <- obsData$y2
  if(!is.null(obsData$y1)) obs <- obs + obsData$y1
  if(!is.null(obsData$y3)) obs <- obs + obsData$y3

  # which sites are occupied?
  temp <- data.frame(cbind(site=dataConstants$site, year=dataConstants$year, t(obs)))
  temp <- melt(temp, id=1:2) %>%
    mutate(value = value>0) %>%
    group_by(site, year, variable) %>%
    summarise(occ = max(value))
  occMatrix <- acast(as.data.frame(temp), variable ~ site ~ year, value.var = "occ", fill = 0)

  occSpSite <- apply(occMatrix,c(1,2),max)

  # Mean transect count: for calculating stats
  if(!is.null(obsData$y3))
    trMean_1 <- trMean <- obsData$y2
  else
    trMean_1 <- trMean <- (obsData$y2 + obsData$y3)/2
  trMean_1[trMean_1 == 0] <- NA

  stats <- data.frame(
    species = dimnames(occMatrix)[[1]],
    naiveOcc = apply(occSpSite, 1, mean),
    reportingRate = ifelse(is.null(obsData$y1), 0.0001, rowMeans(obsData$y1/5)), # per pan trap
    meanCount = rowMeans(trMean), # mean count including zeros
    reportingRate_1 = NA,#mean(with(obsData, y1/5)[occSites,]), # per pan trap. Need to coerce to same shape as above
    meanCount_1 = rowMeans(trMean_1, na.rm=T) # mean count when observed (not on all visits to occupied sites)
  )

  return(list(
    occMatrix = occMatrix,
    stats = stats))
}
