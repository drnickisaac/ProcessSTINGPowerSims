#' summariseData
#'
#' @details Summary statistics on the input data
#' @param obsData dataframe produced by ProcessSimDatFile()
#' @param dataConstants dataframe produced by ProcessSimDatFile()
#' @return summary statistics
#' @import dplyr
#' @import reshape2
#' @export

summariseData <- function(obsData, dataConstants){

  # which sites are occupied?
  temp <- data.frame(cbind(site=dataConstants$site, year=dataConstants$year, t(obsData$y1)))
  temp <- melt(temp, id=1:2) %>%
    mutate(value = value>0) %>%
    group_by(site, year, variable) %>%
    summarise(occ = max(value))
  occMatrix <- acast(as.data.frame(temp), variable ~ site ~ year, value.var = "occ")

  occSpSite <- apply(occMatrix,c(1,2),max)

  return(list(
    occMatrix = occMatrix,
    naiveOcc = apply(occSpSite, 1, mean),
    reportingRate = rowMeans(obsData$y1/5), # per pan trap
    meanCount = rowMeans(obsData$y2), # not counting the second transect
    reportingRate_1 = NULL#mean(with(obsData, y1/5)[occSites,]), # per pan trap. Need to coerce to same shape as above
    #meanCount_1 = NULL#mean(obsData$y2[occSites,])
  ))
}
