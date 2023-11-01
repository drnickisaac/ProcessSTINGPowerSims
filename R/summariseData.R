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
  temp <- data.frame(cbind(site=dataConstants$site, t(obsData$y1)))
  temp <- melt(temp, id=1) %>%
    mutate(value = value>0) %>%
    group_by(site, variable) %>%
    summarise(occ = max(value))
  occSites <- acast(as.data.frame(temp), site ~ variable, value.var = "occ")

  return(list(
    occMatrix = t(occSites),
    naiveOcc = colMeans(occSites),
    reportingRate = colMeans(with(obsData, y1/5)), # per pan trap
    meanCount = colMeans(obsData$y2), # not counting the second transect
    reportingRate_1 = NULL#mean(with(obsData, y1/5)[occSites,]), # per pan trap. Need to coerce to same shape as above
    #meanCount_1 = NULL#mean(obsData$y2[occSites,])
  ))
}
