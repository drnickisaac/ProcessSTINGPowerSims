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

  # melt the occupancy per species:site combo
  mOSS <- melt(occSpSite)
  names(mOSS) <- c("species", "site", "occ")

  ############## Reporting Rate
  RR <- if(is.null(obsData$y1)) 0.000001 else rowMeans(obsData$y1/5)

  RR1 <- if(is.null(obsData$y1)) NA else {
    temp <- melt(as.data.frame(cbind(site=dataConstants$site, t(obsData$y1/5))), id=1)
    names(temp)[c(2,3)] <- c("species", "TrapRate")
    as.data.frame(merge(mOSS, temp, all=TRUE) %>%
                    group_by(species) %>%
                    summarise(mean(trMean[occ == 1])))
  }

  ############### Mean transect count: for calculating stats
  if(is.null(obsData$y3)) {
    trMean <- obsData$y2
  } else trMean <- (obsData$y2 + obsData$y3)/2

  temp <- melt(as.data.frame(cbind(site=dataConstants$site, t(trMean))), id=1)
  names(temp)[2:3] <- c("species", "trMean")

  MC1<- as.data.frame(merge(mOSS, temp, all=TRUE) %>%
                        group_by(species) %>%
                        summarise(mean(trMean[occ == 1])))

  ##############

  stats <- data.frame(
    species = dimnames(occMatrix)[[1]],
    naiveOcc = apply(occSpSite, 1, mean),
    reportingRate = RR, # per pan trap
    meanCount = rowMeans(trMean), # mean count including zeros
    reportingRate_1 = as.numeric(RR1[,2]), # per pan trap. Need to coerce to same shape as above
    meanCount_1 = as.numeric(MC1[,2]) # mean count when observed (not on all visits to occupied sites)
  )
  return(list(
    occMatrix = occMatrix,
    stats = stats))
}
