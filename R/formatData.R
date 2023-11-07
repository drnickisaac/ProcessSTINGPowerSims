#' formatData
#'
#' @details Formats data ready for Nimble
#'
#' @param inData A dataset produced by the simulations
#' @return list of two data frames
#' @import reshape2
#' @export

formatData <- function(inData){

  castDat <- dcast(inData, year + round + siteID + jday + total_pantraps ~ "nsp",
                   value.var = "abundance", fun = length, fill = 0)

  md <- formatMetadata(inData)

  dataConstants <- list(nsp = md$sp_obs,
                        nsite = md$sites,
                        nvisit = nrow(castDat),
                        nyear = md$years,
                        year = castDat$year,
                        site = as.numeric(gsub(castDat$siteID, patt="site_", repl="")),
                        JulDate = castDat$jday,
                        nT = castDat$total_pantraps)

  # extract the observations
  ObsPan <- acast(inData, year + round + siteID ~ species,
                  value.var = "presences_pan", fill = 0)
  trCount1 <- acast(inData, year + round + siteID ~ species,
                    value.var = "obs", fill=0)
  trCount2 <- acast(inData, year + round + siteID ~ species,
                    value.var = "obs2", fill=0)

  # observations have to be transposed because we have coded species as the first dimension
  obsData <- list(y1 = t(ObsPan),
                  y2 = t(trCount1),
                  y3 = t(trCount2))

  return(list(dataConstants = dataConstants,
              obsData = obsData))
}
