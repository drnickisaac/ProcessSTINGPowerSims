#' formatData
#'
#' @details Formats data ready for Nimble
#'
#' @param indata A dataset produced by the simulations
#' @return list of two data frames
#' @import reshape2
#' @export

formatData <- function(indata){

  castDat <- dcast(indata, year + round + siteID + jday + total_pantraps ~ "nsp",
                   value.var = "abundance", fun = length, fill = 0)

  dataConstants <- list(nyear = md$years,
                        nsp = md$sp_obs,
                        nsite = md$sites,
                        nvisit = nrow(castDat),
                        site = as.numeric(gsub(castDat$siteID, patt="site_", repl="")),
                        year = castDat$year,
                        round = castDat$round,
                        JulDate = castDat$jday,
                        nT = castDat$total_pantraps)

  # extract the observations
  ObsPan <- acast(indata, year + round + siteID ~ species,
                  value.var = "presences_pan", fill = 0)
  trCount1 <- acast(indata, year + round + siteID ~ species,
                    value.var = "obs", fill=0)
  trCount2 <- acast(indata, year + round + siteID ~ species,
                    value.var = "obs2", fill=0)

  # observations have to be transposed because we have coded species as the first dimension
  obsData <- list(y1 = t(ObsPan),
                  y2 = t(trCount1),
                  y3 = t(trCount2))

  return(list(dataConstants = dataConstants,
              obsData = obsData))
}
