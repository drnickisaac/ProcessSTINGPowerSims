#' formatData
#'
#' @details Formats data ready for Nimble
#'
#' @param inData A dataset produced by the simulations
#' @param minSite the threshold minimum number of sites for a species to be considered for modelling
#' @return list of two data frames
#' @import reshape2
#' @export

formatData <- function(inData, minSite = 1){

  castDat <- dcast(inData, year + round + siteID + jday + total_pantraps ~ "nsp",
                   value.var = "abundance", fun = length, fill = 0)

  md <- formatMetadata(inData)

  # set the minimum number of sites, as there are some species that were never observed.
  if(minSite < 1) minSite <- 1

  # now restrict the data to species that occur on at least `minSite` sites
  # sum across all three data types: was the species ever observed?
  inData$obs <- rowSums(inData[, c("presences_pan", "obs_transect1", "obs_transect2")], na.rm=TRUE) > 0

  # apparent occupancy matrix across all species:site combos for all data types
  sp_site <- (acast(inData, species~siteID, value.var = "obs", function(x) max(x) > 0, fill = 0))

  sp_n_Site <- rowSums(sp_site)
  sp2incl <- which(sp_n_Site > minSite)
  nExcl <- length(sp_n_Site) - length(sp2incl)
  print(paste('Note:',nExcl,'species out of', length(sp_n_Site), 'have been excluded because they occur on fewer than', minSite, 'sites'))
  print(paste('We proceed to modelling with', length(sp2incl), 'species'))
  md$sp_modelled <- length(sp2incl)

  md$minSite <- minSite

  dataConstants <- list(nsp = md$sp_modelled,
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
                    value.var = "obs_transect1", fill=0)
  trCount2 <- acast(inData, year + round + siteID ~ species,
                    value.var = "obs_transect2", fill=0) # not quite sure why, but this turns NAs into 0s

  # observations have to be transposed because we have coded species as the first dimension
  obsData <- list(y1 = t(ObsPan)[sp2incl,],
                  y2 = t(trCount1)[sp2incl,],
                  y3 = t(trCount2)[sp2incl,])

  return(list(dataConstants = dataConstants,
              obsData = obsData,
              md = md))
}
