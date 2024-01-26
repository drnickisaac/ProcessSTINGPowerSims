#' formatData
#'
#' @details Formats data ready for Nimble
#' @param minSite the threshold minimum number of sites for a species to be considered for modelling
#' @param inData A dataset produced by the simulations
#' @param inclPanTrap should the model include pan trap data?
#' @param incl2ndTransect should the model include data from the second transect walk?
#' @return list of two data frames
#' @import reshape2
#' @export

formatData <- function(inData,
                       incl2ndTransect = TRUE,
                       inclPanTrap = TRUE,
                       minSite = 1){

  castDat <- dcast(inData, year + round + siteID + jday + total_pantraps ~ "nsp",
                   value.var = "abundance", fun = length, fill = 0)

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

  # create metadata object
  md <- formatMetadata(inData)
  md$datastr$sp_n_Site <- data.frame(species = names(sp_n_Site), nSite = as.numeric(sp_n_Site))
  md$settings <- list(sp_modelled = length(sp2incl),
                   minSite = minSite)

  dataConstants <- list(nsp = md$settings["sp_modelled"],
                        nsite = md$simpars$sites,
                        nvisit = nrow(castDat),
                        nyear = md$simpars$years,
                        year = castDat$year,
                        site = as.numeric(gsub(castDat$siteID, patt="site_", repl="")),
                        JulDate = castDat$jday,
                        nT = castDat$total_pantraps)

  # extract the observations and populate the obsData list
  obsData <- list()

  if(inclPanTrap) {
    ObsPan <- acast(inData, year + round + siteID ~ species,
                    value.var = "presences_pan", fill = 0)
    obsData$y1 <- t(ObsPan)[sp2incl,]
  }

  trCount1 <- acast(inData, year + round + siteID ~ species,
                    value.var = "obs_transect1", fill=0)
  obsData$y2 = t(trCount1)[sp2incl,]

  if(incl2ndTransect) {
    trCount2 <- acast(inData, year + round + siteID ~ species,
                      value.var = "obs_transect2", fill=0)
    obsData$y3 <- t(trCount2)[sp2incl,]
  }

  return(list(dataConstants = dataConstants,
              obsData = obsData,
              md = md))
}
