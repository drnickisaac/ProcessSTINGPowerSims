#' formatData
#'
#' @details Formats data ready for Nimble
#'
#' @param inData A dataset produced by the simulations
#' @param minSite the threshold minimum number of sites for a species to be considered for modelling
#' @return list of two data frames
#' @import reshape2
#' @export

formatData <- function(inData, minSite){

  castDat <- dcast(inData, year + round + siteID + jday + total_pantraps ~ "nsp",
                   value.var = "abundance", fun = length, fill = 0)

  md <- formatMetadata(inData)

  # now restrict the data to species that occur on at least `minSite` sites
  if(minSite > 1){
    site_sp <- reshape2::acast(inData, siteID ~ species,
                               value.var = "obs",
                               fun = function(x) length(x) > 0, fill = 0)
    sp_n_Site <- colSums(site_sp)
    sp2incl <- names(sp_n_Site[sp_n_Site > minSite])
    inData <- subset(inData, species %in% sp2incl)
    nExcl <- length(sp_n_Site) - length(sp2incl)
    print(paste('Note:',nExcl,'species out of', length(sp_n_Site), 'have been excluded because they occur on fewer than', minSite, 'sites'))
    print(paste('We proceed to modelling with', length(sp2incl), 'species'))
  }

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
              obsData = obsData,
              md = md))
}
