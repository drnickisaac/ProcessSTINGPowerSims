#' formatData
#'
#' @details Formats data ready for Nimble
#' @param minSite the threshold minimum number of sites for a species to be considered for modelling
#' @param inData A dataset produced by the simulations
#' @param inclPanTrap should the model include pan trap data?
#' @param incl2ndTransect should the model include data from the second transect walk?
#' @param inclPhenology should the model account for seasonal variation?
#' @param minSite the threshold minimum number of sites for a species to be considered for modelling
#' @param maxSite defines a limit on the number of sites in the database
#' @return list of two data frames
#' @import reshape2
#' @export

formatData <- function(inData,
                       incl2ndTransect = TRUE,
                       inclPanTrap = TRUE,
                       inclPhenology = TRUE,
                       minSite = 1, maxSite = 999){

  # first perform some basic checks on the data
  if(any(!paste0("site_",1:attr(inData, "sites")) %in% inData$siteID)){
    missingSites <- setdiff(paste0("site_",1:attr(inData, "sites")), inData$siteID)
    siteNum2remove <- attr(inData, "sites")
    warning("At least one site has no records in the input data")
    # we need to renumber the sites
    # what we'll do is sequentially renumber the
    while(length(missingSites) > 0){
      inData$siteID <- gsub(inData$siteID,
                            pa = paste0("site_", siteNum2remove),
                            repl = missingSites[1])
      print(paste0("renaming site_", siteNum2remove, " as ", missingSites[1]))
      siteNum2remove <- siteNum2remove - 1
      missingSites <- missingSites[-1]
    }
  }
  if(any(!1:attr(inData, "years") %in% inData$year)){
    missingYear <- setdiff(1:attr(inData, "years"), inData$year)
    stop(paste0(missingYear, " has no records in the input data"))
  }
  #### data checks complete

  ### now limit the data to MaxSite whilst preserving the attributes that will be needed later.
  if(maxSite < length(unique(inData$siteID))){
    if(maxSite < 10) {maxSite <- 10}
    print(paste("Subsetting the dataset to", maxSite,"sites"))
    temp <- list(attr(inData, "trend"), attr(inData, "sp_pool"))
    inData <- subset(inData, siteID %in% paste0("site_",1:maxSite))
    attr(inData, "trend") <- temp[[1]]
    attr(inData, "sp_pool") <- temp[[2]]
  }

  castDat <- dcast(inData, year + round + siteID + jday + total_pantraps ~ "nsp",
                   value.var = "abundance", fun = length, fill = 0)

  # set the minimum number of sites, as there are some species that were never observed.
  if(minSite < 1) minSite <- 1

  # now restrict the data to species that occur on at least `minSite` sites
  # create a new variable, obs, that sum across all data types: was the species observed on this visit?
  inData$obs <- inData$obs_transect1 > 0
  if(incl2ndTransect) inData$obs <- apply(cbind(inData$obs, (inData$obs_transect2 > 0)),1,max)
  if(inclPanTrap) inData$obs <- apply(cbind(inData$obs, (inData$presences_pan > 0)),1,max)

  # apparent occupancy matrix across all species:site combos for all data types
  sp_site <- (acast(inData, species~siteID, value.var = "obs", function(x) max(x) > 0, fill = 0))

  sp_n_Site <- rowSums(sp_site)
  sp2incl <- which(sp_n_Site > minSite)
  nExcl <- length(sp_n_Site) - length(sp2incl)
  print(paste('Note:',nExcl,'species out of', length(sp_n_Site), 'have been excluded because they occur on', minSite, 'sites or fewer'))
  if(length(sp2incl) > 0) {
    print(paste('We proceed to modelling with', length(sp2incl), 'species'))
  } else {
    stop(paste0("There are no species with enough sites model"))
  }

  # create metadata object
  md <- formatMetadata(inData,
                       incl2ndTransect=incl2ndTransect,
                       inclPanTrap=inclPanTrap,
                       trueTrend = attr(inData, "trend"),
                       spPool = attr(inData, "sp_pool"))

  md$datastr$sp_n_Site <- data.frame(species = names(sp_n_Site), nSite = as.numeric(sp_n_Site))
  md$settings <- list(sp_modelled = length(sp2incl),
                   minSite = minSite)

  dataConstants <- list(nsp = as.numeric(md$settings["sp_modelled"]),
                        nsite = length(unique(inData$site)),
                        nvisit = nrow(castDat),
                        nyear = md$simpars$years,
                        year = castDat$year,
                        site = as.numeric(gsub(castDat$siteID, patt="site_", repl="")))

  if(inclPhenology){dataConstants$JulDate <- castDat$jday}

  # extract the observations and populate the obsData list
  obsData <- list()

  if(inclPanTrap) {
    dataConstants$nT <- castDat$total_pantraps
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
