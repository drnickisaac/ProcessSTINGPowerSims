#' Format Metadata
#'
#' @param indat A dataset produced by the simulations
#' @param inclPanTrap should the model include pan trap data?
#' @param incl2ndTransect should the model include data from the second transect walk?
#' @return A one line dataframe
#' @export


formatMetadata <- function(indat,
                           incl2ndTransect = TRUE,
                           inclPanTrap = TRUE,
                           trueTrend = NULL,
                           spPool = NULL){
  nTr <- ifelse(incl2ndTransect, as.numeric(attr(indat, "transects")), 1)
  pans <- ifelse(inclPanTrap, attr(indat, "pantraps"), FALSE)

  list(simpars = data.frame(
                trend = trueTrend,
                pantraps = pans,
                transects = nTr,
                rounds = length(unique(indat$round)),
                sites = length(unique(indat$siteID)),
                years = length(unique(indat$year)),
                sp_pool = spPool
                ),
  datastr = list(sp_obs = length(unique(indat$species)))
  )
}
