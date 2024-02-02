#' Format Metadata
#'
#' @param indat A dataset produced by the simulations
#' @param inclPanTrap should the model include pan trap data?
#' @param incl2ndTransect should the model include data from the second transect walk?
#' @return A one line dataframe
#' @export


formatMetadata <- function(indat,
                           incl2ndTransect = TRUE,
                           inclPanTrap = TRUE){
  nTr <- ifelse(incl2ndTransect, as.numeric(attr(indat, "transects")), 1)
  pans <- ifelse(inclPanTrap, attr(indat, "pantraps"), FALSE)

  list(simpars = data.frame(
        trend = attr(indat, "trend"),
        pantraps = pans,
        transects = nTr,
        rounds = attr(indat, "rounds"),
        sites = length(unique(indat$siteID)),
        years = attr(indat, "years"),
        sp_pool = attr(indat, "sp_pool")
      ),
      datastr = list(sp_obs = length(unique(indat$species)))
    )
}
