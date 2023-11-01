#' Format Metadata
#'
#' @param indat A dataset produced by the simulations
#' @return A one line dataframe
#' @export


formatMetadata <- function(indat){
  data.frame(
    trend = attr(indat, "trend"),
    pantraps = attr(indat, "pantraps"),
    transects = as.numeric(attr(indat, "transects")),
    rounds = attr(indat, "rounds"),
    sites = attr(indat, "sites"),
    years = attr(indat, "years"),
    sp_pool = attr(indat, "sp_pool"),
    sp_obs = length(unique(indat$species))
    )
}
