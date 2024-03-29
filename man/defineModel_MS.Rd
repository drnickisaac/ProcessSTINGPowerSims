% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defineModel_MS.R
\name{defineModel_MS}
\alias{defineModel_MS}
\title{defineModel_MS}
\usage{
defineModel_MS(
  incl2ndTransect = TRUE,
  inclPanTrap = TRUE,
  inclPhenology = TRUE,
  inclStateRE = FALSE
)
}
\arguments{
\item{incl2ndTransect}{should the model include data from the second transect walk?}

\item{inclPanTrap}{should the model include pan trap data?}

\item{inclPhenology}{should the model account for seasonal variation?}

\item{inclStateRE}{should there be a site-level random effect in the state model?}
}
\value{
a set of code
}
\description{
defineModel_MS
}
\details{
Defines the Nimble model.
}
