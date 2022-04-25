\name{santa.cruz.audit}
\alias{santa.cruz.audit}
\docType{data}
\title{ Santa Cruz Election Data }
\description{
\code{santa.cruz} and \code{santa.cruz.audit}  hold data from a Santa Cruz County, CA, contest held in
November, 2008, for County Supervisor in the 1st District.  The
competitive candidates were John Leopold and Betty Danner.  According to
the semi-official results provided to us by the Santa Cruz County
Clerk's office, Leopold won with votes on 45\% of the 26,655 ballots.
Danner received the votes on 37\% of the ballots.  The remaining ballots were
undervoted, overvoted, or had votes for minor candidates.


}
\usage{data(santa.cruz.audit)}
\format{
A data frame with 16 observations on the following 4 variables.
  \describe{
    \item{\code{PID}}{Precinct IDs (unique) for all precincts involved
      in race}
    \item{\code{leopold}}{Total number of ballots marked for John Leopold.}
    \item{\code{danner}}{Total number of ballots marked for Betty
      Danner.}
    \item{\code{count}}{The number of times precinct was sampled in the
      PPEB sample taken.}    
  }

}
\details{
  \code{santa.cruz.audit} holds the audit totals for the random sample of
  precincts selected for the audit.  Note the \code{santa.cruz.audit} vote
  counts are larger for some precincts due the missing provisional
  ballot counts in the semi-official results.

}
\source{
Santa Cruz County, CA, Clerk Gail 
Pellerin, and their staffs, which we thank for their generous cooperation and the
onsiderable time and effort they spent counting ballots by hand in order
to collect these data.

}
\seealso{
  \link{santa.cruz}. For an illustration of analyzing this data, see the
  example in \code{\link{trinomial.bound}}.
}
\examples{
data(santa.cruz.audit)
data(santa.cruz)
santa.cruz = elec.data(santa.cruz, C.names=c("leopold","danner"))
trinomial.audit( santa.cruz, santa.cruz.audit )
}
\keyword{datasets}
