\name{santa.cruz}
\alias{santa.cruz}
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
\usage{data(santa.cruz)}
\format{
A data frame with 152 observations on the following 5 variables.
  \describe{
    \item{\code{PID}}{Precinct IDs (unique) for all precincts involved
      in race}
    \item{\code{r}}{Total number of registered voters in the precinct.}
    \item{\code{tot.votes}}{Total number of ballots cast in the precinct.}
    \item{\code{leopold}}{Total number of ballots marked for John Leopold.}
    \item{\code{danner}}{Total number of ballots marked for Betty
      Danner.}
  }

}
\details{
\code{santa.cruz} holds the semi-official results for the race.
  \code{santa.cruz.audit} holds the audit totals for the random sample of
  precincts selected for the audit.  Note the \code{santa.cruz.audit} vote
  counts are larger for some precincts due the missing provisional
  ballot counts in the semi-official results.

}
\seealso{
  \link{santa.cruz.audit}
}
\source{
Santa Cruz County, CA, Clerk Gail 
Pellerin, and their staff.

}
\examples{
data(santa.cruz)
elec.data( santa.cruz, C.names=c("danner","leopold") )
}
\keyword{datasets}
