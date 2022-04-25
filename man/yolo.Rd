\name{yolo}
\alias{yolo}
\docType{data}
\title{ Yolo County, CA Election Data }
\description{
This is for measure W in Yolo County, CA, November 2008.  The file includes
precinct-level reports.
}

\usage{data(yolo)}
\format{
  A data frame with 114 observations on the following 8 variables.
  \describe{
    \item{\code{PID}}{ Unique identifier for the batches of ballots}
    \item{\code{Pct}}{The precinct id of the batch}
    \item{\code{how}}{Vote by mail (VBM) or walk-in (PCT)}
    \item{\code{b}}{Number of votes cast in that unit}
    \item{\code{under}}{Number of undervotes (ballots not voted).}
    \item{\code{over}}{Number of overvotes (where someone marked both
      yes and no).}
    \item{\code{y}}{Reported number of valid ballots marked yes.}
    \item{\code{n}}{Reported number of valid ballots marked no.}
  }
}
\details{
In the actual audit, 6 precincts were selected (see example) and audited
by hand-to-eye count by a group of 4 people cross-checking each other.
One of the 6 batches had underreported the "yes" votes by 1, and one had
overreported the "yes" votes by 1.  There were no other errors.
}
\source{
Yolo County, CA.  Special thanks to Freddie Oakley and Tom Stanionis.
}
\references{
See Stark et al. for papers using this data to illustrate risk-limiting
audits of election data.
}
\examples{
# Make an elec.data object out of precicnt-level results
data(yolo)
yolo = make.Z( yolo, C.names=c("y","n","under","over"), tot.votes.col="b" ) 

# Look at different sample sizes and cuts for setting aside
# small precincts
CAST.calc.opt.cut( yolo, beta=0.75, stages=1, t=5, plot=TRUE )

print( yolo )

# Get details of the audit plan -- expected work, etc.
ap <- CAST.calc.sample( yolo, beta=0.75, stages=1, t=5, small.cut=5 )
print( ap )

# Draw a sample (seed not used for actual audit)
CAST.sample(yolo, ap, seed=12345678)


}
\keyword{datasets}
