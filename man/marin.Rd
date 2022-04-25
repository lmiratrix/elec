\name{marin}
\Rdversion{1.1}
\alias{marin}
\docType{data}
\title{
Marin Measure B Reported Results
}
\description{
These are the reported vote totals from the 2009 election in Marin, CA
for Measure B.
}
\usage{data(marin)}
\format{
  A data frame with 544 observations on the following 5 variables.
  \describe{
    \item{\code{PID}}{Batch ID}
    \item{\code{strata}}{There are two levels, \code{ST-IB}
      \code{ST-VBM} for in-precinct and Vote-by-Mail.}
    \item{\code{tot.votes}}{total ballots cast in the batch.}
    \item{\code{Yes}}{Number recorded for Yes}
    \item{\code{No}}{Number recorded for No}
  }
}
\details{
Note the vote totals for the VBM strata are made up.  The batches are
the ``Decks'', which could not be individually tallied with ease.  The
work-around
was complex.  See the references, below.
}
\source{
  Marin, CA 2009 reported election results.
  
}
\references{
See J. L. Hall, L. W. Miratrix, P. B. Stark, M. Briones, E. Ginnold, F. Oakley, M. Peaden, G. Pellerin, T. Stanionis, 
and T. Webber.  Implementing risk-limiting audits in california.  USENIX EVT/WOTE in press, July 2009.

}
\examples{
data(marin)
marin = elec.data( marin, C.names=c("Yes","No") )

# Hand fixing error bound due to unknown
# vote totals in the VBM decks
marin$V$e.max = maximumMarginBound(marin)
sum( marin$V$e.max )   # 7.128
vbm = marin$V$strata=="ST-VBM"
marin$V[ vbm, "e.max" ] = 2 * marin$V[ vbm, "tot.votes" ] / marin$margin

sum( marin$V$e.max )   # 9.782

}
\keyword{datasets}
