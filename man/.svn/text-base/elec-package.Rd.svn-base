\name{elec-package}
\alias{elec-package}
\alias{elec}
\docType{package}
\title{
Statistical Election Audits Package
}
\description{
 This is a bizzare collection of functions written to do various sorts
 of statistical election audits.  There are also functions to generate
 simulated voting data, and simulated ``truth'' so as to do simulations
 to check charactaristics of these methods.  The package includes two
 data sets consisting of actual reported voting results and, in one case, audit
 data for races held November, 2008, in California.
 
}

\details{
\tabular{ll}{
Package: \tab elec\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2009-01-14\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

This package is a bit of a mess.  Start with the \code{\link{elec.data}} data object
which holds votes, counts, margins, etc.  The simplest functions are the
CAST functions, such as \code{\link{CAST.calc.sample}}, which implement Philip
B. Stark's CAST method of election auditing.  Another collection of
functions, primarily \code{\link{tri.calc.sample}} and \code{\link{trinomial.bound}}, implement the
trinomial bound audit method designed by Luke W. Miratrix and Philip
B. Stark.

The examples primarily use a data set included in the package,
\code{\link{santa.cruz}} and \code{\link{santa.cruz.audit}}, which holds the ballot counts for a Santa
Cruz, CA race that we audited using these methods.  See
\code{\link{trinomial.bound}} for how these data were analyzed.  The
\code{\link{yolo}} data set
holds precinct level counts for a race in Yolo county. 

The package grew out of an earlier, messier (can you believe it?)
package that implemented general routines for election auditing.  Pieces
of this package are used by the aforementioned cleaner methods, but all
the individual functions are still there for specific uses, such as
making different tests.  Start with \link{stark.test}, which has an
index of the pieces in its ``see also'' section.
}
\author{
  Luke W. Miratrix
  
  Maintainer: Luke W. Miratrix <luke@vzvz.org>

}
\references{
  For general papers on election auditing see the list at
  http://www.stat.berkeley.edu/~stark/Vote/index.htm. 
  For the trinomial bound, see Luke W. Miratrix and Philip B. Stark. (2009) Election Audits
  using a Trinomial Bound (in press).
  For an overview of the races and the methods, see 
Joseph Lorenzo Hall1, Philip B. Stark, Luke W. Miratrix, Elaine Ginnold, Freddie 
Oakley, Tom Stanionis, and Gail Pellerin. (2009) Implementing
Risk-Limiting Audits in California (in press).

}

\keyword{ package }

