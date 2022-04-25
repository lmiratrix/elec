\name{CAST}
\alias{CAST}
\alias{test.CAST}
\alias{CAST.calc.sample}
\alias{CAST.sample}
\alias{CAST.audit}

\title{ CAST Functions }
\description{
Collection of functions for planning and evaluating results of a CAST
election audit.  CAST is a system devised by Dr. Philip B., Stark, UC
Berkeley Department of Statistics.

\code{CAST.calc.sample} determines what size SRS sample should be drawn
to have a reasonable chance of certification if the election does not
have substantial error.  It returns an \code{audit.plan}.
\code{CAST.sample} takes the audit.plan and draws a sample to audit.
\code{CAST.audit} takes audit data (presumably from the audit of the
sample drawn in previous step) and analyzes it.  
}
\usage{
CAST.calc.sample(Z, beta = 0.9, stages = 1, t = 3, as.taint=FALSE,
            small.cut = NULL, strata = NULL, drop = NULL,
            method = c("select", "binomial", "hypergeometric"),
            calc.e.max = TRUE, bound.function = maximumMarginBound)
CAST.sample(Z, ns, strata = NULL, 
         seed = NULL, print.trail = FALSE,
         known = "known")
CAST.audit(Z, audit = NULL, plan = NULL, ...)
test.CAST()
}

\arguments{
  \item{Z}{ elec.data object }
  \item{beta}{ the confidence level desired }
  \item{stages}{ number of auditing stages. Each stage will have the
    same confidence level, determined by a function of beta.  A value of
  1 is a single-stage audit.}
  \item{t}{The maximum amount of error, in votes, expected.}
  \item{as.taint}{Boolean value.  TRUE means interpret $t$ as a taint in
  $[0,1]$ by batch (so the threshold error will be batch-specific).  FALSE means interpret $t$ as a proportion of the
               margin or as number of votes (as described above).}
  \item{small.cut}{ Cut-off in votes--any precincts with potential error
    smaller than this value will not be audited and be assumed to be worst
    case error. }
  \item{strata}{ Name of the stratification column of Z.  Not needed if
  audit plan also being passed in case of CAST.sample. NULL means single
  strata. }
  \item{drop}{ Vector of precincts to drop for whatever reasons (such as
    they are already known). }
  \item{method}{ Method of calculation. }
  \item{calc.e.max}{ Should the e.max be taken as given, or
    recalculated? }
  \item{bound.function}{ What function should be used to calculate
    worst-case potential error of precincts.  }
  \item{ns}{ EITHER an audit.plan or a vector of sample sizes for the strata.  Names must
    correspond ot the names of the strata.  If ns is an audit plan, then
  the strata variable should not be passed as well. }
  \item{seed}{ Seed to use--for reproducability. }
  \item{print.trail}{ Print out diagnostics. }
  \item{known}{ The column of known precincts that should thus not be
    selected.  Similar to "drop", above. }
  \item{plan}{ An audit.plan object that the audit was conducted
    under. }
  \item{audit}{ A data.matrix holding the audit data, if the Z object
  does not have one, or if it is desirable to override it.  If both the
  Z object has an audit object and audit is not null, it will use this
  parameter and ignore the one in Z. }
  \item{...}{ Passed to CAST.calc.sample if plan is null and needs to be
  regenerated. }
}

\seealso{
  \code{\link{elec.data}} for a description of the object that holds
  precinct-level vote records.  See \code{\link{tri.calc.sample}} for a
  PPEB auditing method.  See
  \code{\link{CAST.calc.opt.cut}} for calculating optimal cut-offs to
  keep needed sample size low.
Also see \code{\link{sim.race}}, \code{\link{do.audit}}, \code{\link{make.sample}}, and
  \code{\link{make.truth}} for doing simulation studies of this method.
}
\examples{
        ## Make an example cartoon race (from Stark paper)
	Z = make.cartoon()

        ## What should we do?
	samp.info = CAST.calc.sample( Z )
	samp.info

        ## Draw a sample.
	samp = CAST.sample( Z, samp.info$ns )
        samp

        ## Analyze what a CAST audit of santa cruz would entail
        data(santa.cruz)
        Z = elec.data( santa.cruz, C.names=c("leopold","danner") )
        CAST.calc.sample( Z, beta=0.75, stages=1, t=5, small.cut=60)
}
\references{ Philip B. Stark. CAST: Canvass Audits by Sampling and
        Testing. University of California at Berkeley Department of
        Statistics, 2009. URL:
        http://statistics.berkeley.edu/~stark/Preprints/cast09.pdf.
        Also see http://www.stat.berkeley.edu/~stark/Vote/index.htm for other
        relevant information. }
\author{ Luke W. Miratrix }

