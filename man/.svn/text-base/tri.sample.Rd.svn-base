\name{tri.sample}
\alias{tri.sample}
\alias{tri.sample.stats}

\title{ Sample from List of Precincts PPEB }
\description{
  tri.sample selects a sample of precincts PPEB.  Namely, samples n times, with replacement, from
  the precincts proportional to the weights of the precincts.  
}
\usage{
tri.sample(Z, n, seed = NULL, print.trail = FALSE, simplify = TRUE,
          return.precincts = TRUE, PID = "PID", known = "known")
tri.sample.stats(samp)
}

\arguments{
  \item{Z}{ elec.data object }
  \item{n}{ Either a audit.plan.tri object (that contains n) or an
  integer which is the size of the sample }
  \item{seed}{ Seed to use. }
  \item{print.trail}{ Print diagnostics and info on the selection
  process. } 
  \item{simplify}{ If TRUE, return a data frame of unique precincts sampled, with
  counts of how many times they were sampled.  Otherwise return
  repeatedly sampled
  precincts seperately. }
  \item{return.precincts}{ Return the precincts, or just the precint IDs
  }
  \item{PID}{ The name of the column in Z\$V holding unique precinct IDs }
  \item{known}{ Name of column in Z\$V of TRUE/FALSE, where TRUE are
  precincts that are considered ``known'', and thus should not be
  sampled for whatever reason. }
  \item{samp}{ A sample, such as one returned from tri.sample }
}
\details{
   The weights, if passed, are in the ``e.max'' column of Z\$V.
}
\value{
  tri.sample returns a sample of precincts.  tri.sample.stats is a
  utility function returning the total number of unique precincts and
  ballots given a sample.
}
\seealso{
  \code{\link{trinomial.bound}}
  \code{\link{elec.data}}
  \code{\link{tri.calc.sample}}
  \code{\link{audit.plan.tri}}
}
\examples{
data(santa.cruz)
Z = elec.data( santa.cruz, C.names=c("danner","leopold") )
samp = tri.calc.sample( Z, beta=0.75, guess.N = 10, p_d = 0.05,
               swing=10, power=0.9, bound="e.plus" )
tri.sample( Z, samp, seed=541227 )
}


\author{ Luke W. Miratrix }

