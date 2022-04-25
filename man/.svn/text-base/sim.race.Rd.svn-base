\name{sim.race}
\alias{sim.race}

\title{ Simulate CAST audits to assess performance }
\description{
  Simulate a race (using the \code{\link{make.cartoon}} method) and run a CAST audit on that simulation.
  CAST is a system devised by Dr. Philip B., Stark, UC
Berkeley Department of Statistics.
}
\usage{
sim.race(n = 800, beta = 0.75, stages = 2,
         truth.maker = make.truth.opt.bad,
         print.trail = FALSE)
}

\arguments{
  \item{beta}{ the confidence level desired }
  \item{stages}{ number of auditing stages. Each stage will have the
    same confidence level, determined by a function of beta. }
  \item{print.trail}{ Print out diagnostics. }
  \item{n}{ Desired sample size. }
  \item{truth.maker}{ Function to generate "truth" }
}

\value{
  A vector of 3 numbers.  The first is the stage reached.  The second is
  the total number of precincts audited.  The third is 0 if the audit
  failed to certify (i.e. found large error in the final stage), and 1
  if the audit certified the election (did not find large error in the
  final stage).
}

\seealso{
  See \code{\link{CAST}} and 
  \code{\link{CAST.calc.opt.cut}} for methods regarding CAST audits.
Also see \code{\link{do.audit}}, \code{\link{make.sample}}, and
  \code{\link{make.truth}} for doing other simulation studies of this method.
}
\examples{
     ## See how many times the CAST method fails to catch a wrong
     ##  election in 20 trials.
     replicate( 20, sim.race( beta=0.75, stages=2, truth.maker=make.truth.opt.bad) )

     ## Now see how much work the CAST method does for typical elections.
     replicate( 20, sim.race( beta=0.75, stages=2, truth.maker=make.ok.truth) )
}
\references{ See http://www.stat.berkeley.edu/~stark/Vote/index.htm for relevant  information. }
\author{ Luke W. Miratrix }

