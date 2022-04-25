\name{KM.calc.sample}
\Rdversion{1.1}
\alias{KM.calc.sample}
\alias{print.audit.plan.KM}
\title{
  Calculate sample size for KM-audit.
}
\description{
  Calculate the size of a sample needed to certify a correct election if a KM audit is planned.
}
\usage{
KM.calc.sample(Z, beta = 0.75, taint = 0, bound = c("e.plus", "WPM", "passed"))
\method{print}{audit.plan.KM}(x, ...)
}
\arguments{
  \item{Z}{
    elec.data object
  }
  \item{beta}{
Desired level of confidence.  This is 1-risk, where risk is the maximum chance
of not going to a full recount if the results are wrong.  Note that in
Stark's papers, the value of interest is typically risk, denoted $alpha$.
}
  \item{taint}{
  Assumed taint.  Taint is assumed to be the taint for all batches (very conservative).
 If taint=0 then we produce a good baseline.  
  }
  \item{bound}{
Type of bound on the maximum error one could find in a batch.
}
\item{x}{A audit.plan.KM object, such as one returned by
  KM.calc.sample.}
 \item{\dots}{ Unused.}
}
\value{
  A audit.plan.KM object.
}
\author{
  Based on the KM audit by Stark.
}
\seealso{
  KM.audit
}
\examples{
  data(santa.cruz)
  Z = elec.data( santa.cruz, C.names=c("danner","leopold") )
  KM.calc.sample( Z, beta=0.75, taint=0 )
}

