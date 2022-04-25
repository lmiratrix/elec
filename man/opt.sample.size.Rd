\name{opt.sample.size}
\Rdversion{1.1}
\alias{opt.sample.size}

\title{
KM Audit Sample Size Calc
}
\description{
Calc KM Optimal Sample Size
}
\usage{
opt.sample.size(Z, beta = 0.25)
}

\arguments{
  \item{Z}{
elec.data object
}
  \item{beta}{
risk
}
}
\details{
This is how many steps would be needed if no error was found with each
  step.  Obviously a bit idealistic, but still useful.
}
\value{
  Single number of batches to sample.
}
