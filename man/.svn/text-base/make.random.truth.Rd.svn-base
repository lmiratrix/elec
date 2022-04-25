\name{make.random.truth}
\Rdversion{1.1}
\alias{make.random.truth}

\title{
making fake truth for electios
}
\description{
Make a random truth that is with the reported outcome, but has random
error scattered throughout.
}
\usage{
make.random.truth(Z, p_d = 0.1, swing = 10, uniform = TRUE, seed = NULL, PID = "PID")
}

\arguments{
  \item{Z}{
elec.data object.  The original reported results.
}
  \item{p_d}{
chance a batch has error
}
  \item{swing}{
max amount of error in votes.
}
  \item{uniform}{
if yes, then error is from 1 to swing.  If no, then error is swing.
}
  \item{seed}{
random seed to ease replication
}
  \item{PID}{
which column has batch IDs.
}
}
\details{
 Given reported results (Z), make a new data.frame which is the
 truth (that can be 'audited' by looking at relevant precincts).
 
 This is the generic small error generation used in trinomial paper
 and elsewhere as a baseline "normal" mode of operations.



}
\value{
  # Return: elec.data object holding the 'truth'.
  
}
