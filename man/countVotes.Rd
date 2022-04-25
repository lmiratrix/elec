\name{countVotes}
\alias{countVotes}

\title{ countVotes }
\description{
Given a elec.data object, count the votes as reported and determine
winner(s) and loser(s).
}
\usage{
countVotes(Z)
}

\arguments{
  \item{Z}{ the elec.data object. }
}

\value{
  Updated 'Z' matrix with the total votes as components inside it.
}
\examples{
  Z = make.cartoon()
  ## Take away 20 percent of C1's votes.
  Z$V$C1 = Z$V$C1 * 0.8
  ## Count again to find winner.
  Z = countVotes(Z)
  Z
}
\author{ Luke W. Miratrix }


