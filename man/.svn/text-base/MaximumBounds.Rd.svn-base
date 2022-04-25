\name{maximumBounds}
\alias{MaximumBounds}
\alias{maximumMarginBound}
\alias{fractionOfVotesBound}

\title{ Election Audit Error Bound Functions }
\description{
  
Various bounding functions used to bound the maximum amount of error one
could see in a single audit unit.

}
\usage{
maximumMarginBound(Z, votes = NULL)
fractionOfVotesBound(Z, frac = 0.4)
}

\arguments{
  \item{Z}{ The elec.data object. }
  \item{votes}{ The data.frame to compute the maximumMarginBounds for.
  If null, will return all bounds for all precincts in Z.}
  \item{frac}{ Fraction of total votes that could be a winner
  overstatement/loser understatement.  So if the worst-case is a 20\% flip then enter 0.4 }
}
\details{
  maximumMarginBound
    return the maximum margin reduction for each precint by computing
	 all margin reductions between pairs of winners & losers and then
	scaling by that pair's total margin to get a proportion and then
	taking the max of all such proportions (usually will be the last 
	 winner to the closest loser).
         
         fractionOfVotesBound: WPM.  The maximum error of the unit is a fixed
         percentage of the total votes cast in the unit.  Typically the
         20\% WPM is used--meaning a swing of 40\% is the largest error
         possible as 20\% of the votes go from the winner to the loser.
}
\value{
 Vector (of length of precincts) of maximum possible error for 
    each precinct.
}

\author{ Luke W. Miratrix }



