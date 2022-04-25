\name{find.q}
\alias{find.q}

\title{ find.q }
\description{
 Find q, the minimum number of precints with w\_p's greater than given t.stat
   that can hold an entire election shift in them.
   
   This number is behind the SRS methods such as CAST.  If we know how
   many precincts, at minimum, would have to hold substantial error in
   order to have the reported outcome be wrong, we can compute the
   chance of finding at least one such precinct given a SRS draw of size n.
}
\usage{
find.q(V, t.stat, bound.col, M, threshold = 1,
          w_p = weight.function("no.weight"), drop = NULL)
}

\arguments{
  \item{V}{ The data.frame of votes--the subwing of a elec.data object,
  usually. }
  \item{t.stat}{ The worst error found in the audit (weighted, etc.) }
  \item{bound.col}{The name of the column in V to be used for the passed size 
                     (max number of votes, total votes, incl undervotes, etc.) to the error 
                     function. }
  \item{M}{ The margin to close.  Usually 1 for proportional.  Can be
  less if error from other sources is assumed. }
  \item{threshold}{ The total amount of error to pack in the set of tainted precincts}
  \item{w_p}{ The weight function for errors. }
  \item{drop}{ Drop precincts with this column having a "true" value--they are
                     previously audited or otherwise known, and thus can't hold error.
                     Can also pass a logical T/F vector of the length of nrow(V)}
}
\details{
  
  Find the number of precints that need to have "large taint" in order to
   flip the election.  This is, essentially, finding a collection of precints
   such that the max error (e.max) plus the background error (the w\_p-inverse of the
   t.stat) for the rest of the precints is greater than the margin (or 1 if done
   by proportions).
  
  


}
\value{
   integer, number of badly tainted precints needed to hold 'threshold' error
}

\author{ Luke W. Miratrix }

