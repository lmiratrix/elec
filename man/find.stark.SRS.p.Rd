\name{find.stark.SRS.p}
\alias{find.stark.SRS.p}

\title{ find.stark.SRS.p }
\description{
  Find the p-value for a given q, n, and N.  Helper function for a
  simple hypergeometric calculaton--see reports.
  
                                      
}
\usage{
find.stark.SRS.p(N, n, q)
}

\arguments{
  \item{N}{total number of precints }
  \item{n}{ total number of audited precints (must be less than N) }
  \item{q}{  min number of precints that could hold taint to flip election }
}

\value{
  Chance that 1 or more of the q 'bad' things will be seen in a size n
  SRS draw from the N sized bucket.
}


\author{ Luke W. Miratrix }

