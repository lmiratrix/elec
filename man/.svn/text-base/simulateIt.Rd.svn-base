\name{simulateIt}
\Rdversion{1.1}
\alias{simulateIt}

\title{
simulate KM audits
}
\description{
This takes an election and a truth and conducts a KM audit.
}
\usage{
simulateIt(data, M = 50, alpha = 0.25, plot = FALSE, debug = FALSE, return.Ps = FALSE, truncate.Ps = TRUE)
}

\arguments{
  \item{data}{
a data frame, one row per patch, with: tot.votes, e.max, taint
}
  \item{M}{
the maximum number of samples to draw before automatically escalating to 
        a full recount.
        
        
}
  \item{alpha}{
level of risk.
}
  \item{plot}{
plot a chart?
}
  \item{debug}{
    debug diag printed?
    
}
  \item{return.Ps}{
Return the sequence of p-values all the way up to N.



}
  \item{truncate.Ps}{
Return Ps only up to where audit stopped.
}
}
\details{

 Given a list of all precincts and their true taints and their sampling weights (in data, 
 a data.frame), do a sequential audit at the specified alpha.
}
\value{ stopPt - number of draws drawn
         n - number of unique precincts audited
}
