\name{make.opt.packed.bad}
\Rdversion{1.1}
\alias{make.opt.packed.bad}

\title{
make.truth.opt.bad  
}
\description{
Generate a ``truth'' that is optimally bad in the sense of the margin in
error is
packed into as few precints as possible.

}
\usage{
make.opt.packed.bad(Z, max.taint = 1, max.taint.good = max.taint, WPM = FALSE, add.good = 0, add.random = FALSE)
}

\arguments{
  \item{Z}{
elec.data object to make bad truth for.
}
  \item{max.taint}{
max taint for any batch
}
  \item{max.taint.good}{
max taint in good direction for any batch
}
  \item{WPM}{
Use WPM bound on error.
}
  \item{add.good}{
add this amount of margin in good error (i.e. for the winner)
}
  \item{add.random}{
add a random tweak to error
}
}
\details{
 Make an audit data.frame with the error being exactly 1 margin, and packed
 into a small number of precincts (with some potential for binding amount of
 error per precinct).

 Warning: error is not necessarily achievable as the discrete nature of 
 whole votes is disregarded.
}
\value{
Return the vote matrix (a data.frame) with tot.votes, e.max, and taint
computed (NOT the elec data object).

}
