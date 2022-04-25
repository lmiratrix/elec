\name{AuditErrors}
\alias{AuditErrors}
\alias{calc.pairwise.e_p}
\alias{calc.overstatement.e_p}
\alias{compute.audit.errors}

\title{ Functions that Compute Error Levels Given Audit Data }
\description{
Calculate the error amounts for all precincts in Z that were
audited from the audit data, given as overstatement errors for all candidates.
}
\usage{
compute.audit.errors(Z, audit = NULL, calc.e_p = calc.pairwise.e_p, w_p =
                 weight.function("no.weight"), bound.col = "tot.votes",
                 err.override = NULL)
calc.pairwise.e_p(Z, audit=NULL, err.override = NULL)
calc.overstatement.e_p(Z)
}

\arguments{
  \item{Z}{ elec.data object }
  \item{err.override}{ Assume a baserate of this amount of error
    everywhere, ignoring audit data. If non-null, use this as the found error in votes rather than
	                      the actual errors found in the audit. }
  \item{bound.col}{This is the vector (in audit) containing the maximum number
   	                      of votes possible in the various precincts. }
  \item{calc.e_p}{ Calculate e\_p or take as given. }
  \item{w_p}{ The weight function to use to reweight the errors of
    precincts. }
  \item{audit}{ The audit object, if it is not in the Z object, or if
  some other object other than the one in the Z object is desired to be
  considered as the audit object.  Used by the simulation functions to
  generate errors for some fixed amount of error in conjunction with the
  err.override. }

}
\value{
   compute.audit.errors returns a new audit table from Z with two new columns, err and err.weighted, 
	            corresponding
	            to the errors found in each audited precinct before and after the 
	            weight function has been applied to them.
}
\details{
   compute.audit.errors uses the calc functions and the weight functions
   in a 1-2 combination.

   calc.pairwise.e\_p() is often used with an err.override for
   simulation studies and whatnot to see what a fixed vote impact would
   have on taints for trinomial.
}
\note{
  Z must have an audit component, or one must be passed, for this
   function to make sense!  Remember that audit objects have
   overstatements, NOT total votes for candidates.  With err.override
   being set this is less relevant as the actual votes are usually
   ignored.
 }
\seealso{
     See \code{\link{audit.totals.to.OS}} for a utility function that handles processing of
     audit data.
}
     
\author{ Luke W. Miratrix }


