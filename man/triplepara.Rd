\name{triplepara}
\alias{triplepara}
\title{ Internal function }
\description{
This is an internal function used to calculate the loglikelihood of triples.
}
\usage{
triplepara(inode,jnode,nodematrix,nspecies)
}
\arguments{
  \item{inode}{ the decendant node in the triple }
 \item{jnode}{ the ancestral node in the triple }
\item{nodematrix}{the species tree}
\item{nspecies}{the number of species}
}
\details{This function is used to calculate triple likelihoods.
 }
\value{
  The function returns the theta and gamma in a triple.
}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\keyword{IO}
\seealso{ \code{\link{write.subtree}}, \code{\link{read.tree.string}} }
