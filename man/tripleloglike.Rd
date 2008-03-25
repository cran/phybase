\name{tripleloglike}
\alias{tripleloglike}
\title{ Loglikehood of Triples  }
\description{
  The function calculates the loglikelihood for DNA sequences (snip data)
}
\usage{
tripleloglike(sptree,spname,dna)
}
\arguments{
  \item{sptree}{ species tree }
 \item{spname}{ species names }
\item{dna}{  dna sequences }
}
\details{This function is used to calculate the loglikelihood of triples.
 }
\value{
  The function returns the loglikehood of triples.
}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\keyword{IO}
\seealso{ \code{\link{write.subtree}}, \code{\link{read.tree.string}} }
