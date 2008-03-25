\name{genetree.vector}
\alias{genetree.vector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Construct gene tree vectors from multiple loci }
\description{
This function constructs gene tree vectors from gene trees across loci. The gene tree vectors can be used to construct maximum tree by the function \code{\link{maxtree}}.
}

\usage{
genetree.vector(filenames,outputfile)
}
\arguments{
  \item{filenames}{ the gene tree files }
  \item{outputfile}{ the output file  }
}
\value{
The function returns a matrix of gene trees. Each row represents a gene tree vector. The gene tree vector consists of trees from multiple gene tree files.}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\seealso{ \code{\link{maxtree}} }
\references{
Liu, L. and D.K. Pearl. Species trees from gene trees: reconstructing Bayesian posterior distributions of a species phylogeny using estimated gene tree distributions. Systematic Biology, 2007, 56:504-514. 

Edwards, S.V., L. Liu., and D.K. Pearl. High resolution species trees without concatenation. PNAS, 2007, 104:5936-5941.

}
\keyword{ programming }
