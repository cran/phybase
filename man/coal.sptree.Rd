\name{coal.sptree}
\alias{coal.sptree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimating species trees using average coalescence times }
\description{
  For a given set of gene trees, the UPGMA tree is constructed from the distance matrix based on the average coalescence times among taxa.
}
\usage{
coal.sptree(trees, speciesname, nspecies, outgroup=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trees}{ a vector of trees in newick format }
  \item{speciesname}{ species names  }
  \item{nspecies}{ number of species }
  \item{outgroup}{ the node number of the species used to root the tree}
}
\details{
If the gene trees are not clocklike trees, they are first converted to clock trees using function \code{noclock2clock} and then construct a distance matrix in which the entries are twice the coalescence times among species. The distance matrix is used to build an UPGMA tree as the estimate of the species tree. This function is different from \code{steac.sptree} in that \code{steac.sptree} uses nucleotide distances to construct distance matrix.
}
\value{
  The function returns the tree node matrix and the estimate of the species tree.
}
\author{ Liang Liu }
\seealso{ See also to \code{\link{steac.sptree}} }
\examples{
data(rooted.tree)
genetrees<-rooted.tree
spname<-species.name(genetrees[1])
coal.sptree(genetrees,spname,nspecies=4,outgroup=4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
