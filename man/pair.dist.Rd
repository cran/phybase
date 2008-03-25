\name{pair.dist}
\alias{pair.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate all pairwise distances among taxa in the tree}
\description{
  The function computes all pairwise distances among taxa in the tree. 
}
\usage{
pair.dist(nodematrix, nspecies)
}
\arguments{
  \item{nodematrix}{ the tree node matrix }
  \item{nspecies}{the number of taxa in the tree}
}
\value{
  The function returns a distance matrix.}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\examples{
treestr<-"((((H:0.00402,C:0.00402):0.00304,G:0.00705):0.00929,O:0.01635):0.1,W:0.11635);"
nodematrix<-read.tree.nodes(treestr)$nodes
pair.dist(nodematrix,5)
}
\seealso{\code{\link{treedist}}, \code{\link{upgma}}, \code{\link{maxtree}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
