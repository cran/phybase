\name{offspring.nodes}
\alias{offspring.nodes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Find the offspring nodes }
\description{
  The function returns the offspring nodes of \code{inode}. 
}

\usage{
offspring.nodes(inode, nodematrix, nspecies)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inode}{ the node of which the the offspring nodes will be found by the function. }
  \item{nodematrix}{ the tree node matrix. }
  \item{nspecies}{ the number of species. }
}
\value{
  The function returns the offspring nodes of \code{inode}.
}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\examples{
treestr<-"((((H:0.00402,C:0.00402):0.00304,G:0.00707):0.00929,O:0.01635):0.1,W:0.12);"
nodematrix<-read.tree.nodes(treestr)$nodes
offspring.nodes(7,nodematrix,5)
}
\seealso{ \code{\link{offspring.species}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
