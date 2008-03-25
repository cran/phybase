\name{offspring.species}
\alias{offspring.species}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Find the species nodes }
\description{
  The function returns the descendant species of \code{inode}. 
}

\usage{
offspring.species(inode, nodematrix, nspecies)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inode}{ the node. }
  \item{nodematrix}{ the tree node matrix }
  \item{nspecies}{ the number of species }
}
\value{
  This function returns the descendant species of \code{inode}, while the function \code{offspring.nodes} returns all the descendant nodes of \code{inode} including internal nodes in the tree.
}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\examples{
treestr<-"((((H:0.00402,C:0.00402):0.00304,G:0.00707):0.00929,O:0.01635):0.1,W:0.12);"
nodematrix<-read.tree.nodes(treestr)$nodes
offspring.species(7,nodematrix,5)

}
\seealso{ \code{\link{offspring.nodes}}}
\keyword{ programming }
