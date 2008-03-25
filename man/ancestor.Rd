\name{ancestor}
\alias{ancestor}
\title{ Find the ancestral nodes of a node }
\description{
  The function returns the ancestral nodes of inode including inode itself. 
}
\usage{
ancestor(inode, nodematrix)
}
\arguments{
  \item{inode}{ the node number }
  \item{nodematrix}{ the tree node matrix. it must be a rooted tree. }
}
\value{
  The function returns a vector of ancestoral nodes of \code{inode} including \code{inode} itself.
}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\seealso{ \code{\link{mrca.2nodes}}, \code{\link{mrca.nodes}} }
\examples{
treestr<-"((((H:0.00402,C:0.00402):0.00304,G:0.00706):0.00929,O:0.01635):0.1,W:0.11635);"
nodematrix<-read.tree.nodes(treestr)$nodes
ancestor(6,nodematrix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
