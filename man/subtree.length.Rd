\name{subtree.length}
\alias{subtree.length}
\title{ Calculate total branch length of a tree }
\description{
  calculate the total branch length of a sub-tree under \code{inode}. 
}
\usage{
subtree.length(inode, nodes, nspecies)
}
\arguments{
  \item{inode}{ the root node of the sub-tree }
  \item{nodes}{ the tree node matrix }
  \item{nspecies}{ the number of species in the tree }
}
\details{
  The node matrix is the output of the function read.unrooted.nodes or read.rooted.nodes. The function can calculate the total branch length of a tree if inode is set to be the root node. If inode is not the root node, subtree.length calculates the total branch length of a sub-tree.
}
\value{
  The function returns the total branch length of a sub-tree.
}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\keyword{univar}
\seealso{ \code{\link{node.height}} }
\examples{
tree.string<-"(((H:4.2,C:4.2):3.1,G:7.3):6.3,O:13.5);"
nodes<-read.tree.nodes(tree.string)$nodes
subtree.length(6,nodes,4)

}
