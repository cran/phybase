\name{write.subtree}
\alias{write.subtree}
\title{ Write a sub-tree into a string }
\description{
  write a tree or a sub-tree into a string in parenthetical format
}
\usage{
write.subtree(inode, nodes, nspecies,inodeindex)
}
\arguments{
  \item{inode}{ the root node of a sub-tree }
  \item{nodes}{ a tree node matrix }
  \item{nspecies}{ the number of species }
  \item{inodeindex}{ the root node of a sub-tree}
}
\details{
  If inode is the root of the tree, the function will write the whole tree into a string in parenthetical format. If inode is not the root node, the function will write the sub-tree into a string. The function works for both rooted trees and unrooted trees.
}
\value{
  The function returns a tree string in parenthetical format
}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\keyword{IO}
\seealso{ \code{\link{write.tree.string}}, \code{\link{read.tree.nodes}}}
\examples{
data(rooted.tree)
tree<-read.tree.nodes(rooted.tree[1])
tree$nodes
tree$names
write.subtree(7,tree$nodes,length(tree$names),7)
}