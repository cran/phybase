\name{swap.nodes}
\alias{swap.nodes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Swap two nodes }
\description{
  The function swapps two subtrees.
}
\usage{
swap.nodes(inode, jnode, name, nodematrix)
}
\arguments{
  \item{inode}{ the root node of the first subtree }
  \item{jnode}{ the root node of the second subtree }
  \item{name}{ the species names }
  \item{nodematrix}{ the tree node matrix }
}
\value{
  \item{nodes }{the tree node matrix after swapping}
  \item{treestr }{the tree string after swapping}
}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\note{ 
The function is unable to swap two overlapped subtrees.
}
\seealso{ \code{\link{del.node}} }
\examples{
treestr<-"((((H:0.00402,C:0.00402):0.00304,G:0.00707):0.00929,O:0.01635):0.1,W:0.12);"
nodematrix<-read.tree.nodes(treestr)$nodes
spname<-read.tree.nodes(treestr)$names
swap.nodes(1,2,spname,nodematrix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
