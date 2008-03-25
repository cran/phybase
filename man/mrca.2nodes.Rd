\name{mrca.2nodes}
\alias{mrca.2nodes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Find the most recent common ancestor of two nodes }
\description{
  The function can find the most recent common ancestor of two nodes \code{inode} and \code{jnode}
}
\usage{
mrca.2nodes(inode, jnode, nodematrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inode}{ the node \code{inode} }
  \item{jnode}{ the node \code{jnode} }
  \item{nodematrix}{ the tree node matrix }
}
\value{
  \item{anc}{the node number of the most recent common ancestor of \code{inode} and \code{jnode}.}
  \item{dist}{the distance between the two nodes.}
}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\seealso{\code{\link{mrca.nodes}}, \code{\link{coaltime}}, \code{\link{popsize}}}
\examples{
treestr<-"((((H:0.00402,C:0.00402):0.00304,G:0.00707):0.00929,O:0.01635):0.1,W:0.12);"
nodematrix<-read.tree.nodes(treestr)$nodes
mrca.2nodes(1,2,nodematrix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
