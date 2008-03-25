\name{popsize}
\alias{popsize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Population size of the most recent common ancestor of two nodes }
\description{
 This function computes the population size of the most recent common ancestor of two nodes.
}
\usage{
popsize(inode, jnode, nodematrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inode}{ the first node, it could be an internode. }
  \item{jnode}{ the second node, it could be an internode. }
  \item{nodematrix}{ the tree node matrix }
}
\value{
  The function returns the population size of the most recent common ancestor of \code{inode} and \code{jnode}.
}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\seealso{ \code{\link{coaltime}} }
\examples{
treestr<-"((((H:0.00402,C:0.00402#0.035):0.00304,G:0.00706):0.00929,O:0.01635):0.1,W:0.11635);"
nodematrix<-read.tree.nodes(treestr)$nodes
popsize(1,2,nodematrix)
#[1] -9   ##this tree does not have values for population size.

popsize(1,1,nodematrix)
#[1] 0.035       ##the population size for the species C is 0.035
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
