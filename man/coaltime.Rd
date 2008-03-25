\name{coaltime}
\alias{coaltime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Coalescence time of two nodes }
\description{
  The function computes the coalescence time of two nodes.
}
\usage{
coaltime(inode, jnode, nodematrix, nspecies)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inode}{ the first node, it could be an internode. }
  \item{jnode}{ the second node, it could be an internode. }
  \item{nodematrix}{ the tree node matrix }
  \item{nspecies}{ the number of species }
}
\value{
  the function returns the coalescence time of \code{inode} and \code{jnode}.
}
\author{ Liang Liu }
\seealso{ \code{\link{popsize}} }
\examples{
treestr<-"((((H:0.00402,C:0.00402):0.00304,G:0.00706):0.00929,O:0.01635):0.1,W:0.11635);"
taxaname<-species.name(treestr)
nodematrix<-read.tree.nodes(treestr,name=taxaname)$nodes 
coaltime(1,2,nodematrix,5) #the coalescence time of taxa H (1) and C (2).
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
