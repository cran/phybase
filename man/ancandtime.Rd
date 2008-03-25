\name{ancandtime}
\alias{ancandtime}
\title{Get ancestors and their divergence times}
\description{
 This function returns the ancestors of a node and their divergence times. 
}
\usage{
ancandtime(inode, nodematrix, nspecies)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inode}{ a node in the tree. }
  \item{nodematrix}{ the tree matrix. }
  \item{nspecies}{ number of species (taxa) in the tree. }
}
\author{ Liang Liu }
\examples{
treestr<-"((((H:0.00402,C:0.00402):0.00304,G:0.00706):0.00929,O:0.01635):0.1,W:0.11635);"
nodematrix<-read.tree.nodes(treestr)$nodes
inode<-6
ancandtime(inode,nodematrix,nspecies=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
