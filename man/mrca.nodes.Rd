\name{mrca.nodes}
\alias{mrca.nodes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Find the most recent common ancestor of multiple nodes }
\description{
  The function can find the most recent common ancestor of multiple nodes specified in \code{nodevector}
}
\usage{
mrca.nodes(nodevector, nodematrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nodevector}{ a set of nodes }
  \item{nodematrix}{ the tree node matrix }
}
\value{
  The function returns the node number of the most recent common ancestor of the nodes in \code{nodevector}.
}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\seealso{\code{\link{mrca.2nodes}}, \code{\link{coaltime}}, \code{\link{popsize}}}
\examples{
treestr<-"((((H:0.00402,C:0.00402):0.00304,G:0.00707):0.00929,O:0.01635):0.1,W:0.12);"
nodematrix<-read.tree.nodes(treestr)$nodes
mrca.nodes(c(1,2,3),nodematrix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }

