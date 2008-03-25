\name{unroottree}
\alias{unroottree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Unroot a tree }
\description{
  unroot a tree. 
}

\usage{
unroottree(nodematrix)
}
\arguments{
  \item{nodematrix}{ the tree node matrix }
}
\value{
  The function returns an unrooted tree.}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\examples{
treestr<-"((((H:0.00402,C:0.00402):0.00304,G:0.00707):0.00929,O:0.01635):0.1,W:0.12);"
nodematrix<-read.tree.nodes(treestr)$nodes
spname<-read.tree.nodes(treestr)$names
unroottree(nodematrix)
}
\seealso{\code{\link{rootoftree}}, \code{\link{root.tree}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
