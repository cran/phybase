\name{rootoftree}
\alias{rootoftree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Root of a tree }
\description{
This function can be used to find the root of a tree. 
}

\usage{
rootoftree(nodematrix)
}
\arguments{
  \item{nodematrix}{ the tree node matrix }
}
\value{
  The function returns the root of the tree.}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\seealso{\code{\link{rootoftree}}, \code{\link{root.tree}}}
\examples{
treestr<-"((((H:0.00402,C:0.00402):0.00304,G:0.00707):0.00929,O:0.01635):0.1,W:0.12);"
nodematrix<-read.tree.nodes(treestr)$nodes
spname<-read.tree.nodes(treestr)$names
rootoftree(nodematrix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }

