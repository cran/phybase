\name{root.tree}
\alias{root.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Root a tree }
\description{
  Root a tree. 
}

\usage{
root.tree(nodematrix,outgroup)
}
\arguments{
  \item{nodematrix}{ the tree node matrix }
  \item{outgroup}{the node used as outgroup}
}
\value{
  The function returns a rooted tree.}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\seealso{\code{\link{rootoftree}}, \code{\link{is.rootedtree}}}
\examples{
data(unrooted.tree)
nodematrix<-read.tree.nodes(unrooted.tree[1])$nodes
root.tree(nodematrix,23)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
