\name{is.rootedtree}
\alias{is.rootedtree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Is the tree rooted or not}
\description{
This function can test if the tree is rooted.
}
\usage{
is.rootedtree(tree)
}
\arguments{
  \item{tree}{ tree string or tree node matrix }
}
\value{
  The function returns TRUE if the tree is a rooted tree. Otherwise, it returns FALSE.
}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\seealso{\code{\link{is.clock}}}
\examples{
data(unrooted.tree)
nodematrix<-read.tree.nodes(unrooted.tree[1])$nodes
is.rootedtree(nodematrix)

data(rooted.tree)
is.rootedtree(rooted.tree[1])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
