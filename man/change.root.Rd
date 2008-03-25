\name{change.root}
\alias{change.root}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Change tree root }
\description{
  The function changes the tree root.
}
\usage{
change.root(nodematrix, newroot)
}
\arguments{
  \item{nodematrix}{ the tree node matrix }
  \item{newroot}{ the node number of the new root }
}
\details{
  The function always returns an unrooted tree. Use the function \code{link{root.tree}} to root the unrooted tree if you need a rooted tree.
}
\value{
  \item{nodes }{the tree node matrix after changing the tree root }
  \item{rootnode }{the node number of the new root}
}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\seealso{\code{\link{root.tree}}, \code{\link{rootoftree}}}
\examples{
treestr<-"((((H:0.00402,C:0.00402):0.00304,G:0.00707):0.00929,O:0.01635):0.1,W:0.12);"
nodematrix<-read.tree.nodes(treestr)$nodes
change.root(nodematrix,6)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
