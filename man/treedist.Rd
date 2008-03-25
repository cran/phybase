\name{treedist}
\alias{treedist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Distance between two trees }
\description{
  This function calculates the distance between two trees.
}

\usage{
treedist(tree1,tree2)
}
\arguments{
  \item{tree1}{ the first tree node matrix }
  \item{tree2}{ the second tree node matrix }
 }
\value{
  The function returns the distance of two trees.}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\seealso{\code{\link{pair.dist}}, \code{\link{partition.tree}}}
\examples{
treestr1<-"((((H:0.00402,C:0.00402):0.00304,G:0.00706):0.00929,O:0.01635):0.1,W:0.11635);"
treestr2<-"((((H:0.00402,G:0.00402):0.00304,C:0.00706):0.00929,O:0.01635):0.1,W:0.11635);"
name<-species.name(treestr1)
nodematrix1<-read.tree.nodes(treestr1,name)$nodes
nodematrix2<-read.tree.nodes(treestr2,name)$nodes
treedist(nodematrix1,nodematrix2)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
