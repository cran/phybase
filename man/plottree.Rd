\name{plottree}
\alias{plottree}
\title{ Write a tree file }
\description{
  The function plots phylogenetic trees.
}
\usage{
plottree(tree)
}
\arguments{
  \item{tree}{a phylogenetic tree in newrick format}
}
\author{use the function "plot.phylo" in package APE to plot phylogenetic trees.}
\seealso{ \code{\link{write.subtree}}, \code{\link{read.tree.string}} }
\examples{
treestr<-"(((H:4.2,C:4.2):3.1,G:7.3):6.3,O:13.5);"
plottree(treestr)
}
\keyword{IO}
