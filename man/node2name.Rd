\name{node2name}
\alias{node2name}
\title{ Replace node numbers by species names in a tree string }
\description{
This function replaces node numbers in a tree string by species names. 
}
\usage{
node2name(treestr,name="")
}
\arguments{
  \item{treestr}{ a tree string }
  \item{name}{ species names }
}
\value{
  The function returns the tree string with the node numbers replaced by the species names.
}
\author{ Liang Liu }
\seealso{ \code{\link{subtree.length}}, \code{\link{name2node}} }
\keyword{univar}
\examples{
treestr<-"(((1:4.2,2:4.2):3.1,3:7.3):6.3,4:13.5);"
name<-c("H","C", "G", "O")
node2name(treestr,name)

}