\name{postdist.tree}
\alias{postdist.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate posterior probabilities of trees }
\description{
  The function summarize a set of trees by calculating the proportion of each tree in the tree set.
}
\usage{
postdist.tree(trees,name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trees}{ a vector of tree strings }
  \item{name}{ the species names }
}
\value{
  \item{trees }{a vector of tree}
  \item{prob}{the probability associated with each tree in the vector \code{tree}}
}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\seealso{ See Also as \code{\link{read.tree.nodes}}}
\examples{
library(phybase)
tree<-"(((H:0.005 , C:0.005 ) : 0.00025 #.01, G:0.00525):0.00025 #0.01 , O:0.0055) #.01;"
name<-species.name(tree)
nodematrix<-read.tree.nodes(tree,name)$nodes
rootnode<-7
seq<-rep(1,4)
nsim<-100
str<-rep(0,nsim)

for(i in 1:nsim){
	str[i]<-sim.coaltree.sp(rootnode,nodematrix,4,seq,name=name)$gt
}
postdist.tree(str,name)
}
\keyword{ programming}
