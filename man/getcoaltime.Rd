\name{getcoaltime}
\alias{getcoaltime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get coalescence times }
\description{
This function can get gene coalescence times in the species tree.
}
\usage{
getcoaltime(genetree, sptree, ntax, nspecies, species.structure)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genetree}{ a genetree matrix }
  \item{sptree}{ a species tree matrix }
  \item{ntax}{ number of taxa in the gene tree }
  \item{nspecies}{ number of species in the species tree }
  \item{species.structure}{ sequence-species relationship }
}
\value{
The function returns a two-column matrix, the first column is the ancestral node in the species tree, the second column is the gene coalescence time at the corresponding ancestral node in the species tree.
}
\author{ Liang Liu }
\examples{
genetree<-"(((A:1,B:1):3,C:4):2,D:6);"
sptree<-"(((A:0.5,B:0.5):1,C:1.5):1,D:2.5);"
name<-c("A","B","C","D")

genetree<-read.tree.nodes(genetree,name)$nodes
sptree<-read.tree.nodes(sptree,name)$nodes

ntax<-length(name)
nspecies<-length(name)
species.structure<-matrix(0,nrow=nspecies,ncol=ntax)
diag(species.structure)<-1

getcoaltime(genetree,sptree,ntax,nspecies,species.structure)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
