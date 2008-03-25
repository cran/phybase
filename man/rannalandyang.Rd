\name{rannalandyang}
\alias{rannalandyang}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Rannala and Yang's formula }
\description{
  This function calculates the likelihood of a vector of gene trees given the species tree using the Rannala and Yang's formula
}
\usage{
rannalandyang(gtree, stree, taxaname,spname,species.structure)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gtree}{ a collection of gene trees }
  \item{stree}{ a species tree in newick format }
   \item{taxaname}{ the names of taxa }
  \item{spname}{ the names of species }
  \item{species.structure}{ define which sequence belong to which species }
}
\value{
The function returns the log likelihood score.
}
\references{ Rannala, B. and Z. Yang. 2003. Bayes estimation of species divergence times and ancestral population sizes using DNA sequences from multiple loci. Genetics 164: 1645-1656. }
\author{ Liang Liu }
\examples{
gtree<-"(((A:1,B:1):3,C:4):2,D:6);"
stree<-"(((A:0.5,B:0.5):1#0.1,C:1.5):1#0.1,D:2.5)#0.1;"
taxaname<-c("A","B","C","D")
spname<-taxaname
ntax<-length(taxaname)
nspecies<-length(spname)
species.structure<-matrix(0,nrow=nspecies,ncol=ntax)
diag(species.structure)<-1
rannalandyang(gtree,stree,taxaname,spname,species.structure)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
