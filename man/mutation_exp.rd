\name{mutation_exp}
\alias{mutation_exp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate mutation rates for populations in the species tree }
\description{
In the non-clock species tree model (Liu, et.al), the lineages (populations) in the species tree are allowed to have variable mutation rates. This function is used to simulate mutation rates for the non-clock species tree model. There are many other ways to simulate variable mutation rates across populations in the species tree.
}
\usage{
mutation_exp(sptree, root, inode, nspecies,alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sptree}{ the species tree matrix }
  \item{root}{ the root of the species tree }
  \item{inode}{ the root of the species tree }
  \item{nspecies}{ the number of species in the species tree }
 \item{alpha}{the parameter in the gamma distribution used to generate mutation rates.  }
}
\details{
mutation rates are generated from gamma (alpha, alpha/w) where w is the mutation rate of the parent population of the current node. Thus the mean of the mutation rate of the current node equals to the mutation rate of its parent population.
}
\value{
  The function returns a species tree matrix with mutation rates in the last column.
}
\author{ Liang Liu }
\examples{
sptree<-"((((H:0.00402#0.01,C:0.00402#0.01):0.00304#0.01,G:0.00707#0.01):0.00929#0.01,O:0.01635#0.01):0.1#0.01,W:0.12#0.01);"
nodematrix<-read.tree.nodes(sptree)$nodes
mutation_exp(nodematrix, root=9, inode=9, nspecies=5, alpha=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
