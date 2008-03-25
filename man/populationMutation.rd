\name{populationMutation}
\alias{populationMutation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Change branch lengths of a gene tree in the non-clocklike species tree model (internal use only)}
\description{
This function changes branch lengths of a gene tree with the mutation rates in the species tree.
}
\usage{
populationMutation(sptree, spnodedepth, genetree, genenodedepth, speciesmatrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sptree}{ the species tree }
  \item{spnodedepth}{depth of the species tree  }
  \item{genetree}{ a gene tree }
  \item{genenodedepth}{depth of the gene tree }
  \item{speciesmatrix}{tree node matrix of the species tree }
}
\value{
It returns a gene tree.
}

\author{Liang Liu }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
