\name{rank.nodes}
\alias{rank.nodes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Node ranks (internal use only) }
\description{
  The function returns the rank of each node in the tree. 
}
\usage{
rank.nodes(treenode, inode, ntaxa, start, rank)
}
\arguments{
  \item{treenode}{ tree node matrix }
  \item{inode}{ the tree root }
 \item{ntaxa}{the number of taxa in the tree}
 \item{start}{the maximum rank}
 \item{rank}{a dummy vector }
}
\value{
  The function returns a vector of ranks for the nodes in the tree.
}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\seealso{ \code{\link{mrca.2nodes}}, \code{\link{mrca.nodes}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
