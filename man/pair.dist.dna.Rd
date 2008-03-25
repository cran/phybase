\name{pair.dist.dna}
\alias{pair.dist.dna}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate pairwise distances among DNA sequences }
\description{
  Calculate pairwise distances among DNA sequences. The DNA sequences are coded as 1:A, 2:G, 3:C, 4:T.
}
\usage{
pair.dist.dna(sequences, nst = 0)
}
\arguments{
  \item{sequences}{ DNA sequences }
  \item{nst}{ substitution model. 0:no model, 1:JC }
}
\details{
  If nst=0, the distance is equal to the proportion of sites having different nucleotides between two sequences.
}
\value{
  The function returns a distance matrix.
  }
\references{ 
Jukes, TH and Cantor, CR. 1969. Evolution of protein molecules. Pp. 21-123 in H. N. Munro, ed. Mammalian protein metabolism. Academic Press, New York. }

\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\seealso{ \code{\link{upgma}} }
\examples{
tree<-"(((H:0.00402#0.01,C:0.00402#0.01):0.00304#0.01,G:0.00707#0.01):0.00929#0.01,O:0.01635#0.01)#0.01;"
nodematrix<-read.tree.nodes(tree)$nodes
sequences<-sim.dna(nodematrix,10000,model=1)
pair.dist.dna(sequences,nst=1)
}
\keyword{ programming }
