\name{simSeqfromSp}
\alias{simSeqfromSp}
\title{simulate DNA sequences from a species tree }
\description{
  The function simulates sequences from a species tree.
}
\usage{
simSeqfromSp(sptree, spname, ntaxasp, ngene, theta=0, noclock=0, simsequence=1, murate="Dirichlet",alpha=5, seqlength=100, model=1, kappa=2, rate=c(1,1,1,1,1,1), frequency=c(1/4,1/4,1/4,1/4), outfile, format="phylip")}
\arguments{
  \item{sptree}{A species tree which must be a rooted tree. }
  \item{spname}{species names}
  \item{ntaxasp}{a vector of the number of individuals in each species}
\item{ngene}{number of genes}
\item{theta}{population size}
\item{noclock}{0: clocklike species tree 1: nonclocklike species tree}
\item{simsequence}{1: simulate sequences and gene trees, 0: simulate gene trees} 
\item{murate}{distribution of mutation rates}
\item{alpha}{the shape parameter of dirichlet distribution}
\item{seqlength}{the number of nucleotides along the sequences}
\item{model}{substitution model}
\item{kappa}{transition/transversion ratio}
\item{rate}{rates}
\item{frequency}{nucleotide frequency}
\item{outfile}{the full path of the output file}
\item{format}{either "phylip" or "nexus"}
}
\value{
  The function writes sequences into a file.
}
\examples{
#read the species tree from a data file
data(sptree)
outfile<-"out.txt"
spname <- paste("S",1:20,sep="")
outgroup <- "S20"
ntaxasp <- rep(2,length(spname))
ntaxasp[length(spname)]<-1
ngene<-2
seqlength<-100
simSeqfromSp(sptree,spname,ntaxasp,noclock=1,ngene=ngene,seqlength=seqlength,model=1,outfile=outfile)
simSeqfromSp(sptree,spname,ntaxasp,noclock=0,ngene=ngene,simsequence=0,seqlength=seqlength,model=1,outfile=outfile)
}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\keyword{IO}
\references{
Felsenstein, J. The Newick tree format. \url{http://evolution.genetics.washington.edu/phylip/newicktree.html} 
}
\seealso{ \code{\link{write.subtree}}, \code{\link{read.tree.string}} }
