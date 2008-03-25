\name{consense}
\alias{consense}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Consensus tree}
\description{
  The function returns a consensus tree from multiple gene trees. 
}
\usage{
consense(treestr, name,type="freq")
}
\arguments{
  \item{treestr}{ a vector of tree strings }
  \item{name}{the species names}
  \item{type}{if type="freq", the frequency of each clade in the consensus tree is presented at the node of the clade. if type="prop", the proportion of each clade is presented at the node of the clade"}
}
\value{
  The function returns the consensus tree and species names.}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\seealso{\code{\link{maxtree}}, \code{\link{partition.tree}}}
\examples{
treestr<-c("((((H:0.00402,C:0.00402):0.00304,G:0.00707):0.00929,O:0.01635):0.1,W:0.12);","((((H:0.00402,G:0.00402):0.00304,C:0.00707):0.00929,O:0.01635):0.1,W:0.12);","((((H:0.00402,C:0.00402):0.00304,G:0.00707):0.00929,O:0.01635):0.1,W:0.12);")
name<-species.name(treestr[1])
consense(treestr,name)


###unrooted trees
data(unrooted.tree)
name<-paste("S",1:29,sep="")
consense(unrooted.tree,name)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }

