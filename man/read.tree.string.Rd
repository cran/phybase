\name{read.tree.string}
\alias{read.tree.string}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read tree strings from a tree file }
\description{
  This function reads tree strings in the parenthetical format from a tree file. The output of the function is a vector of tree strings that can be converted to a matrix of nodes by the function \code{\link{read.tree.nodes}}.
}
\usage{
read.tree.string(file = "", format="nexus")
}
\arguments{
  \item{file}{ the tree file that contains trees in the parenthetical format.}
  \item{format}{phylip or nexus} 
}
\details{
The function can read NEXUS and PHYLIP tree files. It works for other types of tree files as long as the trees in the tree files are parenthetical trees. This function combining with \code{\link{write.tree.string}} can change the tree file format. 
}
\value{
  \item{tree}{a vector of tree strings.}
  \item{names}{species names.}
  \item{root}{TRUE for rooted trees, FALSE for unrooted trees}
}
\author{ Liang Liu \email{lliu@oeb.harvard.edu} }
\keyword{IO}
\seealso{ \code{\link{write.tree.string}}, \code{\link{read.tree.nodes}}}
\examples{
##read rooted trees in PHYLIP format
cat("(((H:4.2,C:4.2):3.1,G:7.3):6.3,O:13.5);",file = "phylip.tre", sep = "\n")
tree.string<-read.tree.string("phylip.tre")
tree.string

##read unrooted trees in NEXUS format
cat("#NEXUS
[ID: 4045516090]
begin trees;
   translate
       1 WW_7A03_1,
       2 WW_7H06_2,
       3 WW_7H05_1,
       4 WW_N03__5,
       5 WW_Snnr_1,
       6 WW_7P10__1,
       7 WW_7A05_1,
       8 WW_B03__1,
       9 WW_B04_1,
      10 WW_D07_9,
      11 WW_7K01_1,
      12 WW_7K04_1,
      13 WW_7N13_1,
      14 WW_M02_1,
      15 WW_N04_1,
      16 WW_UK6_1,
      17 WW_7A04_1,
      18 Pfuscatus_PF2_4,
      19 Pfuscatus_PF1_1,
      20 Pfuscatus_PF3_2,
      21 PCabietinus_331_1,
      22 PCabietinus_333_6,
      23 PCabietinus_336_1,
      24 PCcollybita_GB_1,
      25 PCtristis_GB_1,
      26 PCbrehmii_GB_1,
      27 Psibilatrix_GB_1,
      28 Pbonelli_GB_1,
      29 PTviridanus_1;
   tree rep.1 = (((((((((((((((16:0.100000,2:0.100000):0.100000,20:0.100000):0.100000,21:0.100000):0.100000,8:0.100000):0.100000,((19:0.100000,9:0.100000):0.100000,13:0.100000):0.100000):0.100000,23:0.100000):0.100000,27:0.100000):0.100000,5:0.100000):0.100000,26:0.100000):0.100000,((22:0.100000,28:0.100000):0.100000,11:0.100000):0.100000):0.100000,(24:0.100000,10:0.100000):0.100000):0.100000,6:0.100000):0.100000,(18:0.100000,((15:0.100000,14:0.100000):0.100000,(25:0.100000,12:0.100000):0.100000):0.100000):0.100000):0.100000,17:0.100000):0.100000,(7:0.100000,(3:0.100000,(1:0.100000,4:0.100000):0.100000):0.100000):0.100000,29:0.100000);
   tree rep.100 = (((((((((((((((16:0.100000,2:0.100000):0.100000,20:0.100000):0.100000,21:0.100000):0.100000,8:0.100000):0.100000,((19:0.100000,9:0.100000):0.100000,13:0.100000):0.100000):0.100000,23:0.100000):0.100000,27:0.100000):0.100000,5:0.100000):0.100000,26:0.100000):0.100000,((22:0.100000,28:0.100000):0.100000,11:0.100000):0.100000):0.100000,(24:0.100000,10:0.100000):0.100000):0.100000,6:0.100000):0.100000,(18:0.100000,((15:0.100000,14:0.100000):0.100000,(25:0.100000,12:0.100000):0.100000):0.100000):0.100000):0.100000,17:0.100000):0.100000,(7:0.100000,(3:0.100000,(1:0.100000,4:0.100000):0.100000):0.100000):0.100000,29:0.100000);
   tree rep.200 = (((((((((((((((16:0.100000,2:0.100000):0.100000,20:0.100000):0.100000,21:0.100000):0.100000,8:0.100000):0.100000,((19:0.100000,9:0.100000):0.100000,13:0.100000):0.100000):0.100000,23:0.100000):0.100000,27:0.100000):0.100000,5:0.100000):0.100000,26:0.100000):0.100000,((22:0.100000,28:0.100000):0.100000,11:0.100000):0.100000):0.100000,(24:0.100000,10:0.100000):0.100000):0.100000,6:0.100000):0.100000,(18:0.100000,((15:0.100000,14:0.100000):0.100000,(25:0.100000,12:0.100000):0.100000):0.100000):0.100000):0.100000,17:0.100000):0.100000,(7:0.100000,(3:0.100000,(1:0.100000,4:0.100000):0.100000):0.100000):0.100000,29:0.100000);
end;",file="tree.nexus")
tree.string<-read.tree.string("tree.nexus")
tree.string
}