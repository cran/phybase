\name{sortmat}
\alias{sortmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Sort a matrix }
\description{
  The function returns a sorted matrix
}
\usage{
sortmat(mat, columns)
}
\arguments{
  \item{mat}{ a matrix }
  \item{columns}{ the columns upon which the matrix is sorted }
}
\value{
  The function returns a sorted matrix.}
\seealso{ \code{\link{del.node}} }
\examples{
mat<-matrix(1:9,ncol=3)
sortmat(mat,1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
