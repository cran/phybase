\name{rdirichlet}
\alias{rdirichlet}
\title{ Generate random numbers from the dirichlet distribution }
\description{
 This function can generate random numbers from a dirichlet distribution.
}
\usage{
rdirichlet(n,a)
}
\arguments{
  \item{n}{ the number of random numbers to be generated }
  \item{a}{ shape parameters of the dirichlet distribution }
}
\value{
  The function returns random numbers from a dirichlet distribution.
}
\author{ Code is taken from Greg's Miscellaneous Functions (gregmisc). His code was based on code posted by Ben Bolker to R-News on Fri Dec 15 2000.}
\keyword{IO}
\examples{
rdirichlet(1,c(3,3,3))
}