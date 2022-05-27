\name{MGGD}
\alias{MGGD}
\title{The Mixture Generalized Gamma Distribution}
\usage{
X <- rMGG(n,lambda,beta,alpha,p)
dMGG(X,lambda,beta,alpha,p)
pMGG(X,lambda,beta,alpha,p)
EM(X,n,lambda,beta,alpha,p)
}
\arguments{
\item{X}{vector of data.}
\item{n}{a number of observations.}
\item{lambda}{a value of the parameter lambda.}
\item{beta}{a value of the parameter beta.}
\item{alpha}{a value of the parameter alpha.}
\item{p}{a value of the parameter p.}
}
\value{
\item{dMGG}{gives the density of a mixture generalized gamma distributionts.}
\item{pMGG}{gives the distribution function of a mixture generalized gamma distribution.}
\item{rMGG}{generates random numbers of a mixture generalized gamma distribution.}
\item{EM}{gives parameter estimation function a mixture  generalized gamma distribution by using the expectation-maximization algorithm (EM).}
}
\description{
Density, distribution function, random generation, and parameter estimation function for the mixture generalized gamma distribution as Suksaengrakcharoen and Bodhisuwan  (2014) introduced.
}
\details{
Suksaengrakcharoen and Bodhisuwan (2014) introduced a mixture of generalized gamma distribution  which this distribution was a combination of the generalized gamma distribution and the length biased generalized gamma distribution. This distribution consisted of four parameters: (1) lambda  is the scale parameter, (2) beta is a shape parameter, (3) alpha is a shape parameter and (4) p is the weight parameter.
}
\examples{
#install.packages("devtools")
library(devtools)
#install_github("wikanda-phaphan/MGGD")
library(MGGD)

? MGGD
lambda <-0.7
beta <-1
alpha <-3
p <- 0.5
n <- 30
X <- rMGG(n,lambda,beta,alpha,p)
X
dMGG(X,lambda,beta,alpha,p)
pMGG(X,lambda,beta,alpha,p)
EM(X,n,lambda,beta,alpha,p)

}
\references{
Suksaengrakcharoen,  S. and Bodhisuwan, W. 2014. A new family of Generalized Gamma distribution and its application.  J. Math. Stat. 10, 2,  211-220.
}
