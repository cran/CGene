\name{CGdisc}
\alias{CGdisc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Causal Gene - Discrete Outcome}
\description{ Tests if a genetic marker is associated with a primary discrete outcome through pathways other than that of an intermediate secondary phenotype with population-based data.}
\usage{
CGdisc(X, L, K, K.family, Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{matrix of marker genotypes. 1 row for each subject, 1 column for each marker}
  \item{L}{matrix of diagnostic criteria for K, 1 row for each subject, 1 column for each variable}
  \item{K}{secondary phenotype, 1 observation for each subject}
  \item{K.family}{family of GLM used to model K}
  \item{Y}{discrete primary outcome}
}

\value{
  The function outputs a single p-value for each marker (column of X), testing whether there is a direct effect of X on Y through pathways other than X to K to Y. The function also outputs an effect size estimate, which is the simple linear regression coefficient of the genetic marker on the primary outcome, after the primary outcome has been adjusted for the effect of the secondary phenotype. This estimate is appropriate for population-based data, assuming no population substructure. 
}
\references{Vansteelandt S. et al: On the Adjustment for Covariates in Genetic Association Studies: A Novel, Simple Principle to Infer Direct Causal Effects. Genetic Epidemiology 2009; 33(5): 394-405.}
\author{Peter Lipman}
\note{No missing values should appear for any variables input into the function. All variables input into the function must be numeric.
}
\examples{
#example
#n<-1000				   #sample size
#X<-rbinom(n,2,.25)			   #marker, additive
#L<-rnorm(n,.1*X,1)			   #diagnostic criteria, X explains 1%
#K<-rnorm(n,.25*L+.2*X,1)	           #secondary phenotype
#Y<-rbinom(n,1,exp(.5*K)/(1+exp(.5*K)))    #primary outcome

#CGdisc(X,L,K,"gaussian",Y)

}

