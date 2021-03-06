\name{CGcoxph}
\alias{CGcoxph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Causal Gene - Continuous Survival Outcome modelled Semi-Parametrically}
\description{ Tests if a genetic marker is associated with a primary survival outcome through pathways other than that of an intermediate secondary phenotype with population-based data. Models the primary survival outcome using a Cox Proportional Hazards model. }
\usage{
CGcoxph(X, L, K, K.family, Y, Y.censured)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{matrix of marker genotypes. 1 row for each subject, 1 column for each marker}
  \item{L}{matrix of diagnostic criteria for K, 1 row for each subject, 1 column for each variable}
  \item{K}{secondary phenotype, 1 observation for each subject}
  \item{K.family}{family of GLM used to model K}
  \item{Y}{continuous primary outcome}
  \item{Y.censured}{censuring mechanism, 0=censured, 1=event occurred, as explained by Surv}
}

\value{
  The function outputs a single p-value for each marker (column of X), testing whether there is a direct effect of X on Y through pathways other than X to K to Y. The function also outputs an effect size estimate, which is the simple linear regression coefficient of the genetic marker on the primary outcome, after the primary outcome has been adjusted for the effect of the secondary phenotype. This estimate is appropriate for population-based data, assuming no population substructure. 
}
\references{Lipman P. et al: Inferring Genetic Causal Effects on Survival Data with Associated Endo-Phenotypes. Genetic Epidemiology 2011 (in press).}
\author{Peter Lipman}
\note{No missing values should appear for any variables input into the function. Uses Package Survival. All variables input into the function must be numeric.
}

\examples{
#example
#n<-1000				   #sample size
#X<-rbinom(n,2,.25)		           #marker, additive
#L<-rnorm(n,.1*X,1)			   #diagnostic criteria
#K<-rnorm(n,4+.25*L+.2*X,1)	           #secondary phenotype
#Y<-rweibull(n,abs(10*K),10)          	   #primary survival phenotype
#censuring<-rbinom(n,1,.95)		   #censuring mechanism
#CGcoxph(X,L,K,"gaussian",Y,censuring)
}

