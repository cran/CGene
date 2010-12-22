\name{CGsurvreg}
\alias{CGsurvreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Causal Gene - Continuous Survival Outcome modelled Parametrically}
\description{ Tests if a genetic marker is associated with a primary survival outcome through pathways other than that of an intermediate secondary phenotype with population-based data. Models the primary survival outcome using a fully parametric survival model.}
\usage{
CGsurvreg(X, L, K, K.family, Y, Y.censured, Y.dist = "weibull", cum.hazard = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{matrix of marker genotypes. 1 row for each subject, 1 column for each marker}
  \item{L}{matrix of diagnostic criteria for K, 1 row for each subject, 1 column for each variable}
  \item{K}{secondary phenotype, 1 observation for each subject}
  \item{K.family}{family of GLM used to model K}
  \item{Y}{continuous primary outcome}
  \item{Y.censured}{censuring mechanism, 0=censured, 1=event occurred, as explained by Surv}
  \item{Y.dist}{family of survreg used to model Y, default is weibull}
  \item{cum.hazard}{user-provided function that takes output of a survreg object and computes the estimated cumulative hazard. must be provided if the Y.dist default (weibull) is not used. See example below}
}

\value{
  The function outputs a single p-value for each marker (column of X), testing whether there is a direct effect of X on Y through pathways other than X to K to Y.
}

\references{Lipman P. et al: Inferring Genetic Causal Effects on Survival Data with Associated Endo-Phenotypes. Genetic Epidemiology 2011 (in press).}
\author{Peter Lipman}
\note{No missing values should appear for any variables input into the function. Uses Package Survival. All variables input into the function must be numeric.
}

\examples{
#example of cum.hazard input:
#cum.hazard.weibull<-function(model,Y){
				#lambda<-1/exp(coef(model)[1])
				#gamma<-model$scale
				#cum.haz<-lambda*(Y^gamma)          
				#return(cum.haz)
				#}
#example:
#n<-1000					 #sample size
#X<-rbinom(n,2,.25)				 #marker, additive
#L<-rnorm(n,.1*X,1)			         #diagnostic criteria
#K<-rnorm(n,4+.25*L+.2*X,1)	                 #secondary phenotype
#Y<-rweibull(n,abs(10*K),10)          	         #primary survival phenotype
#censuring<-rbinom(n,1,.95)			 #censuring mechanism
#CGsurvreg(X,L,K,"gaussian",Y,censuring)
}

