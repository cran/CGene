CGdisc <-
function(X,L,K,K.family,Y){
	
	X<-as.matrix(X)
	L<-as.matrix(L)
	
	p.values<-NULL
	effects<-NULL
	
	for(i in 1:ncol(X)){
		
		X1<-X[,i]
		
		glmY<-glm(Y~K+X1+L,family=poisson)
		
		mu.i<-fitted.values(glmY)
			
		Y.tilde<-Y*exp(-coef(glmY)[2]*K)-mu.i
		eps<-residuals(glmY,type="deviance")
		
		glmK<-suppressWarnings(glm(K~X1+L,family=K.family,weights=mu.i))
		fit.K<-fitted.values(glmK)
		
		effects[i]<-summary(lm(Y.tilde~X1))$coefficients[2,1]

		lambda<-mean((X1-mean(X1))*Y*exp(-coef(glmY)[2]*K)*K)/(fit.K*summary(glmK)$dispersion)
		
		T.tilde<-(X1-mean(X1))*Y.tilde-lambda*(K-fit.K)*eps
		T<-sum(T.tilde)
		sigma<-var(T.tilde)
		test.stat<-T^2/(length(eps)*sigma)
		p.values[i]<-1-pchisq(test.stat,1)
		}
		
	final<-rbind(effects,p.values)
	rownames(final)<-c("effect size","p values")
	return(final)
	
	}

