CGcont <-
function(X,L,K,K.family,Y,Y.family){
	
	p.values<-NULL
	effects<-NULL
	
	X<-as.matrix(X)
	L<-as.matrix(L)
	
	for(i in 1:ncol(X)){
		
		X1<-X[,i]-mean(X[,i],na.rm=T)
		
		glmK<-glm(K~X1+L,family=K.family)
		non.null<-as.integer(rownames(as.matrix(glmK$fitted.values))) #actual values used	
		glmY<-glm(Y~K+X1+L,family=Y.family)
		
		Y.tilde<-Y-mean(Y)-coef(glmY)[2]*(K-mean(K))
		T.i<-X1*Y.tilde	
		T<-sum(T.i)	
			
		eps<-residuals(glmY,type="deviance")
		T.tilde<-T.i[non.null]-mean(X1[non.null]*K[non.null])*(K[non.null]-fitted.values(glmK))*eps[non.null]/summary(glmK)$dispersion
		
		
		effects[i]<-summary(lm(Y.tilde~X1))$coefficients[2,1]
		
		sigma<-var(T.tilde)
		test.stat<-T^2/(length(non.null)*sigma)
		p.values[i]<-1-pchisq(test.stat,1)
		}
		
	final<-rbind(effects,p.values)
	rownames(final)<-c("effect size","p values")
	return(final)
	
	}

