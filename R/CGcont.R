`CGcont` <-
function(X,L,K,K.family,Y,Y.family){
	
	p.values<-NULL
	
	X<-as.matrix(X)
	L<-as.matrix(L)
	
	for(i in 1:ncol(X)){
		
		X1<-X[,i]-mean(X[,i],na.rm=T)
		
		glmK<-glm(K~X1+L,family=K.family)
		non.null<-as.integer(rownames(as.matrix(glmK$fitted.values))) #actual values used	
		glmY<-glm(Y~K+X1+L,family=Y.family)
			
		t.prime<-Y-mean(Y)-coef(glmY)[2]*(K-mean(K))
		T.i<-X1*t.prime	
		T<-sum(T.i)	
			
		eps<-residuals(glmY,type="deviance")
		T.tilda<-T.i[non.null]-mean(X1[non.null]*K[non.null])*(K[non.null]-fitted.values(glmK))*eps[non.null]/summary(glmK)$dispersion
		
		sigma<-var(T.tilda)
		test.stat<-T^2/(length(non.null)*sigma)
		p.values[i]<-1-pchisq(test.stat,1)
		}
	return(p.values)
	}

