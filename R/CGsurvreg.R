CGsurvreg <-
function(X,L,K,K.family,Y,Y.censured,Y.dist="weibull",cum.hazard=NULL){
	
	p.values<-NULL
	effects<-NULL

	X<-as.matrix(X)
	L<-as.matrix(L)
	
	for(i in 1:ncol(X)){
		
		X1<-X[,i]-mean(X[,i],na.rm=T)
		
		glmK<-glm(K~X1+L,family=K.family)
		
		model1<-survreg(Surv(Y,Y.censured)~K+X1+L,dist=Y.dist)
			
			#calculate partial and full residuals
			
		if (is.null(cum.hazard)){
			# calculate partial residual: this code only works for weibull
	
			lambda<-1/exp(coef(model1)[1])
			gamma<-model1$scale
			r.c<-exp((coef(model1)[2]*(K-mean(K))))*lambda*(Y^gamma)         #partial cox-snell residual
			r.m<-Y.censured-r.c				                                 #partial martingale residual
			r.d<-sign(r.m)*sqrt(-2*(r.m+Y.censured*log(Y.censured-r.m)))     #partial deviance residual
			
			# calculate full residual: this code only works for weibull
			r.c.full<-exp(cbind(X1,K,L)%*%as.matrix(coef(model1)[-1]))*lambda*(Y^gamma)
			r.m.full<-Y.censured-r.c.full	              
			eps<-sign(r.m.full)*sqrt(-2*(r.m.full+Y.censured*log(Y.censured-r.m.full)))  
				
			}
			
		if (!is.null(cum.hazard)){
			# user defined way of calculating cum.hazard
			
			cum.haz<-cum.hazard(model1,Y)
				
			#calculate partial residual
			r.c<-exp((coef(model1)[2]*(K-mean(K))))*cum.haz                  #partial cox-snell residual
			r.m<-Y.censured-r.c				                                 #partial martingale residual
			r.d<-sign(r.m)*sqrt(-2*(r.m+Y.censured*log(Y.censured-r.m)))     #partial deviance residual
				
			# calculate full residual
			r.c.full<-exp(cbind(X1,K,L)%*%as.matrix(coef(model1)[-1]))*cum.haz
			r.m.full<-Y.censured-r.c.full	#martingale residual
			eps<-sign(r.m.full)*sqrt(-2*(r.m.full+Y.censured*log(Y.censured-r.m.full)))  
			}
			
        #variance adjustment
		t.prime<-Y-r.d-mean(Y)  #adjusted phenotype
		T.i<-X1*t.prime			#multiply marker score by adjusted phenotype
		T<-sum(T.i)
		T.tilda<-T.i-mean(X1*K)*(K-fitted.values(glmK))*eps/summary(glmK)$dispersion
		sigma<-var(T.tilda)
		test.stat<-T^2/(length(K)*sigma)
		p.values[i]<-1-pchisq(test.stat,1)
		effects[i]<-summary(lm(t.prime~X1))$coefficients[2,1]
		
		}
	
	final<-rbind(effects,p.values)
	rownames(final)<-c("effect size","p values")
	return(final)
	
	}

