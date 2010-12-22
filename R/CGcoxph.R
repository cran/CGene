CGcoxph <-
function(X,L,K,K.family,Y,Y.censured){
	
	p.values<-NULL
	effects<-NULL
	
	X<-as.matrix(X)
	L<-as.matrix(L)
	
	for(i in 1:ncol(X)){
		
		
		X1<-X[,i]-mean(X[,i],na.rm=T)
		
		glmK<-glm(K~X1+L,family=K.family)
		
		model1<-coxph(Surv(Y,Y.censured)~K+X1+L)
			
		#calculate partial (r.d) and full (eps) residuals
		data<-cbind(K,X1,L)
		
		xbase<-matrix(c(mean(K),mean(X1),mean(L)),1,3)
		
		data2<-t(data-matrix(xbase,dim(data)[1],dim(data)[2],byrow=T))
		
		sfit<-survfit(model1,newdata=xbase)
		chaz<--log(sfit$surv)
		
		cum.haz<-c()
		cum.haz[Y.censured==1]<-chaz*exp(model1$coef%*%data2[,Y.censured==1])
		cum.haz[Y.censured==0]<-mean(chaz*exp(model1$coef%*%data2[,Y.censured==1]))
		
		
		r.c<-exp((coef(model1)[2]*(K-mean(K))))*cum.haz                  #partial cox-snell residual
		r.m<-Y.censured-r.c				                                 #partial martingale residual
		r.d<-sign(r.m)*sqrt(-2*(r.m+Y.censured*log(Y.censured-r.m)))     #partial deviance residual
		
		eps<-residuals(model1,type="deviance")	#martingale residuals

					
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

