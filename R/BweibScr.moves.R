




move.RP.BweibScr	<- function(survObj, ini, g){

	n <- survObj$n
	p <- survObj$p
	x <- survObj$x		
	
	gamma.ini		<- ini$gamma.ini	

	y1		<- survObj$y1
	y2		<- survObj$y2
	delta1	<- survObj$delta1
	delta2	<- survObj$delta2	
		
	accept	<- rep(0, p);	
				
	if(g == 1){
		xbeta		<- ini$xbeta1.ini		
		beta.ini	<- ini$beta1.ini
		alpha.ini	<- ini$alpha1.ini
		kappa.ini	<- ini$kappa1.ini	
		case.g		<- survObj$case1_
		y.alpha		<- y1^alpha.ini
		}	

	if(g == 2){
		xbeta		<- ini$xbeta2.ini		
		beta.ini	<- ini$beta2.ini
		alpha.ini	<- ini$alpha2.ini
		kappa.ini	<- ini$kappa2.ini	
		case.g		<- survObj$case01	
		y.alpha		<- y1^alpha.ini				
		}			

	if(g == 3){
		xbeta		<- ini$xbeta3.ini		
		beta.ini	<- ini$beta3.ini
		alpha.ini	<- ini$alpha3.ini
		kappa.ini	<- ini$kappa3.ini	
		case.g		<- survObj$case11	
		y.alpha		<- y2^alpha.ini - y1^alpha.ini
		}
	
#	for(j in 1:p){
		
		j <- sample(1:p, 1)
	
		exp.xbeta	<- exp(xbeta)
		loglh.ini	<- sum(xbeta[case.g]) - sum(kappa.ini * y.alpha * gamma.ini * exp.xbeta)
		
		D1 <- sum(x[case.g, j]) - sum(gamma.ini * kappa.ini * y.alpha * x[,j] * exp.xbeta)
		D2 <- -sum(gamma.ini * kappa.ini * y.alpha * x[,j]^2 * exp.xbeta)
		
		beta.prop.me <- beta.ini[j] - D1/D2
		beta.prop.var <- -2.4^2/D2
								
		beta.prop	<- beta.ini
		beta.prop[j] <- rnorm(1, mean= beta.prop.me, sd = sqrt(beta.prop.var))
		
		xbeta.prop	<- xbeta - x[,j] * beta.ini[j] + x[,j] * beta.prop[j]
		exp.xbeta.prop	<- exp(xbeta.prop)
		loglh.prop	<- sum(xbeta.prop[case.g]) - sum(kappa.ini * y.alpha * gamma.ini * exp.xbeta.prop)
				
		D1.prop <- sum(x[case.g, j]) - sum(gamma.ini * kappa.ini * y.alpha * x[,j] * exp.xbeta.prop)
		D2.prop <- -sum(gamma.ini * kappa.ini * y.alpha * x[,j]^2 * exp.xbeta.prop)
		
		beta.prop.me.ini <- beta.prop[j] - D1.prop/D2.prop
		beta.prop.var.ini <- -2.4^2/D2.prop		
												
		logprop.iniTOprop	<- dnorm(beta.prop[j], mean = beta.prop.me, sd = sqrt(beta.prop.var), log = TRUE)
		logprop.propTOini	<- dnorm(beta.ini[j], mean = beta.prop.me.ini, sd = sqrt(beta.prop.var.ini), log = TRUE)
		
		logR  <- loglh.prop - loglh.ini + logprop.propTOini - logprop.iniTOprop;

		u = log(runif(1)) < logR
	
		if(u == 1){
			beta.ini[j] <- beta.prop[j]
			xbeta	<- xbeta.prop
			}
		
		accept[j]	<- accept[j] + u;			
		
#		}	# end of loop for j
		

	list(beta.ini = beta.ini, accept = accept, xbeta = xbeta)

		

	}

	

	




move.WSC.BweibScr	<- function(survObj, ini, priorPara, g){

	n <- survObj$n
	p <- survObj$p
	x <- survObj$x	

	gamma.ini		<- ini$gamma.ini	

	y1		<- survObj$y1
	y2		<- survObj$y2
	delta1	<- survObj$delta1
	delta2	<- survObj$delta2	
	
	alpha.propVar	<- 0.5
	accept	<- 0	
				
	if(g == 1){
		xbeta		<- ini$xbeta1.ini		
		beta.ini	<- ini$beta1.ini
		alpha.ini	<- ini$alpha1.ini
		kappa.ini	<- ini$kappa1.ini	
		case.g		<- survObj$case1_
		y.alpha		<- y1^alpha.ini
		log.y 		<- log(y1)
		y.alpha.log.y		<- y.alpha * log.y
		y.alpha.log.y.sq	<- y.alpha * log.y^2
		a		<- priorPara$a1		
		b		<- priorPara$b1
		}	

	if(g == 2){
		xbeta		<- ini$xbeta2.ini		
		beta.ini	<- ini$beta2.ini
		alpha.ini	<- ini$alpha2.ini
		kappa.ini	<- ini$kappa2.ini	
		case.g		<- survObj$case01	
		y.alpha		<- y1^alpha.ini
		log.y 		<- log(y1)
		y.alpha.log.y		<- y.alpha * log.y
		y.alpha.log.y.sq	<- y.alpha * log.y^2
		a		<- priorPara$a2		
		b		<- priorPara$b2		
		}			

	if(g == 3){
		xbeta		<- ini$xbeta3.ini		
		beta.ini	<- ini$beta3.ini
		alpha.ini	<- ini$alpha3.ini
		kappa.ini	<- ini$kappa3.ini	
		case.g		<- survObj$case11	
		y.alpha		<- y2^alpha.ini - y1^alpha.ini
		log.y 		<- log(y2)
		y.alpha.log.y		<- y2^alpha.ini*log(y2) - y1^alpha.ini*log(y1)
		y.alpha.log.y.sq	<- y2^alpha.ini*(log(y2))^2 - y1^alpha.ini*(log(y1))^2
		a		<- priorPara$a3		
		b		<- priorPara$b3		
		}

		exp.xbeta	<- exp(xbeta)	
		loglh.ini <- (length(case.g) + a - 1) * log(alpha.ini) + alpha.ini*sum(log.y[case.g])						- sum(kappa.ini * y.alpha * gamma.ini * exp.xbeta) - b*alpha.ini
		
		D1 <- (length(case.g) + a - 1)/alpha.ini - sum(gamma.ini*kappa.ini*exp.xbeta*y.alpha.log.y)				-b + sum(log.y[case.g])
		
		D2 <- -(length(case.g) + a - 1)/alpha.ini^2 - sum(gamma.ini*kappa.ini*exp.xbeta*y.alpha.log.y.sq)
			
		alpha.prop.me <- alpha.ini - D1/D2
		alpha.prop.var <- -2.4^2/D2
		
		alpha.prop <- rtnorm(1, alpha.prop.me, sqrt(alpha.prop.var), lower = 0)
		
			
	if(g == 1 | g == 2){
		y.alpha.prop		<- y1^alpha.prop
		y.alpha.log.y.prop		<- y.alpha.prop * log.y
		y.alpha.log.y.sq.prop	<- y.alpha.prop * log.y^2
		}	
		

	if(g == 3){
		y.alpha.prop		<- y2^alpha.prop - y1^alpha.prop
		y.alpha.log.y.prop		<- y2^alpha.prop*log(y2) - y1^alpha.prop*log(y1)
		y.alpha.log.y.sq.prop	<- y2^alpha.prop*(log(y2))^2 - y1^alpha.prop*(log(y1))^2
		}
		
		
		loglh.prop <- (length(case.g) + a - 1) * log(alpha.prop) + alpha.prop*sum(log.y[case.g])						- sum(kappa.ini * y.alpha.prop * gamma.ini * exp.xbeta) - b*alpha.prop
				
		logprior.prop	<- dgamma(alpha.prop, shape = a, rate = b, log = TRUE)
		logprior.ini	<- dgamma(alpha.ini, shape = a, rate = b, log = TRUE)
		
		
		D1.prop <- (length(case.g) + a - 1)/alpha.prop - sum(gamma.ini*kappa.ini*exp.xbeta*y.alpha.log.y.prop)				-b + sum(log.y[case.g])
		
		D2.prop <- -(length(case.g) + a - 1)/alpha.prop^2 - sum(gamma.ini*kappa.ini*exp.xbeta*y.alpha.log.y.sq.prop)
	
		alpha.prop.me.prop <- alpha.prop - D1.prop/D2.prop
		alpha.prop.var.prop <- -2.4^2/D2.prop
							
		logprop.iniTOprop	<- dtnorm(alpha.prop, alpha.prop.me, sqrt(alpha.prop.var), lower = 0, log = TRUE)
		logprop.propTOini	<- dtnorm(alpha.ini, alpha.prop.me.prop, sqrt(alpha.prop.var.prop), lower = 0, log = TRUE)
				
		logR  <- loglh.prop - loglh.ini + logprior.prop - logprior.ini + logprop.propTOini - logprop.iniTOprop

		u = log(runif(1)) < logR
	
		if(u == 1){alpha.ini <- alpha.prop}
		
		accept	<- accept + u;			
		
	list(alpha.ini = alpha.ini, accept = accept)
		
	}
	
	
	
	







	
	
	
	
	

move.WSH.BweibScr	<- function(survObj, ini, priorPara, g){

	n <- survObj$n
	p <- survObj$p
	x <- survObj$x
	
	gamma.ini		<- ini$gamma.ini	

	y1		<- survObj$y1
	y2		<- survObj$y2
	delta1	<- survObj$delta1
	delta2	<- survObj$delta2	

				
	if(g == 1){
		xbeta		<- ini$xbeta1.ini		
		beta.ini	<- ini$beta1.ini
		alpha.ini	<- ini$alpha1.ini
		kappa.ini	<- ini$kappa1.ini	
		case.g		<- survObj$case1_
		y.alpha		<- y1^alpha.ini
		c		<- priorPara$c1		
		d		<- priorPara$d1					
			
		}	

	if(g == 2){
		xbeta		<- ini$xbeta2.ini		
		beta.ini	<- ini$beta2.ini
		alpha.ini	<- ini$alpha2.ini
		kappa.ini	<- ini$kappa2.ini	
		case.g		<- survObj$case01	
		y.alpha		<- y1^alpha.ini	
		c		<- priorPara$c2		
		d		<- priorPara$d2	
		}			

	if(g == 3){
		xbeta		<- ini$xbeta3.ini		
		beta.ini	<- ini$beta3.ini
		alpha.ini	<- ini$alpha3.ini
		kappa.ini	<- ini$kappa3.ini	
		case.g		<- survObj$case11	
		y.alpha		<- y2^alpha.ini - y1^alpha.ini
		c		<- priorPara$c3		
		d		<- priorPara$d3	
		}
		
	exp.xbeta	<- exp(xbeta)	
	
	kappa.sh	<- length(case.g) + c
	kappa.rate  <- sum(gamma.ini*y.alpha*exp.xbeta) + d
	
	kappa.ini <- rgamma(1, shape = kappa.sh, rate = kappa.rate)
		
	return(kappa.ini)
		
	}
	
	
	




		



move.FP.BweibScr		<- function(survObj, ini){

	n <- survObj$n
	p <- survObj$p
	x <- survObj$x	
		
	beta1.ini	<- ini$beta1.ini
	beta2.ini	<- ini$beta2.ini
	beta3.ini	<- ini$beta3.ini	
	alpha1.ini	<- ini$alpha1.ini		
	alpha2.ini	<- ini$alpha2.ini
	alpha3.ini	<- ini$alpha3.ini
	kappa1.ini	<- ini$kappa1.ini
	kappa2.ini	<- ini$kappa2.ini
	kappa3.ini	<- ini$kappa3.ini	
	theta.ini	<- ini$theta.ini
	xi.ini		<- ini$xi.ini
	
	y1		<- survObj$y1
	y2		<- survObj$y2
	delta1	<- survObj$delta1
	delta2	<- survObj$delta2	
	
	xbeta1		<- ini$xbeta1.ini
	exp.xbeta1	<- exp(xbeta1)	

	xbeta2		<- ini$xbeta2.ini
	exp.xbeta2	<- exp(xbeta2)	
	
	xbeta3		<- ini$xbeta3.ini
	exp.xbeta3	<- exp(xbeta3)	
	
	gamma.sh	<- delta1 + delta2 + 1/theta.ini
	
	gamma.rate	<- kappa1.ini*y1^alpha1.ini*exp.xbeta1 + kappa2.ini*y1^alpha2.ini*exp.xbeta2 					+ kappa3.ini*(y2^alpha3.ini - y1^alpha3.ini)*exp.xbeta3 + 1/theta.ini
	

	
	repeat{
		gamma.ini <- rgamma(n, shape = gamma.sh, rate = gamma.rate)
		if(all(gamma.ini != 0)){break}
		}



	return(gamma.ini);

	

	}		







move.DP.XiGam.BweibScr	<- function(survObj, ini, priorPara){

	n <- survObj$n
	p <- survObj$p
	
	psi		<- priorPara$psi
	omega	<- priorPara$omega	
				
	theta.ini	<- ini$theta.ini
	xi.ini		<- 1/theta.ini
	gamma.ini	<- ini$gamma.ini	

	accept <- 0
	
		D1 <- (psi-1)/xi.ini + n*log(xi.ini) + n - n*digamma(xi.ini) + 				sum(log(gamma.ini)) - (omega + sum(gamma.ini))
		D2 <- -(psi-1)/xi.ini^2 + n/xi.ini - n*trigamma(xi.ini)
			
		xi.prop.me  <-	xi.ini - D1 / D2;
		xi.prop.var <- -2.4^2 / D2;	
		
		if(xi.ini - D1 / D2	<= 0){
			xi.prop.me	<- xi.ini
			}
		
		xi.prop	<- rgamma(1, shape = xi.prop.me^2/xi.prop.var, rate = xi.prop.me/xi.prop.var)
		
		D1.prop <- (psi-1)/xi.prop + n*log(xi.prop) + n - n*digamma(xi.prop) + 				sum(log(gamma.ini)) - (omega + sum(gamma.ini))
		D2.prop <- -(psi-1)/xi.prop^2 + n/xi.prop - n*trigamma(xi.prop)
		
		xi.prop.me.ini  <-	xi.prop - D1.prop  / D2.prop ;
		xi.prop.var.ini <- -2.4^2 / D2.prop ;
		
		if(xi.ini - D1 / D2	<= 0){
			xi.prop.me.ini	<- xi.prop
		}		
								
		theta.prop		<- 1/xi.prop
	
		log.post.ini	<- (psi-1+n*xi.ini)*log(xi.ini) + sum(log(gamma.ini) - gamma.ini)*xi.ini 						- omega * xi.ini - n * lgamma(xi.ini);
		log.post.prop	<- (psi-1+n*xi.prop)*log(xi.prop) + sum(log(gamma.ini) - gamma.ini)*xi.prop 						- omega * xi.prop - n * lgamma(xi.prop);

		log.prop.ini	<- dgamma(xi.ini, shape = xi.prop.me.ini^2/xi.prop.var.ini, 								rate = xi.prop.me.ini/xi.prop.var.ini, log = TRUE)
		log.prop.prop	<- dgamma(xi.prop, shape = xi.prop.me^2/xi.prop.var, 								rate = xi.prop.me/xi.prop.var, log = TRUE)
	
	
	logR 	<- log.post.prop - log.post.ini + log.prop.ini - log.prop.prop
	
	u = log(runif(1)) < logR;
	
	if(u == 1){
		theta.ini	<- theta.prop
		xi.ini		<- xi.prop
		}

	accept	<- accept + u;
	
	list(theta.ini = theta.ini, xi.ini = xi.ini, accept = accept);
	
	}









	