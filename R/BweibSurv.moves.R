


move.RP.BweibSurv	<- function(survObj, ini){
	
	y		<- survObj$y
	delta	<- survObj$delta
	
	n <- survObj$n
	p <- survObj$p
	x <- survObj$x		
	
	accept	<- rep(0, p);	
				
		beta.ini	<- ini$beta.ini
		alpha.ini	<- ini$alpha.ini
		kappa.ini	<- ini$kappa.ini	
		xbeta		<- ini$xbeta.ini		
			
		y.alpha		<- y^alpha.ini
		case.g 		<- which(delta == 1)
	
#	for(j in 1:p){
		
		j <- sample(1:p, 1)
	
		exp.xbeta	<- exp(xbeta)
		loglh.ini	<- sum(xbeta[case.g]) - sum(kappa.ini * y.alpha * exp.xbeta)
		
		D1 <- sum(x[case.g, j]) - sum(kappa.ini * y.alpha * x[,j] * exp.xbeta)
		D2 <- -sum(kappa.ini * y.alpha * x[,j]^2 * exp.xbeta)
		
		beta.prop.me <- beta.ini[j] - D1/D2
		beta.prop.var <- -2.4^2/D2
								
		beta.prop	<- beta.ini
		beta.prop[j] <- rnorm(1, mean= beta.prop.me, sd = sqrt(beta.prop.var))
		
		xbeta.prop	<- xbeta - x[,j] * beta.ini[j] + x[,j] * beta.prop[j]
		exp.xbeta.prop	<- exp(xbeta.prop)
		loglh.prop	<- sum(xbeta.prop[case.g]) - sum(kappa.ini * y.alpha * exp.xbeta.prop)
				
		D1.prop <- sum(x[case.g, j]) - sum(kappa.ini * y.alpha * x[,j] * exp.xbeta.prop)
		D2.prop <- -sum(kappa.ini * y.alpha * x[,j]^2 * exp.xbeta.prop)
		
		beta.prop.me.ini <- beta.prop[j] - D1.prop/D2.prop
		beta.prop.var.ini <- -2.4^2/D2.prop		
												
		logprop.iniTOprop	<- dnorm(beta.prop[j], mean = beta.prop.me, sd = sqrt(beta.prop.var), log = TRUE)
		logprop.propTOini	<- dnorm(beta.ini[j], mean = beta.prop.me.ini, sd = sqrt(beta.prop.var.ini), log = TRUE)
		
		logR  <- loglh.prop - loglh.ini + logprop.propTOini - logprop.iniTOprop;

		u = log(runif(1)) < logR
	
		if(u == 1){
			beta.ini[j] <- beta.prop[j]
			xbeta <- xbeta.prop
			}
		
		accept[j]	<- accept[j] + u;			
		
#		}	# end of loop for j
		
	list(beta.ini = beta.ini, accept = accept, xbeta = xbeta)
		
	}
	
	



move.WSC.BweibSurv	<- function(survObj, ini, priorPara){


	y		<- survObj$y
	delta	<- survObj$delta
	
	n <- survObj$n
	p <- survObj$p	
	
	alpha.propVar	<- 0.5
	accept	<- 0	
				

		beta.ini	<- ini$beta.ini
		xbeta		<- ini$xbeta.ini		
		alpha.ini	<- ini$alpha.ini
		kappa.ini	<- ini$kappa.ini	
		case.g		<- which(delta == 1)
		y.alpha		<- y^alpha.ini
		log.y 		<- log(y)
		y.alpha.log.y		<- y.alpha * log.y
		y.alpha.log.y.sq	<- y.alpha * log.y^2
		a		<- priorPara$a	
		b		<- priorPara$b


		exp.xbeta	<- exp(xbeta)	
		loglh.ini <- (length(case.g) + a - 1) * log(alpha.ini) + alpha.ini*sum(log.y[case.g])						- sum(kappa.ini * y.alpha * exp.xbeta) - b*alpha.ini
		
		D1 <- (length(case.g) + a - 1)/alpha.ini - sum(kappa.ini*exp.xbeta*y.alpha.log.y)				-b + sum(log.y[case.g])
		
		D2 <- -(length(case.g) + a - 1)/alpha.ini^2 -sum(kappa.ini*exp.xbeta*y.alpha.log.y.sq)
			
		alpha.prop.me <- alpha.ini - D1/D2
		alpha.prop.var <- -2.4^2/D2
		
		alpha.prop <- rtnorm(1, alpha.prop.me, sqrt(alpha.prop.var), lower = 0)
		
		y.alpha.prop			<- y^alpha.prop
		y.alpha.log.y.prop		<- y.alpha.prop * log.y
		y.alpha.log.y.sq.prop	<- y.alpha.prop * log.y^2

				
		loglh.prop <- (length(case.g) + a - 1) * log(alpha.prop) + alpha.prop*sum(log.y[case.g])						- sum(kappa.ini * y.alpha.prop * exp.xbeta) - b*alpha.prop
				
		logprior.prop	<- dgamma(alpha.prop, shape = a, rate = b, log = TRUE)
		logprior.ini	<- dgamma(alpha.ini, shape = a, rate = b, log = TRUE)
				
		D1.prop <- (length(case.g) + a - 1)/alpha.prop - sum(kappa.ini*exp.xbeta*y.alpha.log.y.prop)				-b + sum(log.y[case.g])
		
		D2.prop <- -(length(case.g) + a - 1)/alpha.prop^2 - sum(kappa.ini*exp.xbeta*y.alpha.log.y.sq.prop)
	
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
	
	
	
	


	
	
	

move.WSH.BweibSurv	<- function(survObj, ini, priorPara){
	

	y		<- survObj$y
	delta	<- survObj$delta
	
	n <- survObj$n
	p <- survObj$p

	xbeta		<- ini$xbeta.ini
	beta.ini	<- ini$beta.ini
	alpha.ini	<- ini$alpha.ini
	kappa.ini	<- ini$kappa.ini	
	case.g		<- which(delta == 1)
	y.alpha		<- y^alpha.ini
	c		<- priorPara$c	
	d		<- priorPara$d					
			

	exp.xbeta	<- exp(xbeta)	
	
	kappa.sh	<- length(case.g) + c
	kappa.rate  <- sum(y.alpha*exp.xbeta) + d
	
	kappa.ini <- rgamma(1, shape = kappa.sh, rate = kappa.rate)
		
	return(kappa.ini)
		
	}
	
	
	


		




	