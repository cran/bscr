


move.RP.BpeSurv	<- function(survObj, ini){
	
	n <- survObj$n
	p <- survObj$p
	x <- survObj$x		

	beta.ini	<- ini$beta.ini
	lambda		<- ini$lambda
	s			<- ini$s
	J			<- ini$J
	ind.r		<- ini$ind.r
	ind.d		<- ini$ind.d
	d			<- ini$d
	Delta		<- ini$Delta
	xbeta		<- ini$xbeta.ini

	accept <- rep(0, p)

#	updatebeta = sample(c(1:p), mh$numUpdateBeta, replace = FALSE)
	updatebeta = sample(c(1:p), 1, replace = FALSE)
		
	
	for(j in updatebeta){
		
		xbeta.mat					<- matrix(rep(xbeta, J+1), n, J+1)
		xbeta.dum					<- xbeta
		xbeta.dum[xbeta.dum > 700] 	<- 700	
		exp.xbeta					<- exp(xbeta.dum)
		exp.xbeta.mat				<- matrix(rep(exp.xbeta, J+1), n, J+1)
		xbeta.mat					<- matrix(rep(xbeta, J+1), n, J+1)
	
		x.j.mat <- matrix(rep(x[,j],J+1),J+1, nrow = n);
	
		lam.dum	<- lambda;
		lam.dum[lam.dum > 700] <- 700;
	
		D1.1st	<- -exp(lam.dum) * colSums(exp.xbeta.mat * x.j.mat * Delta * ind.r);
		D1.2nd	<- colSums(x.j.mat * ind.d);
		D1		<- sum(D1.1st + D1.2nd)	
			
		D2		<- sum(-exp(lam.dum) * colSums(exp.xbeta.mat * x.j.mat^2 * Delta * ind.r)) 
		
		beta.prop.me  <- beta.ini[j] - D1 / D2
		beta.prop.var <- - 2.4^2 / D2
		
		beta.prop <- beta.ini	
			
#		if(mh$betaUB){
#			if(abs(D1/D2) > mh$betaPropUBval){	
#				beta.prop.me  	<- beta.ini[j] - sign(D1/D2) * mh$betaPropUBval;
#				beta.prop.var	<- mh$betaPropVar
#				}				
#			if(D1 == D2){
#				beta.prop.me  <- beta.ini[j];
#				beta.prop.var	<- mh$betaPropVar
#				}
#			}					
		
		beta.prop[j]		<- rnorm(1, mean = beta.prop.me, sd = sqrt(beta.prop.var))

		
		# Calculating acceptance probability	
	
		xbeta.prop					<- xbeta - x[,j] * beta.ini[j] + x[,j] * beta.prop[j]
		xbeta.mat.prop				<- matrix(rep(xbeta.prop, J+1), n, J+1)
		xbeta.dum.prop				<- xbeta.prop
		xbeta.dum.prop[xbeta.dum.prop > 700] 	<- 700
		exp.xbeta.prop				<- exp(xbeta.dum.prop)
		exp.xbeta.mat.prop			<- matrix(rep(exp.xbeta.prop, J+1), n, J+1)
		xbeta.mat.prop				<- matrix(rep(xbeta.prop, J+1), n, J+1)
		
		D1.1st.prop	<- -exp(lam.dum) * colSums(exp.xbeta.mat.prop * x.j.mat * Delta * ind.r);
		D1.2nd.prop	<- colSums(x.j.mat * ind.d);
		D1.prop		<- sum(D1.1st.prop + D1.2nd.prop)
			
		D2.prop		<- sum(-exp(lam.dum) * colSums(exp.xbeta.mat.prop * x.j.mat^2 * Delta * ind.r)) 
		
		beta.prop.me.ini  <- beta.prop[j] - D1.prop / D2.prop;
		beta.prop.var.ini <- -2.4^2 / D2.prop;	
		
		log.lh.ini.1st	<- -exp(lam.dum) * colSums(exp.xbeta.mat * Delta * ind.r)
		log.lh.ini.2nd	<- colSums(xbeta.mat * ind.d)
		log.lh.ini		<- sum(log.lh.ini.1st + log.lh.ini.2nd)
				
		log.lh.prop.1st	<- -exp(lam.dum) * colSums(exp.xbeta.mat.prop * Delta * ind.r)
		log.lh.prop.2nd	<- colSums(xbeta.mat.prop * ind.d)
		log.lh.prop		<- sum(log.lh.prop.1st + log.lh.prop.2nd)
		
		log.prop.prop  <- dnorm(beta.prop[j] , mean = beta.prop.me, 							sd = sqrt(beta.prop.var), log = TRUE)
		log.prop.ini   <- dnorm(beta.ini[j]  , mean = beta.prop.me.ini, 							sd = sqrt(beta.prop.var.ini), log = TRUE)

		logR  <- log.lh.prop - log.lh.ini + log.prop.ini - log.prop.prop

		u = log(runif(1)) < logR
	
		if(u == 1){
			beta.ini[j] <- beta.prop[j]
			xbeta	<- xbeta.prop
			}

		accept[j] <- accept[j] + u	
		
		}  # end of for loop
	
	
	list(beta.ini = beta.ini, accept = accept, xbeta.ini = xbeta);
	
	} # end of MoveRP







# using joint prior for the prior and conditional for jumping kernel

			
move.BH.BpeSurv	<- function(survObj, ini){
	
	n <- survObj$n
	p <- survObj$p		

	beta.ini	<- ini$beta.ini
	lambda		<- ini$lambda
	s			<- ini$s
	J			<- ini$J
	ind.r		<- ini$ind.r
	d			<- ini$d
	Delta		<- ini$Delta
	mu.lam		<- ini$mu.lam
	sig.sq.lam	<- ini$sig.sq.lam
	Sigma.lam	<- ini$Sigma.lam
	inv.Sigma.lam	<- ini$inv.Sigma.lam
	W			<- ini$W
	Q			<- ini$Q
	xbeta		<- ini$xbeta.ini	
	 
	accept	<- rep(0, J+1);
	
	xbeta.mat					<- matrix(rep(xbeta, J+1), n, J+1)
	xbeta.dum					<- xbeta
	xbeta.dum[xbeta.dum > 700] 	<- 700	
	exp.xbeta					<- exp(xbeta.dum)
							
#	updatelambda = sample(c(1:(J+1)), mh$numUpdateLambda, replace = FALSE)
	updatelambda = sample(c(1:(J+1)), 1, replace = FALSE)
		
	
	for(j in updatelambda){
		
		if((J+1) > 1){
			if(j == 1){		
				nu.lam	<- mu.lam + W[1, 2] * (lambda[2] - mu.lam)
				}
		
			if(j == J+1){
				nu.lam	<- mu.lam + W[J+1, J] * (lambda[J] - mu.lam);
				}
		
			if(j != 1 & j != J+1){
				nu.lam	<- mu.lam + W[j, (j-1)] * (lambda[j-1] - mu.lam) 							+ W[j, (j+1)] * (lambda[j+1] - mu.lam);
				}
			}
		
		if(J+1 == 1){
			nu.lam	<- mu.lam 
			}
				
		lam.dum	<- lambda[j];
	
		if(lam.dum > 700){
			lam.dum <- 700;
			}	
			
		if(J+1 >= 2){
			D.1st	<- -exp(lam.dum) * sum(exp.xbeta * Delta[,j] * ind.r[,j])
			D1		<- D.1st + d[j] - 1/(sig.sq.lam * Q[j,j])*(lambda[j] - nu.lam);
			D2		<- D.1st - 1/(sig.sq.lam * Q[j,j]);
			}
		
		if(J+1 == 1){
			D.1st	<- -exp(lam.dum) * sum(exp.xbeta * Delta[,j] * ind.r[,j])
			D1		<- D.1st + d[j] - 1/(sig.sq.lam * Q)*(lambda[j] - nu.lam);
			D2		<- D.1st - 1/(sig.sq.lam * Q);
			}	

		# proposal value
							
		lambda.prop.me  <- lambda[j] - D1 / D2
		lambda.prop.var <- - 2.4^2 / D2
		
#		if(mh$lambdaUB){
#			if(abs(D1/D2) > mh$lambdaUBval){	
#				lambda.prop.me  <- lambda[j] - sign(D1/D2) * mh$lambdaUBval;
#				lambda.prop.var	<- mh$lambdaPropVar;
#				}			
#			}
	
		lambda.prop 	<- lambda
		lambda.prop[j] 	<- rnorm(1, mean = lambda.prop.me, sd = sqrt(lambda.prop.var))
		
		if((J+1) > 1){
			if(j == 1){		
				nu.lam.prop	<- mu.lam + W[1, 2] * (lambda.prop[2] - mu.lam)
				}
		
			if(j == J+1){
				nu.lam.prop	<- mu.lam + W[J+1, J] * (lambda.prop[J] - mu.lam);
				}
		
			if(j != 1 & j != J+1){
				nu.lam.prop	<- mu.lam + W[j, (j-1)] * (lambda.prop[j-1] - mu.lam) 								+ W[j, (j+1)] * (lambda.prop[j+1] - mu.lam);
				}
			}
		
		if(J+1 == 1){
			nu.lam.prop	<- mu.lam 
			}
				
		lam.dum.prop	<- lambda.prop[j];
	
		if(lam.dum.prop > 700){
			lam.dum.prop <- 700;
			}
			
		if(J+1 >= 2){
			D.1st.prop	<- -exp(lam.dum.prop) * sum(exp.xbeta * Delta[,j] * ind.r[,j])
			D1.prop		<- D.1st.prop + d[j] - 1/(sig.sq.lam * Q[j,j])*							(lambda.prop[j] - nu.lam.prop);
			D2.prop		<- D.1st.prop - 1/(sig.sq.lam * Q[j,j]);
			}
		
		if(J+1 == 1){
			D.1st.prop	<- -exp(lam.dum.prop) * sum(exp.xbeta * Delta[,j] * ind.r[,j])
			D1.prop		<- D.1st.prop + d[j] - 1/(sig.sq.lam * Q)*							(lambda.prop[j] - nu.lam.prop);
			D2.prop		<- D.1st.prop - 1/(sig.sq.lam * Q);
			}	
	
		lambda.prop.me.ini  <- lambda.prop[j] - D1.prop / D2.prop
		lambda.prop.var.ini <- - 2.4^2 / D2.prop
		
		
		# Calculating acceptance probability
				
		loglh.ini	<- D.1st + lambda[j]*d[j]
		loglh.prop	<- D.1st.prop + lambda.prop[j]*d[j]
		
		if(J+1 > 1){	
			logprior.ini	<- dmvnorm(lambda, mean = rep(mu.lam, J+1), 										sigma = sig.sq.lam * Sigma.lam, log = TRUE)
			logprior.prop	<- dmvnorm(lambda.prop, mean = rep(mu.lam, J+1), 										sigma = sig.sq.lam * Sigma.lam, log = TRUE)
			}	
	
		if(J+1 == 1){	
			logprior.ini	<- dnorm(lambda, mean = mu.lam, 									sd = sqrt(sig.sq.lam * Sigma.lam), log = TRUE)
			logprior.prop	<- dnorm(lambda.prop, mean = mu.lam, 									sd = sqrt(sig.sq.lam * Sigma.lam), log = TRUE)
			}

		log.prop.prop    	<- dnorm(lambda.prop[j] , mean = lambda.prop.me, 									sd = sqrt(lambda.prop.var), log = TRUE)
		log.prop.ini     	<- dnorm(lambda[j]  , mean = lambda.prop.me.ini, 									sd = sqrt(lambda.prop.var.ini), log = TRUE)

		logR  <- loglh.prop - loglh.ini + logprior.prop - logprior.ini + 				log.prop.ini - log.prop.prop

		u = log(runif(1)) < logR
	
		if(u == 1){lambda[j] <- lambda.prop[j]}

		accept[j] <- accept[j] + u	
	
		}	# end of for loop	
	
	
	
	list(lambda = lambda, accept = accept, u = u)
	
	
	}	# end of MoveBH			
			



















			
move.SP.BpeSurv	<- function(survObj, ini, priorPara){
	
	n <- survObj$n
	p <- survObj$p	

	lambda			<- ini$lambda
	s				<- ini$s
	J				<- ini$J
	mu.lam			<- ini$mu.lam
	sig.sq.lam		<- ini$sig.sq.lam
	Delta			<- ini$Delta
	Sigma.lam		<- ini$Sigma.lam
	inv.Sigma.lam	<- ini$inv.Sigma.lam	
	a				<- priorPara$a
	b				<- priorPara$b
	
	one.inv.Sigma.lam.lam	<- t(rep(1, J+1)) %*% solve(Sigma.lam) %*% lambda;
	one.inv.Sigma.lam.one	<- t(rep(1, J+1)) %*% solve(Sigma.lam) %*% rep(1, J+1);

	mu.lam.mean	<- as.vector(one.inv.Sigma.lam.lam / one.inv.Sigma.lam.one);
	mu.lam.var	<- as.vector(sig.sq.lam / one.inv.Sigma.lam.one);
	
	mu.lam.sample		<- rnorm(1, mean = mu.lam.mean, sd = sqrt(mu.lam.var));

	mu.minus.lam	<- rep(mu.lam, J+1) - lambda;
	inv.sig.sq.lam.rate = b + 1/2 * as.vector(t(mu.minus.lam) %*% inv.Sigma.lam %*% mu.minus.lam);
	inv.sig.sq.lam	<- rgamma(1, shape = a + (J+1)/2, rate = inv.sig.sq.lam.rate);
	
	sig.sq.lam.sample	<- 1/inv.sig.sq.lam;
	
	list(mu.lam.sample = mu.lam.sample, sig.sq.lam.sample = sig.sq.lam.sample);
	
	}	# end of MovdSP			










move.BI.BpeSurv	<- function(survObj, ini, priorPara){

	n <- survObj$n
	p <- survObj$p	
	
	y		<- survObj$y
	delta	<- survObj$delta
	x		<- survObj$x
	s.propBI 	<- priorPara$s.propBI
	s			<- ini$s	
	J			<- ini$J		
	beta.ini	<- ini$beta.ini		
	lambda		<- ini$lambda		
	sig.sq.lam	<- ini$sig.sq.lam		
	mu.lam		<- ini$mu.lam
	ind.r		<- ini$ind.r
	ind.d		<- ini$ind.d
	d			<- ini$d		
	Delta		<- ini$Delta
	Sigma.lam	<- ini$Sigma.lam
	inv.Sigma.lam	<- ini$inv.Sigma.lam
	W			<- ini$W
	Q			<- ini$Q
	xbeta		<- ini$xbeta.ini	
	
	a		<- priorPara$a
	b		<- priorPara$b
	c.lam	<- priorPara$c.lam
	C		<- priorPara$C
	alpha	<- priorPara$alpha
	delPert <- priorPara$delPert
		
	accept <- 0
	
	
	### Birth ######
	
	s.propBI	<- setdiff(s.propBI, s)
	num.s.prop	<- length(s.propBI)
	
	s.star	<- s.propBI[sample(1:num.s.prop, 1)]
	j.old	<- min(which(s - s.star > 0))	
	
	s.new	<- sort(c(s, s.star))
	J.new	<- length(s.new) - 1
	
	U.pert	<- runif(1, 0.5 - delPert, 0.5 + delPert)

	if(j.old != 1){
		new.lam1	<- lambda[j.old] - (s[j.old] - s.star)/(s[j.old] - s[j.old - 1]) * log((1-U.pert)/U.pert);
		new.lam2	<- lambda[j.old] + (s.star - s[j.old-1])/(s[j.old] - s[j.old - 1]) * log((1-U.pert)/U.pert);
		}

	if(j.old == 1){
		new.lam1	<- lambda[j.old] - (s[j.old] - s.star)/(s[j.old] - 0) * log((1-U.pert)/U.pert);
		new.lam2	<- lambda[j.old] + (s.star - 0)/(s[j.old] - 0) * log((1-U.pert)/U.pert);
		}

	lambda.new	<- append(lambda, c(new.lam1, new.lam2), after = j.old)
	lambda.new	<- lambda.new[-j.old]	
	
	ini.new	<- ini
	
	ini.new$lambda	<- lambda.new
	ini.new$s		<- s.new
	ini.new$J		<- J.new
	
	old <- list()
	old$j.old <- j.old
	old$s.old <- s
	old$J.old <- J
	old$ind.r.old	<- ini$ind.r
	old$ind.d.old	<- ini$ind.d
	old$ind.d.temp.old	<- ini$ind.d.temp	
	old$d.old		<- ini$d
	
	intv.new		<- add.interval.BpeSurv(survObj, ini, priorPara, old, s.new)

	ini.new$ind.r	<- intv.new$ind.r
	ini.new$ind.d	<- intv.new$ind.d
	ini.new$ind.d.temp <- intv.new$ind.d.temp
	ini.new$d		<- intv.new$d

	ini.new$Delta			<- calDelta.BI.BpeSurv(survObj, ini.new, old)
	cal.Sig.new				<- cal.Sigma.lam.BpeSurv(ini.new, priorPara)
	ini.new$Sigma.lam		<- cal.Sig.new$Sigma.lam;
	ini.new$inv.Sigma.lam	<- solve(ini.new$Sigma.lam)
	ini.new$W				<- cal.Sig.new$W;
	ini.new$Q				<- cal.Sig.new$Q;
	
	Delta.new	<- ini.new$Delta
	ind.r.new	<- ini.new$ind.r
	ind.d.new	<- ini.new$ind.d
	ind.d.temp.new	<- ini.new$ind.d.temp	
	d.new		<- ini.new$d
	Sigma.lam.new	<- ini.new$Sigma.lam
	inv.Sigma.lam.new	<- ini.new$inv.Sigma.lam


	# log-likelhood ratio
		
	xbeta.mat					<- matrix(rep(xbeta, J+1), n, J+1)
	xbeta.dum					<- xbeta
	xbeta.dum[xbeta.dum > 700] 	<- 700	
	exp.xbeta					<- exp(xbeta.dum)
	
	lam.dum	<- lambda[j.old];
	
	if(lam.dum > 700){
		lam.dum <- 700;
		}	

	log.lh.ini.1st	<- -exp(lam.dum) * sum(exp.xbeta * Delta[,j.old] * 						ind.r[,j.old]) + lambda[j.old]*d[j.old]
	log.lh.ini		<- 	log.lh.ini.1st + sum(xbeta * ind.d[,j.old]) 
	
	lam.dum.new	<- lambda.new[c(j.old, (j.old + 1))];
	lam.dum.new[lam.dum.new > 700] <- 700
	
	log.lh.new.1st	<- -exp(lam.dum.new) * colSums(exp.xbeta * Delta.new[,c(j.old, j.old + 1)] * 						ind.r.new[,c(j.old, j.old + 1)]) + lambda.new[c(j.old, j.old + 1)]*						d.new[c(j.old, j.old + 1)]
	log.lh.new		<- 	sum(log.lh.new.1st + colSums((matrix(rep(xbeta, 2), n, 2) * 						ind.d.new[,c(j.old, j.old + 1)]))) 
		
	BI.lh.ratio		<- 	log.lh.new - log.lh.ini
	
	
	# log-prior ratio
	
	if(J+1 != 1){

		if(j.old != 1){
			BI.prior.ratio	<- log(alpha/(J+1)) + dmvnorm(lambda.new, mean = rep(mu.lam, J.new+1), 								sigma = sig.sq.lam * Sigma.lam.new, log = TRUE)- 								dmvnorm(lambda, mean = rep(mu.lam, J+1), 									sigma = sig.sq.lam * Sigma.lam, log = TRUE)+ 									log((2*J+3)*(2*J+2)*s[J+1]^(-2)*(s.star - s[j.old-1])*									(s[j.old] - s.star)/(s[j.old] - s[j.old-1]));
			}
		
		if(j.old == 1){
			BI.prior.ratio	<- log(alpha/(J+1)) + dmvnorm(lambda.new, mean = rep(mu.lam, J.new+1), 								sigma = sig.sq.lam * Sigma.lam.new, log = TRUE)- 								dmvnorm(lambda, mean = rep(mu.lam, J+1), 									sigma = sig.sq.lam * Sigma.lam, log = TRUE)+									log((2*J+3)*(2*J+2)*s[J+1]^(-2)*(s.star - 0)*									(s[j.old] - s.star)/(s[j.old] - 0));
			}	
		}
		
	if(J+1 == 1){
		BI.prior.ratio	<- log(alpha/(J+1)) + dmvnorm(lambda.new, mean = rep(mu.lam, J.new+1), 							sigma = sig.sq.lam * Sigma.lam.new, log = TRUE) - 							dnorm(lambda, mean = rep(mu.lam, J+1), sd = sqrt(sig.sq.lam * Sigma.lam), 								log = TRUE) + log((2*J+3)*(2*J+2)*s[J+1]^(-2)*								(s.star - 0)*(s[j.old] - s.star)/(s[j.old] - 0));
		}	

	# log-proposal ratio	

	BI.prop.ratio	<- log(num.s.prop/alpha) - dunif(U.pert, 0.5 - delPert, 0.5 + delPert, log = TRUE);
	
	# log-Jacobian	
	
	BI.Jacob		<- log(1/U.pert/(1-U.pert))
	
	# Acceptance probability	
	
	logR.BI	<- BI.lh.ratio + BI.prior.ratio + BI.prop.ratio + BI.Jacob
	
	u.BI	= log(runif(1)) < logR.BI

	if(u.BI == 1){
		ini	<- ini.new
		}
		accept <- accept + u.BI

	list(ini = ini, accept = accept);
	
	}	# end of BI

		



			
move.DI.BpeSurv	<- function(survObj, ini, priorPara){

	n <- survObj$n
	p <- survObj$p		
	y		<- survObj$y
	delta	<- survObj$delta
	x		<- survObj$x

	s.propBI 	<- priorPara$s.propBI		
	s			<- ini$s
	J			<- ini$J		
	beta.ini	<- ini$beta.ini		
	lambda		<- ini$lambda		
	sig.sq.lam	<- ini$sig.sq.lam		
	mu.lam		<- ini$mu.lam
	ind.r		<- ini$ind.r
	ind.d		<- ini$ind.d
	ind.d.temp	<- ini$ind.d.temp	
	d			<- ini$d		
	Delta		<- ini$Delta
	Sigma.lam	<- ini$Sigma.lam
	inv.Sigma.lam	<- ini$inv.Sigma.lam
	W			<- ini$W
	Q			<- ini$Q
	xbeta		<- ini$xbeta.ini	
	
	
	a		<- priorPara$a
	b		<- priorPara$b
	c.lam	<- priorPara$c.lam
	C		<- priorPara$C
	alpha	<- priorPara$alpha
	delPert <- priorPara$delPert
	
	accept <- 0
	
	s.propBI	<- setdiff(s.propBI, s)
	num.s.prop	<- length(s.propBI)
		
		
	# new time split and new lambda
	
	j.old	<- sample(c(1:(J)), 1)
	
	s.new	<- s[-j.old]
	J.new	<- length(s.new)-1
	
	U.pert	<- 1/(exp(lambda[j.old+1] - lambda[j.old]) + 1);
	
	if(j.old != 1){
		new.lam	<- ((s[j.old] - s[j.old-1]) * lambda[j.old] + 					(s[j.old+1] - s[j.old]) * lambda[j.old + 1]) / (s[j.old + 1] - s[j.old - 1]);
		}
		
	if(j.old == 1){
		
		new.lam	<- ((s[j.old] - 0) * lambda[j.old] + 					(s[j.old+1] - s[j.old]) * lambda[j.old + 1]) / (s[j.old + 1] - 0);
		}		
		
	lambda.new 			<- lambda[-j.old]
	lambda.new[j.old]	<- new.lam


	ini.new	<- ini
	
	ini.new$lambda	<- lambda.new
	ini.new$s		<- s.new
	ini.new$J		<- J.new
	
	old <- list()
	old$j.old <- j.old
	old$s.old <- s
	old$J.old <- J
	old$ind.r.old <- ini$ind.r
	old$ind.d.old <- ini$ind.d
	old$ind.d.temp.old <- ini$ind.d.temp	
	old$d.old <- ini$d

	intv.new		<- remove.interval.BpeSurv(survObj, ini, priorPara, old, s.new)
	ini.new$ind.r	<- intv.new$ind.r
	ini.new$ind.d	<- intv.new$ind.d
	ini.new$ind.d.temp	<- intv.new$ind.d.temp	
	ini.new$d		<- intv.new$d	
		
	ini.new$Delta			<- cal.Del.DI.BpeSurv(survObj, ini.new, old)
	cal.Sig.new				<- cal.Sigma.lam.BpeSurv(ini.new, priorPara)
	ini.new$Sigma.lam		<- cal.Sig.new$Sigma.lam;
	ini.new$inv.Sigma.lam	<- solve(ini.new$Sigma.lam)
	ini.new$W				<- cal.Sig.new$W;
	ini.new$Q				<- cal.Sig.new$Q;
		
	Delta.new	<- ini.new$Delta
	ind.r.new	<- ini.new$ind.r
	ind.d.new	<- ini.new$ind.d
	d.new		<- ini.new$d
	Sigma.lam.new	<- ini.new$Sigma.lam
	
	
	# log-likelhood ratio
	
	xbeta.mat					<- matrix(rep(xbeta, J+1), n, J+1)
	xbeta.dum					<- xbeta
	xbeta.dum[xbeta.dum > 700] 	<- 700	
	exp.xbeta					<- exp(xbeta.dum)
	
	lam.dum	<- lambda[c(j.old, j.old + 1)];
	lam.dum[lam.dum> 700] <- 700

	log.lh.ini.1st	<- -exp(lam.dum) * 									colSums(exp.xbeta * Delta[,c(j.old, j.old + 1)] * ind.r[,c(j.old, j.old + 1)]) + 						lambda[c(j.old, j.old + 1)]*d[c(j.old, j.old + 1)]
	log.lh.ini		<- 	sum(log.lh.ini.1st + 							colSums((matrix(rep(xbeta, 2), n, 2) * ind.d[,c(j.old, j.old + 1)]))) 

	lam.dum.new	<- lambda.new[j.old];
	
	if(lam.dum.new > 700){
		lam.dum.new <- 700;
		}	

	log.lh.new.1st	<- -exp(lam.dum.new) * sum(exp.xbeta * Delta.new[,j.old] * ind.r.new[,j.old]) + lambda.new[j.old]*d.new[j.old]
	log.lh.new		<- 	log.lh.new.1st + sum(xbeta * ind.d.new[,j.old]) 
		
	DI.lh.ratio		<- 	log.lh.new - log.lh.ini
	
	
	# log-prior ratio	
	
	if(J+1 != 2){
		if(j.old != 1){
			DI.prior.ratio	<- log((J)/alpha) + dmvnorm(lambda.new, mean = rep(mu.lam, J), sigma = sig.sq.lam * Sigma.lam.new, log = TRUE) -							 dmvnorm(lambda, mean = rep(mu.lam, J+1), sigma = sig.sq.lam * Sigma.lam, log = TRUE) + 							log(1/(2*J+1)/(2*J)*s[J+1]^2*(s[j.old+1] - s[j.old-1])/(s[j.old] - s[j.old - 1])/(s[j.old+1] - s[j.old]));
			}
		
		if(j.old == 1){
			DI.prior.ratio	<- log((J)/alpha) + dmvnorm(lambda.new, mean = rep(mu.lam, J), sigma = sig.sq.lam * Sigma.lam.new, log = TRUE) -							dmvnorm(lambda, mean = rep(mu.lam, J+1), sigma = sig.sq.lam * Sigma.lam, log = TRUE) + 							log(1/(2*J+1)/(2*J)*s[J+1]^2*(s[j.old+1] - 0)/(s[j.old] - 0)/(s[j.old+1] - s[j.old]));
			}		
		}	

	if(J+1 == 2){
		DI.prior.ratio	<- log((J)/alpha) + dnorm(lambda.new, mean = mu.lam, sd = sqrt(sig.sq.lam * Sigma.lam.new), log = TRUE) - 						dmvnorm(lambda, mean = rep(mu.lam, J+1), sigma = sig.sq.lam * Sigma.lam, log = TRUE) + 						log(1/(2*J+1)/(2*J)*s[J+1]^2*(s[j.old+1] - 0)/(s[j.old] - 0)/(s[j.old+1] - s[j.old]));
		}	


	# log-proposal ratio	

	DI.prop.ratio	<- 	log(alpha/num.s.prop) + dunif(U.pert, 0.5 - delPert, 0.5 + delPert, log = TRUE);


	# log-Jacobian

	DI.Jacob		<- log(U.pert*(1-U.pert))
	
	
	# Acceptance probability	
	
	logR.DI	<- DI.lh.ratio + DI.prior.ratio + DI.prop.ratio + DI.Jacob
	
	u.DI	= log(runif(1)) < logR.DI		

	if(u.DI == 1){

		ini	<- ini.new			

		accept	<- accept + u.DI

		}
		
	list(ini = ini, accept = accept);
	
	}	# end of MoveDI		
			










		
	