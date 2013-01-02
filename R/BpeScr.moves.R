move.RP.BpeScr <- function(survObj, ini, g){
	
	n <- survObj$n 
	p <- survObj$p
	x <- survObj$x
	
	mh <- list()

	mh$betaUB = FALSE
	mh$betaPropUBval			<- 1
	mh$betaPropVar				<- 0.01	
			
	gamma.ini	<- ini$gamma.ini	
	
	if(g == 1){
		xbeta		<- ini$xbeta1.ini
		beta.ini	<- ini$beta1.ini
		lambda		<- ini$lambda1
		s			<- ini$s1
		J			<- ini$J1
		ind.r		<- ini$ind.r1.sc1
		ind.d		<- ini$ind.d1.sc1
		d			<- ini$d1.sc1
		Delta		<- ini$Delta1
		}	

	if(g == 2){
		xbeta		<- ini$xbeta2.ini		
		beta.ini	<- ini$beta2.ini
		lambda		<- ini$lambda2
		s			<- ini$s2
		J			<- ini$J2
		ind.r		<- ini$ind.r1.sc2
		ind.d		<- ini$ind.d2.sc2
		d			<- ini$d2.sc2
		Delta		<- ini$Delta2
		}
		
	if(g == 3){
		xbeta		<- ini$xbeta3.ini		
		beta.ini	<- ini$beta3.ini
		lambda		<- ini$lambda3
		s			<- ini$s3
		J			<- ini$J3
		ind.r		<- ini$ind.r2.sc3
		ind.d		<- ini$ind.d3.sc3
		d			<- ini$d3.sc3
		Delta		<- ini$Delta3
		}		
						
	accept	<- rep(0, p);
	
	
#	updatebeta = sample(c(1:p), mh$numUpdateBeta, replace = FALSE)
	updatebeta = sample(c(1:p), 1, replace = FALSE)

	
	for(j in updatebeta){	
		
		# For Gaussian jumping distribution
				
		xbeta.mat		<- matrix(rep(xbeta, J+1), n, J+1);
		xbeta.dum		<- xbeta;
		xbeta.dum[xbeta.dum > 700] <- 700;	
		exp.xbeta		<- exp(xbeta.dum);
		exp.xbeta.mat	<- matrix(rep(exp.xbeta, J+1), n, J+1);
		xbeta.mat		<- matrix(rep(xbeta, J+1), n, J+1);
		
		gamma.mat		<- matrix(rep(gamma.ini, J+1), n, J+1);
		x.j.mat <- matrix(rep(x[,j],J+1),J+1, nrow = n);
		lambda[lambda > 700] <- 700;
	
		D1.1st <- -exp(lambda) * colSums(Delta * gamma.mat * exp.xbeta.mat * x.j.mat * ind.r);
		D1.2nd <- colSums(x.j.mat * ind.d);
		D1	<- sum(D1.1st + D1.2nd);	
		D2	<- sum(-exp(lambda) * colSums(Delta * gamma.mat * exp.xbeta.mat * x.j.mat^2 * ind.r));
	
		# proposal value
	
		beta.prop.me  <- beta.ini[j] - D1 / D2;
		beta.prop.var <- -2.4^2 / D2;
		
		beta.prop		<- beta.ini;
		
		if(mh$betaUB){
			if(abs(D1/D2) > mh$betaPropUBval){	
				beta.prop.me  	<- beta.ini[j] - sign(D1/D2) * mh$betaPropUBval;
				beta.prop.var	<- mh$betaPropVar
				}				
			if(D1 == D2){
				beta.prop.me  <- beta.ini[j];
				beta.prop.var	<- mh$betaPropVar
				}
			}
	
		beta.prop[j]	<- rnorm(1, mean = beta.prop.me, sd = sqrt(beta.prop.var));
		
		# Calculating acceptance probability

		xbeta.prop	<- xbeta - x[,j] * beta.ini[j] + x[,j] * beta.prop[j]
		
		xbeta.prop			<- as.vector(x %*% beta.prop)
		xbeta.mat.prop		<- matrix(rep(xbeta.prop, J+1), n, J+1);
		xbeta.dum.prop		<- xbeta.prop;
		xbeta.dum.prop[xbeta.dum.prop > 700] <- 700;
		exp.xbeta.prop		<- exp(xbeta.dum.prop);
		exp.xbeta.mat.prop	<- matrix(rep(exp.xbeta.prop, J+1), n, J+1);
		xbeta.mat.prop		<- matrix(rep(xbeta.prop, J+1), n, J+1);
		
		D1.1st.prop <- -exp(lambda) * colSums(Delta * gamma.mat * exp.xbeta.mat.prop * x.j.mat * ind.r);
		D1.2nd.prop <- colSums(x.j.mat * ind.d);
		D1.prop		<- sum(D1.1st.prop + D1.2nd.prop);
		D2.prop		<- sum(-exp(lambda) * colSums(Delta * gamma.mat * exp.xbeta.mat.prop * x.j.mat^2 * ind.r));
		
		beta.prop.me.ini  <- beta.prop[j] - D1.prop / D2.prop;
		beta.prop.var.ini <- -2.4^2 / D2.prop;	
			
		loglh.ini.1st	<- -exp(lambda) * colSums(Delta * gamma.mat * exp.xbeta.mat * ind.r);
		loglh.ini.2nd	<- colSums(xbeta.mat*ind.d);
		loglh.ini		<- sum(loglh.ini.1st + loglh.ini.2nd);
		
		loglh.prop.1st	<- -exp(lambda) * colSums(Delta * gamma.mat * exp.xbeta.mat.prop * ind.r);
		loglh.prop.2nd	<- colSums(xbeta.mat.prop*ind.d);
		loglh.prop		<- sum(loglh.prop.1st + loglh.prop.2nd);	 		
		logprop.prop    	<- dnorm(beta.prop[j] , mean = beta.prop.me.ini, sd = sqrt(beta.prop.var.ini), log = TRUE);
		logprop.ini     	<- dnorm(beta.ini[j]  , mean = beta.prop.me, sd = sqrt(beta.prop.var), log = TRUE);

		logR  <- loglh.prop - loglh.ini + logprop.ini - logprop.prop;

		u = log(runif(1)) < logR
	
		if(u == 1){
			beta.ini[j] <- beta.prop[j]
			xbeta	<- xbeta.prop
			}
		
		accept[j]	<- accept[j] + u;

		} #end of 'for' loop
		
	list(beta.ini = beta.ini, accept = accept, xbeta.ini = xbeta)
		
	}







move.BH.BpeScr <-
function(survObj, ini, g){
	
	n <- survObj$n
	p <- survObj$p		

	mh <- list()

	mh$lambdaUB = TRUE
	mh$lambdaUBval				<- 1
	mh$lambdaPropVar			<- 1
	
	gamma.ini	<- ini$gamma.ini	
	
	if(g == 1){
		xbeta		<- ini$xbeta1.ini
		beta.ini	<- ini$beta1.ini
		lambda		<- ini$lambda1
		s			<- ini$s1
		J			<- ini$J1
		ind.r		<- ini$ind.r1.sc1
		d			<- ini$d1.sc1
		Delta		<- ini$Delta1
		mu.lam		<- ini$mu.lam1
		sigma.sq.lam<- ini$sig.sq.lam1
		Sigma.lam	<- ini$Sigma.lam1
		W			<- ini$W1
		Q			<- ini$Q1 
		}		
	
	if(g == 2){
		xbeta		<- ini$xbeta2.ini		
		beta.ini	<- ini$beta2.ini
		lambda		<- ini$lambda2
		s			<- ini$s2
		J			<- ini$J2
		ind.r		<- ini$ind.r1.sc2
		d			<- ini$d2.sc2
		Delta		<- ini$Delta2
		mu.lam		<- ini$mu.lam2
		sigma.sq.lam<- ini$sig.sq.lam2
		Sigma.lam	<- ini$Sigma.lam2
		W			<- ini$W2
		Q			<- ini$Q2 
		}
		
	
	if(g == 3){
		xbeta		<- ini$xbeta3.ini		
		beta.ini	<- ini$beta3.ini
		lambda		<- ini$lambda3
		s			<- ini$s3
		J			<- ini$J3
		ind.r		<- ini$ind.r2.sc3
		d			<- ini$d3.sc3
		Delta		<- ini$Delta3
		mu.lam		<- ini$mu.lam3
		sigma.sq.lam<- ini$sig.sq.lam3
		Sigma.lam	<- ini$Sigma.lam3
		W			<- ini$W3
		Q			<- ini$Q3 
		}		
		
	inv.Sigma.lam	<- solve(Sigma.lam)		
					
	accept	<- rep(0, J+1);

	xbeta.mat	<- matrix(rep(xbeta, J+1), n, J+1)
	xbeta.dum	<- xbeta;
	xbeta.dum[xbeta.dum > 700] <- 700;	
	exp.xbeta	<- exp(xbeta.dum)	
	
	#j = sample(c(1:(J+1)), mh$numUpdateLambda, replace = FALSE)
	
	for(j in 1:(J+1)){
		
		# For Gaussian jumping distribution
				
		sum.exp.xbeta	<- sum(gamma.ini * exp.xbeta * Delta[,j] * ind.r[,j]);
						
		if(J+1 > 1){
			if(j == 1){		
				nu.lam	<- mu.lam + W[1, 2] * (lambda[2] - mu.lam)
				}
		
			if(j == J+1){
				nu.lam	<- mu.lam + W[J+1, (J)] * (lambda[J] - mu.lam);
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
			D1	<- -exp(lam.dum) * sum.exp.xbeta + d[j] - 1/(sigma.sq.lam * Q[j,j])*(lambda[j] - nu.lam);
			D2	<- -exp(lam.dum) * sum.exp.xbeta - 1/(sigma.sq.lam * Q[j,j]);
			}
		
		if(J+1 == 1){
			D1	<- -exp(lam.dum) * sum.exp.xbeta + d[j] - 1/(sigma.sq.lam * Q)*(lambda[j] - nu.lam);
			D2	<- -exp(lam.dum) * sum.exp.xbeta - 1/(sigma.sq.lam * Q);
			}
					
		# proposal value
	
		lambda.prop.me  <- lambda[j] - D1 / D2;
		lambda.prop.var <- -2.4^2 / D2;	

		
		
		if(mh$lambdaUB){
			if(abs(D1/D2) > mh$lambdaUBval){	
				lambda.prop.me  <- lambda[j] - sign(D1/D2) * mh$lambdaUBval;
				lambda.prop.var	<- mh$lambdaPropVar;
				}			
			}

		lambda.prop			<- rnorm(1, mean = lambda.prop.me, sd = sqrt(lambda.prop.var));
		lambda.propVec		<- lambda
		lambda.propVec[j]	<- lambda.prop		
		
		if(J+1 > 1){
			if(j == 1){		
				nu.lam.prop	<- mu.lam + W[1, 2] * (lambda.propVec[2] - mu.lam)
				}
		
			if(j == J+1){
				nu.lam.prop	<- mu.lam + W[J+1, (J)] * (lambda.propVec[J] - mu.lam);
				}
		
			if(j != 1 & j != J+1){
				nu.lam.prop	<- mu.lam + W[j, (j-1)] * (lambda.propVec[j-1] - mu.lam) 							+ W[j, (j+1)] * (lambda.propVec[j+1] - mu.lam);
				}
			}
		
		if(J+1 == 1){
			nu.lam.prop	<- mu.lam 
			}
				
		lam.dum.prop	<- lambda.prop;
		
		if(lam.dum.prop > 700){lam.dum.prop <- 700};
			
		if(J+1 >= 2){
			D1.prop	<- -exp(lam.dum.prop) * sum.exp.xbeta + d[j] - 1/(sigma.sq.lam * Q[j,j])*(lambda.propVec[j] - nu.lam.prop);
			D2.prop	<- -exp(lam.dum.prop) * sum.exp.xbeta - 1/(sigma.sq.lam * Q[j,j]);
			}
		
		if(J+1 == 1){
			D1.prop	<- -exp(lam.dum.prop) * sum.exp.xbeta + d[j] - 1/(sigma.sq.lam * Q)*(lambda.propVec[j] - nu.lam.prop);
			D2.prop	<- -exp(lam.dum.prop) * sum.exp.xbeta - 1/(sigma.sq.lam * Q);
			}
				

		
		lambda.prop.me.ini  <- lambda.prop - D1.prop / D2.prop;
		lambda.prop.var.ini <- -2.4^2 / D2.prop;
						
		# Calculating acceptance probability
		
		loglh.ini	<- -exp(lam.dum) * sum.exp.xbeta + lambda[j] * d[j] 
		loglh.prop	<- -exp(lam.dum.prop) * sum.exp.xbeta + lambda.prop * d[j]
	
		if(J+1 > 1){	
			logprior.ini	<- dmvnorm(lambda, mean = rep(mu.lam, J+1), sigma = sigma.sq.lam * Sigma.lam, log = TRUE)
			logprior.prop	<- dmvnorm(lambda.propVec, mean = rep(mu.lam, J+1), sigma = sigma.sq.lam * Sigma.lam, log = TRUE)
			}	
	
		if(J+1 == 1){	
			logprior.ini	<- dnorm(lambda[j], mean = mu.lam, sd = sqrt(sigma.sq.lam * Sigma.lam), log = TRUE)
			logprior.prop	<- dnorm(lambda.prop, mean = mu.lam, sd = sqrt(sigma.sq.lam * Sigma.lam), log = TRUE)
			}
					
		logprop.ini			<- dnorm(lambda[j], mean = lambda.prop.me.ini, sd = sqrt(lambda.prop.var.ini), log = TRUE);
		logprop.prop		<- dnorm(lambda.prop, mean = lambda.prop.me, sd = sqrt(lambda.prop.var), log = TRUE);
		
		logR  <- loglh.prop - loglh.ini + logprior.prop - logprior.ini + logprop.ini - logprop.prop;

		u = log(runif(1)) < logR;
	
		if(u == 1){lambda[j] <- lambda.prop}

		accept[j] <- accept[j] + u;
		
		} #end of the 'for' loop
		
	list(lambda = lambda, accept = accept)
	
	}




move.SP.BpeScr <-
function(survObj, ini, priorPara, g){
	
	n <- survObj$n
	p <- survObj$p
	
	if(g == 1){
		lambda	<- ini$lambda1
		s		<- ini$s1
		J		<- ini$J1
		mu.lam	<- ini$mu.lam1
		sigma.sq.lam	<- ini$sig.sq.lam1
		Delta			<- ini$Delta1
		Sigma.lam		<- ini$Sigma.lam1
		a		<- priorPara$a1
		b		<- priorPara$b1
		}
	
	if(g == 2){
		lambda	<- ini$lambda2
		s		<- ini$s2
		J		<- ini$J2
		mu.lam	<- ini$mu.lam2
		sigma.sq.lam	<- ini$sig.sq.lam2
		Delta			<- ini$Delta2
		Sigma.lam		<- ini$Sigma.lam2
		a		<- priorPara$a2
		b		<- priorPara$b2
		}
		
	if(g == 3){
		lambda	<- ini$lambda3
		s		<- ini$s3
		J		<- ini$J3
		mu.lam	<- ini$mu.lam3
		sigma.sq.lam	<- ini$sig.sq.lam3
		Delta			<- ini$Delta3
		Sigma.lam		<- ini$Sigma.lam3
		a		<- priorPara$a3
		b		<- priorPara$b3
		}
				
			
	inv.Sigma.lam	<- solve(Sigma.lam);
		
	one.inv.Sigma.lam.lam	<- t(rep(1, J+1)) %*% solve(Sigma.lam) %*% lambda;
	one.inv.Sigma.lam.one	<- t(rep(1, J+1)) %*% solve(Sigma.lam) %*% rep(1, J+1);
	
	mu.lam.mean	<- as.vector(one.inv.Sigma.lam.lam / one.inv.Sigma.lam.one);
	mu.lam.var	<- as.vector(sigma.sq.lam / one.inv.Sigma.lam.one);
	
	mu.lam.sample		<- rnorm(1, mean = mu.lam.mean, sd = sqrt(mu.lam.var));
	
	mu.minus.lam	<- rep(mu.lam.sample, J+1) - lambda;
	inv.sig.sq.lam.rate = b + 1/2 * as.vector(t(mu.minus.lam) %*% inv.Sigma.lam %*% mu.minus.lam);
	inv.sig.sq.lam	<- rgamma(1, shape = a + (J+1)/2, rate = inv.sig.sq.lam.rate);
	
	sigma.sq.lam.sample	<- 1/inv.sig.sq.lam;	
	
	list(mu.lam.sample = mu.lam.sample, sigma.sq.lam.sample = sigma.sq.lam.sample);
	
	}





move.BI.BpeScr <-
function(survObj, ini, priorPara, g){
		
	n		<- length(survObj$y1)
	p		<- ncol(survObj$x)
	
	y1		<- survObj$y1
	y2		<- survObj$y2	
	delta1	<- survObj$delta1
	delta2	<- survObj$delta2	
	x		<- survObj$x

	if(g == 1){
		xbeta		<- ini$xbeta1.ini
		s.propBI 	<- priorPara$s.propBI1		
		s			<- ini$s1	
		J			<- ini$J1		
		beta.ini	<- ini$beta1.ini		
		lambda		<- ini$lambda1		
		sig.sq.lam	<- ini$sig.sq.lam1		
		mu.lam		<- ini$mu.lam1
		ind.r		<- ini$ind.r1.sc1
		ind.d		<- ini$ind.d1.sc1
		ind.d.temp	<- ini$ind.d.temp.sc1		
		d			<- ini$d1.sc1		
		Delta		<- ini$Delta1
		Sigma.lam	<- ini$Sigma.lam1
		inv.Sigma.lam	<- ini$inv.Sigma.lam1
		W			<- ini$W1
		Q			<- ini$Q1
		a		<- priorPara$a1
		b		<- priorPara$b1
		c.lam	<- priorPara$c.lam1
		C		<- priorPara$C1
		alpha	<- priorPara$alpha1
		delPert <- priorPara$delPert1
		}	
		
	if(g == 2){
		xbeta		<- ini$xbeta2.ini		
		s.propBI 	<- priorPara$s.propBI2		
		s			<- ini$s2	
		J			<- ini$J2		
		beta.ini	<- ini$beta2.ini		
		lambda		<- ini$lambda2		
		sig.sq.lam	<- ini$sig.sq.lam2		
		mu.lam		<- ini$mu.lam2
		ind.r		<- ini$ind.r1.sc2
		ind.d		<- ini$ind.d2.sc2
		ind.d.temp	<- ini$ind.d.temp.sc2		
		d			<- ini$d2.sc2		
		Delta		<- ini$Delta2
		Sigma.lam	<- ini$Sigma.lam2
		inv.Sigma.lam	<- ini$inv.Sigma.lam2
		W			<- ini$W2
		Q			<- ini$Q2
		a		<- priorPara$a2
		b		<- priorPara$b2
		c.lam	<- priorPara$c.lam2
		C		<- priorPara$C2
		alpha	<- priorPara$alpha2
		delPert <- priorPara$delPert2
		}	

	if(g == 3){
		xbeta		<- ini$xbeta3.ini		
		s.propBI 	<- priorPara$s.propBI3		
		s			<- ini$s3	
		J			<- ini$J3		
		beta.ini	<- ini$beta3.ini		
		lambda		<- ini$lambda3		
		sig.sq.lam	<- ini$sig.sq.lam3		
		mu.lam		<- ini$mu.lam3
		ind.r		<- ini$ind.r2.sc3
		ind.d		<- ini$ind.d3.sc3
		ind.d.temp	<- ini$ind.d.temp.sc3		
		d			<- ini$d3.sc3		
		Delta		<- ini$Delta3
		Sigma.lam	<- ini$Sigma.lam3
		inv.Sigma.lam	<- ini$inv.Sigma.lam3
		W			<- ini$W3
		Q			<- ini$Q3
		a		<- priorPara$a3
		b		<- priorPara$b3
		c.lam	<- priorPara$c.lam3
		C		<- priorPara$C3
		alpha	<- priorPara$alpha3
		delPert <- priorPara$delPert3
		}	
		
	theta.ini	<- ini$theta.ini
	gamma.ini	<- ini$gamma.ini
	accept = 0
	
	
	### Birth ######
	
	s.propBI	<- setdiff(s.propBI, s)
	num.s.prop	<- length(s.propBI)
	
	s.star	<- s.propBI[sample(1:num.s.prop, 1)]
	j.old	<- min(which(s - s.star > 0))	
	
	#s.star	<- runif(1, 0, max(s))
	#j.old	<- min(which(s - s.star > 0))		
	
	s.new	<- sort(c(s, s.star))
	J.new	<- length(s.new) - 1
	
	U.pert	<- runif(1, 0.5 - delPert, 0.5 + delPert)

	if(j.old != 1){
		new.lam1	<- lambda[j.old] - (s[j.old] - s.star)/(s[j.old] - s[j.old - 1]) * log((1-U.pert)/U.pert);
		new.lam2	<- lambda[j.old] + (s.star - s[j.old-1])/(s[j.old] - s[j.old - 1]) * log((1-U.pert)/U.pert);
		}

	if(j.old == 1){s
		new.lam1	<- lambda[j.old] - (s[j.old] - s.star)/(s[j.old] - 0) * log((1-U.pert)/U.pert);
		new.lam2	<- lambda[j.old] + (s.star - 0)/(s[j.old] - 0) * log((1-U.pert)/U.pert);
		}

	lambda.new	<- append(lambda, c(new.lam1, new.lam2), after = j.old)
	lambda.new	<- lambda.new[-j.old]
		
	ini.new	<- ini
	
	if(g == 1){
		ini.new$lambda1	<- lambda.new
		ini.new$s1		<- s.new
		ini.new$J1		<- J.new				

		old <- list()
		old$j.old <- j.old
		old$s.old <- s
		old$J.old <- J
		old$ind.r1.old <- ini$ind.r1.sc1
		old$ind.r2.old <- ini$ind.r2.sc1
		old$ind.d1.old <- ini$ind.d1.sc1
		old$ind.d2.old <- ini$ind.d2.sc1
		old$ind.d3.old <- ini$ind.d3.sc1
		old$ind.d.temp.old <- ini$ind.d.temp.sc1
		old$d1.old <- ini$d1.sc1
		old$d2.old <- ini$d2.sc1
		old$d3.old <- NULL
		
		intv.new		<- add.interval.BpeScr(survObj, ini, priorPara, old, s.new, g = 1)
				
		ini.new$ind.r1.sc1	<- intv.new$ind.r1;
		ini.new$ind.d1.sc1	<- intv.new$ind.d1;
		ini.new$ind.d2.sc1	<- intv.new$ind.d2;
		ini.new$ind.d.temp.sc1	<- intv.new$ind.d.temp
		ini.new$d1.sc1		<- intv.new$d1;
		ini.new$d2.sc1		<- intv.new$d2;		
		
		ini.new$Delta1 <- cal.Del1.2.BI.BpeScr(survObj, ini.new, old, g = 1)
		
		cal.Sig1.new			<- cal.Sigma.lam.BpeScr(ini.new, priorPara, g = 1);
		ini.new$Sigma.lam1		<- cal.Sig1.new$Sigma.lam;
		ini.new$inv.Sigma.lam1	<- solve(ini.new$Sigma.lam1)
		ini.new$W1				<- cal.Sig1.new$W;
		ini.new$Q1				<- cal.Sig1.new$Q;
		
		Delta.new	<- ini.new$Delta1
		ind.r.new	<- ini.new$ind.r1.sc1
		ind.d.new	<- ini.new$ind.d1.sc1		
		d.new		<- ini.new$d1.sc1
		Sigma.lam.new	<- ini.new$Sigma.lam1
		}	
	
	if(g == 2){
		ini.new$lambda2	<- lambda.new
		ini.new$s2		<- s.new
		ini.new$J2		<- J.new				

		old <- list()
		old$j.old <- j.old
		old$s.old <- s
		old$J.old <- J
		old$ind.r1.old <- ini$ind.r1.sc2
		old$ind.r2.old <- ini$ind.r2.sc2
		old$ind.d1.old <- ini$ind.d1.sc2
		old$ind.d2.old <- ini$ind.d2.sc2
		old$ind.d3.old <- ini$ind.d3.sc2
		old$ind.d.temp.old <- ini$ind.d.temp.sc2
		old$d1.old <- ini$d1.sc2
		old$d2.old <- ini$d2.sc2
		old$d3.old <- ini$d3.sc2
		
		intv.new		<- add.interval.BpeScr(survObj, ini, priorPara, old, s.new, g = 2)
	
		ini.new$ind.r1.sc2	<- intv.new$ind.r1;
		ini.new$ind.d1.sc2	<- intv.new$ind.d1;
		ini.new$ind.d2.sc2	<- intv.new$ind.d2;
		ini.new$ind.d.temp.sc2	<- intv.new$ind.d.temp
		ini.new$d1.sc2		<- intv.new$d1;
		ini.new$d2.sc2		<- intv.new$d2;		

		ini.new$Delta2			<- cal.Del1.2.BI.BpeScr(survObj, ini.new, old, g = 2)
		cal.Sig2.new			<- cal.Sigma.lam.BpeScr(ini.new, priorPara, g = 2);
		ini.new$Sigma.lam2		<- cal.Sig2.new$Sigma.lam;
		ini.new$inv.Sigma.lam2	<- solve(ini.new$Sigma.lam2)
		ini.new$W2				<- cal.Sig2.new$W;
		ini.new$Q2				<- cal.Sig2.new$Q;
		
		Delta.new	<- ini.new$Delta2
		ind.r.new	<- ini.new$ind.r1.sc2
		ind.d.new	<- ini.new$ind.d2.sc2		
		d.new		<- ini.new$d2.sc2
		Sigma.lam.new	<- ini.new$Sigma.lam2
		}	
		
	if(g == 3){
		ini.new$lambda3	<- lambda.new
		ini.new$s3		<- s.new
		ini.new$J3		<- J.new				

		old <- list()
		old$j.old <- j.old
		old$s.old <- s
		old$J.old <- J
		old$ind.r1.old <- ini$ind.r1.sc3
		old$ind.r2.old <- ini$ind.r2.sc3
		old$ind.d1.old <- ini$ind.d1.sc3
		old$ind.d2.old <- ini$ind.d2.sc3
		old$ind.d3.old <- ini$ind.d3.sc3
		old$ind.d.temp.old <- ini$ind.d.temp.sc3
		old$d1.old <- ini$d1.sc3
		old$d2.old <- ini$d2.sc3		
		old$d3.old <- ini$d3.sc3
		
		intv.new		<- add.interval.BpeScr(survObj, ini, priorPara, old, s.new, g = 3)
		
		ini.new$ind.r2.sc3	<- intv.new$ind.r2;
		ini.new$ind.d1.sc3	<- intv.new$ind.d1;
		ini.new$ind.d3.sc3	<- intv.new$ind.d3;
		ini.new$ind.d.temp.sc3	<- intv.new$ind.d.temp
		ini.new$d1.sc3		<- intv.new$d1;
		ini.new$d3.sc3		<- intv.new$d3;	
		
		ini.new$Delta3			<- cal.Del3.BpeScr(survObj, ini.new, g = 3);
		
		#ini.new$Delta3			<- cal.Del3.BI.BpeScr(survObj, ini.new, old, g = 3)
		
		cal.Sig3.new			<- cal.Sigma.lam.BpeScr(ini.new, priorPara, g = 3);
		ini.new$Sigma.lam3		<- cal.Sig3.new$Sigma.lam;
		ini.new$inv.Sigma.lam3	<- solve(ini.new$Sigma.lam3)
		ini.new$W3				<- cal.Sig3.new$W;
		ini.new$Q3				<- cal.Sig3.new$Q;
		
		Delta.new	<- ini.new$Delta3
		ind.r.new	<- ini.new$ind.r2.sc3
		ind.d.new	<- ini.new$ind.d3.sc3		
		d.new		<- ini.new$d3.sc3
		Sigma.lam.new	<- ini.new$Sigma.lam3
		}
				
	# log-likelhood ratio
	
		xbeta.dum.ini		<- xbeta;
		xbeta.dum.ini[xbeta.dum.ini > 700] <- 700;
		exp.xbeta.ini		<- exp(xbeta.dum.ini);
		
		loglh.ini.1st	<- -exp(lambda[j.old]) * sum(Delta[,j.old] * gamma.ini * exp.xbeta.ini * ind.r[,j.old])
		loglh.ini.2nd	<- lambda[j.old] * d[j.old] 
		loglh.ini		<- loglh.ini.1st + loglh.ini.2nd
			
		loglh.new.1st	<- -exp(lambda.new[c(j.old, j.old+1)])*colSums(Delta.new[,c(j.old, j.old+1)] * 							gamma.ini * exp.xbeta.ini * ind.r.new[,c(j.old, j.old+1)])
		loglh.new.2nd	<- lambda.new[c(j.old, j.old+1)] * d.new[c(j.old, j.old+1)]
		loglh.new		<- sum(loglh.new.1st + loglh.new.2nd)
			
		BI.lh.ratio	<- loglh.new - loglh.ini
		
		
	# log-prior ratio
	
	if(J+1 != 1){
		if(j.old != 1){
			BI.prior.ratio	<- log(alpha/(J+1)) + dmvnorm(lambda.new, mean = rep(mu.lam, J.new+1), 								sigma = sig.sq.lam * Sigma.lam.new, log = TRUE)- 								dmvnorm(lambda, mean = rep(mu.lam, J+1), sigma = sig.sq.lam * Sigma.lam, log = TRUE)+								log((2*J+3)*(2*J+2)*s[J+1]^(-2)*(s.star - s[j.old-1])*(s[j.old] - s.star)/								(s[j.old] - s[j.old-1]));
			}
		
		if(j.old == 1){
			BI.prior.ratio	<- log(alpha/(J+1)) + dmvnorm(lambda.new, mean = rep(mu.lam, J.new+1), 								sigma = sig.sq.lam * Sigma.lam.new, log = TRUE)- 								dmvnorm(lambda, mean = rep(mu.lam, J+1), sigma = sig.sq.lam * Sigma.lam, log = TRUE)+								log((2*J+3)*(2*J+2)*s[J+1]^(-2)*(s.star - 0)*(s[j.old] - s.star)/(s[j.old] - 0));
			}	
		}
		
	if(J+1 == 1){
		BI.prior.ratio	<- log(alpha/(J+1)) + dmvnorm(lambda.new, mean = rep(mu.lam, J.new+1), 							sigma = sig.sq.lam * Sigma.lam.new, log = TRUE) - 							dnorm(lambda, mean = rep(mu.lam, J+1), sd = sqrt(sig.sq.lam * Sigma.lam), log = TRUE) + 							log((2*J+3)*(2*J+2)*s[J+1]^(-2)*(s.star - 0)*(s[j.old] - s.star)/(s[j.old] - 0));
		}
		
	# log-proposal ratio			
	
	BI.prop.ratio	<- log(num.s.prop/alpha) - dunif(U.pert, 0.5 - delPert, 0.5 + delPert, log = TRUE);

	#BI.prop.ratio	<- log(max(s)/alpha) - dunif(U.pert, 0.5 - delPert, 0.5 + delPert, log = TRUE);
	
	
	# log-Jacobian	
	
	BI.Jacob		<- log(1/U.pert/(1-U.pert))
	
	# Acceptance probability
	
	logR.BI	<- BI.lh.ratio + BI.prior.ratio + BI.prop.ratio + BI.Jacob;
	
	u.BI	= log(runif(1)) < logR.BI;
	
	if(u.BI == 1){
		ini	<- ini.new			
		accept	<- accept + u.BI;
		}
		
	list(ini = ini, accept = accept);
	
	}





move.DI.BpeScr <-
function(survObj, ini, priorPara, g){
	
	n		<- length(survObj$y1)
	p		<- ncol(survObj$x)
	
	y1		<- survObj$y1
	y2		<- survObj$y2	
	delta1	<- survObj$delta1
	delta2	<- survObj$delta2	
	x		<- survObj$x

	if(g == 1){	
		xbeta		<- ini$xbeta1.ini		
		s.propBI 	<- priorPara$s.propBI1		
		s			<- ini$s1	
		J			<- ini$J1		
		beta.ini	<- ini$beta1.ini		
		lambda		<- ini$lambda1		
		sig.sq.lam	<- ini$sig.sq.lam1		
		mu.lam		<- ini$mu.lam1
		ind.r		<- ini$ind.r1.sc1
		ind.d		<- ini$ind.d1.sc1
		ind.d.temp	<- ini$ind.d.temp.sc1		
		d			<- ini$d1.sc1		
		Delta		<- ini$Delta1
		Sigma.lam	<- ini$Sigma.lam1
		inv.Sigma.lam	<- ini$inv.Sigma.lam1
		W			<- ini$W1
		Q			<- ini$Q1
		a		<- priorPara$a1
		b		<- priorPara$b1
		c.lam	<- priorPara$c.lam1
		C		<- priorPara$C1
		alpha	<- priorPara$alpha1
		delPert <- priorPara$delPert1
		}	
		
	if(g == 2){
		xbeta		<- ini$xbeta2.ini		
		s.propBI 	<- priorPara$s.propBI2		
		s			<- ini$s2	
		J			<- ini$J2		
		beta.ini	<- ini$beta2.ini		
		lambda		<- ini$lambda2		
		sig.sq.lam	<- ini$sig.sq.lam2		
		mu.lam		<- ini$mu.lam2
		ind.r		<- ini$ind.r1.sc2
		ind.d		<- ini$ind.d2.sc2
		ind.d.temp	<- ini$ind.d.temp.sc2		
		d			<- ini$d2.sc2		
		Delta		<- ini$Delta2
		Sigma.lam	<- ini$Sigma.lam2
		inv.Sigma.lam	<- ini$inv.Sigma.lam2
		W			<- ini$W2
		Q			<- ini$Q2
		a		<- priorPara$a2
		b		<- priorPara$b2
		c.lam	<- priorPara$c.lam2
		C		<- priorPara$C2
		alpha	<- priorPara$alpha2
		delPert <- priorPara$delPert2
		}	

	if(g == 3){
		xbeta		<- ini$xbeta3.ini		
		s.propBI 	<- priorPara$s.propBI3
		s			<- ini$s3	
		J			<- ini$J3		
		beta.ini	<- ini$beta3.ini		
		lambda		<- ini$lambda3		
		sig.sq.lam	<- ini$sig.sq.lam3		
		mu.lam		<- ini$mu.lam3
		ind.r		<- ini$ind.r2.sc3
		ind.d		<- ini$ind.d3.sc3
		ind.d.temp	<- ini$ind.d.temp.sc3		
		d			<- ini$d3.sc3		
		Delta		<- ini$Delta3
		Sigma.lam	<- ini$Sigma.lam3
		inv.Sigma.lam	<- ini$inv.Sigma.lam3
		W			<- ini$W3
		Q			<- ini$Q3
		a		<- priorPara$a3
		b		<- priorPara$b3
		c.lam	<- priorPara$c.lam3
		C		<- priorPara$C3
		alpha	<- priorPara$alpha3
		delPert <- priorPara$delPert3
		}	
		
	theta.ini	<- ini$theta.ini
	gamma.ini	<- ini$gamma.ini
	accept = 0
	
	s.propBI	<- setdiff(s.propBI, s)
	num.s.prop	<- length(s.propBI)	

	# new time split and new lambda

	j.old	<- sample(c(1:J), 1);
	s.new	<- s[-j.old];
	J.new	<- length(s.new) - 1;
	
	U.pert	<- 1/(exp(lambda[j.old+1] - lambda[j.old]) + 1);
	
	if(j.old != 1){
		new.lam	<- ((s[j.old] - s[j.old-1]) * lambda[j.old] + 					(s[j.old+1] - s[j.old]) * lambda[j.old + 1]) / (s[j.old + 1] - s[j.old - 1]);
		}
		
	if(j.old == 1){
		new.lam	<- ((s[j.old] - 0) * lambda[j.old] + 					(s[j.old+1] - s[j.old]) * lambda[j.old + 1]) / (s[j.old + 1] - 0);
		}		
		
	lambda.new <- lambda[-j.old];
	lambda.new[j.old]	<- new.lam;


	ini.new	<- ini
	
	if(g == 1){
		ini.new$lambda1	<- lambda.new
		ini.new$s1		<- s.new
		ini.new$J1		<- J.new				

		old <- list()
		old$j.old <- j.old
		old$s.old <- s
		old$J.old <- J
		old$ind.r1.old <- ini$ind.r1.sc1
		old$ind.r2.old <- ini$ind.r2.sc1
		old$ind.d1.old <- ini$ind.d1.sc1
		old$ind.d2.old <- ini$ind.d2.sc1
		old$ind.d3.old <- ini$ind.d3.sc1
		old$ind.d.temp.old <- ini$ind.d.temp.sc1
		old$d1.old <- ini$d1.sc1
		old$d2.old <- ini$d2.sc1
		old$d3.old <- NULL
		
		intv.new		<- remove.interval.BpeScr(survObj, ini, priorPara, old, s.new, g = 1)
		ini.new$ind.r1.sc1	<- intv.new$ind.r1;
		ini.new$ind.d1.sc1	<- intv.new$ind.d1;
		ini.new$ind.d2.sc1	<- intv.new$ind.d2;
		ini.new$ind.d.temp.sc1	<- intv.new$ind.d.temp
		ini.new$d1.sc1		<- intv.new$d1;
		ini.new$d2.sc1		<- intv.new$d2;		
				
		ini.new$Delta1			<- cal.Del.DI.BpeScr(survObj, ini.new, old, g = 1)
		
		cal.Sig1.new			<- cal.Sigma.lam.BpeScr(ini.new, priorPara, g = 1);
		ini.new$Sigma.lam1		<- cal.Sig1.new$Sigma.lam;
		ini.new$inv.Sigma.lam1	<- solve(ini.new$Sigma.lam1)
		ini.new$W1				<- cal.Sig1.new$W;
		ini.new$Q1				<- cal.Sig1.new$Q;
		
		Delta.new	<- ini.new$Delta1
		ind.r.new	<- ini.new$ind.r1.sc1
		ind.d.new	<- ini.new$ind.d1.sc1		
		d.new		<- ini.new$d1.sc1
		Sigma.lam.new	<- ini.new$Sigma.lam1
		}	
	
	if(g == 2){
		ini.new$lambda2	<- lambda.new
		ini.new$s2		<- s.new
		ini.new$J2		<- J.new				

		old <- list()
		old$j.old <- j.old
		old$s.old <- s
		old$J.old <- J
		old$ind.r1.old <- ini$ind.r1.sc2
		old$ind.r2.old <- ini$ind.r2.sc2
		old$ind.d1.old <- ini$ind.d1.sc2
		old$ind.d2.old <- ini$ind.d2.sc2
		old$ind.d3.old <- ini$ind.d3.sc2
		old$ind.d.temp.old <- ini$ind.d.temp.sc2
		old$d1.old <- ini$d1.sc2
		old$d2.old <- ini$d2.sc2
		old$d3.old <- ini$d3.sc2
		
		intv.new		<- remove.interval.BpeScr(survObj, ini, priorPara, old, s.new, g = 2)
		ini.new$ind.r1.sc2	<- intv.new$ind.r1;
		ini.new$ind.d1.sc2	<- intv.new$ind.d1;
		ini.new$ind.d2.sc2	<- intv.new$ind.d2;
		ini.new$ind.d.temp.sc2	<- intv.new$ind.d.temp
		ini.new$d1.sc2		<- intv.new$d1;
		ini.new$d2.sc2		<- intv.new$d2;		

		ini.new$Delta2			<- cal.Del.DI.BpeScr(survObj, ini.new, old, g = 2)
		
		cal.Sig2.new			<- cal.Sigma.lam.BpeScr(ini.new, priorPara, g = 2);
		ini.new$Sigma.lam2		<- cal.Sig2.new$Sigma.lam;
		ini.new$inv.Sigma.lam2	<- solve(ini.new$Sigma.lam2)
		ini.new$W2				<- cal.Sig2.new$W;
		ini.new$Q2				<- cal.Sig2.new$Q;
		
		Delta.new	<- ini.new$Delta2
		ind.r.new	<- ini.new$ind.r1.sc2
		ind.d.new	<- ini.new$ind.d2.sc2		
		d.new		<- ini.new$d2.sc2
		Sigma.lam.new	<- ini.new$Sigma.lam2
		}	
		
	if(g == 3){
		ini.new$lambda3	<- lambda.new
		ini.new$s3		<- s.new
		ini.new$J3		<- J.new				
		
		old <- list()
		old$j.old <- j.old
		old$s.old <- s
		old$J.old <- J
		old$ind.r1.old <- ini$ind.r1.sc3
		old$ind.r2.old <- ini$ind.r2.sc3
		old$ind.d1.old <- ini$ind.d1.sc3
		old$ind.d2.old <- ini$ind.d2.sc3
		old$ind.d3.old <- ini$ind.d3.sc3
		old$ind.d.temp.old <- ini$ind.d.temp.sc3
		old$d1.old <- ini$d1.sc3
		old$d2.old <- ini$d2.sc3
		old$d3.old <- ini$d3.sc3
		
		intv.new		<- remove.interval.BpeScr(survObj, ini, priorPara, old, s.new, g = 3)
		
		ini.new$ind.r2.sc3	<- intv.new$ind.r2;
		ini.new$ind.d1.sc3	<- intv.new$ind.d1;
		ini.new$ind.d3.sc3	<- intv.new$ind.d3;
		ini.new$ind.d.temp.sc3	<- intv.new$ind.d.temp
		ini.new$d1.sc3		<- intv.new$d1;
		ini.new$d3.sc3		<- intv.new$d3;	
		
		ini.new$Delta3			<- cal.Del.DI.BpeScr(survObj, ini.new, old, g = 3)
		
		cal.Sig3.new			<- cal.Sigma.lam.BpeScr(ini.new, priorPara, g = 3);
		ini.new$Sigma.lam3		<- cal.Sig3.new$Sigma.lam;
		ini.new$inv.Sigma.lam3	<- solve(ini.new$Sigma.lam3)
		ini.new$W3				<- cal.Sig3.new$W;
		ini.new$Q3				<- cal.Sig3.new$Q;
		
		Delta.new	<- ini.new$Delta3
		ind.r.new	<- ini.new$ind.r2.sc3
		ind.d.new	<- ini.new$ind.d3.sc3		
		d.new		<- ini.new$d3.sc3
		Sigma.lam.new	<- ini.new$Sigma.lam3
		}

	# log-likelhood ratio
	
		xbeta.dum.ini		<- xbeta;
		xbeta.dum.ini[xbeta.dum.ini > 700] <- 700;
		exp.xbeta.ini		<- exp(xbeta.dum.ini);
			
		loglh.ini.1st	<- -exp(lambda[c(j.old, j.old+1)]) * colSums(Delta[,c(j.old, j.old+1)] * 							gamma.ini * exp.xbeta.ini * ind.r[,c(j.old, j.old+1)])
		loglh.ini.2nd	<- lambda[c(j.old, j.old+1)] * d[c(j.old, j.old+1)] 
		loglh.ini		<- sum(loglh.ini.1st + loglh.ini.2nd)
			
		loglh.new.1st	<- -exp(lambda.new[j.old])*sum(Delta.new[,c(j.old)]*gamma.ini * exp.xbeta.ini * ind.r.new[,j.old])
		loglh.new.2nd	<- lambda.new[j.old] * d.new[j.old]
		loglh.new		<- loglh.new.1st + loglh.new.2nd
			
		DI.lh.ratio	<- loglh.new - loglh.ini
				
	
	# log-prior ratio
		
	if(J+1 != 2){
		if(j.old != 1){
			DI.prior.ratio	<- log((J)/alpha) + dmvnorm(lambda.new, mean = rep(mu.lam, J.new + 1), 								sigma = sig.sq.lam * Sigma.lam.new, log = TRUE) -							 	dmvnorm(lambda, mean = rep(mu.lam, J+1), sigma = sig.sq.lam * Sigma.lam, log = TRUE) + 								log(1/(2*J+1)/(2*J)*s[J+1]^2*(s[j.old+1] - s[j.old-1])/(s[j.old] - s[j.old - 1])/								(s[j.old+1] - s[j.old]));
			}
		
		if(j.old == 1){
			DI.prior.ratio	<- log((J)/alpha) + dmvnorm(lambda.new, mean = rep(mu.lam, J.new + 1), 								sigma = sig.sq.lam * Sigma.lam.new, log = TRUE) -								dmvnorm(lambda, mean = rep(mu.lam, J+1), sigma = sig.sq.lam * Sigma.lam, log = TRUE) + 								log(1/(2*J+1)/(2*J)*s[J+1]^2*(s[j.old+1] - 0)/(s[j.old] - 0)/(s[j.old+1] - s[j.old]));
			}		
		}	

	if(J+1 == 2){
		DI.prior.ratio	<- log((J)/alpha) + dnorm(lambda.new, mean = mu.lam, 							sd = sqrt(sig.sq.lam * Sigma.lam.new), log = TRUE) - 							dmvnorm(lambda, mean = rep(mu.lam, J+1), sigma = sig.sq.lam * Sigma.lam, log = TRUE) + 							log(1/(2*J+1)/(2*J)*s[J+1]^2*(s[j.old+1] - 0)/(s[j.old] - 0)/(s[j.old+1] - s[j.old]));
		}

				
	# log-proposal ratio			
	
	DI.prop.ratio	<- log(alpha/num.s.prop) + dunif(0.5, 0.5 - delPert, 0.5 + delPert, log = TRUE);
	
	#DI.prop.ratio	<- log(alpha/max(s)) + dunif(0.5, 0.5 - delPert, 0.5 + delPert, log = TRUE);
	
	
	# log-Jacobian	
	
	DI.Jacob		<- log(U.pert*(1-U.pert));
		
		
	# Acceptance probability
	
	logR.DI	<- DI.lh.ratio + DI.prior.ratio + DI.prop.ratio + DI.Jacob;
	
	u.DI	= log(runif(1)) < logR.DI;
	
	if(u.DI == 1){
		ini	<- ini.new			
		accept	<- accept + u.DI
		}
	
	list(ini = ini, accept = accept);
		
	}





move.DP.XiGam.BpeScr <-
function(survObj, ini, priorPara){

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





move.FP.BpeScr <-
function(survObj, ini){
	
		n <- survObj$n
		p <- survObj$p
			
		s1			<- ini$s1
		s2			<- ini$s2
		s3			<- ini$s3		
		J1			<- ini$J1
		J2			<- ini$J2
		J3			<- ini$J3		
		beta1.ini	<- ini$beta1.ini
		beta2.ini	<- ini$beta2.ini
		beta3.ini	<- ini$beta3.ini		
		lambda1		<- ini$lambda1
		lambda2		<- ini$lambda2
		lambda3		<- ini$lambda3		
		sig.sq.lam1	<- ini$sig.sq.lam1
		sig.sq.lam2	<- ini$sig.sq.lam2
		sig.sq.lam3	<- ini$sig.sq.lam3		
		mu.lam1		<- ini$mu.lam1
		mu.lam2		<- ini$mu.lam2
		mu.lam3		<- ini$mu.lam3	
		theta.ini	<- ini$theta.ini
		gamma.ini	<- ini$gamma.ini
		Delta1		<- ini$Delta1
		Delta2		<- ini$Delta2
		Delta3		<- ini$Delta3
		xbeta1		<- ini$xbeta1.ini
		xbeta2		<- ini$xbeta2.ini
		xbeta3		<- ini$xbeta3.ini		
				
		case00		<- survObj$case00
		case10		<- survObj$case10
		case01		<- survObj$case01
		case11		<- survObj$case11
								
	xbeta1.mat		<- matrix(rep(xbeta1, J1+1), n, J1+1);
	xbeta1.dum		<- xbeta1;
	xbeta1.dum[xbeta1.dum > 700] <- 700;	
	exp.xbeta1		<- exp(xbeta1.dum);
	exp.xbeta1.mat	<- matrix(rep(exp.xbeta1, J1+1), n, J1+1);
	
	xbeta2.mat		<- matrix(rep(xbeta2, J2+1), n, J2+1);
	xbeta2.dum		<- xbeta2;
	xbeta2.dum[xbeta2.dum > 700] <- 700;	
	exp.xbeta2		<- exp(xbeta2.dum);
	exp.xbeta2.mat	<- matrix(rep(exp.xbeta2, J2+1), n, J2+1);
	
	xbeta3.mat		<- matrix(rep(xbeta3, J3+1), n, J3+1);
	xbeta3.dum		<- xbeta3;
	xbeta3.dum[xbeta3.dum > 700] <- 700;	
	exp.xbeta3		<- exp(xbeta3.dum);
	exp.xbeta3.mat	<- matrix(rep(exp.xbeta3, J3+1), n, J3+1);
	
	if(length(case00) != 0){		
		gam.rate00.1st	<- matrix(rep(exp(lambda1), length(case00)), ncol = J1+1, byrow = TRUE) * Delta1[case00,] * exp.xbeta1.mat[case00,];
		gam.rate00.2nd	<- matrix(rep(exp(lambda2), length(case00)), ncol = J2+1, byrow = TRUE) * Delta2[case00,] * exp.xbeta2.mat[case00,];
		gam.rate00	<- rowSums(gam.rate00.1st) + rowSums(gam.rate00.2nd) + 1/theta.ini;
		gamma.ini[case00]	<- rgamma(length(case00), shape = 1/theta.ini, rate = gam.rate00);
		}

	if(length(case10) != 0){
		gam.rate10.1st	<- matrix(rep(exp(lambda1), length(case10)), ncol = J1+1, byrow = TRUE) * Delta1[case10,] * exp.xbeta1.mat[case10,];
		gam.rate10.2nd	<- matrix(rep(exp(lambda2), length(case10)), ncol = J2+1, byrow = TRUE) * Delta2[case10,] * exp.xbeta2.mat[case10,];
		gam.rate10.3rd	<- matrix(rep(exp(lambda3), length(case10)), ncol = J3+1, byrow = TRUE) * Delta3[case10,] * exp.xbeta3.mat[case10,];
		gam.rate10	<- rowSums(gam.rate10.1st) + rowSums(gam.rate10.2nd) + rowSums(gam.rate10.3rd) + 1/theta.ini;
		gamma.ini[case10]	<- rgamma(length(case10), shape = 1/theta.ini + 1, rate = gam.rate10);
		}
		
	if(length(case01) != 0){
		gam.rate01.1st	<- matrix(rep(exp(lambda1), length(case01)), ncol = J1+1, byrow = TRUE) * Delta1[case01,] * exp.xbeta1.mat[case01,];
		gam.rate01.2nd	<- matrix(rep(exp(lambda2), length(case01)), ncol = J2+1, byrow = TRUE) * Delta2[case01,] * exp.xbeta2.mat[case01,];
		gam.rate01	<- rowSums(gam.rate01.1st) + rowSums(gam.rate01.2nd) + 1/theta.ini;
		gamma.ini[case01]	<- rgamma(length(case01), shape = 1/theta.ini + 1, rate = gam.rate01); 
		}
	
	if(length(case11) != 0){
		gam.rate11.1st	<- matrix(rep(exp(lambda1), length(case11)), ncol = J1+1, byrow = TRUE) * Delta1[case11,] * exp.xbeta1.mat[case11,];
		gam.rate11.2nd	<- matrix(rep(exp(lambda2), length(case11)), ncol = J2+1, byrow = TRUE) * Delta2[case11,] * exp.xbeta2.mat[case11,];
		gam.rate11.3rd	<- matrix(rep(exp(lambda3), length(case11)), ncol = J3+1, byrow = TRUE) * Delta3[case11,] * exp.xbeta3.mat[case11,];
		gam.rate11	<- rowSums(gam.rate11.1st) + rowSums(gam.rate11.2nd) + rowSums(gam.rate11.3rd) + 1/theta.ini;
		gamma.ini[case11]	<- rgamma(length(case11), shape = 1/theta.ini + 2, rate = gam.rate11);
		}

	return(gamma.ini);
	
	}

