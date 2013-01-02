
BpeSurv	<- function(survObj, priorPara, initial, num.reps, thin = 1, chain = 1, save=1000, RJ = TRUE){
	
	ini	<- list()

	y		<- survObj$y
	delta	<- survObj$delta
	x		<- survObj$x
	
	n <- survObj$n	<- nrow(survObj$x)
	p <- survObj$p	<- ncol(survObj$x)
	
	smax	<- priorPara$smax
	
	survObj$case0 <- which(delta == 0)
	survObj$case1 <- which(delta == 1)	
	
	survObj$case0yleq <- which(delta == 0 & y <= smax)
	survObj$case0ygeq <- which(delta == 0 & y > smax)
	survObj$case1yleq <- which(delta == 1 & y <= smax)
	survObj$case1ygeq <- which(delta == 1 & y > smax)
			
	ini$s			<- initial$s
	ini$J			<- initial$J
	ini$beta.ini	<- initial$beta.ini
	ini$lambda		<- initial$lambda
	ini$sig.sq.lam	<- initial$sig.sq.lam
	ini$mu.lam		<- initial$mu.lam
	
	ini$xbeta.ini	<- as.vector(x%*%ini$beta.ini)

	c.lam	<- priorPara$c.lam
	C		<- priorPara$C
	a		<- priorPara$a
	b		<- priorPara$b
	alpha	<- priorPara$alpha
	psi		<- priorPara$psi
	omega	<- priorPara$omega
	delPert <- priorPara$delPert		

	intv	<- setting.interval.BpeSurv(survObj, ini, priorPara);
	ini$ind.d 		<- intv$ind.d;
	ini$ind.d.temp 	<- intv$ind.d.temp	 
	ini$ind.r 		<- intv$ind.r; 
	ini$d 			<- intv$d;		
	
	ini$Delta	<- calDelta.BpeSurv(survObj, ini)
		
	cal.Sig			<- cal.Sigma.lam.BpeSurv(ini, priorPara)
	
	ini$Sigma.lam 		<- cal.Sig$Sigma.lam;
	ini$inv.Sigma.lam	<- solve(ini$Sigma.lam);
	ini$W				<- cal.Sig$W;	
	ini$Q				<- cal.Sig$Q;
	
	
	# for posterior samples
	
	mcmcOutcome	<- list()
	
	mcmcOutcome$priorPara	<- priorPara
	mcmcOutcome$initial		<- initial	
	mcmcOutcome$num.reps	<- num.reps
	mcmcOutcome$thin		<- thin
	mcmcOutcome$save		<- save
	mcmcOutcome$RJ			<- RJ
	mcmcOutcome$chain		<- chain	
	
	mcmcOutcome$move.ind		<- NULL;
	
	mcmcOutcome$beta.p			<- initial$beta.ini
	mcmcOutcome$lambda.p		<- as.list(c(rep(NA, num.reps)));
	mcmcOutcome$s.p				<- as.list(c(rep(NA, num.reps)));
	mcmcOutcome$J.p				<- initial$J;
	mcmcOutcome$s.p[[1]]		<- initial$s;
	mcmcOutcome$lambda.p[[1]]	<- initial$lambda;

	mcmcOutcome$log.jpost 		<- NULL
	mcmcOutcome$log.like		<- NULL
	mcmcOutcome$mu.lam.p		<- NULL 
	mcmcOutcome$sig.sq.lam.p	<- NULL
	mcmcOutcome$accept.beta		<- c(rep(0, p))
	mcmcOutcome$accept.lambda	<- as.list(c(rep(0, num.reps/thin)))
	mcmcOutcome$accept.bi		<- 0
	mcmcOutcome$accept.di		<- 0


	# for move probability
	
	Jmax 				<- length(priorPara$s.propBI) - 1
		
	pB   <- pmin(1, alpha/c(1:(Jmax)+1))
	pD   <- pmin(1, (c(1:Jmax))/alpha)
	rho.lam  <- min(C / (pB + pD))


	dir.create('mcmcOutcome')		

	### Start MCMC sampling
	

	for(M in 1:num.reps){
	
		if(M %% 100 == 0){
			cat("Chain :", chain, "Scan :", M,fill=TRUE)
			fsh()
			} 
		
		if(RJ){
			if(ini$J < Jmax){	
				p.BI	<- rho.lam * min(1, alpha/(ini$J+1))
				p.DI	<- rho.lam * min(1, (ini$J)/alpha)
				}
			
			if(ini$J >= Jmax){
				p.BI	<- 0
				p.DI	<- rho.lam * 2
				}
	
			p.RP	<- p.BH	<- p.SP	<- (1 - p.BI - p.DI) / 3
			move.inx <- sample(c("RP", "BH", "SP", "BI", "DI"), size = 1, 						prob = c(p.RP, p.BH, p.SP, p.BI, p.DI));
			}
			
		if(!RJ){
			p.BI	<- 0
			p.DI	<- 0
			p.RP	<- p.BH	<- p.SP	<- (1 - p.BI - p.DI) / 3
			move.inx <- sample(c("RP", "BH", "SP"), size = 1, 								prob = c(p.RP, p.BH, p.SP));
			}		

		mcmcOutcome$move.ind	<- c(mcmcOutcome$move.ind, move.inx);


		# Updating Regression coefficients
			
		if(move.inx == "RP"){		# update beta
			sample.beta			<- move.RP.BpeSurv(survObj, ini);
			ini$beta.ini		<- sample.beta$beta.ini;
			ini$xbeta.ini		<- sample.beta$xbeta.ini;
			mcmcOutcome$accept.beta		<- mcmcOutcome$accept.beta + sample.beta$accept;
			}		
				
	
		# update the log-baseline hazard	
					
		if(move.inx == "BH"){		# update lambda
			mcmcOutcome$accept.lambda[[M]]	<- rep(0, ini$J+1)
			sample.lambda		<- move.BH.BpeSurv(survObj, ini);
			ini$lambda			<- sample.lambda$lambda;
			mcmcOutcome$accept.lambda[[M]]	<- sample.lambda$accept;
			}


		# Updating second stage survival components: mu.lam and sigma.sq.lam
	
		if(move.inx == "SP"){		# update mu.lam and sigma.sq.lam
			sample.sp				<- move.SP.BpeSurv(survObj, ini, priorPara);
			ini$mu.lam				<- sample.sp$mu.lam.sample;
			ini$sig.sq.lam			<- sample.sp$sig.sq.lam.sample;
			}

		# Creating a new time split: Birth move	

		if(move.inx == "BI"){ 			# for s and lambda
			sample.bi	<- move.BI.BpeSurv(survObj, ini, priorPara);
			if(sample.bi$accept == 1){
				ini	<- sample.bi$ini			
				mcmcOutcome$accept.bi	<- mcmcOutcome$accept.bi + sample.bi$accept;
				}
			}

						
		# Removing a time split: Reverse move
	
		if(move.inx == "DI"){ 			# for s and lambda
			sample.di	<- move.DI.BpeSurv(survObj, ini, priorPara);
			if(sample.di$accept == 1){
				ini	<- sample.di$ini			
				mcmcOutcome$accept.di		<- mcmcOutcome$accept.di + sample.di$accept;
				}
			}			

	

	if((M %% thin) == 0){
				
      # storage of posterior samples from MCMC	
      
		mcmcOutcome$ini	<- ini
      		
		mcmcOutcome$beta.p					<- rbind(mcmcOutcome$beta.p, ini$beta.ini, deparse.level = 0)
		mcmcOutcome$lambda.p[[M/thin+1]] 	<- ini$lambda;
		mcmcOutcome$J.p						<- c(mcmcOutcome$J.p, ini$J)
		mcmcOutcome$s.p[[M/thin+1]]			<- ini$s
		mcmcOutcome$mu.lam.p		<- c(mcmcOutcome$mu.lam.p, ini$mu.lam)
		mcmcOutcome$sig.sq.lam.p	<- c(mcmcOutcome$sig.sq.lam.p, ini$sig.sq.lam)
		}
	

	# for monitoring the posterior samples	

	if(M %% save == 0 | M == num.reps){
		
		betaMed		<- apply(mcmcOutcome$beta.p[(M/thin/2):(M/thin),], 2, mean, na.rm = T)
		beta0.975	<- apply(mcmcOutcome$beta.p[(M/thin/2):(M/thin),], 2, quantile, prob = 0.975, na.rm = T)
		beta0.025	<- apply(mcmcOutcome$beta.p[(M/thin/2):(M/thin),], 2, quantile, prob = 0.025, na.rm = T)
		
				
		#if(M <= 100000){

		#	pdf(file = paste("mcmcOutcome/profJpost.ch",chain,".pdf", sep = ""), height = 10, width = 20)
	
		#	par(mfrow = c(1,2))

		#	plot(mcmcOutcome$log.jpost, type = "l", main = "log(joint posterior)", xlab = "scan", ylab = "log.joint.posterior"); 
		#	plot(mcmcOutcome$log.like, type = "l", main = "log(likelihood)", xlab = "scan", ylab = "log.likelihood");

		#	dev.off()
		
		#	}
		
		
		# for beta)


		pdf(file = paste("mcmcOutcome/ch",chain,"regPara.pdf", sep = ""), height = 10, width = 20)
		

		par(mfrow = c(3,6))		
		
		for(i in 1:p){
			plot(mcmcOutcome$beta.p[(M/thin/2):(M/thin),i], type = "l", 			main = colnames(x)[i], xlab = "scan", ylab = "reg. coef.");
#			abline(h = coef(fit.cox)[i], col = "blue", lty = 4)
			abline(h = betaMed[i], col = "red")
			abline(h = beta0.975[i], lty = 2, col = "red")
			abline(h = beta0.025[i], lty = 2, col = "red")
			}
		
		dev.off()



		pdf(file = paste("mcmcOutcome/otherPara.ch",chain,".pdf", sep = ""), height = 10, width = 30)

		
		par(mfrow = c(1, 3))
		plot(log(mcmcOutcome$sig.sq.lam.p[(length(mcmcOutcome$sig.sq.lam.p)/2):length(mcmcOutcome$sig.sq.lam.p)]), 			type = "l", main = "log(sigmaL)", xlab = "scan", ylab = "log.sigmaL");
		plot(mcmcOutcome$mu.lam.p[(length(mcmcOutcome$mu.lam.p)/2):length(mcmcOutcome$mu.lam.p)], 				type = "l", main = "muL", xlab = "scan", ylab = "muL");
		plot(mcmcOutcome$J.p[(M/thin/2):(M/thin)], type = "l", main = "#intervals", xlab = "scan", ylab = "J");

		dev.off()
		
		}



	if(M %% save == 0 | M == num.reps){


		ss	<- seq(1, max(initial$s) , 1)
		JJ	<- length(ss) - 1;

		lambda.fin	<- matrix(NA, M/thin/2, JJ+1);
		ind.lam		<- rep(NA, JJ+1);


		for(i in (M/thin/2 + 1):(M/thin)){
	
			ss.n	<- length(mcmcOutcome$s.p[[i]]);
			ind.lam[(ss <= mcmcOutcome$s.p[[i]][1]) & (ss > 0)] <- 1;
	
			for(j in 2:ss.n){
				ind.lam[(ss <= mcmcOutcome$s.p[[i]][j]) & (ss > mcmcOutcome$s.p[[i]][j-1])] <- j;
				}

			lambda.fin[i - (M/thin/2), ] <- mcmcOutcome$lambda.p[[i]][ind.lam];
	
			}

		lambda.mean		<- apply(lambda.fin, 2, mean, na.rm = T);
		lambda.sd		<- apply(lambda.fin, 2, sd, na.rm = T);
		lambda.med		<- apply(lambda.fin, 2, median, na.rm = T);
		lambda0.975		<- apply(lambda.fin, 2, quantile, prob = 0.975, na.rm = T)
		lambda0.025		<- apply(lambda.fin, 2, quantile, prob = 0.025, na.rm = T)
		lb	<- lambda.mean - 1.96 * lambda.sd
		ub	<- lambda.mean + 1.96 * lambda.sd
		



		pdf(file = paste("mcmcOutcome/ch", chain, "baseHaz.pdf", sep = ""), height = 10, width = 20)

		
		par(mfrow = c(3, 3))


		plot(ss, lambda0.975, type = "s", col = "blue", lty = 2, ylim = c(min(lambda0.025), max(lambda0.975)), 			main = "Estimate of the log-baseline hazard", xlab = "time", ylab = "log-baseline hazard")
		lines(ss, lambda.mean, type = "s", col = "red")
		lines(ss, lb, type = "s", col = "red", lty = 2)
		lines(ss, ub, type = "s", col = "red", lty = 2);
		lines(ss, lambda.med, type = "s", col = "blue")
		lines(ss, lambda0.025, type = "s", col = "blue", lty = 2)
		lines(ss, lambda0.975, type = "s", col = "blue", lty = 2)
		# lines(supsmu(hazt, log(haz)), type = "l")
		#lines(c(0, sSim), c(logBH1, logBH1[length(logBH1)]), type = "s", col = "green")

		for(i in 1:8){
			plot(lambda.fin[,floor(max(initial$s)*i/8)], type = "l", 				main = paste("Log-baseline hazard at t = ", floor(max(initial$s)*i/8)))
			}
		dev.off()
		
		}


	# store the mcmc outcomes

	if(M %% save == 0 | M == num.reps){	
			if(!RJ){
				save(mcmcOutcome, file = paste("mcmcOutcome/BpeSurvIntv",initial$J + 1,"clam",c.lam,					"a", a, "b", b, "ch", chain, ".RData", sep = ""))
				}
			if(RJ){
				save(mcmcOutcome, file = paste("mcmcOutcome/BpeSurvRjAlp",alpha, "pert", delPert, "clam",c.lam,					"a", a, "b", b, "ch", chain, ".RData", sep = ""))
				}				
			}
	
	
	}  # end of MCMC sampling

	return(mcmcOutcome)

	
	} # end of BpeSurv function




