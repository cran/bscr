		

BweibSurv	<- function(survObj, priorPara, initial, num.reps, thin = 1, chain = 1, save = 1000){


			
	ini	<- list()

	y		<- survObj$y
	delta	<- survObj$delta
	x		<- survObj$x
	
	n <- survObj$n	<- nrow(survObj$x)
	p <- survObj$p	<- ncol(survObj$x)
	
	ini$beta.ini	<- initial$beta.ini		
	ini$alpha.ini	<- initial$alpha.ini
	ini$kappa.ini	<- initial$kappa.ini
	
	ini$xbeta.ini	<- as.vector(x%*%ini$beta.ini)
							
	a		<- priorPara$a
	b		<- priorPara$b
	c		<- priorPara$c
	d		<- priorPara$d
								
	
	# for posterior samples
	
	mcmcOutcome	<- list()
	
	mcmcOutcome$priorPara	<- priorPara
	mcmcOutcome$initial		<- initial
	mcmcOutcome$num.reps	<- num.reps
	mcmcOutcome$thin		<- thin
	mcmcOutcome$save		<- save
	mcmcOutcome$chain		<- chain		
		
	mcmcOutcome$move.ind	<- NULL;
		
	mcmcOutcome$beta.p			<- initial$beta.ini

	mcmcOutcome$alpha.p 		<- initial$alpha
				
	mcmcOutcome$kappa.p 		<- initial$kappa
		
	mcmcOutcome$log.jpost 		<- log.like <- NULL;
					
	mcmcOutcome$accept.beta		<- c(rep(0, p))
	
	mcmcOutcome$accept.alpha		<- 0
		
	
	dir.create('mcmcOutcome')
	
		
	### Start MCMC sampling
	

	for(M in 1:num.reps){
		
		cat("Chain :", chain, "Scan :", M,fill=TRUE);
		fsh();
				
		
		# Updating Regression coefficients

			sample.beta		<- move.RP.BweibSurv(survObj, ini);
			ini$beta.ini		<- sample.beta$beta.ini;
			ini$xbeta.ini		<- sample.beta$xbeta;
						
			mcmcOutcome$accept.beta		<- mcmcOutcome$accept.beta + sample.beta$accept;

						
		# Updating Weibull scale parameter : alpha
		
			sample.alpha		<- move.WSC.BweibSurv(survObj, ini, priorPara);
			ini$alpha.ini		<- sample.alpha$alpha.ini;
			mcmcOutcome$accept.alpha		<- mcmcOutcome$accept.alpha + sample.alpha$accept;
						

						
			
			
			
		# Updating Weibull shape parameter : kappa
	
			sample.kappa		<- move.WSH.BweibSurv(survObj, ini, priorPara);
			ini$kappa.ini		<- sample.kappa

			


		

	# Profiling the log-likelihood and log-joint posterior

	if((M %% thin) == 0){
		
      
      # storage of posterior samples from MCMC	
      
		mcmcOutcome$ini	<- ini      
      
      	mcmcOutcome$beta.p	<- rbind(mcmcOutcome$beta.p, ini$beta.ini, deparse.level = 0)
			
		mcmcOutcome$alpha.p			<- c(mcmcOutcome$alpha.p, ini$alpha.ini);

		mcmcOutcome$kappa.p			<- c(mcmcOutcome$kappa.p, ini$kappa.ini);
								
		}	
	

		

		# for monitoring the posterior samples	
		
	
	
	if(M %% save == 0 | M == num.reps){
		
		beta.Med		<- apply(mcmcOutcome$beta.p[(M/thin/2):(M/thin),], 2, mean, na.rm = T)
		beta.0.975		<- apply(mcmcOutcome$beta.p[(M/thin/2):(M/thin),], 2, quantile, prob = 0.975, na.rm = T)
		beta.0.025		<- apply(mcmcOutcome$beta.p[(M/thin/2):(M/thin),], 2, quantile, prob = 0.025, na.rm = T)

		
		# for beta

		pdf(file = paste("mcmcOutcome/ch",chain,"regPara.pdf", sep = ""), height = 10, width = 20)

		par(mfrow = c(3,6))		
		
		for(i in 1:p){
			plot(mcmcOutcome$beta.p[(M/thin/2):(M/thin),i], type = "l", 			main = colnames(x)[i], xlab = "scan", ylab = "reg. coef.");
#			abline(h = coef(fit.cox)[i], col = "blue", lty = 4)
			abline(h = beta.Med[i], col = "red")
			abline(h = beta.0.975[i], lty = 2, col = "red")
			abline(h = beta.0.025[i], lty = 2, col = "red")
			}
		
		dev.off()
		



		# hazard function
		
		 # alpha, kappa, baseline hazard

		alpha.med	<- median(mcmcOutcome$alpha.p[(length(mcmcOutcome$alpha.p)/2):				length(mcmcOutcome$alpha.p)])
		kappa.med	<-  median(mcmcOutcome$kappa.p[(length(mcmcOutcome$kappa.p)/2):				length(mcmcOutcome$kappa.p)])
				
		pdf(file = paste("mcmcOutcome/ch", chain, "basehaz.pdf", sep = ""), height = 10, width = 30)
		
		par(mfrow = c(1, 3))
		
		plot(mcmcOutcome$alpha.p[(length(mcmcOutcome$alpha.p)/2):length(mcmcOutcome$alpha.p)], 			type = "l", main = expression(alpha), xlab = "scan", ylab = expression(alpha));
		abline(h = alpha.med, col = "red")
							
		plot(mcmcOutcome$kappa.p[(length(mcmcOutcome$kappa.p)/2):length(mcmcOutcome$kappa.p)], 			type = "l", main = expression(kappa), xlab = "scan", ylab = expression(kappa));
		abline(h = kappa.med, col = "red")
				
		plot(0:max(y[delta == 1]),  log(alpha.med*kappa.med* c(0:max(y[delta == 1]))^(alpha.med - 1)), type = "l", 				main = "log-baseline hazard fucntion", xlab = "scan", ylab = "log-baselne hazard", col = "red");
#		lines(supsmu(hazt, log(haz)), type = "l")

				
		dev.off()
		
		}
		




		# store the mcmc outcomes

		if(M %% save == 0 | M == num.reps){
			
			save(mcmcOutcome, file = paste("mcmcOutcome/BwbUni.a", a, "b", b, "c", c, "d", d, "ch", chain, ".RData", sep = ""))
		
			}
	
		} # end of the 'for' loop for MCMC sampling
	

	return(mcmcOutcome)
	
	} # end of BpeScr function
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		
	