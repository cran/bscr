		

BweibScr	<- function(survObj, priorPara, initial, num.reps, thin = 1, chain = 1, save = 1000){
		
	ini	<- list()

	y1		<- survObj$y1
	y2		<- survObj$y2	
	delta1	<- survObj$delta1
	delta2	<- survObj$delta2	
	x		<- survObj$x
	
	n <- survObj$n	<- nrow(survObj$x)
	p <- survObj$p	<- ncol(survObj$x)

	survObj$case1_	<- which(delta1 == 1);
	survObj$case00	<- which(delta1 == 0 & delta2 == 0);
	survObj$case10	<- which(delta1 == 1 & delta2 == 0);
	survObj$case01	<- which(delta1 == 0 & delta2 == 1);
	survObj$case11	<- which(delta1 == 1 & delta2 == 1);

	
	ini$beta1.ini	<- initial$beta1.ini
	ini$beta2.ini	<- initial$beta2.ini
	ini$beta3.ini	<- initial$beta3.ini		
	ini$alpha1.ini	<- initial$alpha1.ini
	ini$alpha2.ini	<- initial$alpha2.ini
	ini$alpha3.ini	<- initial$alpha3.ini
	ini$kappa1.ini	<- initial$kappa1.ini
	ini$kappa2.ini	<- initial$kappa2.ini
	ini$kappa3.ini	<- initial$kappa3.ini
	ini$theta.ini	<- initial$theta.ini
	ini$xi.ini		<- 1/initial$theta.ini
	ini$gamma.ini	<- initial$gamma.ini	

	ini$xbeta1.ini	<- as.vector(x%*%ini$beta1.ini)
	ini$xbeta2.ini	<- as.vector(x%*%ini$beta2.ini)
	ini$xbeta3.ini	<- as.vector(x%*%ini$beta3.ini)
								
	a1		<- priorPara$a1
	a2		<- priorPara$a2
	a3		<- priorPara$a3		
	b1		<- priorPara$b1
	b2		<- priorPara$b2
	b3		<- priorPara$b3
	c1		<- priorPara$c1
	c2		<- priorPara$c2
	c3		<- priorPara$c3		
	d1		<- priorPara$d1
	d2		<- priorPara$d2
	d3		<- priorPara$d3			
	psi		<- priorPara$psi
	omega	<- priorPara$omega

								
	
	# for posterior samples
	
	mcmcOutcome	<- list()
	
	mcmcOutcome$priorPara	<- priorPara
	mcmcOutcome$initial		<- initial
	mcmcOutcome$num.reps	<- num.reps
	mcmcOutcome$thin		<- thin
	mcmcOutcome$save		<- save
	mcmcOutcome$chain		<- chain		
		
	mcmcOutcome$move.ind	<- NULL;
		
	mcmcOutcome$beta1.p			<- initial$beta1.ini
	mcmcOutcome$beta2.p			<- initial$beta2.ini
	mcmcOutcome$beta3.p			<- initial$beta3.ini

	mcmcOutcome$alpha1.p 		<- initial$alpha1
	mcmcOutcome$alpha2.p 		<- initial$alpha2
	mcmcOutcome$alpha3.p 		<- initial$alpha3
				
	mcmcOutcome$kappa1.p 		<- initial$kappa1
	mcmcOutcome$kappa2.p 		<- initial$kappa2
	mcmcOutcome$kappa3.p 		<- initial$kappa3
			
	mcmcOutcome$gamma.p			<- initial$gamma.ini[1:10]
	mcmcOutcome$theta.p			<- initial$theta.ini
	mcmcOutcome$xi.p			<- initial$xi.ini
		
	mcmcOutcome$log.jpost 		<- log.like <- NULL;
					
	mcmcOutcome$accept.beta1		<- c(rep(0, p))
	mcmcOutcome$accept.beta2		<- c(rep(0, p))
	mcmcOutcome$accept.beta3		<- c(rep(0, p))
	
	mcmcOutcome$accept.alpha1		<- 0
	mcmcOutcome$accept.alpha2		<- 0
	mcmcOutcome$accept.alpha3		<- 0
						
	mcmcOutcome$accept.theta	<- 0
	mcmcOutcome$accept.xi		<- 0		
	
	dir.create('mcmcOutcome')
	
		
	### Start MCMC sampling
	

	for(M in 1:num.reps){
		
		cat("Chain :", chain, "Scan :", M,fill=TRUE);
		fsh();
		

				
		# move probability
		
		p.DP <- p.WSC1 <- p.WSC2 <- p.WSC3 <- 0.10
		p.RP1 <- p.RP2 <- p.RP3 <- p.WSH1 <- p.WSH2 <- p.WSH3 <- p.FP <- (1 - p.DP - p.WSC1 - p.WSC2 - p.WSC3)/7
						 
	
		move.inx	<- sample(c("RP1", "RP2", "RP3", "WSH1", "WSH2", "WSH3", "WSC1", "WSC2", "WSC3", "FP", "DP"),							 size=1, prob=c(p.RP1, p.RP2, p.RP3, p.WSH1, p.WSH2, p.WSH3 , p.WSC1, p.WSC2, p.WSC3, p.FP, p.DP));

			
		mcmcOutcome$move.ind	<- c(mcmcOutcome$move.ind, move.inx);

		
		
		# Updating Regression coefficients

		if(move.inx == "RP1"){		# update beta1
			sample.beta1		<- move.RP.BweibScr(survObj, ini, g = 1);
			ini$beta1.ini		<- sample.beta1$beta.ini;
			ini$xbeta1.ini		<- sample.beta1$xbeta;
						
			mcmcOutcome$accept.beta1		<- mcmcOutcome$accept.beta1 + sample.beta1$accept;
			}

		if(move.inx == "RP2"){		# update beta2
			sample.beta2		<- move.RP.BweibScr(survObj, ini, g = 2);
			ini$beta2.ini		<- sample.beta2$beta.ini;
			ini$xbeta2.ini		<- sample.beta2$xbeta;
						
			mcmcOutcome$accept.beta2		<- mcmcOutcome$accept.beta2 + sample.beta2$accept;
			}
			
		if(move.inx == "RP3"){		# update beta3
			sample.beta3		<- move.RP.BweibScr(survObj, ini, g = 3);
			ini$beta3.ini		<- sample.beta3$beta.ini;
			ini$xbeta3.ini		<- sample.beta3$xbeta;
			mcmcOutcome$accept.beta3		<- mcmcOutcome$accept.beta3 + sample.beta3$accept;
			}
			
						
		# Updating Weibull scale parameter : alpha
		
		if(move.inx == "WSC1"){		# update alpha1
			sample.alpha1		<- move.WSC.BweibScr(survObj, ini, priorPara, g = 1);
			ini$alpha1.ini		<- sample.alpha1$alpha.ini;
			mcmcOutcome$accept.alpha1		<- mcmcOutcome$accept.alpha1 + sample.alpha1$accept;
			}			
			
		if(move.inx == "WSC2"){		# update alpha2
			sample.alpha2		<- move.WSC.BweibScr(survObj, ini, priorPara, g = 2);
			ini$alpha2.ini		<- sample.alpha2$alpha.ini;
			mcmcOutcome$accept.alpha2		<- mcmcOutcome$accept.alpha2 + sample.alpha2$accept;
			}			
			
		if(move.inx == "WSC3"){		# update alpha3
			sample.alpha3		<- move.WSC.BweibScr(survObj, ini, priorPara, g = 3);
			ini$alpha3.ini		<- sample.alpha3$alpha.ini;
			mcmcOutcome$accept.alpha3		<- mcmcOutcome$accept.alpha3 + sample.alpha3$accept;
			}
						
			
			
			
		# Updating Weibull shape parameter : kappa
	
		if(move.inx == "WSH1"){		# update kappa1
			sample.kappa1		<- move.WSH.BweibScr(survObj, ini, priorPara, g = 1);
			ini$kappa1.ini		<- sample.kappa1
			}	
			
		if(move.inx == "WSH2"){		# update kappa2
			sample.kappa2		<- move.WSH.BweibScr(survObj, ini, priorPara, g = 2);
			ini$kappa2.ini		<- sample.kappa2
			}
						
		if(move.inx == "WSH3"){		# update kappa3
			sample.kappa3		<- move.WSH.BweibScr(survObj, ini, priorPara, g = 3);
			ini$kappa3.ini		<- sample.kappa3
			}
											

		# Updating the shared frailty 
	
		if(move.inx == "FP"){    #update gamma
			sample.fp		<- move.FP.BweibScr(survObj, ini)
			ini$gamma.ini	<- sample.fp;
			}
			
			
		# Updating the frailty parameter (dependence parameter)

		if(move.inx == "DP"){    #update theta
			sample.dp		<- move.DP.XiGam.BweibScr(survObj, ini, priorPara);
			ini$theta.ini	<- sample.dp$theta.ini;
			ini$xi.ini		<- sample.dp$xi.ini;
			mcmcOutcome$accept.theta	<- mcmcOutcome$accept.theta + sample.dp$accept;
			mcmcOutcome$accept.xi		<- mcmcOutcome$accept.xi + sample.dp$accept;
			}
			
		


	if((M %% thin) == 0){
			
      
      # storage of posterior samples from MCMC	
      
		mcmcOutcome$ini	<- ini      
      
      	mcmcOutcome$beta1.p	<- rbind(mcmcOutcome$beta1.p, ini$beta1.ini, deparse.level = 0)
		mcmcOutcome$beta2.p	<- rbind(mcmcOutcome$beta2.p, ini$beta2.ini, deparse.level = 0)
		mcmcOutcome$beta3.p	<- rbind(mcmcOutcome$beta3.p, ini$beta3.ini, deparse.level = 0)
		
		
		mcmcOutcome$alpha1.p			<- c(mcmcOutcome$alpha1.p, ini$alpha1.ini);
		mcmcOutcome$alpha2.p			<- c(mcmcOutcome$alpha2.p, ini$alpha2.ini);
		mcmcOutcome$alpha3.p			<- c(mcmcOutcome$alpha3.p, ini$alpha3.ini);
				
		mcmcOutcome$kappa1.p			<- c(mcmcOutcome$kappa1.p, ini$kappa1.ini);
		mcmcOutcome$kappa2.p			<- c(mcmcOutcome$kappa2.p, ini$kappa2.ini);
		mcmcOutcome$kappa3.p			<- c(mcmcOutcome$kappa3.p, ini$kappa3.ini);
								
		
		mcmcOutcome$gamma.p			<- rbind(mcmcOutcome$gamma.p, ini$gamma.ini[1:10], deparse.level = 0)
		write.table(t(ini$gamma.ini[priorPara$indGamma]), file = paste("mcmcOutcome/ch",chain,"gammaP.txt", sep = ""), 						append = T, col.names = F, row.names = F)
				
		
				
		
		
		mcmcOutcome$theta.p			<- c(mcmcOutcome$theta.p, ini$theta.ini);
		mcmcOutcome$xi.p			<- c(mcmcOutcome$xi.p, ini$xi.ini);
					
		}	
	

		

		# for monitoring the posterior samples	
		
	
	
	if(M %% save == 0 | M == num.reps){
		
		beta1.Med		<- apply(mcmcOutcome$beta1.p[(M/thin/2):(M/thin),], 2, mean, na.rm = T)
		beta1.0.975		<- apply(mcmcOutcome$beta1.p[(M/thin/2):(M/thin),], 2, quantile, prob = 0.975, na.rm = T)
		beta1.0.025		<- apply(mcmcOutcome$beta1.p[(M/thin/2):(M/thin),], 2, quantile, prob = 0.025, na.rm = T)
		
		beta2.Med		<- apply(mcmcOutcome$beta2.p[(M/thin/2):(M/thin),], 2, mean, na.rm = T)
		beta2.0.975		<- apply(mcmcOutcome$beta2.p[(M/thin/2):(M/thin),], 2, quantile, prob = 0.975, na.rm = T)
		beta2.0.025		<- apply(mcmcOutcome$beta2.p[(M/thin/2):(M/thin),], 2, quantile, prob = 0.025, na.rm = T)
		
		beta3.Med		<- apply(mcmcOutcome$beta3.p[(M/thin/2):(M/thin),], 2, mean, na.rm = T)
		beta3.0.975		<- apply(mcmcOutcome$beta3.p[(M/thin/2):(M/thin),], 2, quantile, prob = 0.975, na.rm = T)
		beta3.0.025		<- apply(mcmcOutcome$beta3.p[(M/thin/2):(M/thin),], 2, quantile, prob = 0.025, na.rm = T)
		

		
		# for beta1

		pdf(file = paste("mcmcOutcome/ch",chain,"regPara1.pdf", sep = ""), height = 10, width = 20)

		par(mfrow = c(3,6))		
		
		for(i in 1:p){
			plot(mcmcOutcome$beta1.p[(M/thin/2):(M/thin),i], type = "l", 			main = colnames(x)[i], xlab = "scan", ylab = "reg. coef.");
			abline(h = beta1.Med[i], col = "red")
			abline(h = beta1.0.975[i], lty = 2, col = "red")
			abline(h = beta1.0.025[i], lty = 2, col = "red")
			}
		
		dev.off()
		
	
		# for beta2

		pdf(file = paste("mcmcOutcome/ch",chain,"regPara2.pdf", sep = ""), height = 10, width = 20)

		par(mfrow = c(3,6))		
		
		for(i in 1:p){
			plot(mcmcOutcome$beta2.p[(M/thin/2):(M/thin),i], type = "l", 			main = colnames(x)[i], xlab = "scan", ylab = "reg. coef.")
			abline(h = beta2.Med[i], col = "red")
			abline(h = beta2.0.975[i], lty = 2, col = "red")
			abline(h = beta2.0.025[i], lty = 2, col = "red")
			}
		
		dev.off()	
						
		# for beta3

		pdf(file = paste("mcmcOutcome/ch",chain,"regPara3.pdf", sep = ""), height = 10, width = 20)

		par(mfrow = c(3,6))		
		
		for(i in 1:p){
			plot(mcmcOutcome$beta3.p[(M/thin/2):(M/thin),i], type = "l", 			main = colnames(x)[i], xlab = "scan", ylab = "reg. coef.");
			abline(h = beta3.Med[i], col = "red")
			abline(h = beta3.0.975[i], lty = 2, col = "red")
			abline(h = beta3.0.025[i], lty = 2, col = "red")
			}
		
		dev.off()			







		# hazard function
		
		 # alpha1, kappa1, baseline hazard 1

		alpha1.med	<- median(mcmcOutcome$alpha1.p[(length(mcmcOutcome$alpha1.p)/2):				length(mcmcOutcome$alpha1.p)])
		kappa1.med	<-  median(mcmcOutcome$kappa1.p[(length(mcmcOutcome$kappa1.p)/2):				length(mcmcOutcome$kappa1.p)])
				
		pdf(file = paste("mcmcOutcome/ch", chain, "basehaz1.pdf", sep = ""), height = 10, width = 30)
		
		par(mfrow = c(1, 3))
		
		plot(mcmcOutcome$alpha1.p[(length(mcmcOutcome$alpha1.p)/2):length(mcmcOutcome$alpha1.p)], 			type = "l", main = expression(alpha[1]), xlab = "scan", ylab = expression(alpha[1]));
		abline(h = alpha1.med, col = "red")
							
		plot(mcmcOutcome$kappa1.p[(length(mcmcOutcome$kappa1.p)/2):length(mcmcOutcome$kappa1.p)], 			type = "l", main = expression(kappa[1]), xlab = "scan", ylab = expression(kappa[1]));
		abline(h = kappa1.med, col = "red")
				
		plot(0:max(y1[delta1 == 1]),  log(alpha1.med*kappa1.med* c(0:max(y1[delta1 == 1]))^(alpha1.med - 1)), type = "l", 				main = "baseline hazard fucntion 1", xlab = "scan", ylab = "baselne hazard 1", col = "red");

				
		dev.off()
		
	
	
		 # alpha2, kappa2, baseline hazard 2
	
		
		alpha2.med	<-  median(mcmcOutcome$alpha2.p[(length(mcmcOutcome$alpha2.p)/2):				length(mcmcOutcome$alpha2.p)])
		kappa2.med	<-  median(mcmcOutcome$kappa2.p[(length(mcmcOutcome$kappa2.p)/2):				length(mcmcOutcome$kappa2.p)])
				
		pdf(file = paste("mcmcOutcome/ch", chain, "basehaz2.pdf", sep = ""), height = 10, width = 30)
		
		par(mfrow = c(1, 3))
		
		plot(mcmcOutcome$alpha2.p[(length(mcmcOutcome$alpha2.p)/2):length(mcmcOutcome$alpha2.p)], 			type = "l", main = expression(alpha[2]), xlab = "scan", ylab = expression(alpha[2]));
		abline(h = alpha2.med, col = "red")
							
		plot(mcmcOutcome$kappa2.p[(length(mcmcOutcome$kappa2.p)/2):length(mcmcOutcome$kappa2.p)], 			type = "l", main = expression(kappa[2]), xlab = "scan", ylab = expression(kappa[2]));
		abline(h = kappa2.med, col = "red")
						
		plot(0:max(y2[delta1 == 0 & delta2 == 1]),  				log(alpha2.med*kappa2.med* c(0:max(y2[delta1 == 0 & delta2 == 1]))^(alpha2.med - 1)), type = "l", 				main = "baseline hazard fucntion 2", xlab = "scan", ylab = "baselne hazard 2", col = "red");

				
		dev.off()
		



		 # alpha3, kappa3, baseline hazard 3

		
		alpha3.med	<-  median(mcmcOutcome$alpha3.p[(length(mcmcOutcome$alpha3.p)/2):				length(mcmcOutcome$alpha3.p)])
		kappa3.med	<-  median(mcmcOutcome$kappa3.p[(length(mcmcOutcome$kappa3.p)/2):				length(mcmcOutcome$kappa3.p)])
				
		pdf(file = paste("mcmcOutcome/ch", chain, "basehaz3.pdf", sep = ""), height = 10, width = 30)
		
		par(mfrow = c(1, 3))
		
		plot(mcmcOutcome$alpha3.p[(length(mcmcOutcome$alpha3.p)/2):length(mcmcOutcome$alpha3.p)], 			type = "l", main = expression(alpha[3]), xlab = "scan", ylab = expression(alpha[3]));
		abline(h = alpha3.med, col = "red")
							
		plot(mcmcOutcome$kappa3.p[(length(mcmcOutcome$kappa3.p)/2):length(mcmcOutcome$kappa3.p)], 			type = "l", main = expression(kappa[3]), xlab = "scan", ylab = expression(kappa[3]));
		abline(h = kappa3.med, col = "red")
				

		plot(0:max(y2[delta1 == 1 & delta2 == 1]),  				log(alpha3.med*kappa3.med* c(0:max(y2[delta1 == 1 & delta2 == 1]))^(alpha3.med - 1)), type = "l", 				main = "baseline hazard fucntion 3", xlab = "scan", ylab = "baselne hazard 3", col = "red");

				
		dev.off()
		

		
			
		# dependency parameter and frailty parameters
			
		pdf(file = paste("mcmcOutcome/ch", chain, "DepFrailtyPara.pdf", sep = ""), height = 10, width = 20)
		
		par(mfrow = c(2, 3))
		plot(mcmcOutcome$theta.p[(M/thin/2):(M/thin)], type = "l", main = "theta", xlab = "scan", ylab = "theta");
		abline(h = median(mcmcOutcome$theta.p[(M/thin/2):(M/thin)]), col = "red")
				
		plot(log(mcmcOutcome$gamma.p[(M/thin/2):(M/thin), 1]), type = "l", main = "log-gamma1", xlab = "scan", ylab = "log-gamma1");
		plot(log(mcmcOutcome$gamma.p[(M/thin/2):(M/thin), 2]), type = "l", main = "log-gamma2", xlab = "scan", ylab = "log-gamma2");
		plot(log(mcmcOutcome$gamma.p[(M/thin/2):(M/thin), 3]), type = "l", main = "log-gamma3", xlab = "scan", ylab = "log-gamma3");
		plot(log(mcmcOutcome$gamma.p[(M/thin/2):(M/thin), 4]), type = "l", main = "log-gamma4", xlab = "scan", ylab = "log-gamma4");
		plot(log(mcmcOutcome$gamma.p[(M/thin/2):(M/thin), 5]), type = "l", main = "log-gamma5", xlab = "scan", ylab = "log-gamma5");
										
		dev.off()
		
		}
		




		# store the mcmc outcomes

		if(M %% save == 0 | M == num.reps){
			
			save(mcmcOutcome, file = paste("mcmcOutcome/BwbScr.a", a1, "b", b1, "c", c1, "d", d1, 				"psi", priorPara$psi, "omega", priorPara$omega, "ch", chain, ".RData", sep = ""))
		
			}
	
		} # end of the 'for' loop for MCMC sampling
	

	return(mcmcOutcome)
	
	} # end of BpeScr function
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		
	