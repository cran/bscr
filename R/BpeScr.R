BpeScr <-
function(survObj, priorPara, initial, num.reps, thin = 1, chain = 1, save = 1000, RJ = TRUE){
			
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


	s1max 	<- priorPara$s1max
	s2max 	<- priorPara$s2max
	s3max 	<- priorPara$s3max

	initial$J1				<- length(initial$s1)-1
	initial$J2				<- length(initial$s2)-1
	initial$J3				<- length(initial$s3)-1
			
	survObj$case1_y1leq.s1 	<- which(delta1 == 1 & y1 <= s1max)
	survObj$case01y2leq.s1 	<- which(delta1 == 0 & delta2 == 1 & y2 <= s1max)
	survObj$case00yleq.s1	<- which(delta1 == 0 & delta2 == 0 & y1 <= s1max)
	survObj$case1_or01or00y1y2leq.s1 <- which( (delta1 == 1 & y1 <= s1max) 												| (delta1 == 0 & delta2 == 1 & y2 <= s1max)												| (delta1 == 0 & delta2 == 0 & y1 <= s1max))
	survObj$case1_or01or00y1y2gre.s1 <- which( (delta1 == 1 & y1 > s1max) 												| (delta1 == 0 & delta2 == 1 & y2 > s1max)												| (delta1 == 0 & delta2 == 0 & y1 > s1max))

	survObj$case1_y1leq.s2 	<- which(delta1 == 1 & y1 <= s2max)
	survObj$case01y2leq.s2 	<- which(delta1 == 0 & delta2 == 1 & y2 <= s2max)
	survObj$case00yleq.s2	<- which(delta1 == 0 & delta2 == 0 & y1 <= s2max)
	survObj$case1_or01or00y1y2leq.s2 <- which( (delta1 == 1 & y1 <= s2max) 												| (delta1 == 0 & delta2 == 1 & y2 <= s2max)												| (delta1 == 0 & delta2 == 0 & y1 <= s2max))
	survObj$case1_or01or00y1y2gre.s2 <- which( (delta1 == 1 & y1 > s2max) 												| (delta1 == 0 & delta2 == 1 & y2 > s2max)												| (delta1 == 0 & delta2 == 0 & y1 > s2max))

	survObj$case1_y1leq.s3 	<- which(delta1 == 1 & y1 <= s3max)
	survObj$case01y2leq.s3 	<- which(delta1 == 0 & delta2 == 1 & y2 <= s3max)
	survObj$case00yleq.s3	<- which(delta1 == 0 & delta2 == 0 & y1 <= s3max)
	survObj$case1_or01or00y1y2leq.s3 <- which( (delta1 == 1 & y1 <= s3max) 												| (delta1 == 0 & delta2 == 1 & y2 <= s3max)												| (delta1 == 0 & delta2 == 0 & y1 <= s3max))
	survObj$case1_or01or00y1y2gre.s3 <- which( (delta1 == 1 & y1 > s3max) 												| (delta1 == 0 & delta2 == 1 & y2 > s3max)												| (delta1 == 0 & delta2 == 0 & y1 > s3max))
	
	survObj$case11y2leq.s3	<- which(delta1 == 1 & delta2 == 1 & y2 <= s3max)
	survObj$case10y2leq.s3	<- which(delta1 == 1 & delta2 == 0 & y2 <= s3max)
	survObj$case11or10y2leq.s3	<- which((delta1 == 1 & delta2 == 1 & y2 <= s3max)										|(delta1 == 1 & delta2 == 0 & y2 <= s3max))
	survObj$case11or10y2gre.s3	<- which((delta1 == 1 & delta2 == 1 & y2 > s3max)										|(delta1 == 1 & delta2 == 0 & y2 > s3max))
		
		
	
			
	ini$s1			<- initial$s1
	ini$s2			<- initial$s2
	ini$s3			<- initial$s3		
	ini$J1			<- initial$J1
	ini$J2			<- initial$J2
	ini$J3			<- initial$J3		
	ini$beta1.ini	<- initial$beta1.ini
	ini$beta2.ini	<- initial$beta2.ini
	ini$beta3.ini	<- initial$beta3.ini		
	ini$lambda1		<- initial$lambda1
	ini$lambda2		<- initial$lambda2
	ini$lambda3		<- initial$lambda3		
	ini$sig.sq.lam1	<- initial$sig.sq.lam1
	ini$sig.sq.lam2	<- initial$sig.sq.lam2
	ini$sig.sq.lam3	<- initial$sig.sq.lam3		
	ini$mu.lam1		<- initial$mu.lam1
	ini$mu.lam2		<- initial$mu.lam2
	ini$mu.lam3		<- initial$mu.lam3	
	ini$theta.ini	<- initial$theta.ini
	ini$xi.ini		<- 1/initial$theta.ini
	ini$gamma.ini	<- initial$gamma.ini
	ini$xbeta1.ini	<- as.vector(x %*% ini$beta1.ini)
	ini$xbeta2.ini	<- as.vector(x %*% ini$beta2.ini)
	ini$xbeta3.ini	<- as.vector(x %*% ini$beta3.ini)
	
	
			
	c.lam1	<- priorPara$c.lam1
	c.lam2	<- priorPara$c.lam2
	c.lam3	<- priorPara$c.lam3		
	C1		<- priorPara$C1
	C2		<- priorPara$C2
	C3		<- priorPara$C3		
	a1		<- priorPara$a1
	a2		<- priorPara$a2
	a3		<- priorPara$a3		
	b1		<- priorPara$b1
	b2		<- priorPara$b2
	b3		<- priorPara$b3		
	alpha1	<- priorPara$alpha1
	alpha2	<- priorPara$alpha2
	alpha3	<- priorPara$alpha3		
	psi		<- priorPara$psi
	omega	<- priorPara$omega
	delPert1 <- priorPara$delPert1
	delPert2 <- priorPara$delPert2
	delPert3 <- priorPara$delPert3

								

	intv1	<- setting.interval.BpeScr(survObj, ini, priorPara, g = 1);
	intv2	<- setting.interval.BpeScr(survObj, ini, priorPara, g = 2);
	intv3	<- setting.interval.BpeScr(survObj, ini, priorPara, g = 3);
		
	ini$ind.r1.sc1	<- intv1$ind.r1;
	ini$ind.r1.sc2	<- intv2$ind.r1;	
	ini$ind.r2.sc3	<- intv3$ind.r2;
	
	ini$ind.d1.sc1	<- intv1$ind.d1;
	ini$ind.d1.sc2	<- intv2$ind.d1;
	ini$ind.d1.sc3	<- intv3$ind.d1;
	ini$ind.d2.sc1	<- intv1$ind.d2;
	ini$ind.d2.sc2	<- intv2$ind.d2;
	ini$ind.d3.sc3	<- intv3$ind.d3;	
	
	ini$ind.d.temp.sc1	<- intv1$ind.d.temp
	ini$ind.d.temp.sc2	<- intv2$ind.d.temp		
	ini$ind.d.temp.sc3	<- intv3$ind.d.temp
		
	
	ini$d1.sc1	<- intv1$d1;
	ini$d1.sc2	<- intv2$d1;
	ini$d1.sc3	<- intv3$d1;	
	
	ini$d2.sc1	<- intv1$d2;
	ini$d2.sc2	<- intv2$d2;
	ini$d3.sc3	<- intv3$d3;
			
	ini$Delta1	<- cal.Del1.2.BpeScr(survObj, ini, g = 1);
	ini$Delta2	<- cal.Del1.2.BpeScr(survObj, ini, g = 2);
	ini$Delta3  <- cal.Del3.BpeScr(survObj, ini, g = 3);
	
	cal.Sig1	<- cal.Sigma.lam.BpeScr(ini, priorPara, g = 1);
	cal.Sig2	<- cal.Sigma.lam.BpeScr(ini, priorPara, g = 2);
	cal.Sig3	<- cal.Sigma.lam.BpeScr(ini, priorPara, g = 3);
	
	ini$Sigma.lam1	<- cal.Sig1$Sigma.lam;
	ini$Sigma.lam2	<- cal.Sig2$Sigma.lam;
	ini$Sigma.lam3	<- cal.Sig3$Sigma.lam;
	
	ini$inv.Sigma.lam1	<- solve(ini$Sigma.lam1)
	ini$inv.Sigma.lam2	<- solve(ini$Sigma.lam2)
	ini$inv.Sigma.lam3	<- solve(ini$Sigma.lam3)
		
	ini$W1			<- cal.Sig1$W;
	ini$W2			<- cal.Sig2$W;
	ini$W3			<- cal.Sig3$W;
	
	ini$Q1			<- cal.Sig1$Q;
	ini$Q2			<- cal.Sig2$Q;
	ini$Q3			<- cal.Sig3$Q;			
							
	
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
		
	mcmcOutcome$beta1.p			<- initial$beta1.ini
	mcmcOutcome$beta2.p			<- initial$beta2.ini
	mcmcOutcome$beta3.p			<- initial$beta3.ini
			
	mcmcOutcome$lambda1.p		<- as.list(c(rep(NA, num.reps/thin)));
	mcmcOutcome$lambda2.p		<- as.list(c(rep(NA, num.reps/thin)));
	mcmcOutcome$lambda3.p		<- as.list(c(rep(NA, num.reps/thin)));
			
	mcmcOutcome$s1.p			<- as.list(c(rep(NA, num.reps/thin)));
	mcmcOutcome$s2.p			<- as.list(c(rep(NA, num.reps/thin)));
	mcmcOutcome$s3.p			<- as.list(c(rep(NA, num.reps/thin)));
			
	mcmcOutcome$J1.p			<- initial$J1;
	mcmcOutcome$J2.p			<- initial$J2;
	mcmcOutcome$J3.p			<- initial$J3;	
		
	mcmcOutcome$s1.p[[1]]		<- initial$s1;
	mcmcOutcome$s2.p[[1]]		<- initial$s2;
	mcmcOutcome$s3.p[[1]]		<- initial$s3;
			
	mcmcOutcome$lambda1.p[[1]]	<- initial$lambda1;
	mcmcOutcome$lambda2.p[[1]]	<- initial$lambda2;
	mcmcOutcome$lambda3.p[[1]]	<- initial$lambda3;
			
	mcmcOutcome$log.jpost 		<- log.like <- NULL;
	
	mcmcOutcome$mu.lam1.p 		<- initial$mu.lam1
	mcmcOutcome$mu.lam2.p 		<- initial$mu.lam2
	mcmcOutcome$mu.lam3.p 		<- initial$mu.lam3
		
	mcmcOutcome$sig.sq.lam1.p 	<- initial$sig.sq.lam1
	mcmcOutcome$sig.sq.lam2.p 	<- initial$sig.sq.lam2
	mcmcOutcome$sig.sq.lam3.p 	<- initial$sig.sq.lam3
	
	mcmcOutcome$gamma.p			<- initial$gamma.ini[1:5]
	mcmcOutcome$theta.p			<- initial$theta.ini
	mcmcOutcome$xi.p			<- initial$xi.ini
		
	mcmcOutcome$accept.beta1		<- c(rep(0, p))
	mcmcOutcome$accept.beta2		<- c(rep(0, p))
	mcmcOutcome$accept.beta3		<- c(rep(0, p))
	
	mcmcOutcome$accept.lambda1		<- as.list(c(rep(0, num.reps/thin)))
	mcmcOutcome$accept.lambda2		<- as.list(c(rep(0, num.reps/thin)))
	mcmcOutcome$accept.lambda3		<- as.list(c(rep(0, num.reps/thin)))
				
	mcmcOutcome$accept.theta	<- 0
	mcmcOutcome$accept.xi		<- 0		
	
	mcmcOutcome$accept.bi1		<- 0
	mcmcOutcome$accept.bi2		<- 0
	mcmcOutcome$accept.bi3		<- 0
			
	mcmcOutcome$accept.di1		<- 0
	mcmcOutcome$accept.di2		<- 0
	mcmcOutcome$accept.di3		<- 0
		
	accept.lambda1	<- as.list(c(rep(NA, num.reps/thin)))
	accept.lambda2	<- as.list(c(rep(NA, num.reps/thin)))
	accept.lambda3	<- as.list(c(rep(NA, num.reps/thin)))
	

	# for move probability	


	J1max <- length(priorPara$s.propBI1) - 1
	J2max <- length(priorPara$s.propBI2) - 1
	J3max <- length(priorPara$s.propBI3) - 1
			
	pB1			<- pmin(1, alpha1/c(1:J1max))
	pD1			<- pmin(1, (c(1:J1max)-1)/alpha1)
	rho.lam1	<- min(C1 / (pB1 + pD1))
	
	pB2			<- pmin(1, alpha2/c(1:J2max))
	pD2			<- pmin(1, (c(1:J2max)-1)/alpha2)
	rho.lam2	<- min(C2 / (pB2 + pD2))	
	
	pB3			<- pmin(1, alpha3/c(1:J3max))
	pD3			<- pmin(1, (c(1:J3max)-1)/alpha3)
	rho.lam3	<- min(C3 / (pB3 + pD3))	
	
	dir.create('mcmcOutcome')
	
		
	### Start MCMC sampling
	

	for(M in 1:num.reps){
		
		if(M %% 100 == 0){
			cat("Chain :", chain, "Scan :", M,fill=TRUE)
			fsh()
			}
	
	
		# move probability with RJ MCMC	
		
		if(RJ){

			if(ini$J1 < J1max){
				p.BI1	<- rho.lam1 * min(1, alpha1/(ini$J1+1))
				p.DI1	<- rho.lam1 * min(1, (ini$J1)/alpha1)
				}
					
			if(ini$J2 < J2max){
				p.BI2	<- rho.lam2 * min(1, alpha2/(ini$J2+1))
				p.DI2	<- rho.lam2 * min(1, (ini$J2)/alpha2)
				}
		
			if(ini$J3 < J3max){
				p.BI3	<- rho.lam3 * min(1, alpha3/(ini$J3+1))
				p.DI3	<- rho.lam3 * min(1, (ini$J3)/alpha3)
				}			

			if(ini$J1 >= J1max){
				p.BI1	<- 0
				p.DI1	<- rho.lam1 * 2
				}
					
			if(ini$J2 >= J2max){
				p.BI2	<- 0
				p.DI2	<- rho.lam2 * 2
				}
		
			if(ini$J3 >= J3max){
				p.BI3	<- 0
				p.DI3	<- rho.lam3 * 2
				}		
					

			#p.RP1 <- p.BH1 <- p.SP1 <- (1 - p.BI1 - p.DI1 - p.BI2 - p.DI2 - p.BI3 - p.DI3) / 11;
			#p.RP2 <- p.BH2 <- p.SP2 <- (1 - p.BI1 - p.DI1 - p.BI2 - p.DI2 - p.BI3 - p.DI3) / 11;
			#p.RP3 <- p.BH3 <- p.SP3 <- (1 - p.BI1 - p.DI1 - p.BI2 - p.DI2 - p.BI3 - p.DI3) / 11;
			#p.FP  <- p.DP			<- (1 - p.BI1 - p.DI1 - p.BI2 - p.DI2 - p.BI3 - p.DI3) / 11;
									
			#move.inx	<- sample(c("RP1", "RP2", "RP3", "BH1", "BH2", "BH3", 								"SP1", "SP2", "SP3", "FP", "DP", "BI1", "BI2", "BI3", "DI1", "DI2", "DI3"), 								size=1, prob=c(p.RP1, p.RP2, p.RP3, p.BH1, p.BH2, p.BH3, 								p.SP1, p.SP2, p.SP3, p.FP, p.DP, p.BI1, p.BI2, p.BI3, p.DI1, p.DI2, p.DI3));
			
			
			
			p.DP	<- 0.2
			p.RP1 	<- p.BH1 <- p.SP1 <- (1 - p.BI1 - p.DI1 - p.BI2 - p.DI2 - p.BI3 - p.DI3 - p.DP) / 10;
			p.RP2 	<- p.BH2 <- p.SP2 <- (1 - p.BI1 - p.DI1 - p.BI2 - p.DI2 - p.BI3 - p.DI3 - p.DP) / 10;
			p.RP3 	<- p.BH3 <- p.SP3 <- (1 - p.BI1 - p.DI1 - p.BI2 - p.DI2 - p.BI3 - p.DI3 - p.DP) / 10;
			p.FP  			<- (1 - p.BI1 - p.DI1 - p.BI2 - p.DI2 - p.BI3 - p.DI3 - p.DP) / 10;
			
									
			move.inx	<- sample(c("RP1", "RP2", "RP3", "BH1", "BH2", "BH3", 								"SP1", "SP2", "SP3", "FP", "DP", "BI1", "BI2", "BI3", "DI1", "DI2", "DI3"), 								size=1, prob=c(p.RP1, p.RP2, p.RP3, p.BH1, p.BH2, p.BH3, 								p.SP1, p.SP2, p.SP3, p.FP, p.DP, p.BI1, p.BI2, p.BI3, p.DI1, p.DI2, p.DI3));
						
			} 			


		# move probability without RJ MCMC
		
		if(!RJ){

			p.BI1 <- p.BI2 <- p.BI3 <- p.DI1 <- p.DI2 <- p.DI3 <- 0
		 	 	
			p.RP1 <- p.BH1 <- p.SP1 <- (1 - p.BI1 - p.DI1 - p.BI2 - p.DI2 - p.BI3 - p.DI3) / 11;
			p.RP2 <- p.BH2 <- p.SP2 <- (1 - p.BI1 - p.DI1 - p.BI2 - p.DI2 - p.BI3 - p.DI3) / 11;
			p.RP3 <- p.BH3 <- p.SP3 <- (1 - p.BI1 - p.DI1 - p.BI2 - p.DI2 - p.BI3 - p.DI3) / 11;
			p.FP  <- p.DP			<- (1 - p.BI1 - p.DI1 - p.BI2 - p.DI2 - p.BI3 - p.DI3) / 11;
									
			move.inx	<- sample(c("RP1", "RP2", "RP3", "BH1", "BH2", "BH3", "SP1", "SP2", "SP3", "FP", "DP"),							 size=1, prob=c(p.RP1, p.RP2, p.RP3, p.BH1, p.BH2, p.BH3, p.SP1, p.SP2, p.SP3, 								p.FP, p.DP));
						
			} 			
			
		mcmcOutcome$move.ind	<- c(mcmcOutcome$move.ind, move.inx);

		
			
		
		# Updating Regression coefficients

		if(move.inx == "RP1"){		# update beta1
			sample.beta1		<- move.RP.BpeScr(survObj, ini, g = 1)
			ini$xbeta1.ini		<- sample.beta1$xbeta.ini
			ini$beta1.ini		<- sample.beta1$beta.ini;
			mcmcOutcome$accept.beta1		<- mcmcOutcome$accept.beta1 + sample.beta1$accept;
			}

		if(move.inx == "RP2"){		# update beta2
			sample.beta2		<- move.RP.BpeScr(survObj, ini, g = 2)
			ini$xbeta2.ini		<- sample.beta2$xbeta.ini
			ini$beta2.ini		<- sample.beta2$beta.ini;
			mcmcOutcome$accept.beta2	<- mcmcOutcome$accept.beta2 + sample.beta2$accept;
			}
		
		if(move.inx == "RP3"){		# update beta3
			sample.beta3		<- move.RP.BpeScr(survObj, ini, g = 3)
			ini$xbeta3.ini		<- sample.beta3$xbeta.ini
			ini$beta3.ini		<- sample.beta3$beta.ini;
			mcmcOutcome$accept.beta3	<- mcmcOutcome$accept.beta3 + sample.beta3$accept;
			}	
				

		# Updating log-baseline hazards
		
		if(move.inx == "BH1"){		# update lambda1
			mcmcOutcome$accept.lambda1[[M]]	<- rep(0, ini$J1+1)
			sample.lambda1		<- move.BH.BpeScr(survObj, ini, g = 1)
			ini$lambda1			<- sample.lambda1$lambda;
			mcmcOutcome$accept.lambda1[[M]]	<- sample.lambda1$accept;
			}

		if(move.inx == "BH2"){		# update lambda2
			mcmcOutcome$accept.lambda2[[M]]	<- rep(0, ini$J2+1)
			sample.lambda2		<- move.BH.BpeScr(survObj, ini, g = 2)
			ini$lambda2			<- sample.lambda2$lambda;
			mcmcOutcome$accept.lambda2[[M]]	<- sample.lambda2$accept;
			}
			
		if(move.inx == "BH3"){		# update lambda3
			mcmcOutcome$accept.lambda3[[M]]	<- rep(0, ini$J3+1)
			sample.lambda3		<- move.BH.BpeScr(survObj, ini, g = 3)
			ini$lambda3			<- sample.lambda3$lambda;
			mcmcOutcome$accept.lambda3[[M]]	<- sample.lambda3$accept;
			}	
								
			
		# Updating second stage survival components: mu.lam and sigma.sq.lam
	
		if(move.inx == "SP1"){		# update mu.lam1 and sigma.sq.lam1
			sample.sp1				<- move.SP.BpeScr(survObj, ini, priorPara, g = 1)
			ini$mu.lam1				<- sample.sp1$mu.lam.sample;
			ini$sig.sq.lam1			<- sample.sp1$sigma.sq.lam.sample;
			}
			
		if(move.inx == "SP2"){		# update mu.lam2 and sigma.sq.lam2
			sample.sp2				<- move.SP.BpeScr(survObj, ini, priorPara, g = 2)
			ini$mu.lam2				<- sample.sp2$mu.lam.sample;
			ini$sig.sq.lam2			<- sample.sp2$sigma.sq.lam.sample;
			}	
			
		if(move.inx == "SP3"){		# update mu.lam3 and sigma.sq.lam3
			sample.sp3				<- move.SP.BpeScr(survObj, ini, priorPara, g = 3)
			ini$mu.lam3				<- sample.sp3$mu.lam.sample;
			ini$sig.sq.lam3			<- sample.sp3$sigma.sq.lam.sample;
			}
												

		# Updating the shared frailty 
	
		if(move.inx == "FP"){    #update gamma
			sample.fp		<- move.FP.BpeScr(survObj, ini)
			ini$gamma.ini	<- sample.fp;
			}
			
			
		# Updating the frailty parameter (dependence parameter)

		if(move.inx == "DP"){    #update theta
			sample.dp		<- move.DP.XiGam.BpeScr(survObj, ini, priorPara)
			ini$theta.ini	<- sample.dp$theta.ini;
			ini$xi.ini		<- sample.dp$xi.ini;
			mcmcOutcome$accept.theta	<- mcmcOutcome$accept.theta + sample.dp$accept;
			mcmcOutcome$accept.xi		<- mcmcOutcome$accept.xi + sample.dp$accept;
			}
			


		# Creating a new time split: Birth move
	
		if(move.inx == "BI1"){ 			# for s_1
			sample.bi1	<- move.BI.BpeScr(survObj, ini, priorPara, g = 1);
			if(sample.bi1$accept == 1){
				ini	<- sample.bi1$ini			
				mcmcOutcome$accept.bi1		<- mcmcOutcome$accept.bi1 + sample.bi1$accept;
				}
			}
			
		if(move.inx == "BI2"){ 			# for s_2
			sample.bi2	<- move.BI.BpeScr(survObj, ini, priorPara, g = 2); 
			if(sample.bi2$accept == 1){
				ini	<- sample.bi2$ini			
				mcmcOutcome$accept.bi2		<- mcmcOutcome$accept.bi2 + sample.bi2$accept;
				}
			}		
					
		if(move.inx == "BI3"){ 			# for s_3
			sample.bi3	<- move.BI.BpeScr(survObj, ini, priorPara, g = 3); 
			if(sample.bi3$accept == 1){
				ini	<- sample.bi3$ini			
				mcmcOutcome$accept.bi3		<- mcmcOutcome$accept.bi3 + sample.bi3$accept;
				}
			}
			
							
		# Removing a time split: Reverse move
	
		if(move.inx == "DI1"){ 			# for s_1
			sample.di1	<- move.DI.BpeScr(survObj, ini, priorPara, g = 1);
			if(sample.di1$accept == 1){
				ini	<- sample.di1$ini			
				mcmcOutcome$accept.di1		<- mcmcOutcome$accept.di1 + sample.di1$accept;
				}
			}
			
		if(move.inx == "DI2"){ 			# for s_2
			sample.di2	<- move.DI.BpeScr(survObj, ini, priorPara, g = 2);
			if(sample.di2$accept == 1){
				ini	<- sample.di2$ini			
				mcmcOutcome$accept.di2		<- mcmcOutcome$accept.di2 + sample.di2$accept;
				}
			}
			
		if(move.inx == "DI3"){ 			# for s_3
			sample.di3	<- move.DI.BpeScr(survObj, ini, priorPara, g = 3);
			if(sample.di3$accept == 1){
				ini	<- sample.di3$ini			
				mcmcOutcome$accept.di3		<- mcmcOutcome$accept.di3 + sample.di3$accept;
				}
			}	
		


	if((M %% thin) == 0){
			
      
      # storage of posterior samples from MCMC	
      
		mcmcOutcome$ini	<- ini      
      
      	mcmcOutcome$beta1.p	<- rbind(mcmcOutcome$beta1.p, ini$beta1.ini, deparse.level = 0)
		mcmcOutcome$beta2.p	<- rbind(mcmcOutcome$beta2.p, ini$beta2.ini, deparse.level = 0)
		mcmcOutcome$beta3.p	<- rbind(mcmcOutcome$beta3.p, ini$beta3.ini, deparse.level = 0)
		
		mcmcOutcome$lambda1.p[[M/thin+1]] <- ini$lambda1;
		mcmcOutcome$lambda2.p[[M/thin+1]] <- ini$lambda2;
		mcmcOutcome$lambda3.p[[M/thin+1]] <- ini$lambda3;
		
		mcmcOutcome$mu.lam1.p			<- c(mcmcOutcome$mu.lam1.p, ini$mu.lam1);
		mcmcOutcome$mu.lam2.p			<- c(mcmcOutcome$mu.lam2.p, ini$mu.lam2);
		mcmcOutcome$mu.lam3.p			<- c(mcmcOutcome$mu.lam3.p, ini$mu.lam3);
		mcmcOutcome$sig.sq.lam1.p		<- c(mcmcOutcome$sig.sq.lam1.p, ini$sig.sq.lam1);
		mcmcOutcome$sig.sq.lam2.p		<- c(mcmcOutcome$sig.sq.lam2.p, ini$sig.sq.lam2);
		mcmcOutcome$sig.sq.lam3.p		<- c(mcmcOutcome$sig.sq.lam3.p, ini$sig.sq.lam3);

		mcmcOutcome$gamma.p			<- rbind(mcmcOutcome$gamma.p, ini$gamma.ini[1:5], deparse.level = 0)
		
		write.table(t(ini$gamma.ini[priorPara$indGamma]), file = paste("mcmcOutcome/ch",chain,"gammaP.txt", sep = ""), 						append = T, col.names = F, row.names = F)
				
		
		
		mcmcOutcome$theta.p			<- c(mcmcOutcome$theta.p, ini$theta.ini);
		mcmcOutcome$xi.p			<- c(mcmcOutcome$xi.p, ini$xi.ini);
				
			
		mcmcOutcome$s1.p[[M/thin+1]]		<- ini$s1;
		mcmcOutcome$s2.p[[M/thin+1]]		<- ini$s2;
		mcmcOutcome$s3.p[[M/thin+1]]		<- ini$s3;
			
		mcmcOutcome$J1.p			<- c(mcmcOutcome$J1.p, ini$J1);
		mcmcOutcome$J2.p			<- c(mcmcOutcome$J2.p, ini$J2);
		mcmcOutcome$J3.p			<- c(mcmcOutcome$J3.p, ini$J3);
		
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
			plot(mcmcOutcome$beta2.p[(M/thin/2):(M/thin),i], type = "l", 			main = colnames(x)[i], xlab = "scan", ylab = "reg. coef.");
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


		# second stage survival components and the number of intervals
			
		pdf(file = paste("mcmcOutcome/ch", chain, "otherPara1.pdf", sep = ""), height = 10, width = 30)
		
		par(mfrow = c(1, 3))
		
		plot(log(mcmcOutcome$sig.sq.lam1.p[(length(mcmcOutcome$sig.sq.lam1.p)/2):				length(mcmcOutcome$sig.sq.lam1.p)]), type = "l", 				main = "log(sigmaL1)", xlab = "scan", ylab = "log.sigmaL1");
		abline(h = median(log(mcmcOutcome$sig.sq.lam1.p[(length(mcmcOutcome$sig.sq.lam1.p)/2):				length(mcmcOutcome$sig.sq.lam1.p)])), col = "red")
							
		plot(mcmcOutcome$mu.lam1.p[(length(mcmcOutcome$mu.lam1.p)/2):				length(mcmcOutcome$mu.lam1.p)], type = "l", main = "muL1", xlab = "scan", ylab = "muL1");
		abline(h = median(mcmcOutcome$mu.lam1.p[(length(mcmcOutcome$mu.lam1.p)/2):				length(mcmcOutcome$mu.lam1.p)]), col = "red")
		
		plot(mcmcOutcome$J1.p[(M/thin/2):(M/thin)], type = "l", 				main = "#intervals", xlab = "scan", ylab = "J1");
		abline(h = median(mcmcOutcome$J1.p[(M/thin/2):(M/thin)]), col = "red")
				
		dev.off()
		
		
		
		pdf(file = paste("mcmcOutcome/ch", chain, "otherPara2.pdf", sep = ""), height = 10, width = 30)
		
		par(mfrow = c(1, 3))
		plot(log(mcmcOutcome$sig.sq.lam2.p[(length(mcmcOutcome$sig.sq.lam2.p)/2):				length(mcmcOutcome$sig.sq.lam2.p)]), type = "l", 				main = "log(sigmaL2)", xlab = "scan", ylab = "log.sigmaL2");
		abline(h = median(log(mcmcOutcome$sig.sq.lam2.p[(length(mcmcOutcome$sig.sq.lam2.p)/2):				length(mcmcOutcome$sig.sq.lam2.p)])), col = "red")
				
		plot(mcmcOutcome$mu.lam2.p[(length(mcmcOutcome$mu.lam2.p)/2):				length(mcmcOutcome$mu.lam2.p)], type = "l", 				main = "muL2", xlab = "scan", ylab = "muL2");
		abline(h = median(mcmcOutcome$mu.lam2.p[(length(mcmcOutcome$mu.lam2.p)/2):				length(mcmcOutcome$mu.lam2.p)]), col = "red")
					
		plot(mcmcOutcome$J2.p[(M/thin/2):(M/thin)], type = "l", 				main = "#intervals", xlab = "scan", ylab = "J2");
		abline(h = median(mcmcOutcome$J2.p[(M/thin/2):(M/thin)]), col = "red")
		
		dev.off()
		
		
				
		pdf(file = paste("mcmcOutcome/ch", chain, "otherPara3.pdf", sep = ""), height = 10, width = 30)
		
		par(mfrow = c(1, 3))
		
		plot(log(mcmcOutcome$sig.sq.lam3.p[(length(mcmcOutcome$sig.sq.lam3.p)/2):				length(mcmcOutcome$sig.sq.lam3.p)]), type = "l", 				main = "log(sigmaL3)", xlab = "scan", ylab = "log.sigmaL3");
		abline(h = median(log(mcmcOutcome$sig.sq.lam3.p[(length(mcmcOutcome$sig.sq.lam3.p)/2):				length(mcmcOutcome$sig.sq.lam3.p)])), col = "red")
				
		plot(mcmcOutcome$mu.lam3.p[(length(mcmcOutcome$mu.lam3.p)/2):				length(mcmcOutcome$mu.lam3.p)], type = "l", 				main = "muL3", xlab = "scan", ylab = "muL3");
		abline(h = median(mcmcOutcome$mu.lam3.p[(length(mcmcOutcome$mu.lam3.p)/2):				length(mcmcOutcome$mu.lam3.p)]), col = "red")
		
		plot(mcmcOutcome$J3.p[(M/thin/2):(M/thin)], type = "l", 				main = "#intervals", xlab = "scan", ylab = "J3");
		
		abline(h = median(mcmcOutcome$J3.p[(M/thin/2):(M/thin)]), col = "red")
				
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
		


		

	if(M %% save == 0 | M == num.reps){


		ss1	<- seq(1, priorPara$s1max, 1)
		ss2	<- seq(1, priorPara$s2max, 1)
		ss3	<- seq(1, priorPara$s3max, 1)
						
		JJ1	<- length(ss1)-1;
		JJ2	<- length(ss2)-1;
		JJ3	<- length(ss3)-1;				

		lambda1.fin	<- matrix(NA, M/thin/2, JJ1+1);
		lambda2.fin	<- matrix(NA, M/thin/2, JJ2+1);
		lambda3.fin	<- matrix(NA, M/thin/2, JJ3+1);
						
		ind.lam1		<- rep(NA, JJ1+1);
		ind.lam2		<- rep(NA, JJ2+1);
		ind.lam3		<- rep(NA, JJ3+1);
		
		for(i in (M/thin/2 + 1):(M/thin)){
	
			ss1.n	<- length(mcmcOutcome$s1.p[[i]]);
			ss2.n	<- length(mcmcOutcome$s2.p[[i]]);
			ss3.n	<- length(mcmcOutcome$s3.p[[i]]);
									
			ind.lam1[(ss1 <= mcmcOutcome$s1.p[[i]][1]) & (ss1 > 0)] <- 1;
			ind.lam2[(ss2 <= mcmcOutcome$s2.p[[i]][1]) & (ss2 > 0)] <- 1;
			ind.lam3[(ss3 <= mcmcOutcome$s3.p[[i]][1]) & (ss3 > 0)] <- 1;
									
			for(j in 2:ss1.n){
				ind.lam1[(ss1 <= mcmcOutcome$s1.p[[i]][j]) & (ss1 > mcmcOutcome$s1.p[[i]][j-1])] <- j;
				}
			for(j in 2:ss2.n){
				ind.lam2[(ss2 <= mcmcOutcome$s2.p[[i]][j]) & (ss2 > mcmcOutcome$s2.p[[i]][j-1])] <- j;
				}
			for(j in 2:ss3.n){
				ind.lam3[(ss3 <= mcmcOutcome$s3.p[[i]][j]) & (ss3 > mcmcOutcome$s3.p[[i]][j-1])] <- j;
				}	

			lambda1.fin[i - (M/thin/2), ] <- mcmcOutcome$lambda1.p[[i]][ind.lam1];
			lambda2.fin[i - (M/thin/2), ] <- mcmcOutcome$lambda2.p[[i]][ind.lam2];
			lambda3.fin[i - (M/thin/2), ] <- mcmcOutcome$lambda3.p[[i]][ind.lam3];
			}


		lambda1.mean		<- apply(lambda1.fin, 2, mean, na.rm = T);
		lambda2.mean		<- apply(lambda2.fin, 2, mean, na.rm = T);
		lambda3.mean		<- apply(lambda3.fin, 2, mean, na.rm = T);
						
		lambda1.sd		<- apply(lambda1.fin, 2, sd, na.rm = T);
		lambda2.sd		<- apply(lambda2.fin, 2, sd, na.rm = T);
		lambda3.sd		<- apply(lambda3.fin, 2, sd, na.rm = T);
						
		lambda1.med		<- apply(lambda1.fin, 2, median, na.rm = T);
		lambda2.med		<- apply(lambda2.fin, 2, median, na.rm = T);
		lambda3.med		<- apply(lambda3.fin, 2, median, na.rm = T);
						
		lambda1.0.975		<- apply(lambda1.fin, 2, quantile, prob = 0.975, na.rm = T)
		lambda1.0.025		<- apply(lambda1.fin, 2, quantile, prob = 0.025, na.rm = T)
		lambda2.0.975		<- apply(lambda2.fin, 2, quantile, prob = 0.975, na.rm = T)
		lambda2.0.025		<- apply(lambda2.fin, 2, quantile, prob = 0.025, na.rm = T)
		lambda3.0.975		<- apply(lambda3.fin, 2, quantile, prob = 0.975, na.rm = T)
		lambda3.0.025		<- apply(lambda3.fin, 2, quantile, prob = 0.025, na.rm = T)
						
		lb1	<- lambda1.mean - 1.96 * lambda1.sd
		ub1	<- lambda1.mean + 1.96 * lambda1.sd
		lb2	<- lambda2.mean - 1.96 * lambda2.sd
		ub2	<- lambda2.mean + 1.96 * lambda2.sd
		lb3	<- lambda3.mean - 1.96 * lambda3.sd
		ub3	<- lambda3.mean + 1.96 * lambda3.sd

		pdf(file = paste("mcmcOutcome/ch", chain, "baseHaz1.pdf", sep = ""), height = 10, width = 20)

		par(mfrow = c(3, 3))


		plot(ss1, lambda1.0.975, type = "s", col = "blue", lty = 2, ylim = c(min(lambda1.0.025), max(lambda1.0.975)), 			main = "Estimate of the log-baseline hazard1", xlab = "time", ylab = "log-baseline hazard")
		lines(ss1, lambda1.mean, type = "s", col = "red")
		lines(ss1, lb1, type = "s", col = "red", lty = 2)
		lines(ss1, ub1, type = "s", col = "red", lty = 2);
		lines(ss1, lambda1.med, type = "s", col = "blue")
		lines(ss1, lambda1.0.025, type = "s", col = "blue", lty = 2)
		lines(ss1, lambda1.0.975, type = "s", col = "blue", lty = 2)

		for(i in 1:8){
			plot(lambda1.fin[,floor(priorPara$s1max*i/8)], type = "l", 				main = paste("Log-baseline hazard at t = ", floor(priorPara$s1max*i/8)))
			}
			

		dev.off()
		
		
		pdf(file = paste("mcmcOutcome/ch", chain, "baseHaz2.pdf", sep = ""), height = 10, width = 20)

		par(mfrow = c(3, 3))

		plot(ss2, lambda2.0.975, type = "s", col = "blue", lty = 2, ylim = c(min(lambda2.0.025), max(lambda2.0.975)), 			main = "Estimate of the log-baseline hazard2", xlab = "time", ylab = "log-baseline hazard")
		lines(ss2, lambda2.mean, type = "s", col = "red")
		lines(ss2, lb2, type = "s", col = "red", lty = 2)
		lines(ss2, ub2, type = "s", col = "red", lty = 2);
		lines(ss2, lambda2.med, type = "s", col = "blue")
		lines(ss2, lambda2.0.025, type = "s", col = "blue", lty = 2)
				
			
		for(i in 1:8){
			plot(lambda2.fin[,floor(priorPara$s2max*i/8)], type = "l", 				main = paste("Log-baseline hazard at t = ", floor(priorPara$s2max*i/8)))
			}
			
		dev.off()
				

		pdf(file = paste("mcmcOutcome/ch", chain, "baseHaz3.pdf", sep = ""), height = 10, width = 20)

		par(mfrow = c(3, 3))

		plot(ss3, lambda3.0.975, type = "s", col = "blue", lty = 2, ylim = c(min(lambda3.0.025), max(lambda3.0.975)), 			main = "Estimate of the log-baseline hazard3", xlab = "time", ylab = "log-baseline hazard")
		lines(ss3, lambda3.mean, type = "s", col = "red")
		lines(ss3, lb3, type = "s", col = "red", lty = 2)
		lines(ss3, ub3, type = "s", col = "red", lty = 2);
		lines(ss3, lambda3.med, type = "s", col = "blue")
		lines(ss3, lambda3.0.025, type = "s", col = "blue", lty = 2)
			
		
		for(i in 1:8){
			plot(lambda3.fin[,floor(priorPara$s3max*i/8)], type = "l", 				main = paste("Log-baseline hazard at t = ", floor(priorPara$s3max*i/8)))
			}
			
		dev.off()
		
		
		pdf(file = paste("mcmcOutcome/ch", chain, "baseHazAll.pdf", sep = ""), height = 10, width = 30)

		par(mfrow = c(1, 3))

		plot(ss1, lambda1.0.975, type = "s", col = "blue", lty = 2, ylim = c(min(lambda1.0.025), max(lambda1.0.975)), 			main = "Estimate of the log-baseline hazard1", xlab = "time", ylab = "log-baseline hazard")
		lines(ss1, lambda1.mean, type = "s", col = "red")
		lines(ss1, lb1, type = "s", col = "red", lty = 2)
		lines(ss1, ub1, type = "s", col = "red", lty = 2);
		lines(ss1, lambda1.med, type = "s", col = "blue")
		lines(ss1, lambda1.0.025, type = "s", col = "blue", lty = 2)
		lines(ss1, lambda1.0.975, type = "s", col = "blue", lty = 2)

		plot(ss2, lambda2.0.975, type = "s", col = "blue", lty = 2, ylim = c(min(lambda2.0.025), max(lambda2.0.975)), 			main = "Estimate of the log-baseline hazard2", xlab = "time", ylab = "log-baseline hazard")
		lines(ss2, lambda2.mean, type = "s", col = "red")
		lines(ss2, lb2, type = "s", col = "red", lty = 2)
		lines(ss2, ub2, type = "s", col = "red", lty = 2);
		lines(ss2, lambda2.med, type = "s", col = "blue")
		lines(ss2, lambda2.0.025, type = "s", col = "blue", lty = 2)

		
		plot(ss3, lambda3.0.975, type = "s", col = "blue", lty = 2, ylim = c(min(lambda3.0.025), max(lambda3.0.975)), 			main = "Estimate of the log-baseline hazard3", xlab = "time", ylab = "log-baseline hazard")
		lines(ss3, lambda3.mean, type = "s", col = "red")
		lines(ss3, lb3, type = "s", col = "red", lty = 2)
		lines(ss3, ub3, type = "s", col = "red", lty = 2);
		lines(ss3, lambda3.med, type = "s", col = "blue")
		lines(ss3, lambda3.0.025, type = "s", col = "blue", lty = 2)

			
		dev.off()		
		
		
				
		}

		# store the mcmc outcomes

		if(M %% save == 0 | M == num.reps){
			
			if(!RJ){
				save(mcmcOutcome, file = paste("mcmcOutcome/BpeScrIntv",initial$J1 + 1,"clam",c.lam1,					"a", a1, "b", b1, "psi", priorPara$psi, "omega", priorPara$omega, "ch", chain, ".RData", sep = ""))
				}
			if(RJ){
				save(mcmcOutcome, file = paste("mcmcOutcome/BpeScrRjAlp",alpha1, "pert", delPert1, "clam",c.lam1,					"a", a1, "b", b1, "psi", priorPara$psi, "omega", priorPara$omega, "ch", chain, ".RData", sep = ""))
				}			
			}
	
		} # end of the 'for' loop for MCMC sampling
	

	return(mcmcOutcome)
	
	} # end of BpeScr function

