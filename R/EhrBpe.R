

EhrBpe <- function(x, beta2.p, beta3.p, lambda2.p, lambda3.p, s2.p, s3.p, s){
	
	num.sample <- dim(beta2.p)[1]

	J	<- length(s) - 1
	
	lambda2.fin	<- matrix(NA, num.sample, J+1)
	lambda3.fin	<- matrix(NA, num.sample, J+1)
	
	ind2.lam	<- rep(NA, J+1)
	ind3.lam	<- rep(NA, J+1)
	
	ehrSubj		<- matrix(NA, num.sample, J+1)
	
	for(i in 1:num.sample){
		
		if(i %% 100 == 0){
			cat("i", i,fill=TRUE)
			fsh()
			}
		
		s2.n	<- length(s2.p[[i]])
		s3.n	<- length(s3.p[[i]])
				
		ind2.lam[(s <= s2.p[[i]][1]) & (s > 0)] <- 1
		ind3.lam[(s <= s3.p[[i]][1]) & (s > 0)] <- 1

		for(j in 2:s2.n){
			ind2.lam[(s <= s2.p[[i]][j]) & (s > s2.p[[i]][j-1])] <- j
			}
		for(j in 2:s3.n){
			ind3.lam[(s <= s3.p[[i]][j]) & (s > s3.p[[i]][j-1])] <- j
			}	
		
		lambda2.fin[i, ] <- lambda2.p[[i]][ind2.lam]
		lambda3.fin[i, ] <- lambda3.p[[i]][ind3.lam]
		
		bhr	<- exp(lambda3.fin[i, ] - lambda2.fin[i, ]) 
		coefDif <- (beta3.p[i, ] - beta2.p[i, ])
		
		ehrSubj[i,]	<- bhr * exp(as.vector(x %*% coefDif))
				
					
		} # end of 'for' loop
		
		return(ehrSubj)
			
	}



