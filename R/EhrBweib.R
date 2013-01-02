

EhrBweib <- function(x, beta2.p, beta3.p, alpha2.p, alpha3.p, kappa2.p, kappa3.p, s){
	
	num.sample <- dim(beta2.p)[1]

	J	<- length(s) - 1
		
	ehrSubj		<- matrix(NA, num.sample, J+1)
	
	for(i in 1:num.sample){
		
		if(i %% 100 == 0){
			cat("i", i,fill=TRUE)
			fsh()
			}
		
		bhr	<- (alpha3.p[i]*kappa3.p[i]*s^(alpha3.p[i] - 1))/(alpha2.p[i]*kappa2.p[i]*s^(alpha2.p[i] - 1))
		coefDif <- (beta3.p[i, ] - beta2.p[i, ])
		
		ehrSubj[i,]	<- bhr * exp(as.vector(x %*% coefDif))
								
		} 
		
	return(ehrSubj)
			
	}



