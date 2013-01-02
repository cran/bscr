

setting.interval.BpeSurv	<- function(survObj, ini, priorPara){
		
	y		<- survObj$y	
	delta	<- survObj$delta
	n		<- survObj$n
	p		<- survObj$p	

	s		<- ini$s
	J		<- ini$J


	case0	<- survObj$case0
	case1	<- survObj$case1	
	
	case0yleq	<- survObj$case0yleq
	case0ygeq	<- survObj$case0ygeq
	case1yleq	<- survObj$case1yleq
	case1ygeq	<- survObj$case1ygeq


	# Define the indicator matrices for risk sets and failure sets 

	ind.d <- ind.d.temp <- ind.r <- matrix(0, n, J+1)

	for(i in case1yleq){
		d.mat.ind	<- min(which(s - y[i] >=0))
		ind.d[i, d.mat.ind]	<- 1
		ind.d.temp[i, d.mat.ind]	<- 1		
		ind.r[i,1:d.mat.ind] <- 1
		}
		
	
#	if(J > 0){	
		for(i in case0yleq){
			cen.j <- min(which(s - y[i] >=0))
			ind.d.temp[i, cen.j]	<- 1
			ind.r[i, 1:cen.j] <- 1
			}		
#		}
			
	
	if(length(union(case1ygeq, case0ygeq)) > 0){
		ind.r[union(case1ygeq, case0ygeq),]	<- 1
		}

			
#	ind.r[,1]	<- 1		
		
	d	<- colSums(ind.d)
	
	list(ind.r = ind.r, ind.d = ind.d, ind.d.temp = ind.d.temp, d = d)
	}









calDelta.BpeSurv	<- function(survObj, ini){

	y		<- survObj$y
	delta	<- survObj$delta
	s		<- ini$s
	J		<- ini$J	
	ind.r	<- ini$ind.r
	ind.d	<- ini$ind.d			
	ind.d.temp	<- ini$ind.d.temp	
	
	n		<- survObj$n
	p		<- survObj$p		
		
	Delta	<- matrix(0, n, J+1);
	
	case0	<- survObj$case0
	case1	<- survObj$case1	
	
	case0yleq	<- survObj$case0yleq
	case0ygeq	<- survObj$case0ygeq
	case1yleq	<- survObj$case1yleq
	case1ygeq	<- survObj$case1ygeq	
	
	for(i in union(case1yleq, case0yleq)){
		
		y.j	<- min(which(s - y[i] >=0))
				
		if(y.j != 1){
			Delta[i,y.j] <- y[i] - s[y.j - 1]
			}

		if(y.j == 1){
			Delta[i, y.j] <- y[i] - 0
			}			
		}
		
		
	val	<- (matrix(rep(diff(c(0, s)), n), n, J+1, byrow = T) * (ind.r - ind.d.temp) + Delta) * ind.r;
	
	return(val)
	
	
	}







cal.Sigma.lam.BpeSurv <- function(ini, priorPara){
	
	s	<- ini$s
	J	<- ini$J
	c.lam	<- priorPara$c.lam
	
	Del.s	<- diff(c(0, s))
	
	if(J+1 >= 3){
	
		W	<- Q	<- matrix(0, J+1, J+1);
	
		W[1, 2]	<- c.lam * (Del.s[1] + Del.s[2]) / (2*Del.s[1] + Del.s[2]);
		W[J+1, J]	<- c.lam * (Del.s[J] + Del.s[J+1]) / (Del.s[J] + 2*Del.s[J+1]);
		Q[1, 1]		<- 2 / (2*Del.s[1] + Del.s[2]);
		Q[J+1, J+1]		<- 2 / (Del.s[J] + 2*Del.s[J+1]);
	
		for(i in 2:(J)){
			Q[i, i]		<- 2 / (Del.s[i-1] + 2*Del.s[i] + Del.s[i+1]);
			W[i, i-1]	<- c.lam * (Del.s[i-1] + Del.s[i]) / (Del.s[i-1] + 2*Del.s[i] + Del.s[i+1]);
			W[i, i+1]	<- c.lam * (Del.s[i] + Del.s[i+1]) / (Del.s[i-1] + 2*Del.s[i] + Del.s[i+1]);			}

		Sigma.lam	<- solve(diag(1, J+1) - W) %*% Q;
		
		}	
		
	if(J+1 == 2){
	
		W	<- Q	<- matrix(0, J+1, J+1);
		W[1, 2]	<- c.lam * (Del.s[1] + Del.s[2]) / (2*Del.s[1] + Del.s[2]);
		W[J+1, J]	<- c.lam * (Del.s[J] + Del.s[J+1]) / (Del.s[J] + 2*Del.s[J+1]);
		Q[1, 1]		<- 2 / (2*Del.s[1] + Del.s[2]);
		Q[J+1, J+1]		<- 2 / (Del.s[J] + 2*Del.s[J+1]);
		
		Sigma.lam	<- solve(diag(1, J+1) - W) %*% Q;
		
		}	
	
	if(J+1 == 1){
	
		W	<- 0
		Q	<- 2 / (2*Del.s[1]);
		
		Sigma.lam	<- Q
		
		}	

	list(Sigma.lam = Sigma.lam, W = W, Q = Q);
	
	}	

















add.interval.BpeSurv	<- function(survObj, ini, priorPara, old, s.new){

 	y		<- survObj$y
 	delta	<- survObj$delta
	n		<- survObj$n
	p		<- survObj$p	
 	
 	
 	case0	<- survObj$case0
 	case1	<- survObj$case1	
 	
 	case0yleq	<- survObj$case0yleq
 	case0ygeq	<- survObj$case0ygeq
 	case1yleq	<- survObj$case1yleq
 	case1ygeq	<- survObj$case1ygeq
 	
 	s.max 	<- priorPara$smax	
 	
 
 	j.old 	<- old$j.old
 	s.old	<- old$s.old
 	J.old   <- old$J.old
 	ind.r.old <- old$ind.r.old
 	ind.d.old <- old$ind.d.old
 	ind.d.temp.old <- old$ind.d.temp.old
 	d.old	<- old$d.old
 
 		
 	J		<- J.old
 	ind.r <- cbind(ind.r.old, rep(0, n))
 	ind.d <- cbind(ind.d.old, rep(0, n))
  	ind.d.temp <- cbind(ind.d.temp.old, rep(0, n))
  	d	<- c(d.old, 0)
 
 								
 	if(j.old != (J + 1)){
 		ind.r[,(j.old+1):(J+2)] <- ind.r[,c(J+2, (j.old+1):(J+1))]
 		ind.d[,(j.old+1):(J+2)] <- ind.d[,c(J+2, (j.old+1):(J+1))]
  		ind.d.temp[,(j.old+1):(J+2)] <- ind.d.temp[,c(J+2, (j.old+1):(J+1))]
  		d[(j.old+1):(J+2)] <- d[c(J+2, (j.old+1):(J+1))]
 		}
 		
 	
 	# Define the indicator matrices for risk sets and failure sets 
 		
 	for(i in case1yleq){
 		if(ind.d.old[i,j.old] == 1 & y[i] > s.new[j.old]){
 			ind.d[i,c(j.old, j.old+1)] <- c(0, 1)
  			ind.d.temp[i,c(j.old, j.old+1)] <- c(0, 1)
  			}
  			
 		if(y[i] > s.new[j.old]){
 			ind.r[i,c(j.old, j.old+1)] <- c(1, 1)
 			}				
 		}	
 		
 	for(i in case0yleq){
 		if(ind.d.temp.old[i,j.old] == 1 & y[i] > s.new[j.old]){
  			ind.d.temp[i,c(j.old, j.old+1)] <- c(0, 1)
  			}
  			
 		if(y[i] > s.new[j.old]){
 			ind.r[i,c(j.old, j.old+1)] <- c(1, 1)
 			}				
 		}	 		
 	
	if(length(union(case1ygeq, case0ygeq)) > 0){
		ind.r[union(case1ygeq, case0ygeq),c(j.old, j.old+1)] <- 1
		} 	
 	
 				
 	d[c(j.old, j.old+1)]	<- colSums(ind.d[,c(j.old, j.old+1)]);
 		
 	
 	list(ind.r = ind.r, ind.d = ind.d, ind.d.temp = ind.d.temp, d = d);
 	
 	}
 	
 	
 	
 	 		
 	
 	
 	



calDelta.BI.BpeSurv	<- function(survObj, ini.new, old){

 	y		<- survObj$y
 	delta	<- survObj$delta
	n		<- survObj$n
	p		<- survObj$p	
 	
 	
 	case0	<- survObj$case0
 	case1	<- survObj$case1	
 	
 	case0yleq	<- survObj$case0yleq
 	case0ygeq	<- survObj$case0ygeq
 	case1yleq	<- survObj$case1yleq
 	case1ygeq	<- survObj$case1ygeq
 	
 	s	<- ini.new$s
	J	<- ini.new$J
	ind.r	<- ini.new$ind.r
	ind.d	<- ini.new$ind.d
	ind.d.temp	<- ini.new$ind.d.temp
			
 	j.old 	<- old$j.old
 	s.old	<- old$s.old
 	J.old   <- old$J.old
 	ind.r.old <- old$ind.r.old
 	ind.d.old <- old$ind.d.old
 	ind.d.temp.old <- old$ind.d.temp.old 	
 	d.old	<- old$d.old
		
	Delta.old <- ini.new$Delta
		
	Delta <- cbind(Delta.old, rep(0, n))
		
	diff.s <- diff(c(0, s))
	diff.s.old <- diff(c(0, s.old))
				
	if(j.old != (J.old + 1)){
		Delta[,(j.old+1):(J.old+2)] <- Delta[,c(J.old+2, (j.old+1):(J.old+1))]
		}		
			
#	for(i in case0){
#		Delta[i,] <- diff.s
#		}	

	for(i in intersect(union(case1yleq, case0yleq), which(Delta.old[,j.old] != 0))){
		if(j.old == 1){
			Delta[i, j.old] <- max(0, min(y[i], s[j.old]) - 0)
			}
		if(j.old != 1){
			Delta[i, j.old] <- max(0, min(y[i], s[j.old]) - s[j.old-1])
			}			
		Delta[i, j.old+1] <- max(0, min(y[i], s[j.old+1]) - s[j.old])
		}
		
		
	if(length(union(case1ygeq, case0ygeq)) > 0){
		for(i in union(case1ygeq, case0ygeq)){
			Delta[i, c(j.old, j.old +1)] <- diff.s[c(j.old, j.old +1)]
			}
		} 
		
				
	return(Delta);
	
	}	
	
	
	
	 	
 	
 	
 	
 	
 	
 
 
remove.interval.BpeSurv	<- function(survObj, ini, priorPara, old, s.new){

	y		<- survObj$y
 	delta	<- survObj$delta
	n		<- survObj$n
	p		<- survObj$p	
 	
 	
 	case0	<- survObj$case0
 	case1	<- survObj$case1	
 	
 	case0yleq	<- survObj$case0yleq
 	case0ygeq	<- survObj$case0ygeq
 	case1yleq	<- survObj$case1yleq
 	case1ygeq	<- survObj$case1ygeq
 	
	s.max 	<- priorPara$smax	
 	
 
 	j.old 	<- old$j.old
 	s.old	<- old$s.old
 	J.old   <- old$J.old
 	ind.r.old <- old$ind.r.old
 	ind.d.old <- old$ind.d.old
 	ind.d.temp.old <- old$ind.d.temp.old 	
 	d.old	<- old$d.old
		
	J	<- J.old
		
	if(J > 1){
		ind.r <- ind.r.old[,-j.old]
		ind.d <- ind.d.old[,-j.old]
		ind.d.temp <- ind.d.temp.old[,-j.old]	
			
		d	<- d.old[-j.old]
		ind.r[,j.old]	<- as.numeric(ind.r.old[,j.old] | ind.r.old[,j.old+1])
		ind.d[,j.old]	<- as.numeric(ind.d.old[,j.old] | ind.d.old[,j.old+1])
		ind.d.temp[,j.old]	<- as.numeric(ind.d.temp.old[,j.old] | ind.d.temp.old[,j.old+1])
		d[j.old]		<- sum(d.old[c(j.old, j.old + 1)])
		}
			
	if(J == 1){
		ind.r <- matrix(1, nrow = n, ncol = 1)
		ind.d <- matrix(NA, nrow = n, ncol = 1)
		ind.d.temp <- matrix(NA, nrow = n, ncol = 1)
				
		d	<- d.old[-j.old]

		ind.d[,j.old] <- as.numeric(ind.d.old[,j.old] | ind.d.old[,j.old+1])
		ind.d.temp[,j.old] <- as.numeric(ind.d.temp.old[,j.old] | ind.d.temp.old[,j.old+1])
		d[j.old]	<- sum(d.old[c(j.old, j.old + 1)])

		}	
						
	list(ind.r = ind.r, ind.d = ind.d, ind.d.temp = ind.d.temp, d = d);
	
	}
	

	
		

	
cal.Del.DI.BpeSurv	<- function(survObj, ini.new, old){

	y		<- survObj$y
 	delta	<- survObj$delta
	n		<- survObj$n
	p		<- survObj$p	
 	
 	
 	case0	<- survObj$case0
 	case1	<- survObj$case1	
 	
 	case0yleq	<- survObj$case0yleq
 	case0ygeq	<- survObj$case0ygeq
 	case1yleq	<- survObj$case1yleq
 	case1ygeq	<- survObj$case1ygeq
	
	Delta.old <- ini.new$Delta
	J.new	<- ini.new$J

		

	j.old <- old$j.old
	
	if(J.new != 0){
		Delta <- Delta.old[,-j.old]
		Delta[,j.old] <- rowSums(Delta.old[,c(j.old, j.old+1)])
		}		
			
	if(J.new == 0){
		Delta <- as.matrix(Delta.old[,-j.old], nrow = n)
		Delta[,j.old] <- rowSums(Delta.old[,c(j.old, j.old+1)])
		}	
	

	return(Delta);
	
	}	
	
		
	