setting.interval.BpeScr <- function(survObj, ini, priorPara, g){

	y1		<- survObj$y1
	y2		<- survObj$y2	
	delta1	<- survObj$delta1
	delta2	<- survObj$delta2
	n		<- survObj$n

	case00	<- survObj$case00
	
	ind.r1	<- ind.r2	<- NULL
	
	if(g == 1){
		s	<- ini$s1
		J	<- ini$J1
		s.max <- priorPara$s1max
		
		case1_y1leq 			<- survObj$case1_y1leq.s1
		case01y2leq				<- survObj$case01y2leq.s1
		case00yleq				<- survObj$case00yleq.s1
		case1_or01or00y1y2leq 	<- survObj$case1_or01or00y1y2leq.s1
		case1_or01or00y1y2gre 	<- survObj$case1_or01or00y1y2gre.s1
	
		# Define the indicator matrices for risk sets and failure sets 

		ind.d1	<-	matrix(0, n, J+1);
		ind.d2	<-	matrix(0, n, J+1);
		ind.d3	<-	matrix(0, n, J+1);
		ind.r1	<-	matrix(0, n, J+1);
		ind.d.temp	<-	matrix(0, n, J+1);		
			
		
		for(i in case1_y1leq){
			d1.mat.ind	<- min(which(s - y1[i] >=0))
			ind.d1[i, d1.mat.ind]	<- 1
			ind.d.temp[i, d1.mat.ind]	<- 1
			ind.r1[i, 1:d1.mat.ind]	<- 1	
			}	
					
		for(i in case01y2leq){
			d2.mat.ind	<- min(which(s - y2[i] >=0))
			ind.d2[i, d2.mat.ind]	<- 1
			ind.d.temp[i, d2.mat.ind]	<- 1
			ind.r1[i, 1:d2.mat.ind]	<- 1
			}
			
		for(i in case00yleq){
			cen.j	<- min(which(s - y1[i] >=0))
			ind.d.temp[i, cen.j]	<- 1
			ind.r1[i, 1:cen.j]		<- 1
			}			
	
			
	if(length(case1_or01or00y1y2gre) > 0){
		ind.r1[case1_or01or00y1y2gre,]	<- 1
		}
				
		d1	<- colSums(ind.d1);
		d2	<- colSums(ind.d2);	
		d3	<- colSums(ind.d3);			
		}
		
	
	if(g == 2){		
		s	<- ini$s2
		J	<- ini$J2
		s.max <- priorPara$s2max
		
		case1_y1leq 			<- survObj$case1_y1leq.s2
		case01y2leq				<- survObj$case01y2leq.s2
		case00yleq				<- survObj$case00yleq.s2
		case1_or01or00y1y2leq 	<- survObj$case1_or01or00y1y2leq.s2
		case1_or01or00y1y2gre 	<- survObj$case1_or01or00y1y2gre.s2
				
		# Define the indicator matrices for risk sets and failure sets 

		ind.d1	<-	matrix(0, n, J+1);
		ind.d2	<-	matrix(0, n, J+1);
		ind.d3	<-	matrix(0, n, J+1);	
		ind.r1	<-	matrix(0, n, J+1);
		ind.d.temp	<-	matrix(0, n, J+1);	
				
		for(i in case1_y1leq){
			d1.mat.ind	<- min(which(s - y1[i] >=0))
			ind.d1[i, d1.mat.ind]	<- 1
			ind.d.temp[i, d1.mat.ind]	<- 1
			ind.r1[i, 1:d1.mat.ind]	<- 1	
			}	
					
		for(i in case01y2leq){
			d2.mat.ind	<- min(which(s - y2[i] >=0))
			ind.d2[i, d2.mat.ind]	<- 1
			ind.d.temp[i, d2.mat.ind]	<- 1
			ind.r1[i, 1:d2.mat.ind]	<- 1
			}
			
		for(i in case00yleq){
			cen.j	<- min(which(s - y1[i] >=0))
			ind.d.temp[i, cen.j]	<- 1
			ind.r1[i, 1:cen.j]		<- 1
			}			
	
			
		if(length(case1_or01or00y1y2gre) > 0){
			ind.r1[case1_or01or00y1y2gre,]	<- 1
			}
				
		d1	<- colSums(ind.d1);
		d2	<- colSums(ind.d2);	
		d3	<- colSums(ind.d3);			
					
		}
	
	if(g == 3){
		s	<- ini$s3
		J	<- ini$J3	
		s.max <- priorPara$s3max		
		
		case1_y1leq		<- survObj$case1_y1leq.s3
		case01y2leq		<- survObj$case01y2leq.s3
		case00yleq		<- survObj$case00yleq.s3
		case1_or01or00y1y2leq	<- survObj$case1_or01or00y1y2leq.s3
		case1_or01or00y1y2gre	<- survObj$case1_or01or00y1y2gre.s3
		
		case11y2leq		<- survObj$case11y2leq.s3
		case10y2leq		<- survObj$case10y2leq.s3
		case11or10y2leq	<- survObj$case11or10y2leq.s3
		case11or10y2gre	<- survObj$case11or10y2gre.s3
								
		# Define the indicator matrices for risk sets and failure sets 

		ind.d1	<-	matrix(0, n, J+1);
		ind.d2	<-	matrix(0, n, J+1);
		ind.d3	<-	matrix(0, n, J+1);	
		ind.r2	<-	matrix(0, n, J+1);
		ind.d.temp	<-	matrix(0, n, J+1);	
		
		for(i in case1_y1leq){
			d1.mat.ind	<- min(which(s - y1[i] >=0))
			ind.d1[i, d1.mat.ind]	<- 1
			}
	
		for(i in case11y2leq){
			d3.mat.ind	<- min(which(s - y2[i] >=0))
			ind.d3[i, d3.mat.ind]	<- 1
			ind.d.temp[i, d3.mat.ind]	<- 1
			ind.r2[i, which(ind.d1[i,] == 1):d3.mat.ind] <- 1
			}

		for(i in case10y2leq){
			cen.j	<- min(which(s - y2[i] >=0))
			ind.d.temp[i, cen.j]	<- 1
			ind.r2[i, which(ind.d1[i,] == 1):cen.j] <- 1
			}	
			
		if(length(case11or10y2gre) > 0){
			for(i in case11or10y2gre){
				ind.r2[i,which(ind.d1[i,] == 1):(J+1)] <- 1
				}
			}			
					

			d1	<- colSums(ind.d1);
			d2	<- colSums(ind.d2);	
			d3	<- colSums(ind.d3);		
		}				
				
	
	list(ind.r1 = ind.r1, ind.r2 = ind.r2, ind.d1 = ind.d1, ind.d2 = ind.d2, ind.d3 = ind.d3,		ind.d.temp = ind.d.temp, d1 = d1, d2 = d2, d3 = d3);
	
	}






add.interval.BpeScr <-
function(survObj, ini, priorPara, old, s.new, g = 1){

	y1		<- survObj$y1
	y2		<- survObj$y2	
	delta1	<- survObj$delta1
	delta2	<- survObj$delta2
	n		<- survObj$n

	case00	<- survObj$case00
	
	if(g == 1){
		j.old 	<- old$j.old
		s.old	<- old$s.old
		J.old   <- old$J.old
		ind.r1.old <- old$ind.r1.old
		ind.r2.old <- old$ind.r2.old
		ind.d1.old <- old$ind.d1.old
		ind.d2.old <- old$ind.d2.old
		ind.d3.old <- old$ind.d3.old
		ind.d.temp.old <- old$ind.d.temp.old
			
		d1.old	<- old$d1.old
		d2.old	<- old$d2.old
		
		J	<- J.old
		ind.r1 <- cbind(ind.r1.old, rep(0, n))
		ind.r2 <- NULL
		ind.d1 <- cbind(ind.d1.old, rep(0, n))
		ind.d2 <- cbind(ind.d2.old, rep(0, n))
		ind.d3 <- NULL
		ind.d.temp <- cbind(ind.d.temp.old, rep(0, n))
							
		d1	<- c(d1.old, 0)
		d2	<- c(d2.old, 0)
		d3	<- NULL						
								
		if(j.old != (J + 1)){
			ind.r1[,(j.old+1):(J+2)] <- ind.r1[,c(J+2, (j.old+1):(J+1))]
			ind.d1[,(j.old+1):(J+2)] <- ind.d1[,c(J+2, (j.old+1):(J+1))]
			ind.d2[,(j.old+1):(J+2)] <- ind.d2[,c(J+2, (j.old+1):(J+1))]
			ind.d.temp[,(j.old+1):(J+2)] <- ind.d.temp[,c(J+2, (j.old+1):(J+1))]
			d1[(j.old+1):(J+2)] <- d1[c(J+2, (j.old+1):(J+1))]
			d2[(j.old+1):(J+2)] <- d2[c(J+2, (j.old+1):(J+1))]
			}
		
		s.max <- priorPara$s1max
				
		case1_y1leq 			<- survObj$case1_y1leq.s1
		case01y2leq				<- survObj$case01y2leq.s1
		case00yleq				<- survObj$case00yleq.s1
		case1_or01or00y1y2leq 	<- survObj$case1_or01or00y1y2leq.s1
		case1_or01or00y1y2gre 	<- survObj$case1_or01or00y1y2gre.s1
				
	
		# Define the indicator matrices for risk sets and failure sets 
				
		
		for(i in case1_y1leq){
			if(ind.d.temp.old[i,j.old] == 1 & y1[i] > s.new[j.old]){
				ind.d1[i,c(j.old, j.old+1)] <- c(0, 1)
				ind.d.temp[i,c(j.old, j.old+1)] <- c(0, 1)
				}
			if(y1[i] > s.new[j.old]){
				ind.r1[i,c(j.old, j.old+1)] <- c(1, 1)
				}				
			}	
					
		for(i in case01y2leq){
			if(ind.d.temp.old[i,j.old] == 1 & y2[i] > s.new[j.old]){
				ind.d2[i,c(j.old, j.old+1)] <- c(0, 1)
				ind.d.temp[i,c(j.old, j.old+1)] <- c(0, 1)
				}					
			if(y2[i] > s.new[j.old]){
				ind.r1[i,c(j.old, j.old+1)] <- c(1, 1)
				}							
			}
			
		for(i in case00yleq){	
			if(ind.d.temp.old[i,j.old] == 1 & y1[i] > s.new[j.old]){
				ind.d.temp[i,c(j.old, j.old+1)] <- c(0, 1)
				}							
			if(y1[i] > s.new[j.old]){
				ind.r1[i,c(j.old, j.old+1)] <- c(1, 1)
				}							
			}			
			
		
		ind.r1[case1_or01or00y1y2gre, j.old+1] <- 1
				
		d1[c(j.old, j.old+1)]	<- colSums(ind.d1[,c(j.old, j.old+1)]);
		d2[c(j.old, j.old+1)]	<- colSums(ind.d2[,c(j.old, j.old+1)]);
		}
		
	
	if(g == 2){
		j.old 	<- old$j.old
		s.old	<- old$s.old
		J.old   <- old$J.old
		ind.r1.old <- old$ind.r1.old
		ind.r2.old <- old$ind.r2.old
		ind.d1.old <- old$ind.d1.old
		ind.d2.old <- old$ind.d2.old
		ind.d3.old <- old$ind.d3.old
		ind.d.temp.old <- old$ind.d.temp.old
				
		d1.old	<- old$d1.old
		d2.old	<- old$d2.old
		
		J	<- J.old
		ind.r1 <- cbind(ind.r1.old, rep(0, n))
		ind.r2 <- NULL		
		ind.d1 <- cbind(ind.d1.old, rep(0, n))
		ind.d2 <- cbind(ind.d2.old, rep(0, n))
		ind.d3 <- NULL
		ind.d.temp <- cbind(ind.d.temp.old, rep(0, n))
		
		d1	<- c(d1.old, 0)
		d2	<- c(d2.old, 0)
		d3	<- NULL
		
										
								
		if(j.old != (J + 1)){
			ind.r1[,(j.old+1):(J+2)] <- ind.r1[,c(J+2, (j.old+1):(J+1))]
			ind.d1[,(j.old+1):(J+2)] <- ind.d1[,c(J+2, (j.old+1):(J+1))]
			ind.d2[,(j.old+1):(J+2)] <- ind.d2[,c(J+2, (j.old+1):(J+1))]
			ind.d.temp[,(j.old+1):(J+2)] <- ind.d.temp[,c(J+2, (j.old+1):(J+1))]
			d1[(j.old+1):(J+2)] <- d1[c(J+2, (j.old+1):(J+1))]
			d2[(j.old+1):(J+2)] <- d2[c(J+2, (j.old+1):(J+1))]
			}

		s.max <- priorPara$s2max
		
		case1_y1leq 			<- survObj$case1_y1leq.s2
		case01y2leq				<- survObj$case01y2leq.s2
		case00yleq				<- survObj$case00yleq.s2
		case1_or01or00y1y2leq 	<- survObj$case1_or01or00y1y2leq.s2
		case1_or01or00y1y2gre 	<- survObj$case1_or01or00y1y2gre.s2
				
		# Define the indicator matrices for risk sets and failure sets 
				
		
		for(i in case1_y1leq){
			if(ind.d.temp.old[i,j.old] == 1 & y1[i] > s.new[j.old]){
				ind.d1[i,c(j.old, j.old+1)] <- c(0, 1)
				ind.d.temp[i,c(j.old, j.old+1)] <- c(0, 1)
				}
			if(y1[i] > s.new[j.old]){
				ind.r1[i,c(j.old, j.old+1)] <- c(1, 1)
				}				
			}	
					
		for(i in case01y2leq){
			if(ind.d.temp.old[i,j.old] == 1 & y2[i] > s.new[j.old]){
				ind.d2[i,c(j.old, j.old+1)] <- c(0, 1)
				ind.d.temp[i,c(j.old, j.old+1)] <- c(0, 1)
				}					
			if(y2[i] > s.new[j.old]){
				ind.r1[i,c(j.old, j.old+1)] <- c(1, 1)
				}							
			}
			
		for(i in case00yleq){	
			if(ind.d.temp.old[i,j.old] == 1 & y1[i] > s.new[j.old]){
				ind.d.temp[i,c(j.old, j.old+1)] <- c(0, 1)
				}							
			if(y1[i] > s.new[j.old]){
				ind.r1[i,c(j.old, j.old+1)] <- c(1, 1)
				}							
			}			
			
		
		ind.r1[case1_or01or00y1y2gre, j.old+1] <- 1
				
		d1[c(j.old, j.old+1)]	<- colSums(ind.d1[,c(j.old, j.old+1)]);
		d2[c(j.old, j.old+1)]	<- colSums(ind.d2[,c(j.old, j.old+1)]);
								
		}
	
	if(g == 3){
		j.old 	<- old$j.old
		s.old	<- old$s.old
		J.old   <- old$J.old
		ind.r1.old <- old$ind.r1.old
		ind.r2.old <- old$ind.r2.old
		ind.d1.old <- old$ind.d1.old
		ind.d2.old <- old$ind.d2.old
		ind.d3.old <- old$ind.d3.old
		ind.d.temp.old <- old$ind.d.temp.old
				
		d1.old	<- old$d1.old
		d3.old	<- old$d3.old
		
		J	<- J.old
		ind.r1 <- NULL
		ind.r2 <- cbind(ind.r2.old, rep(0, n))
		ind.d1 <- cbind(ind.d1.old, rep(0, n))
		ind.d2 <- NULL
		ind.d3 <- cbind(ind.d3.old, rep(0, n))
		ind.d.temp <- cbind(ind.d.temp.old, rep(0, n))
				
		d1	<- c(d1.old, 0)
		d2	<- NULL
		d3	<- c(d3.old, 0)
								
		if(j.old != (J + 1)){
			ind.r2[,(j.old+1):(J+2)] <- ind.r2[,c(J+2, (j.old+1):(J+1))]
			ind.d1[,(j.old+1):(J+2)] <- ind.d1[,c(J+2, (j.old+1):(J+1))]
			ind.d3[,(j.old+1):(J+2)] <- ind.d3[,c(J+2, (j.old+1):(J+1))]
			ind.d.temp[,(j.old+1):(J+2)] <- ind.d.temp[,c(J+2, (j.old+1):(J+1))]
			d1[(j.old+1):(J+2)] <- d1[c(J+2, (j.old+1):(J+1))]
			d3[(j.old+1):(J+2)] <- d3[c(J+2, (j.old+1):(J+1))]
			}

	
		s.max <- priorPara$s3max		
		
		case1_y1leq		<- survObj$case1_y1leq.s3
		case01y2leq		<- survObj$case01y2leq.s3
		case00yleq		<- survObj$case00yleq.s3
		case1_or01or00y1y2leq	<- survObj$case1_or01or00y1y2leq.s3
		case1_or01or00y1y2gre	<- survObj$case1_or01or00y1y2gre.s3
		
		case11y2leq		<- survObj$case11y2leq.s3
		case10y2leq		<- survObj$case10y2leq.s3
		case11or10y2leq	<- survObj$case11or10y2leq.s3
		case11or10y2gre	<- survObj$case11or10y2gre.s3
						
		# Define the indicator matrices for risk sets and failure sets 
		
		for(i in case1_y1leq){
			if(ind.d1.old[i,j.old] == 1 & y1[i] > s.new[j.old]){
				ind.d1[i,c(j.old, j.old+1)] <- c(0, 1)
				}
			}
	
				
		for(i in case10y2leq){
			if(ind.r2.old[i, j.old] == 1 & y1[i] <= s.new[j.old] & y2[i] > s.new[j.old]){
				ind.r2[i,c(j.old, j.old+1)] <- c(1, 1)
				}	
			if(ind.r2.old[i, j.old] == 1 & y1[i] > s.new[j.old]){
				ind.r2[i,c(j.old, j.old+1)] <- c(0, 1)
				}
			if(ind.d.temp.old[i,j.old] == 1 & y2[i] > s.new[j.old]){
				ind.d.temp[i,c(j.old, j.old+1)] <- c(0, 1)
				}				
			}
				
		for(i in case11y2leq){
			if(ind.r2.old[i, j.old] == 1 & y1[i] <= s.new[j.old] & y2[i] > s.new[j.old]){
				ind.r2[i,c(j.old, j.old+1)] <- c(1, 1)
				}	
			if(ind.r2.old[i, j.old] == 1 & y1[i] > s.new[j.old]){
				ind.r2[i, j.old] <- 0
				ind.r2[i, j.old+1] <- 1
				}
			if(ind.d3.old[i,j.old] == 1 & y2[i] > s.new[j.old]){
				ind.d3[i,c(j.old, j.old+1)] <- c(0, 1)
				ind.d.temp[i,c(j.old, j.old+1)] <- c(0, 1)
				}				
			}		
	
		if(length(case11or10y2gre) > 0){
			for(i in case11or10y2gre){
				if(ind.r2.old[i, j.old] == 1 & y1[i] <= s.new[j.old]){
					ind.r2[i,c(j.old, j.old+1)] <- c(1, 1)
					}
				if(ind.r2.old[i, j.old] == 1 & y1[i] > s.new[j.old]){
					ind.r2[i,c(j.old, j.old+1)] <- c(0, 1)
					}
				}
			}	
				
		d1[c(j.old, j.old+1)]	<- colSums(ind.d1[,c(j.old, j.old+1)]);
		d3[c(j.old, j.old+1)]	<- colSums(ind.d3[,c(j.old, j.old+1)]);

		
		}				
				
	
	list(ind.r1 = ind.r1, ind.r2 = ind.r2, ind.d1 = ind.d1, ind.d2 = ind.d2, ind.d3 = ind.d3, 		ind.d.temp = ind.d.temp, d1 = d1, d2 = d2, d3 = d3);
	
	}






remove.interval.BpeScr <-
function(survObj, ini, priorPara, old, s.new, g = 1){

	y1		<- survObj$y1
	y2		<- survObj$y2	
	delta1	<- survObj$delta1
	delta2	<- survObj$delta2
	n		<- survObj$n

	case00	<- survObj$case00
	
	if(g == 1){
		j.old 	<- old$j.old
		s.old	<- old$s.old
		J.old   <- old$J.old
		ind.r1.old <- old$ind.r1.old
		ind.r2.old <- old$ind.r2.old
		ind.d1.old <- old$ind.d1.old
		ind.d2.old <- old$ind.d2.old
		ind.d3.old <- old$ind.d3.old
		ind.d.temp.old	<- old$ind.d.temp.old
		d1.old	<- old$d1.old
		d2.old	<- old$d2.old
		
		J	<- J.old
		
		if(J > 1){
			ind.r1 <- ind.r1.old[,-j.old]
			ind.r2 <- NULL
			ind.d1 <- ind.d1.old[,-j.old]
			ind.d2 <- ind.d2.old[,-j.old]
			ind.d3 <- NULL
			ind.d.temp <- ind.d.temp.old[,-j.old]
									
			d1	<- d1.old[-j.old]
			d2	<- d2.old[-j.old]
			d3	<- NULL		
			ind.r1[,j.old] <- as.numeric(ind.r1.old[,j.old] | ind.r1.old[,j.old+1])
			ind.d1[,j.old] <- as.numeric(ind.d1.old[,j.old] | ind.d1.old[,j.old+1])
			ind.d2[,j.old] <- as.numeric(ind.d2.old[,j.old] | ind.d2.old[,j.old+1])
			ind.d.temp[,j.old] <- as.numeric(ind.d.temp.old[,j.old] | ind.d.temp.old[,j.old+1])
						
			d1[j.old]	<- sum(d1.old[c(j.old, j.old + 1)])
			d2[j.old]	<- sum(d2.old[c(j.old, j.old + 1)])
			}
			
		if(J == 1){
			ind.r1 <- matrix(1, nrow = n, ncol = 1)
			ind.r2 <- NULL
			ind.d1 <- matrix(NA, nrow = n, ncol = 1)
			ind.d2 <- matrix(NA, nrow = n, ncol = 1)
			ind.d3 <- NULL
			ind.d.temp <- matrix(NA, nrow = n, ncol = 1)
									
			d1	<- d1.old[-j.old]
			d2	<- d2.old[-j.old]
			d3	<- NULL		
			ind.d1[,j.old] <- as.numeric(ind.d1.old[,j.old] | ind.d1.old[,j.old+1])
			ind.d2[,j.old] <- as.numeric(ind.d2.old[,j.old] | ind.d2.old[,j.old+1])
			ind.d.temp[,j.old] <- as.numeric(ind.d.temp.old[,j.old] | ind.d.temp.old[,j.old+1])
						
			d1[j.old]	<- sum(d1.old[c(j.old, j.old + 1)])
			d2[j.old]	<- sum(d2.old[c(j.old, j.old + 1)])
			}	
					
		}
		
	
	if(g == 2){
		j.old 	<- old$j.old
		s.old	<- old$s.old
		J.old   <- old$J.old
		ind.r1.old <- old$ind.r1.old
		ind.r2.old <- old$ind.r2.old
		ind.d1.old <- old$ind.d1.old
		ind.d2.old <- old$ind.d2.old
		ind.d3.old <- old$ind.d3.old
		ind.d.temp.old	<- old$ind.d.temp.old
		d1.old	<- old$d1.old
		d2.old	<- old$d2.old
		
		J	<- J.old
		
		if(J > 1){
			ind.r1 <- ind.r1.old[,-j.old]
			ind.r2 <- NULL
			ind.d1 <- ind.d1.old[,-j.old]
			ind.d2 <- ind.d2.old[,-j.old]
			ind.d3 <- NULL
			ind.d.temp <- ind.d.temp.old[,-j.old]
			
			d1	<- d1.old[-j.old]
			d2	<- d2.old[-j.old]
			d3	<- NULL		
			
			ind.r1[,j.old] <- as.numeric(ind.r1.old[,j.old] | ind.r1.old[,j.old+1])
			ind.d1[,j.old] <- as.numeric(ind.d1.old[,j.old] | ind.d1.old[,j.old+1])
			ind.d2[,j.old] <- as.numeric(ind.d2.old[,j.old] | ind.d2.old[,j.old+1])
			ind.d.temp[,j.old] <- as.numeric(ind.d.temp.old[,j.old] | ind.d.temp.old[,j.old+1])
						
			d1[j.old]	<- sum(d1.old[c(j.old, j.old + 1)])
			d2[j.old]	<- sum(d2.old[c(j.old, j.old + 1)])
			}
			
		if(J == 1){
			ind.r1 <- matrix(1, nrow = n, ncol = 1)
			ind.r2 <- NULL
			ind.d1 <- matrix(NA, nrow = n, ncol = 1)
			ind.d2 <- matrix(NA, nrow = n, ncol = 1)
			ind.d3 <- NULL
			ind.d.temp <- matrix(NA, nrow = n, ncol = 1)
									
			d1	<- d1.old[-j.old]
			d2	<- d2.old[-j.old]
			d3	<- NULL		
			ind.d1[,j.old] <- as.numeric(ind.d1.old[,j.old] | ind.d1.old[,j.old+1])
			ind.d2[,j.old] <- as.numeric(ind.d2.old[,j.old] | ind.d2.old[,j.old+1])
			ind.d.temp[,j.old] <- as.numeric(ind.d.temp.old[,j.old] | ind.d.temp.old[,j.old+1])
						
			d1[j.old]	<- sum(d1.old[c(j.old, j.old + 1)])
			d2[j.old]	<- sum(d2.old[c(j.old, j.old + 1)])
			}				
							
		}
	
	if(g == 3){
		j.old 	<- old$j.old
		s.old	<- old$s.old
		J.old   <- old$J.old
		ind.r1.old <- old$ind.r1.old
		ind.r2.old <- old$ind.r2.old
		ind.d1.old <- old$ind.d1.old
		ind.d2.old <- old$ind.d2.old
		ind.d3.old <- old$ind.d3.old
		ind.d.temp.old	<- old$ind.d.temp.old
		d1.old	<- old$d1.old
		d3.old	<- old$d3.old

		J	<- J.old

		if(J > 1){
			ind.r1 <- NULL
			ind.r2 <- ind.r2.old[,-j.old]
			ind.d1 <- ind.d1.old[,-j.old]
			ind.d2 <- NULL
			ind.d3 <- ind.d3.old[,-j.old]
			ind.d.temp <- ind.d.temp.old[,-j.old]
			d1	<- d1.old[-j.old]
			d2	<- NULL
			d3	<- d3.old[-j.old]
		
			ind.r2[,j.old] <- as.numeric(ind.r2.old[,j.old] | ind.r2.old[,j.old+1])
			ind.d1[,j.old] <- as.numeric(ind.d1.old[,j.old] | ind.d1.old[,j.old+1])
			ind.d3[,j.old] <- as.numeric(ind.d3.old[,j.old] | ind.d3.old[,j.old+1])
			ind.d.temp[,j.old] <- as.numeric(ind.d.temp.old[,j.old] | ind.d.temp.old[,j.old+1])
			d1[j.old]	<- sum(d1.old[c(j.old, j.old + 1)])
			d3[j.old]	<- sum(d3.old[c(j.old, j.old + 1)])
			}
			
		if(J == 1){
			ind.r1 <- NULL
			ind.r2 <- matrix(NA, nrow = n, ncol = 1)
			ind.d1 <- matrix(NA, nrow = n, ncol = 1)
			ind.d2 <- NULL
			ind.d3 <- matrix(NA, nrow = n, ncol = 1)
			ind.d.temp <- matrix(NA, nrow = n, ncol = 1)
			d1	<- d1.old[-j.old]
			d2	<- NULL
			d3	<- d3.old[-j.old]
		
			ind.r2[,j.old] <- as.numeric(ind.r2.old[,j.old] | ind.r2.old[,j.old+1])
			ind.d1[,j.old] <- as.numeric(ind.d1.old[,j.old] | ind.d1.old[,j.old+1])
			ind.d3[,j.old] <- as.numeric(ind.d3.old[,j.old] | ind.d3.old[,j.old+1])
			ind.d.temp[,j.old] <- as.numeric(ind.d.temp.old[,j.old] | ind.d.temp.old[,j.old+1])
			d1[j.old]	<- sum(d1.old[c(j.old, j.old + 1)])
			d3[j.old]	<- sum(d3.old[c(j.old, j.old + 1)])
			}	


				
		}				
				
	
	list(ind.r1 = ind.r1, ind.r2 = ind.r2, ind.d1 = ind.d1, ind.d2 = ind.d2, ind.d3 = ind.d3, 		ind.d.temp = ind.d.temp, d1 = d1, d2 = d2, d3 = d3);
	
	}







cal.Del1.2.BpeScr <-
function(survObj, ini, g){
	
	n <- survObj$n
	p <- survObj$p	
	
	if(g == 1){
		s	<- ini$s1
		J	<- ini$J1
		ind.r1	<- ini$ind.r1.sc1
		ind.d1	<- ini$ind.d1.sc1
		ind.d2	<- ini$ind.d2.sc1
		ind.d.temp	<- ini$ind.d.temp.sc1
		case1_y1leq 			<- survObj$case1_y1leq.s1
		case01y2leq				<- survObj$case01y2leq.s1
		case00yleq				<- survObj$case00yleq.s1
		case1_or01or00y1y2leq 	<- survObj$case1_or01or00y1y2leq.s1
		case1_or01or00y1y2gre 	<- survObj$case1_or01or00y1y2gre.s1
		}	
	
	if(g == 2){
		s	<- ini$s2
		J	<- ini$J2
		ind.r1	<- ini$ind.r1.sc2
		ind.d1	<- ini$ind.d1.sc2
		ind.d2	<- ini$ind.d2.sc2
		ind.d.temp	<- ini$ind.d.temp.sc2		
		case1_y1leq 			<- survObj$case1_y1leq.s2
		case01y2leq				<- survObj$case01y2leq.s2
		case00yleq				<- survObj$case00yleq.s2
		case1_or01or00y1y2leq 	<- survObj$case1_or01or00y1y2leq.s2
		case1_or01or00y1y2gre 	<- survObj$case1_or01or00y1y2gre.s2
		}	
	
	y1		<- survObj$y1
	y2		<- survObj$y2	
	delta1	<- survObj$delta1
	delta2	<- survObj$delta2	
	
	case1_	<- survObj$case1_
	case01	<- survObj$case01
			
	Delta	<- matrix(0, n, J+1)
	
	for(i in case1_y1leq){
				
		y.j	<- which(ind.d.temp[i,] == 1)
				
		if(y.j != 1){
			Delta[i,y.j] <- y1[i] - s[y.j - 1]
			}

		if(y.j == 1){
			Delta[i, y.j] <- y1[i] - 0
			}			
		}	
	
	for(i in case01y2leq){
				
		y.j	<- which(ind.d.temp[i,] == 1)
				
		if(y.j != 1){
			Delta[i,y.j] <- y2[i] - s[y.j - 1]
			}

		if(y.j == 1){
			Delta[i, y.j] <- y2[i] - 0
			}			
		}
		
	for(i in case00yleq){
				
		y.j	<- which(ind.d.temp[i,] == 1)
				
		if(y.j != 1){
			Delta[i,y.j] <- y1[i] - s[y.j - 1]
			}

		if(y.j == 1){
			Delta[i, y.j] <- y1[i] - 0
			}			
		}		
		

		
	val	<- (matrix(rep(diff(c(0, s)), n), n, J+1, byrow = T) * (ind.r1 - ind.d.temp) + Delta) * ind.r1
		
	return(val);
	
	}

cal.Del3.BpeScr <-
function(survObj, ini, g = 3){

	n <- survObj$n
	p <- survObj$p	

	s	<- ini$s3
	J	<- ini$J3
	ind.r2	<- ini$ind.r2.sc3
	ind.d1	<- ini$ind.d1.sc3
	ind.d3	<- ini$ind.d3.sc3
	ind.d.temp	<- ini$ind.d.temp.sc3
	s.max	<- max(s)	

	case1_y1leq		<- survObj$case1_y1leq.s3
	case01y2leq		<- survObj$case01y2leq.s3
	case00yleq		<- survObj$case00yleq.s3
	case1_or01or00y1y2leq	<- survObj$case1_or01or00y1y2leq.s3
	case1_or01or00y1y2gre	<- survObj$case1_or01or00y1y2gre.s3
		
	case11y2leq		<- survObj$case11y2leq.s3
	case10y2leq		<- survObj$case10y2leq.s3
	case11or10y2leq	<- survObj$case11or10y2leq.s3
	case11or10y2gre	<- survObj$case11or10y2gre.s3
			
	y1		<- survObj$y1
	y2		<- survObj$y2	
	delta1	<- survObj$delta1
	delta2	<- survObj$delta2	
	
	case10	<- survObj$case10
	case11	<- survObj$case11
		
	Delta	<- matrix(0, n, J+1);
	
	diff.s <- diff(c(0, s))
	
	for(i in case11or10y2leq){
		start.j <- which(ind.d1[i,] == 1)
		end.j	<- which(ind.d.temp[i,] == 1)
		
		if(start.j == end.j){	# two failures are observed from the same interval
			Delta[i, start.j] <- y2[i] - y1[i];
			}
			
		if(start.j != end.j){	# two failures are NOT observed from the same interval
			Delta[i, start.j:end.j] <- diff.s[start.j:end.j]
			Delta[i, start.j] <- s[start.j] - y1[i]
			Delta[i, end.j] <- y2[i] - s[end.j - 1]
			}			
		}	
	
	if(length(case11or10y2gre) > 0){
		for(i in case11or10y2gre){
			start.j <- which(ind.d1[i,] == 1)
			Delta[i,start.j:(J+1)] <- diff.s[start.j:(J+1)]
			Delta[i, start.j] <- s[start.j] - y1[i]
			}
		}

	return(Delta);

	}




cal.Del1.2.BI.BpeScr <-
function(survObj, ini.new, old, g){

	y1		<- survObj$y1
	y2		<- survObj$y2	
	delta1	<- survObj$delta1
	delta2	<- survObj$delta2
	n		<- survObj$n
	case1_	<- survObj$case1_
	case01	<- survObj$case01	
	case00	<- survObj$case00	
	
	if(g == 1){		
		s	<- ini.new$s1
		J	<- ini.new$J1
		ind.r1	<- ini.new$ind.r1.sc1
		ind.d1	<- ini.new$ind.d1.sc1
		ind.d2	<- ini.new$ind.d2.sc1
		ind.d.temp <- ini.new$ind.d.temp.sc1
		
		case1_y1leq 			<- survObj$case1_y1leq.s1
		case01y2leq				<- survObj$case01y2leq.s1
		case00yleq				<- survObj$case00yleq.s1
		case1_or01or00y1y2leq 	<- survObj$case1_or01or00y1y2leq.s1
		case1_or01or00y1y2gre 	<- survObj$case1_or01or00y1y2gre.s1
		
		j.old 	<- old$j.old
		
		s.old	<- old$s.old
		J.old   <- old$J.old
		ind.r1.old <- old$ind.r1.old
		ind.d1.old <- old$ind.d1.old
		ind.d2.old <- old$ind.d2.old	
		ind.d.temp.old <- old$ind.d.temp.old
				
		Delta.old <- ini.new$Delta1
		
		}
		

	
	if(g == 2){
		s	<- ini.new$s2
		J	<- ini.new$J2
		ind.r1	<- ini.new$ind.r1.sc2
		ind.d1	<- ini.new$ind.d1.sc2
		ind.d2	<- ini.new$ind.d2.sc2
		
		case1_y1leq 			<- survObj$case1_y1leq.s2
		case01y2leq				<- survObj$case01y2leq.s2
		case00yleq				<- survObj$case00yleq.s2
		case1_or01or00y1y2leq 	<- survObj$case1_or01or00y1y2leq.s2
		case1_or01or00y1y2gre 	<- survObj$case1_or01or00y1y2gre.s2
		
		j.old 	<- old$j.old
		
		s.old	<- old$s.old
		J.old   <- old$J.old
		ind.r1.old <- old$ind.r1.old
		ind.d1.old <- old$ind.d1.old
		ind.d2.old <- old$ind.d2.old	
		ind.d.temp.old <- old$ind.d.temp.old
				
		
		Delta.old <- ini.new$Delta2
		
		}		
				
		Delta <- cbind(Delta.old, rep(0, n))
		
		diff.s <- diff(c(0, s))
		diff.s.old <- diff(c(0, s.old))
				
		if(j.old != (J.old + 1)){
			Delta[,(j.old+1):(J.old+2)] <- Delta[,c(J.old+2, (j.old+1):(J.old+1))]
			}		
		
		for(i in case1_or01or00y1y2gre){
			Delta[i,] <- diff.s
			}
		
		for(i in intersect(case1_or01or00y1y2leq, which(Delta.old[,j.old] != 0))){
			if(j.old == 1){
				Delta[i, j.old] <- max(0, min(y1[i], s[j.old]) - 0)
				}
			if(j.old != 1){
				Delta[i, j.old] <- max(0, min(y1[i], s[j.old]) - s[j.old-1])
				}			
			Delta[i, j.old+1] <- max(0, min(y1[i], s[j.old+1]) - s[j.old])
			}
							

	return(Delta);
	
	}





cal.Del3.BI.BpeScr <-
function(survObj, ini.new, old, g = 3){

	n <- survObj$n
	p <- survObj$p		

	s	<- ini.new$s3
	J	<- ini.new$J3
	ind.r2	<- ini.new$ind.r2.sc3
	ind.d1	<- ini.new$ind.d1.sc3
	ind.d3	<- ini.new$ind.d3.sc3
	s.max	<- max(s)	

	case1_y1leq			<- survObj$case1_y1leq.s3
	case01y2leq			<- survObj$case01y2leq.s3
	case11y2leq			<- survObj$case11y2leq.s3
	case1_or01y1y2leq	<- survObj$case1_or01y1y2leq.s3
	case1_or01y1y2gre	<- survObj$case1_or01y1y2gre.s3
	case11y1leqy2gre	<- survObj$case11y1leqy2gre.s3
	case10y1leq			<- survObj$case10y1leq.s3
			
	y1		<- survObj$y1
	y2		<- survObj$y2	
	delta1	<- survObj$delta1
	delta2	<- survObj$delta2	
	
	case10	<- survObj$case10
	case11	<- survObj$case11
	
	j.old 	<- old$j.old
		
	s.old	<- old$s.old
	J.old   <- old$J.old
	ind.r1.old <- old$ind.r1.old
	ind.d1.old <- old$ind.d1.old
	ind.d2.old <- old$ind.d2.old	
		
	Delta.old <- ini.new$Delta3	
	
	Delta <- cbind(Delta.old, rep(0, n))
		
	diff.s <- diff(c(0, s))
	diff.s.old <- diff(c(0, s.old))	
	
	if(j.old != (J.old + 1)){
		Delta[,(j.old+1):(J.old+2)] <- Delta[,c(J.old+2, (j.old+1):(J.old+1))]
		}			

	if(j.old  != J.old + 1){
		for(i in intersect(union(case11y2leq, case10y1leq), which(Delta.old[,j.old] != 0 | Delta.old[,j.old + 1] != 0))){
			if(j.old == 1){
				Delta[i, j.old] <- max(0, min(y2[i], s[j.old]) - max(y1[i], 0))
				}
			if(j.old != 1){
				Delta[i, j.old] <- max(0, min(y2[i], s[j.old]) - max(y1[i], s[j.old-1]))
				}
			
			Delta[i, j.old+1] <- max(0, min(y2[i], s[j.old+1]) - max(y1[i], s[j.old]))
			}		
		}
	

	if(j.old  == J.old + 1){
		for(i in intersect(case11y2leq , which(Delta.old[,j.old] != 0))){
			if(j.old == 1){
				Delta[i, j.old] <- max(0, min(y2[i], s[j.old]) - max(y1[i], 0))
				}
			if(j.old != 1){
				Delta[i, j.old] <- max(0, min(y2[i], s[j.old]) - max(y1[i], s[j.old-1]))
				}
			}		
		}				
	
	for(i in case11y1leqy2gre){
		chk1 	<- which(ind.d1[i, ] == 1);
		Delta[i, chk1:(J+1)] <- diff.s[chk1:(J+1)]
		Delta[i, chk1] <- s[chk1] - y1[i];		
		}	
	


	return(Delta);

	}






cal.Del.DI.BpeScr <-
function(survObj, ini.new, old, g){

	y1		<- survObj$y1
	y2		<- survObj$y2	
	delta1	<- survObj$delta1
	delta2	<- survObj$delta2
	n		<- survObj$n
	
	
	if(g == 1){		
		Delta.old <- ini.new$Delta1
		J.new	<- ini.new$J1
		}
		
	if(g == 2){		
		Delta.old <- ini.new$Delta2
		J.new	<- ini.new$J2
		}
		
	if(g == 3){		
		Delta.old <- ini.new$Delta3
		J.new	<- ini.new$J3		
		}				
		
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





















cal.Sigma.lam.BpeScr <-
function(ini, priorPara, g){
	
	if(g == 1){
		s <- ini$s1
		J <- ini$J1		
		c.lam <- priorPara$c.lam1		
		}
		
	if(g == 2){
		s <- ini$s2
		J <- ini$J2		
		c.lam <- priorPara$c.lam2		
		}
		
	if(g == 3){
		s <- ini$s3
		J <- ini$J3		
		c.lam <- priorPara$c.lam3
		}				
		
	Del.s	<- diff(c(0, s))
	
	if(J+1 >= 3){	
		W	<- Q	<- matrix(0, J+1, J+1);
	
		W[1, 2]	<- c.lam * (Del.s[1] + Del.s[2]) / (2*Del.s[1] + Del.s[2]);
		W[J+1, (J)]	<- c.lam * (Del.s[J] + Del.s[J+1]) / (Del.s[J] + 2*Del.s[J+1]);
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
		W[J+1, (J)]	<- c.lam * (Del.s[J] + Del.s[J+1]) / (Del.s[J] + 2*Del.s[J+1]);
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









