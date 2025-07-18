################### Script: 2fit_ML2.R ##################################

#####################################################
######## Step 4: iterative model estimation #########
#####################################################
## In this step, we used all previously generated results to estimate the coefficients of the disaggregation model

## We set the stop criteria
### 'n.dif' equals the norm between the vector of coefficients in two successive iterations
### low 'n.dif' values indicate that the model has converged
while(it < maxit){
	it <- it + 1

	## At each iteration, we print its number and calculate the necessary matrices (using 'out.fn')
	### At each step, everything is saved to help debug in case of problems
	print(paste0('iteration #', it))
	ll.t1 <- out.fn(b=b0)

	mtcs[[it]] <- ll.t1

	## Then, we optimize the log-likelihood function considering all the matrices generated
	### In order to improve the chance of finding the true minima, we do it multiple times
	### Here I decided to do it only 3 times due to the high number of coarse pixels, but this varies in each case
	### We then take the result with lower log-likelihood
	print(paste0('starting optimization: ', Sys.time()))	
	l <- lapply(1:10, function(i){
		print(i)
		set.seed(i)
		jj = rnorm(np)
		jj[1] = 10 # Attempt to generate feasible starting values

		optim(jj, fn=fitCpp, gr=grfitCpp, Qty=ll.t1[[1]], R=ll.t1[[2]],
			RtR = crossprod(ll.t1[[2]]), r=ll.t1[[3]], N=ll.t1[[4]], b0=ll.t1[[5]],
			S1=S1, S2=S2, S3=S3, S4=S4, S5=S5, S6=S6, S7=S7, S8=S8, S9=S9, S10=S10,
			sldiagA=ll.t1[[6]], it=it, method='BFGS')
	})
	vs <- sapply(l, function(i) i$value)
	print(vs)
	l <- l[[which.min(vs)]]

	### If the likelihood increases compared to the previous iteration, then the model is no longer improving and we can stop.
	ml[it] <- as.numeric(l$value)
	if(it > 2 & !res){
		print(ml[it])
		if(ml[it] > ml[it-1]){
			print('end')
			break
		}
	} else {
		res = FALSE
	}

	mtcs[[it]][[7]] <- l$par
	saveRDS(mtcs, 'output/1mtcs.rds')
	if(it > 1){
		saveRDS(mtcs[[it]], 'output/1mtcs_sol.rds')
		mtcs[[it-1]] <- NA
		gc()
	}

	## After model optimization, we can calculate:
	### 'bnew': the new updated vector of coefficients
	### 'gg': the norm of the gradient at the minimum. Ideally, this should be close to zero, indicating a successful optimization
	RtR <- crossprod(ll.t1[[2]])
	bnew <- hbCpp(l$par, Qty=ll.t1[[1]], R=ll.t1[[2]], RtR = RtR, r=ll.t1[[3]], N=ll.t1[[4]], b0=ll.t1[[5]],
		S1=S1, S2=S2, S3=S3, S4=S4, S5=S5, S6=S6, S7=S7, S8=S8, S9=S9, S10=S10)
	gg <- crossprod(grfitCpp(l$par, Qty=ll.t1[[1]], R=ll.t1[[2]], RtR = RtR, r=ll.t1[[3]], N=ll.t1[[4]], b0=ll.t1[[5]],
		S1=S1, S2=S2, S3=S3, S4=S4, S5=S5, S6=S6, S7=S7, S8=S8, S9=S9, S10=S10, sldiagA=ll.t1[[6]], it=it))

	## We also calculate the effective degrees of freedom (edf) for the fit
	st <- l$par
	rho <- exp(st)
	Sl <- rho[2]*S1 + rho[3]*S2 + rho[4]*S3 + rho[5]*S4 + rho[6]*S5 + rho[7]*S6 + rho[8]*S7 + rho[9]*S8 + rho[10]*S9 + rho[11]*S10
	V <- RtR/rho[1] + Sl
	eigenV <- eigen(V)
	sqrtinv <- sqrt(1/eigenV$values)
	hat.cov <- crossprod(sqrtinv * t(eigenV$vectors))
	## Total edf (i.e., 'edf.t') would be given by the formula below
	### edf.t <- sum(hat.cov * RtR/rho[1]) # Equivalent to: sum(diag(hat.cov %*% RtR/rho[1]))
	## However, we split per variable (i.e., 'edf.mult') to facilitate early model diagnostics
	H <- rowSums(hat.cov * t(RtR/rho[1])) # Equivalent to: diag(hat.cov %*% RtR/rho[1])
	edf.t <- sum(H) # Total edf
	edf.mult <- tot.mult <- NA * dims
	for(di in seq_along(dims)){
		if(di == 1){
			edf.mult[di] <- sum(H[2:dims[di]])
			tot.mult[di] <- length(2:dims[di])
		} else {
			edf.mult[di] <- sum(H[(dims[di-1] + 1):dims[di]])
			tot.mult[di] <- length((dims[di-1] + 1):dims[di])
		}
	}
	
	## Then, we calculate the Akaike and Bayesian information criteria
	aic[it] <- 2 * ml[it] - 2 * edf.t
	bic[it] <- 2 * ml[it] - log(ll.t1[[4]]) * edf.t

	## Then, we update the variables before starting the next iteration, and print some relevant information
	n.dif <- crossprod(b0 - bnew)
	b0 <- bnew
	print(paste0('n.dif: ', round(n.dif, 4)))
	print(paste0('norm gr @ minimum: ', gg))
	print(paste0('ll: ', round(ml[it], 4)))
	print(paste0('AIC: ', round(aic[it], 2)))
	print(paste0('BIC: ', round(bic[it], 4)))
	print(paste0('rhos: ', paste0(round(st, 4), collapse = ' ')))
	print(paste0('s2: ', round(exp(st[1]), 8)))
	print(paste0('edf: ', paste(sprintf("%.2f/%.0f", edf.mult, tot.mult), collapse = ', ')))

	rm(ll.t1); gc()
}
