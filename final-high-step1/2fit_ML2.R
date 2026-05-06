################### Script: 2fit_ML2.R ##################################

#####################################################
######## Step 4: iterative model estimation #########
#####################################################
## In this step, we used all previously generated results to estimate the coefficients of the disaggregation model

tries <- 0
while(it < maxit){
	# (re)starts the cluster
	cl <- start_cl(nodes)
	it <- it + 1

	## At each iteration, we print its number and calculate the necessary matrices (using 'out.fn')
	### At each step, everything is saved to help debug in case of problems
	print(paste0('iteration #', it))
	ll.t1 <- out.fn(b=b0)
	stop_cl(cl)

	## Then, we optimize the log-likelihood function considering all the matrices generated
	### In order to improve the chance of finding the true minima, we do it multiple times
	### Here I decided to do it only 3 times due to the high number of coarse pixels, but this varies in each case
	### We then take the result with lower log-likelihood
	print(paste0('starting optimization: ', Sys.time()))
	RtR = crossprod(ll.t1[[2]])
	l <- lapply(1:10, function(i){
		set.seed(i)
		jj = rnorm(np)
		jj[1] = 10 # Attempt to generate feasible starting values

		tryCatch({
			optim(jj, fn=fitCpp, gr=grfitCpp, Qty=ll.t1[[1]], R=ll.t1[[2]],
				RtR = RtR, r=ll.t1[[3]], N=ll.t1[[4]], b0=ll.t1[[5]],
				S1=S1, S2=S2, S3=S3, S4=S4, S5=S5, S6=S6, S7=S7, S8=S8, S9=S9, S10=S10, S11=S11, S12=S12, S13=S13, S14=S14, S15=S15, S16=S16,
				sldiagA=ll.t1[[6]], it=it, method='BFGS', control = list(maxit = 1000))
		}, error = function(e){
			data.frame(value = NA)
		})
	})
	vs <- sapply(l, function(i) i$value)
	if(all(is.na(vs))){
		# In that case, we try to reduce the step
		print('the optimization failed')
		if(tries < 3){
			tries <- tries + 1
			print(paste0('retrying this iteration with a smaller step: ', tries, '/3'))

			ref[it] <- NA
			ml[it] <- NA
			it = it - 1
			
			# When tries = 0, b0 = b0.cur + (b0.new - b0.cur)
			## Here, we try a smaller step, with: b0 = b0.cur + 2^(-tries) * (b0.new - b0.cur)
			b0 <- b0.cur + 2^(-tries) * (b0.new - b0.cur)
			next
		} else {
			print('end')
			break
		}
	}
	
	print(vs)
	l <- l[[which.min(vs)]]

	## After model optimization, we can calculate:
	### 'bnew': the new updated vector of coefficients
	### 'gg': the norm of the gradient at the minimum. Ideally, this should be close to zero, indicating a successful optimization
	RtR <- crossprod(ll.t1[[2]])
	bnew <- hbCpp(l$par, Qty=ll.t1[[1]], R=ll.t1[[2]],
				RtR = RtR, r=ll.t1[[3]], N=ll.t1[[4]], b0=ll.t1[[5]],
				S1=S1, S2=S2, S3=S3, S4=S4, S5=S5, S6=S6, S7=S7, S8=S8, S9=S9, S10=S10, S11=S11, S12=S12, S13=S13, S14=S14, S15=S15, S16=S16)
	n.dif <- crossprod(b0 - bnew)
	gg <- crossprod(grfitCpp(l$par, Qty=ll.t1[[1]], R=ll.t1[[2]],
				RtR = RtR, r=ll.t1[[3]], N=ll.t1[[4]], b0=ll.t1[[5]],
				S1=S1, S2=S2, S3=S3, S4=S4, S5=S5, S6=S6, S7=S7, S8=S8, S9=S9, S10=S10, S11=S11, S12=S12, S13=S13, S14=S14, S15=S15, S16=S16, sldiagA=ll.t1[[6]], it=it))

	### If the likelihood increases compared to the previous iteration, then the model is no longer improving and we can stop.
	ref.l <- lapply(1:10, function(i){
		set.seed(i)
		jj = rnorm(np)
		jj[1] = 10 # Attempt to generate feasible starting values

		tryCatch({
			optim(jj, fn=ll_at_min, Qty=ll.t1[[1]], R=ll.t1[[2]],
				RtR = RtR, r=ll.t1[[3]], N=ll.t1[[4]], b0=ll.t1[[5]],
				S1=S1, S2=S2, S3=S3, S4=S4, S5=S5, S6=S6, S7=S7, S8=S8, S9=S9, S10=S10, S11=S11, S12=S12, S13=S13, S14=S14, S15=S15, S16=S16,
				sldiagA=ll.t1[[6]], method='BFGS', control = list(maxit = 1000))
		}, error = function(e){
			data.frame(value = NA)
		})
	})
	ref[it] <- min(sapply(ref.l, function(i) i$value), na.rm = TRUE)	
#	ml[it] <- as.numeric(l$value - 0.5 * 100/2^(it*pen-1) * n.dif) # We remove the artificial it penalty
	ml[it] <- as.numeric(l$value - 0.5 * 100/2^(it-1) * n.dif) # We remove the artificial it penalty
	if(it > 1){
		if(!is.na(ref[it-1]) & exists('b0.cur')){
			if(ref[it] > ref[it-1]){
				# In that case, we try to reduce the step
				print(paste0('reference ll: ', round(ref[it], 4), ' (higher than current)'))
				if(tries < 3){
					tries <- tries + 1
					print(paste0('retrying this iteration with a smaller step: ', tries, '/3'))

					ref[it] <- NA
					ml[it] <- NA
					it = it - 1
					
					# When tries = 0, b0 = b0.cur + (b0.new - b0.cur)
					## Here, we try a smaller step, with: b0 = b0.cur + 2^(-tries) * (b0.new - b0.cur)
					b0 <- b0.cur + 2^(-tries) * (b0.new - b0.cur)
					next
				} else {
					print('end')
					break
				}
			} else {
				tries <- 0
				if(abs((ref[it]-ref[it-1])/ref[it-1]) < 1e-4){
					print('convergence')
					break
				}
			}
		}
	}

	# Only when successful, export the matrices
	saveRDS(attr(ll.t1[[1]], 'matrices'), 'output/1gen_mtc-mu.rds')
	conc.mult <- NA * c(NA, dims) # The estimated concurvity before removing the matrices
	conc.mult[1] <- conc_fn(X = attr(ll.t1[[1]], 'matrices')[[2]], idx = 1:npar)
	for(di in seq_along(dims)){
		if(di == 1){
			idx = (npar+1):dims[di]
		} else {
			idx = (dims[di-1]+1):dims[di]
		}
		conc.mult[di+1] <- conc_fn(X = attr(ll.t1[[1]], 'matrices')[[2]], idx)
	}
	attr(ll.t1[[1]], 'matrices') <- NULL
	ll.t1[[7]] <- l$par
	if(it > 1){
		saveRDS(ll.t1, 'output/1mtcs_sol.rds')
	}

	## We also calculate the effective degrees of freedom (edf) for the fit
	st <- l$par
	rho <- exp(st)
	Sl <- rho[2]*S1 + rho[3]*S2 + rho[4]*S3 + rho[5]*S4 + rho[6]*S5 + rho[7]*S6 + rho[8]*S7 + rho[9]*S8 + rho[10]*S9 + rho[11]*S10 + rho[12]*S11 + rho[13]*S12 + rho[14]*S13 + rho[15]*S14 + rho[16]*S15 + rho[17]*S16
	V <- RtR/rho[1] + Sl
	eigenV <- eigen(V)
	sqrtinv <- sqrt(1/eigenV$values)
	hat.cov <- crossprod(sqrtinv * t(eigenV$vectors))
	## Total edf (i.e., 'edf.t') would be given by the formula below
	### edf.t <- sum(hat.cov * RtR/rho[1]) # Equivalent to: sum(diag(hat.cov %*% RtR/rho[1]))
	## However, we split per variable (i.e., 'edf.mult') to facilitate early model diagnostics
	H <- rowSums(hat.cov * t(RtR/rho[1])) # Equivalent to: diag(hat.cov %*% RtR/rho[1])
	edf.t <- sum(H) # Total edf
	edf.mult <- tot.mult <- NA * c(NA, dims)
	edf.mult[1] <- sum(H[1:npar])
	tot.mult[1] <- length(1:npar)
	for(di in seq_along(dims)){
		if(di == 1){
			edf.mult[di+1] <- sum(H[(npar+1):dims[di]])
			tot.mult[di+1] <- length((npar+1):dims[di])
		} else {
			edf.mult[di+1] <- sum(H[(dims[di-1] + 1):dims[di]])
			tot.mult[di+1] <- length((dims[di-1] + 1):dims[di])
		}
	}

	## The gradient during optimization
	grd <- as.numeric(t(ll.t1[[2]]) %*% ll.t1[[1]]/rho[1] - 2 * Sl %*% b0)

	## The n.dif split per smoother
	ndif.mult <- NA * c(NA, dims)
	ndif.mult[1] <- crossprod((b0 - bnew)[1:npar])
	for(di in seq_along(dims)){
		if(di == 1){
			idx = (npar+1):dims[di]
		} else {
			idx = (dims[di-1]+1):dims[di]
		}
		ndif.mult[di+1] <- crossprod((b0 - bnew)[idx])
	}

	## Then, we calculate the Akaike and Bayesian information criteria
	aic[it] <- 2 * ml[it] - 2 * edf.t
	bic[it] <- 2 * ml[it] - log(ll.t1[[4]]) * edf.t

	## Then, we update the variables before starting the next iteration, and print some relevant information
	b0.cur <- b0
	b0.new <- bnew
	b0 <- bnew
	print(paste0('gradient norm: ', round(crossprod(grd), 4)))
	print(paste0('n.dif: ', round(n.dif, 4)))
	print(paste0('n.dif per term: ', paste(sprintf("%.1f", ndif.mult), collapse = ', ')))
	print(paste0('estim. concurvity: ', paste(sprintf("%.2f", conc.mult), collapse = ', ')))
	print(paste0('reference ll: ', round(ref[it], 4)))
	print(paste0('optimized ll [without it penalty]: ', round(ml[it], 4)))
	print(paste0('optimized ll [with it penalty]: ', round(ml[it] + 0.5 * 100/2^(it-1) * n.dif, 4)))
	print(paste0('norm gr @ minimum: ', gg))
	print(paste0('AIC: ', round(aic[it], 2)))
	print(paste0('BIC: ', round(bic[it], 4)))
	print(paste0('rhos: ', paste0(round(st, 4), collapse = ' ')))
	print(paste0('s2: ', round(exp(st[1]), 8)))
	print(paste0('edf: ', paste(sprintf("%.2f/%.0f", edf.mult, tot.mult), collapse = ', ')))

	rm(ll.t1); gc()
}
