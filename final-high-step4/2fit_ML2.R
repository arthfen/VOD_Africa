################### Script: 2fit_ML2.R ##################################

#####################################################
######## Step 4: iterative model estimation #########
#####################################################
## In this step, we used all previously generated results to estimate the coefficients of the disaggregation model

dta <- 500
fit_at_min <- function(rho, ref, v0, Sv1, Sv2, Sv3, Sv4, Sv5, Sv6, Sv7, g0_var, B0_var){
	rho <- exp(rho)
	if(min(rho) < 1e-10 | max(rho) > 1e10) return(Inf)
	Slv = rho[1]*Sv1 + rho[2]*Sv2 + rho[3]*Sv3 + rho[4]*Sv4 + rho[5]*Sv5 + rho[6]*Sv6 + rho[7]*Sv7

#	ldetSlv <- determinant(Slv[-1,-1], logarithm = TRUE)$modulus
	g0v = g0_var + Slv %*% v0;
	B0v = Bk + Slv
	tryCatch({
		hatv = v0 - solve(B0v, g0v)
	}, error = function(e){
		return(Inf)
	})
	
	logLik = ref + 1/2 * t(v0) %*% Slv %*% v0 + 1/2 * t(g0v) %*% (hatv - v0) + 1/4 * t(hatv - v0) %*% B0v %*% (hatv - v0)
	
	attr(logLik, 'pk') <- hatv - v0
	attr(logLik, 'Slv') <- Slv
	logLik
}
grfit_at_min <- function(x, ...) grad(fit_at_min, x, ..., method = 'simple')

ref_loglik <- function(diagAVA){
	diagMMt <- diagAVA
	ql <- diagMMt > 0
	diagL <- 0 * diagMMt
	diagL[ql] <- sqrt(1/diagMMt[ql])
	diagL <- as.numeric(diagL)
	y. <- diagL * mtcs.t1[[1]]
	X. <- diagL * mtcs.t1[[2]]
	X.qr <- qr(X., LAPACK = TRUE)
	Qty. <- qr.qty(X.qr, y.)
	R <- qr.R(X.qr, complete = TRUE)[,sort.list(X.qr$pivot)]
	r <- crossprod(y.) - crossprod(Qty.)
	RtR = crossprod(R)
	sldiagA <- sum(log(diagMMt))
	ref.0 <- ll_at_min(st, Qty=Qty., R=R, RtR=RtR, r=r, N=length(mtcs.t1[[1]]), b0=b0,
		S1=S1, S2=S2, S3=S3, S4=S4, S5=S5, S6=S6, S7=S7, S8=S8, S9=S9, S10=S10, S11=S11, S12=S12, S13=S13, S14=S14, S15=S15, S16=S16, sldiagA=sldiagA)

	return(ref.0)
}

## We set the stop criteria
### 'n.dif' equals the norm between the vector of coefficients in two successive iterations
### low 'n.dif' values indicate that the model has converged
while(it < maxit){
	it <- it + 1

	## At each iteration, we print its number and calculate the necessary matrices (using 'out.fn')
	### At each step, everything is saved to help debug in case of problems
	print(paste0('iteration #', it))
	if(it == 1){ # If it is the first iteration, we have to run out.fn
		# (re)starts the cluster
		cl <- start_cl(nodes)
		vk = list(matrix(v0, ncol = 1))
		ll.t2 <- out.fn(b=b0, v=vk)
		stop_cl(cl)

		print(paste0('starting optimization: ', Sys.time()))	

		ref.0 <- ref_loglik(1/mtcs.t1[[3]])
		print(paste0('previous ll [approx., without any penalties]: ', round(ref.0, 4)))
	} # Otherwise, ll.t2 is already obtained from the previous iteration

	Y_ <- mtcs.t1[[1]]
	rho <- exp(st)
	Sl <- rho[2]*S1 + rho[3]*S2 + rho[4]*S3 + rho[5]*S4 + rho[6]*S5 + rho[7]*S6 + rho[8]*S7 + rho[9]*S8 + rho[10]*S9 + rho[11]*S10 + rho[12]*S11 + rho[13]*S12 + rho[14]*S13 + rho[15]*S14 + rho[16]*S15 + rho[17]*S16

	# The function to be minimized is : f(v) = ll_at_min(v, ...) + 1/2 * v' %*% Slv %*% v
	# The gradient below is exact for ll_at_min(v, ...)
	nk <- ncol(ll.t2[[2]][[1]])
	invAVA <- as.numeric(1/ll.t2[[1]][[1]])
	Vinv <- solve(crossprod(mtcs.t1[[2]] * sqrt(invAVA))/rho[1] + Sl, tol = 1e-30)
	g0_var <- as.numeric(t(Y_) %*% (- invAVA^2 * ll.t2[[2]][[1]] * Y_)/rho[1]) # gradient of the squared error
	g0_var <- g0_var + colSums(ll.t2[[2]][[1]] * invAVA) # gradient of the log determinant (likelihood)
	T.m <- mtcs.t1[[2]] %*% Vinv
	d.v <- rowSums(T.m * mtcs.t1[[2]])
	w.v <- (invAVA^2) * d.v / rho[1]
	g0_var <- g0_var - c(w.v %*% ll.t2[[2]][[1]])
	g0_var <- g0_var * 1/2 # because we actually return 0.5 * logLik in ll_at_min

	# Hessian BFGS
	if(it == 1){
		if(!inBk){
			dg <- as.numeric(crossprod(g0_var))/dta
			Bk <- diag(dg + 0 * g0_var) # non-inverted Hessian
		}
	} else if(it == 2){
		Bk <- diag(as.numeric((t(yk) %*% yk)/(t(yk) %*% sk)) + 0 * g0_var) # non-inverted Hessian
	} else {
		Bk <- Bk - (Bk %*% sk %*% t(sk) %*% Bk)/as.numeric(t(sk) %*% Bk %*% sk) + (yk %*% t(yk))/as.numeric(t(yk) %*% sk)
	}

	# Optimizes the smoothing parameters
	ref_np[it] <- ref_loglik(ll.t2[[1]][[1]]) # This is the current function evaluation (without sigma smoother penalty)
	l <- lapply(1:5, function(i){
		set.seed(i)
		jj = rnorm(np_var)

		optim(jj, fn=fit_at_min, gr=grfit_at_min, ref=ref_np[it], v0=v0, Sv1=Sv1, Sv2=Sv2, Sv3=Sv3, Sv4=Sv4, Sv5=Sv5, Sv6=Sv6, Sv7=Sv7, g0_var=g0_var, B0_var=Bk, method = 'BFGS')
	})
	vs <- sapply(l, function(i) i$value)
	print(vs)
	l <- l[[which.min(vs)]]
	lmin <- fit_at_min(l$par, ref=ref_np[it], v0=v0, Sv1=Sv1, Sv2=Sv2, Sv3=Sv3, Sv4=Sv4, Sv5=Sv5, Sv6=Sv6, Sv7=Sv7, g0_var=g0_var, B0_var=Bk)
	pk <- attr(lmin, 'pk')

	Slv <- attr(lmin, 'Slv')
	ref[it] <- ref_np[it] + 1/2 * t(v0) %*% Slv %*% v0 # The current function evaluation with sigma smoother penalty

	# Now, we call out.fn again to perform the line search with Wolfe conditions:
	c1 <- 1e-4
	c2 <- 0.9
	aps <- c(0.25, 0.5, 1, 2, 4, 8, 16)
	vnew <- lapply(seq_along(aps), function(j){
		v0 + aps[j] * pk
	})

	cl <- start_cl(nodes)
	aps.t2 <- out.fn(b=b0, v=vnew)
	stop_cl(cl)

	grds <- lapply(1:length(aps), function(j){ # We calculate the gradients at the candidates
		nk <- ncol(aps.t2[[2]][[j]])
		invAVA <- as.numeric(1/aps.t2[[1]][[j]])
		Vinv <- solve(crossprod(mtcs.t1[[2]] * sqrt(invAVA))/rho[1] + Sl, tol = 1e-30)
		j.g0_var <- as.numeric(t(Y_) %*% (- invAVA^2 * aps.t2[[2]][[j]] * Y_)/rho[1]) # gradient of the squared error
		j.g0_var <- j.g0_var + colSums(aps.t2[[2]][[j]] * invAVA) # gradient of the log determinant (likelihood)
		T.m <- mtcs.t1[[2]] %*% Vinv
		d.v <- rowSums(T.m * mtcs.t1[[2]])
		w.v <- (invAVA^2) * d.v / rho[1]
		j.g0_var <- j.g0_var - c(w.v %*% aps.t2[[2]][[j]])
		j.g0_var <- j.g0_var * 1/2 # because we actually return 0.5 * logLik in ll_at_min

		return(j.g0_var)
	})
	wlf <- lapply(1:length(aps), function(j){ # Here, we calculate the Wolfe conditions
		# This looks at logLik
		# f.new <- ref_loglik(aps.t2[[1]][[j]])
		# cond1 <- (f.new <= ref_np[it] + c1 * aps[j] * (t(g0_var) %*% pk))
		# cond2 <- (t(grds[[j]]) %*% pk >= c2 * (t(g0_var) %*% pk))

		# This looks at logLik + 1/2 * v'Sv
		f.new <- as.numeric(ref_loglik(aps.t2[[1]][[j]]) + 1/2 * t(vnew[[j]]) %*% Slv %*% vnew[[j]])
		cond1 <- (f.new <= ref[it] + c1 * aps[j] * (t(g0_var + Slv %*% v0) %*% pk))
		cond2 <- (t(grds[[j]] + Slv %*% vnew[[j]]) %*% pk >= c2 * (t(g0_var + Slv %*% v0) %*% pk))

		return(c(f.new, cond1 + cond2))
	})
	wlf <- do.call(rbind, wlf)

	## Then, we pick the highest alpha satisfying the conditions
	if(sum(wlf[,2] == 2) == 0){ # In this case, no one satisfies the conditions, we stop
		print('Wolfe conditions not satisfied, end')
		break
	} else { # In this case, at least one satisfies so we pick the lowest function evaluation
		idw <- which(wlf[,2] == 2)
		nw <- which(wlf[,1] == min(wlf[idw, 1]))
	}
	alpha <- aps[nw]
	vnew <- vnew[[nw]]
	
	sk <- vnew - v0
	yk <- grds[[nw]] - g0_var

	### If the likelihood increases compared to the previous iteration, then the model is no longer improving and we can stop.
	ml[it] <- as.numeric(l$value)
	if(it > 2){
		if(!is.na(ref[it-1])){
			if(ref[it] > ref[it-1]){
				print('end')
				break
			} else {
				if(abs((ref[it]-ref[it-1])/ref[it-1]) < 1e-6){
					print('convergence')
					break
				}
			}
		}
	}

	ll.t2[[4]] <- l$par
	ll.t2[[5]] <- Bk
	ll.t2[[6]] <- v0
	if(it > 1){
		saveRDS(ll.t2, 'output/1mtcs_sol-var.rds')
	}

	## Then, we update the variables before starting the next iteration, and print some relevant information
	n.dif <- crossprod(v0 - vnew)
	v0 <- vnew
	st.var <- l$par
	ll.t2 <- list(list(aps.t2[[1]][[nw]]), list(aps.t2[[2]][[nw]])) # We don't need to run ll.t2 again in the next iteration

	print(paste0('gradient norm: ', round(crossprod(g0_var), 4)))
	print(paste0('alpha: ', round(alpha, 4)))
	print(paste0('n.dif: ', round(n.dif, 4)))
	print(paste0('rhos: ', paste0(round(st.var, 4), collapse = ' ')))
	print(paste0('reference ll [no smoother penalty]: ', round(ref_np[it], 4)))
	print(paste0('reference ll [approx., with smoother penalty]: ', round(ref[it], 4)))
	print(paste0('optimized ll [approx., with smoother penalty]: ', round(ml[it], 4)))
}
