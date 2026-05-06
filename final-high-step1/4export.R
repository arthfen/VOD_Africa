################### Script: 4export.R ##################################
## The following lines load the necessary libraries
library(disag.vAfrica)
library(numDeriv)
library(collapse)
library(terra)
library(mgcv)
library(snow)
library(Rcpp)

npar <- 2 # number of parametric terms: intercept, fire

nodes <- 14
lowres <- FALSE
yrmin <- 2010
yrmax <- 2022
id.min <- yrmin - 2000
id.max <- yrmax - 2000

#####################################################
######## Step 1: setting up the environment #########
#####################################################
## In this step, we just define some necessary functions and loads the coarse variable

if(lowres){
	nr <- nrow(rast('../large_domain_allpixels/input_low/VOD/cellNumber-high.tif')) #lowresolution
	blks <- 2000 #lowresolution
} else {
	nr <- nrow(rast('../large_domain_allpixels/input/VOD/cellNumber-high.tif'))
	blks <- nr
}

#auxiliary functions
source('utils/PredictMat_new.R')
g <- function(x) exp(x)
dg <- function(x) exp(x)
conc_fn <- function(X, idx){
	Xt <- X[,idx, drop = FALSE]
	Xo <- X[,-idx, drop = FALSE]

	#Q <- qr.Q(qr(Xo))
	#Xt_hat <- Q %*% (t(Q) %*% Xt)
	qr.Xo <- qr(Xo)
	Xt_hat <- qr.fitted(qr.Xo, Xt)

	return(sum(Xt_hat^2)/sum(Xt^2))
}
recvOneData_fair <- function(cl){ # this is a modified version of snow::recvOneData
	socklist <- lapply(cl, function(x) x$con)
	repeat {
		ready <- socketSelect(socklist)
		if (length(ready) > 0) break;
	}
	l <- which(ready)
	n <- ifelse(length(l) > 1, sample(l)[1], l)
	list(node = n, value = unserialize(socklist[[n]]))
}
bssize <- function(nr, nb) {
	size <- ceiling(nr/nb)
	nb <- ceiling(nr / size)

	row <- (0:(nb-1))*size + 1
	nrows <- rep(size, length(row))
	dif <- nb * size - nr
	nrows[length(nrows)] = nrows[length(nrows)] - dif
	return(list(row=row, nrows=nrows, n=nb))
}
cellFromRow <- function(object, rownr) {
	cols <- rep(1:ncol(object), times=length(rownr))
	rows <- rep(rownr, each=ncol(object))
	cellFromRowCol(object, rows, cols)
}
coordinates <- function(object) {
	xyFromCell(object, cell=1:ncell(object))
}

### Loads the tree cover fraction and NA fraction for later removal
ds <- readRDS('../large_domain_allpixels/input/TreeCover/TreeCover.rds')
as <- readRDS('../large_domain_allpixels/input/Boundary/Boundary.rds')
ys <- readRDS('../large_domain_allpixels/input/VOD/blockLatitude.rds')
pn <- readRDS('../large_domain_allpixels/input/percentNceo.rds'); pn = pn[pn[,2] < 0.9,]
rmv <- c(66302, 69006, 70210, 69914, 69019, 78564, 85395) # Some isolated grid cells
rmv <- c(rmv, 39693, 39394, 63765, 39694, 38196) # Some very influential points according to the Cook's distance

## Loads the target coarse variables Y to a single dataframe 'df.Y'
for(j in id.min:id.max){
	k <- readRDS(paste0('../large_domain_allpixels/input/VOD/1df_20', j, '.rds'))

	k$treecover <- ds[match(k[,1], ds[,2]), 1]
	k$boundary <- as[match(k[,1], as[,1]), 2]
	k$latitude <- ys[match(k[,1], ys[,1]), 2]

	nak <- which(k[,1] %in% c(pn[,1], rmv))
	if(length(nak) > 0) k = k[-nak,,drop = FALSE]

	## Below is a trick to treat block data separately in time by assigning a unique ID. ID bl = 2010, ID bl+1e5 = 2011, ID bl + 2*1e5 = 2012 etc
	k[,1] <- k[,1] + (j-10) * 1e5
	assign(paste0('df.Y', j), k)
}
df.Y <- lapply(id.min:id.max, function(i){
	get(paste0('df.Y', i))
})
df.Y <- do.call(rbind, df.Y)

## Removes some cells
na <- which(df.Y$treecover == 0 | df.Y$rmean <= 0 | df.Y$latitude > 17)
df.Y <- df.Y[-na,]

## Sets up a training/validation dataset
set.seed(1)
ids.all <- unique(df.Y$blocks %% 1e5)
ids.tra <- ids.all[sample(1:length(ids.all), round(0.75 * length(ids.all), 0))] # Training
ids.val <- setdiff(ids.all, ids.tra) # Validation

ids.tra.all <- intersect(df.Y$blocks, do.call(c, lapply(id.min:id.max, function(j) ids.tra + (j-10) * 1e5)))
ids.val.all <- setdiff(df.Y$blocks, ids.tra.all)
df.Y.tra <- df.Y[which(df.Y$blocks %in% ids.tra.all),]
df.Y.val <- df.Y[which(df.Y$blocks %in% ids.val.all),]
df.Y <- df.Y.tra


## Filters all the cells with full information
df_ids.ct <- table(df.Y$blocks %% 1e5)
df_ids.ct <- as.numeric(names(df_ids.ct)[which(df_ids.ct == 13)])
df.Y_filter <- df.Y[which(df.Y$blocks %% 1e5 %in% df_ids.ct), 'blocks']
saveRDS(df.Y_filter, 'utils/dfY_blocks_filtered.rds')

# Reference smoothers for the log-likelihood calculation
sm = readRDS('output/1arable_sm1.rds')
nsm = length(sm)
for(i in 1:nsm){
	sm[[i]]$X <- NULL
	assign(paste0('sm', i), sm[[i]])
}
bd <- sapply(1:nsm, function(i) get(paste0('sm', i))[[1]]$df)
dims <- npar + cumsum(bd)
bd <- npar + sum(bd)

mt <- 1
idx <- npar
for(i in 1:nsm){
	nl <- get(paste0('sm', i))
	for(j in 1:length(nl[[1]]$S)){
		v <- matrix(0, bd, bd)
		v[idx + (1:nl[[1]]$df), idx + (1:nl[[1]]$df)] <- nl[[1]]$S[[j]]

		assign(paste0('S', mt), v)
		mt <- mt + 1
	}
	idx <- idx + nl[[1]]$df
}
np <- mt # Number of smoothing parameters


## Splits the raster in chunks and randomize to avoid overload
bss <- bssize(nr, blks)
print(paste0('bss: ', bss$n))
d <- bss$n

# The two results
mtcs = readRDS('output/1mtcs_sol.rds')
wt <- 1+0*1/df.Y.tra$rsd^2 # Lower weight for the less reliable observations
wt <- wt/mean(wt) # This avoids modifying N in the log-likelihood
st = mtcs[[7]]
b0 = mtcs[[5]]

cts <- readRDS('output/1gen_mtc-mu.rds')
ll.t2 = readRDS('../final-high-step4/output/1mtcs_sol-var.rds')
diagL <- as.numeric(sqrt(1/ll.t2[[1]][[1]]))
Y_ <- cts[[1]]
X_ <- cts[[2]]
y. <- diagL * sqrt(wt) * Y_
X. <- diagL * sqrt(wt) * X_
X.qr <- qr(X., LAPACK = TRUE)
Qty. <- qr.qty(X.qr, y.)
R <- qr.R(X.qr, complete = TRUE)[,sort.list(X.qr$pivot)]
rho <- exp(st)
Sl <- rho[2]*S1 + rho[3]*S2 + rho[4]*S3 + rho[5]*S4 + rho[6]*S5 + rho[7]*S6 + rho[8]*S7 + rho[9]*S8 + rho[10]*S9 + rho[11]*S10 + rho[12]*S11 + rho[13]*S12 + rho[14]*S13 + rho[15]*S14 + rho[16]*S15 + rho[17]*S16
V <- crossprod(R)/rho[1] + Sl

set.seed(1)
eigenV <- eigen(V)
sqrtinv <- sqrt(1/eigenV$values)
hat.cov <- crossprod(sqrtinv * t(eigenV$vectors))
Cv <- chol(hat.cov, pivot = TRUE)
Cv <- Cv[,order(attr(Cv, 'pivot'))]
n.rep = 500 # Number of simulations, increase for more precise confidence intervals
nb <- dim(hat.cov)[1]
br <- t(Cv) %*% matrix(rnorm(n.rep*nb), nb, n.rep)
br <- as.numeric(b0) + br

## 'ag.m' is the function that calculates a chunk of the aggregated matrices based on two arguments (i = raster line; b = vector of coefficients)
### The initial vector of coefficients is a random guess, which is iteratively updated
ag.m <- function(i, br, b0){
	print(paste0('received, i = ', i, ' : ', Sys.time()))

	## Stores these two for easier access later
#	vb1 = b[[1]]
#	df.Y1 = df.Y[,1]

	bl <- readValues(block, row = bss$row[i], nrows = bss$nrow[i])
#	if(all(!bl %in% (df.Y1 %% 1e5))){
#		return(list(NULL, NULL, NULL, NULL, NULL))
#	}

	## Then, a dataframe containing their values is created
	xy <- xyFromCell(block, cellFromRow(block, rownr = bss$row[i]:(bss$row[i] + bss$nrows[i] - 1)))
	v.cci <- readValues(cci, row = bss$row[i], nrows = bss$nrow[i])
	v.map <- readValues(map, row = bss$row[i], nrows = bss$nrow[i])
	v.snd <- readValues(snd, row = bss$row[i], nrows = bss$nrow[i])
	v.soc <- readValues(soc, row = bss$row[i], nrows = bss$nrow[i])
	v.rbound <- readValues(rbound, row = bss$row[i], nrows = bss$nrow[i])

	na = (v.map == 0 | v.rbound == 0)
	if(all(na)) return(list(NULL, NULL, NULL, NULL, NULL))

	v.cci <- log(v.cci) # warning: zeroes become NA, but we want to remove them
	v.cci[!is.finite(v.cci)] <- NA

	tmp <- rep(list(NA), nrgs) ## Predictors
	lcc <- rep(list(NA), nrgs) ## Land cover, treated separately for being a time-series
	tmax <- rep(list(NA), nrgs) ## Maximum temperature, treated separately for the same reason
	mcwd <- rep(list(NA), nrgs) ## MCWD, treated separately for the same reason
	vpd <- rep(list(NA), nrgs) ## VPD, treated separately for the same reason

	for(k in seq_along(rgs)){
		# tmp
		j = rgs[k]
		yr <- 2000 + j
		v.yfi <- readValues(r.yfi[[k]], row = bss$row[i], nrows = bss$nrow[i])
		tmp[[k]] <- cbind(xy, (j - 10) * 1e5 + bl, v.map, v.cci, v.yfi, v.snd, v.soc, yr, v.rbound) # (j-10)*1e5 is a trick to have an unique id for each combination (coarse area, year)

		# tmax
		lcc[[k]] = cbind(t0 = readValues(r.lcc0[[k]], row = bss$row[i], nrows = bss$nrow[i]),
			t1 = readValues(r.lcc1[[k]], row = bss$row[i], nrows = bss$nrow[i]),
			t2 = readValues(r.lcc2[[k]], row = bss$row[i], nrows = bss$nrow[i]))

		# tmax
		tmax[[k]] = cbind(t0 = readValues(r.tem0[[k]], row = bss$row[i], nrows = bss$nrow[i]),
			t1 = readValues(r.tem1[[k]], row = bss$row[i], nrows = bss$nrow[i]),
			t2 = readValues(r.tem2[[k]], row = bss$row[i], nrows = bss$nrow[i]),
			t3 = readValues(r.tem3[[k]], row = bss$row[i], nrows = bss$nrow[i]))
		tmax[[k]][,'t0'] = tmax[[k]][,'t0'] - tmax[[k]][,'t1']
		tmax[[k]][,'t1'] = tmax[[k]][,'t1'] - tmax[[k]][,'t2']
		tmax[[k]][,'t2'] = tmax[[k]][,'t2'] - tmax[[k]][,'t3']
		tmax[[k]] = tmax[[k]][,1:3]
		
		# mcwd
		mcwd[[k]] <- cbind(t0 = readValues(r.mcw0[[k]], row = bss$row[i], nrows = bss$nrow[i]),
			t1 = readValues(r.mcw1[[k]], row = bss$row[i], nrows = bss$nrow[i]),
			t2 = readValues(r.mcw2[[k]], row = bss$row[i], nrows = bss$nrow[i]),
			t3 = readValues(r.mcw3[[k]], row = bss$row[i], nrows = bss$nrow[i]))
		mcwd[[k]][,'t0'] = mcwd[[k]][,'t0'] - mcwd[[k]][,'t1']
		mcwd[[k]][,'t1'] = mcwd[[k]][,'t1'] - mcwd[[k]][,'t2']
		mcwd[[k]][,'t2'] = mcwd[[k]][,'t2'] - mcwd[[k]][,'t3']
		mcwd[[k]] = mcwd[[k]][,1:3]

		# vpd
		vpd[[k]] <- cbind(t0 = readValues(r.vpd0[[k]], row = bss$row[i], nrows = bss$nrow[i]),
			t1 = readValues(r.vpd1[[k]], row = bss$row[i], nrows = bss$nrow[i]),
			t2 = readValues(r.vpd2[[k]], row = bss$row[i], nrows = bss$nrow[i]),
			t3 = readValues(r.vpd3[[k]], row = bss$row[i], nrows = bss$nrow[i]))
		vpd[[k]][,'t0'] = vpd[[k]][,'t0'] - vpd[[k]][,'t1']
		vpd[[k]][,'t1'] = vpd[[k]][,'t1'] - vpd[[k]][,'t2']
		vpd[[k]][,'t2'] = vpd[[k]][,'t2'] - vpd[[k]][,'t3']
		vpd[[k]] = vpd[[k]][,1:3]
	}
	tmp <- do.call(rbind, tmp)
	lcc <- do.call(rbind, lcc)
	tmax <- do.call(rbind, tmax)
	mcwd <- do.call(rbind, mcwd)
	vpd <- do.call(rbind, vpd)
	colnames(tmp) <- c('x1', 'x2', 'block', 'map',  'cci', 'yfi', 'snd', 'soc', 'year', 'rbound')
	outval <- list(NA * tmp[,1], NA * tmp[,1])

#	## The pixels without fire are convert to NoData
#	tmp[!is.na(tmp[,'yfi']) & tmp[,'yfi'] >= 19, 'yfi'] <- 19 # Some rare (around 2%) values affect the estimation, so we treat them separately.
	tmp[is.na(tmp[,'yfi']), 'yfi'] <- 999 # 999 is used for degradation as NoData

	## The pixels with desert are treated as NA
	lcc[lcc == 10 & !is.na(lcc)] <- NA

	## The groups of lcc for the yfi smoother
	tmp_lcc <- lcc[,'t0']
	lcc_yfi <- NA * tmp_lcc
	lcc_yfi[tmp_lcc == 6] <- 1
	lcc_yfi[tmp_lcc == 7] <- 2
	lcc_yfi[tmp_lcc %in% c(1, 8, 9)] <- 3
	lcc_yfi[tmp_lcc %in% c(2, 3, 4, 5)] <- 4

	## The pixels that are not in the boundaries map are removed
#	nak <- (tmp[,'map'] == 0 & !is.na(tmp[,'map'])) | (tmp[,'rbound'] == 0 & !is.na(tmp[,'rbound'])) | !(tmp[,'block'] %in% df.Y1)
	nak <- (tmp[,'map'] == 0 & !is.na(tmp[,'map'])) | (tmp[,'rbound'] == 0 & !is.na(tmp[,'rbound']))
	tmp[nak, 'block'] <- NA
	na <- (!complete.cases(tmp, tmax, mcwd, lcc))
	if(all(na)) return(list(NULL, NULL, NULL, NULL, NULL))

	## This part is a little bit tricky. The coarse-level average is calculated based on the AGC of forest pixels + non-forest pixels
	### Since we assume that AGC(bare) = 0, we have to count how many non-forest pixels we have in order to add to the calculation
	### I call this "extra zeros", reliable pixels that are removed from the model, but not from the average calculation
#	nak <- (tmp[,'map'] != 0 & !is.na(tmp[,'map'])) & (is.na(lcc_yfi) | is.na(tmp[,'cci']))
#	vd <- !is.na(tmp[, 'block'])
#	s0 <- fsum(as.integer(nak[vd]), g = tmp[vd,'block'], use.g.names = TRUE)
#	s0 <- matrix(c(as.numeric(names(s0)), s0), ncol = 2)
#	s0 <- s0[(s0[,2] > 0),, drop=FALSE]

	## I convert the dataframe to a list (for practical reasons)
	if(any(na)){
		kp <- !na
		tmp <- list(map = tmp[,'map'][kp], cci = tmp[,'cci'][kp],
			x1 = tmp[,'x1'][kp], x2 = tmp[,'x2'][kp], block = tmp[,'block'][kp],
			yfi = tmp[,'yfi'][kp], snd = tmp[,'snd'][kp],
			soc = tmp[,'soc'][kp], lcfi = lcc_yfi[kp], year = tmp[,'year'][kp],
			tmax = tmax[kp,,drop=FALSE], mcwd = mcwd[kp,,drop=FALSE],
			vpd = vpd[kp,,drop=FALSE], lcc = lcc[kp,,drop=FALSE], 
			L = matrix(rep(tmp[,'year'][kp], 3), nrow = sum(kp), ncol = 3))
	} else {
		tmp <- list(map = tmp[,'map'], cci = tmp[,'cci'], 
			x1 = tmp[,'x1'], x2 = tmp[,'x2'], block = tmp[,'block'],
			yfi = tmp[,'yfi'], snd = tmp[,'snd'],
			soc = tmp[,'soc'], lcfi = lcc_yfi, year = tmp[,'year'],
			tmax = tmax, mcwd = mcwd, vpd = vpd, lcc = lcc,
			L = matrix(rep(tmp[,'year'], 3), nrow = nrow(tmax), ncol = 3))
	}
	for(i in 2:3) tmp$L[,i] <- tmp$L[,(i-1)] - 1

	## We save the ordering of groups in the aggregation function, which allows us identify which row belongs to each coarse pixel
	blks.id <- tmp$block
	groups <- sort(unique(blks.id))

	## Then, we evaluate the smoothers at the pixel values extracted, and create the model matrix at the fine scale ('tmpgam')
	### Based on this model matrix, we derive the matrices at the coarse scale, that allow us to recover the coefficients while enforcing the disaggregation constraints

	tmpgam.X1 <- mgcv:::Predict.matrix3(sm1[[1]], tmp)$X
	tmpgam.X1 <- tpmm_mat_noblock(tmpgam.X1, qrQ1[,-1])

	tmpgam.X2 <- FastPredictMat3(sm2[[1]], tmp, qrQ2)

	tmpgam.X3 <- FastPredictMat3(sm3[[1]], tmp, qrQ3)

	tmpgam.X4 <- FastPredictMat3(sm4[[1]], tmp, qrQ4)

	tmpgam.X5 <- FastPredictMat3(sm5[[1]], tmp, qrQ5)

	# sm4 has the NoData encoded as 999, so it's handled separately
	fi_occ <- as.numeric(!(tmp$yfi == 999))
	tmp$lcfi <- tmp$lcfi
	tmpgam.X6 <- mgcv:::Predict.matrix3(sm6[[1]], tmp)$X
	tmpgam.X6 <- tpmm_mat_noblock(tmpgam.X6, qrQ6[,-1])
	tmpgam.X6 <- tmpgam.X6 * fi_occ

	tmpgam.X7 <- mgcv:::Predict.matrix3(sm7[[1]], tmp)$X
	tmpgam.X7 <- tpmm_mat_noblock(tmpgam.X7, qrQ7[,2:ncol(tmpgam.X7)])

	tmpgam.X8 <- mgcv:::Predict.matrix3(sm8[[1]], tmp)$X
	tmpgam.X8 <- tpmm_mat_noblock(tmpgam.X8, qrQ8[,2:ncol(tmpgam.X8)])

	tmpgam.X9 <- mgcv:::Predict.matrix3(sm9[[1]], tmp)$X
	tmpgam.X9 <- tpmm_mat_noblock(tmpgam.X9, qrQ9[,2:ncol(tmpgam.X9)])

	tmpgam <- cbind(1, fi_occ, tmpgam.X1, tmpgam.X2, tmpgam.X3, tmpgam.X4, tmpgam.X5, tmpgam.X6, tmpgam.X7, tmpgam.X8, tmpgam.X9)

	z.med = tpmm_vec_noblock(tmpgam, b0)
	z.sd = tpmm_vec_q_noblock(tmpgam, br, type = 4)

	if(length(na) > 0){
		outval[[1]][kp] <- g(z.med)
		outval[[2]][kp] <- z.sd
	} else {
		outval[[1]] <- g(z.med)
		outval[[2]] <- z.sd
	}
	rm(tmp)
	# Reminder: what we have is y ~ lognormal(logmean = log(outval[[1]]), logsd = outval[[2]])

	return(outval)
}
environment(ag.m) <- .GlobalEnv # Prevents slow serialization: https://stackoverflow.com/a/59117823

save(list = c('bss', 'nsm', 'npar', 'df.Y', 'ag.m', 'g', 'dg', 'cellFromRow', 'coordinates', 'lowres', 'id.min', 'id.max', paste0('sm', 1:nsm)),
	file = 'output/workers.rdata')
start_cl <- function(nodes){
	## Creates the parallel processing cluster of 'nodes' nodes
	cl <- parallel::makeCluster(nodes, type = 'PSOCK')
	clusterEvalQ(cl, load('output/workers.rdata'))
	clusterApply(cl, seq_along(cl), function(i) sink(paste0('output/worker_', i)))
	clusterEvalQ(cl, {
		library(disag.vAfrica)
		library(numDeriv)
		library(collapse)
		library(terra)
		library(mgcv)
		library(Rcpp)
		source('utils/PredictMat_new.R')

		### TODO: manual
		## This block loads the QR decompositions
		qrQ1 <- qr.Q(attr(sm1[[1]], 'qrc'), complete = TRUE)
		qrQ2 <- qr.Q(attr(sm2[[1]], 'qrc'), complete = TRUE)
		qrQ3 <- qr.Q(attr(sm3[[1]], 'qrc'), complete = TRUE)
		qrQ4 <- qr.Q(attr(sm4[[1]], 'qrc'), complete = TRUE)
		qrQ5 <- qr.Q(attr(sm5[[1]], 'qrc'), complete = TRUE)
		qrQ6 <- qr.Q(attr(sm6[[1]], 'qrc'), complete = TRUE)
		qrQ7 <- qr.Q(attr(sm7[[1]], 'qrc'), complete = TRUE)
		qrQ8 <- qr.Q(attr(sm8[[1]], 'qrc'), complete = TRUE)
		qrQ9 <- qr.Q(attr(sm9[[1]], 'qrc'), complete = TRUE)

		## This block loads all static spatial files
		if(lowres){
			block <- rast('../large_domain_allpixels/input_low/VOD/cellNumber-high.tif') #lowresolution
			cci <- rast('../large_domain_allpixels/input_low/NCEO/AGB_map_2017v0-aligned.tif') #lowresolution
			map <- rast('../large_domain_allpixels/input_low/boundaries_high.tif') #lowresolution
			snd <- rast('../large_domain_allpixels/input_low/isdaAfrica/SoilSand.tif') #lowresolution
			soc <- rast('../large_domain_allpixels/input_low/isdaAfrica/SOC.tif') #lowresolution
			rbound <- rast('../large_domain_allpixels/input_low/roughBoundary.tif') #lowresolution
		} else {
			block <- rast('../large_domain_allpixels/input/VOD/cellNumber-high.tif')
			cci <- rast('../large_domain_allpixels/input/NCEO/AGB_map_2017v0-aligned.tif')
			map <- rast('../large_domain_allpixels/input/boundaries_high.tif')
			snd <- rast('../large_domain_allpixels/input/isdaAfrica/SoilSand.tif')
			soc <- rast('../large_domain_allpixels/input/isdaAfrica/SOC.tif')
			rbound <- rast('../large_domain_allpixels/input/roughBoundary.tif') #lowresolution
		}
		readStart(block)
		readStart(cci)
		readStart(map)
		readStart(snd)
		readStart(soc)
		readStart(rbound)

		## This block loads all dynamic spatial files
		rgs <- id.min:id.max
		nrgs <- length(rgs)
		r.yfi <- rep(list(NA), nrgs)
		r.lcc0 <- rep(list(NA), nrgs)
		r.lcc1 <- rep(list(NA), nrgs)
		r.lcc2 <- rep(list(NA), nrgs)
		r.lcc3 <- rep(list(NA), nrgs)
		r.tem0 <- rep(list(NA), nrgs)
		r.tem1 <- rep(list(NA), nrgs)
		r.tem2 <- rep(list(NA), nrgs)
		r.tem3 <- rep(list(NA), nrgs)
		r.mcw0 <- rep(list(NA), nrgs)
		r.mcw1 <- rep(list(NA), nrgs)
		r.mcw2 <- rep(list(NA), nrgs)
		r.mcw3 <- rep(list(NA), nrgs)
		r.vpd0 <- rep(list(NA), nrgs)
		r.vpd1 <- rep(list(NA), nrgs)
		r.vpd2 <- rep(list(NA), nrgs)
		r.vpd3 <- rep(list(NA), nrgs)
		for(k in seq_along(rgs)){
			j = rgs[k]
			y = 2000 + j

			if(lowres){
				r.yfi[[k]] <- rast(paste0('../large_domain_allpixels/input_low/GABAM/YearsAfterFire_20', j, '.tif')) #lowresolution

				r.lcc0[[k]] <- rast(paste0('../large_domain_allpixels/input_low/LandCover/LandCover', y-0, '.tif')) #lowresolution
				r.lcc1[[k]] <- rast(paste0('../large_domain_allpixels/input_low/LandCover/LandCover', y-1, '.tif')) #lowresolution
				r.lcc2[[k]] <- rast(paste0('../large_domain_allpixels/input_low/LandCover/LandCover', y-2, '.tif')) #lowresolution
				r.lcc3[[k]] <- rast(paste0('../large_domain_allpixels/input_low/LandCover/LandCover', y-3, '.tif')) #lowresolution

				r.tem0[[k]] <- rast(paste0('../large_domain_allpixels/input_low/Temperature/MaxTemp_', y-0, '.tif'))[[1]] #lowresolution
				r.tem1[[k]] <- rast(paste0('../large_domain_allpixels/input_low/Temperature/MaxTemp_', y-1, '.tif'))[[1]] #lowresolution
				r.tem2[[k]] <- rast(paste0('../large_domain_allpixels/input_low/Temperature/MaxTemp_', y-2, '.tif'))[[1]] #lowresolution
				r.tem3[[k]] <- rast(paste0('../large_domain_allpixels/input_low/Temperature/MaxTemp_', y-3, '.tif'))[[1]] #lowresolution

				r.mcw0[[k]] <- rast(paste0('../large_domain_allpixels/input_low/TerraClimate/Max-TerraClimate_def_', y-0, '.tif')) #lowresolution
				r.mcw1[[k]] <- rast(paste0('../large_domain_allpixels/input_low/TerraClimate/Max-TerraClimate_def_', y-1, '.tif')) #lowresolution
				r.mcw2[[k]] <- rast(paste0('../large_domain_allpixels/input_low/TerraClimate/Max-TerraClimate_def_', y-2, '.tif')) #lowresolution
				r.mcw3[[k]] <- rast(paste0('../large_domain_allpixels/input_low/TerraClimate/Max-TerraClimate_def_', y-3, '.tif')) #lowresolution

				r.vpd0[[k]] <- rast(paste0('../large_domain_allpixels/input_low/TerraClimate/Max-TerraClimate_vpd_', y-0, '.tif')) #lowresolution
				r.vpd1[[k]] <- rast(paste0('../large_domain_allpixels/input_low/TerraClimate/Max-TerraClimate_vpd_', y-1, '.tif')) #lowresolution
				r.vpd2[[k]] <- rast(paste0('../large_domain_allpixels/input_low/TerraClimate/Max-TerraClimate_vpd_', y-2, '.tif')) #lowresolution
				r.vpd3[[k]] <- rast(paste0('../large_domain_allpixels/input_low/TerraClimate/Max-TerraClimate_vpd_', y-3, '.tif')) #lowresolution
			} else {
				r.yfi[[k]] <- rast(paste0('../large_domain_allpixels/input/GABAM/YearsAfterFire_20', j, '.tif'))

				r.lcc0[[k]] <- rast(paste0('../large_domain_allpixels/input/LandCover/LandCover', y-0, '.tif'))
				r.lcc1[[k]] <- rast(paste0('../large_domain_allpixels/input/LandCover/LandCover', y-1, '.tif'))
				r.lcc2[[k]] <- rast(paste0('../large_domain_allpixels/input/LandCover/LandCover', y-2, '.tif'))
				r.lcc3[[k]] <- rast(paste0('../large_domain_allpixels/input/LandCover/LandCover', y-3, '.tif'))

				r.tem0[[k]] <- rast(paste0('../large_domain_allpixels/input/Temperature/MaxTemp_', y-0, '.tif'))[[1]]
				r.tem1[[k]] <- rast(paste0('../large_domain_allpixels/input/Temperature/MaxTemp_', y-1, '.tif'))[[1]]
				r.tem2[[k]] <- rast(paste0('../large_domain_allpixels/input/Temperature/MaxTemp_', y-2, '.tif'))[[1]]
				r.tem3[[k]] <- rast(paste0('../large_domain_allpixels/input/Temperature/MaxTemp_', y-3, '.tif'))[[1]]

				r.mcw0[[k]] <- rast(paste0('../large_domain_allpixels/input/TerraClimate/Max-TerraClimate_def_', y-0, '.tif'))
				r.mcw1[[k]] <- rast(paste0('../large_domain_allpixels/input/TerraClimate/Max-TerraClimate_def_', y-1, '.tif'))
				r.mcw2[[k]] <- rast(paste0('../large_domain_allpixels/input/TerraClimate/Max-TerraClimate_def_', y-2, '.tif'))
				r.mcw3[[k]] <- rast(paste0('../large_domain_allpixels/input/TerraClimate/Max-TerraClimate_def_', y-3, '.tif'))

				r.vpd0[[k]] <- rast(paste0('../large_domain_allpixels/input/TerraClimate/Max-TerraClimate_vpd_', y-0, '.tif'))
				r.vpd1[[k]] <- rast(paste0('../large_domain_allpixels/input/TerraClimate/Max-TerraClimate_vpd_', y-1, '.tif'))
				r.vpd2[[k]] <- rast(paste0('../large_domain_allpixels/input/TerraClimate/Max-TerraClimate_vpd_', y-2, '.tif'))
				r.vpd3[[k]] <- rast(paste0('../large_domain_allpixels/input/TerraClimate/Max-TerraClimate_vpd_', y-3, '.tif'))
			}

			readStart(r.yfi[[k]])

			readStart(r.lcc0[[k]])
			readStart(r.lcc1[[k]])
			readStart(r.lcc2[[k]])
			readStart(r.lcc3[[k]])

			readStart(r.tem0[[k]])
			readStart(r.tem1[[k]])
			readStart(r.tem2[[k]])
			readStart(r.tem3[[k]])

			readStart(r.mcw0[[k]])
			readStart(r.mcw1[[k]])
			readStart(r.mcw2[[k]])
			readStart(r.mcw3[[k]])

			readStart(r.vpd0[[k]])
			readStart(r.vpd1[[k]])
			readStart(r.vpd2[[k]])
			readStart(r.vpd3[[k]])
		}

		return(TRUE)
	})
	return(cl)
}
stop_cl <- function(cl){
	clusterEvalQ(cl, {
		readStop(block)
		readStop(cci)
		readStop(map)
		readStop(snd)
		readStop(soc)
		readStop(rbound)
		for(k in seq_along(rgs)){
			readStop(r.yfi[[k]])

			readStop(r.lcc0[[k]])
			readStop(r.lcc1[[k]])
			readStop(r.lcc2[[k]])
			readStop(r.lcc3[[k]])

			readStop(r.tem0[[k]])
			readStop(r.tem1[[k]])
			readStop(r.tem2[[k]])
			readStop(r.tem3[[k]])

			readStop(r.mcw0[[k]])
			readStop(r.mcw1[[k]])
			readStop(r.mcw2[[k]])
			readStop(r.mcw3[[k]])

			readStop(r.vpd0[[k]])
			readStop(r.vpd1[[k]])
			readStop(r.vpd2[[k]])
			readStop(r.vpd3[[k]])
		}
		return(TRUE)
	})
	parallel::stopCluster(cl)
	return()
}

exporta <- function(){
	for(i in 1:nodes){
		sendCall(cl[[i]], ag.m, list(i = i, br = br, b0 = b0), tag = i)
	}

	for (i in 1:bss$n) {
		d <- recvOneData_fair(cl)
		if (!d$value$success) {
			saveRDS(d, 'erro.rds')
			cat('erro no numero: ', d$value$tag, '\n'); flush.console();
			stop('cluster error')
		}

		ni <- nodes + i
		if (ni <= bss$n){
			sendCall(cl[[d$node]], ag.m, list(i = ni, br = br, b0 = b0), tag = ni)
		}
		b <- d$value$tag
		print(paste0('recebido: ', b))

		saveRDS(d$value$value, paste0('output/tifs-', ifelse(lowres, 'low', 'high'), '/v_', b))
		rm(d)
	}
}
cl <- start_cl(nodes)
exporta()
stop_cl(cl)
