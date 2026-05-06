################### Script: 4export.R ##################################
## The following lines load the necessary libraries
library(disag.vAfrica)
library(numDeriv)
library(collapse)
library(terra)
library(mgcv)
library(snow)
library(Rcpp)

npar <- 2
npar_v <- 1
np_var <- 7
nodes <- 15
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
	mtcs.t1 <- readRDS('../final-high-step1/output-low/1gen_mtc-mu.rds')
} else {
	mtcs.t1 <- readRDS('../final-high-step1/output/1gen_mtc-mu.rds')
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

if(lowres){
	nr <- nrow(rast('../large_domain_allpixels/input_low/VOD/cellNumber-high.tif')) #lowresolution
	blks <- 2000 #lowresolution
} else {
	nr <- nrow(rast('../large_domain_allpixels/input/VOD/cellNumber-high.tif'))
	blks <- 45000
}

#auxiliary functions
source('utils/PredictMat_new.R')
g <- function(x) exp(x)
dg <- function(x) exp(x)

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

#auxiliary functions
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

# Reference smoothers for the log-likelihood calculation
if(lowres){
	sm = readRDS('../final-high-step1/output-low/1arable_sm1.rds')
} else {
	sm = readRDS('../final-high-step1/output/1arable_sm1.rds')
}
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

#reference smoothers for the variance
if(lowres){
	f1_mtcs_sol <- readRDS('../final-high-step1/output-low/1mtcs_sol.rds')
} else {
	f1_mtcs_sol <- readRDS('../final-high-step1/output/1mtcs_sol.rds')
}
st = f1_mtcs_sol[[7]]
rho = exp(st)

if(lowres){
	tmp.sm <- readRDS('../final-high-step1/output-low/1tmp_sm.rds')
} else {
	tmp.sm <- readRDS('../final-high-step1/output/1tmp_sm.rds')
}
tmp.sm <- na.omit(tmp.sm)
tmp.sm$mcwd <- tmp.sm$mcwd[,1]
tmp.sm$vpd <- tmp.sm$vpd[,1]
tmp.sm$tmax <- tmp.sm$tmax[,1]
sm_var1 <- smoothCon(te(x1, x2, bs = c('cs', 'cs'), m = 2, k = c(8, 8)), data = tmp.sm, absorb.cons = TRUE) # s(long, lat)
sm_var1[[1]]$X <- sm_var1[[1]]$margin[[1]]$X <- sm_var1[[1]]$margin[[2]]$X <- NULL
sm_var2 <- smoothCon(s(mcwd, bs = 'cs', m = 2, k = 12), data = tmp.sm, absorb.cons = TRUE) # s(mcwd)
sm_var3 <- smoothCon(s(vpd, bs = 'cs', m = 2, k = 12), data = tmp.sm, absorb.cons = TRUE) # s(mcwd)
sm_var4 <- smoothCon(s(tmax, bs = 'cs', m = 2, k = 12), data = tmp.sm, absorb.cons = TRUE) # s(tmax)
sm_var5 <- smoothCon(s(cci, bs = 'cs', m = 2, k = 12), data = tmp.sm, absorb.cons = TRUE) # s(cci)
sm_var6 <- smoothCon(s(lcc, bs = 'cs', m = 2, k = 9), data = tmp.sm, absorb.cons = TRUE) # s(lcc)

bd_var <- c(sm_var1[[1]]$df, sm_var2[[1]]$df, sm_var3[[1]]$df, sm_var4[[1]]$df, sm_var5[[1]]$df, sm_var6[[1]]$df)
dims_var <- npar_v + cumsum(bd_var)
bd_var <- npar_v + sum(bd_var)
Sv1 <- Sv2 <- Sv3 <- Sv4 <- Sv5 <- Sv6 <- Sv7 <- matrix(0, bd_var, bd_var); idx <- npar_v
Sv1[idx + (1:sm_var1[[1]]$df), idx + (1:sm_var1[[1]]$df)] <- sm_var1[[1]]$S[[1]]
Sv2[idx + (1:sm_var1[[1]]$df), idx + (1:sm_var1[[1]]$df)] <- sm_var1[[1]]$S[[2]]; idx = idx + sm_var1[[1]]$df
Sv3[idx + (1:sm_var2[[1]]$df), idx + (1:sm_var2[[1]]$df)] <- sm_var2[[1]]$S[[1]]; idx = idx + sm_var2[[1]]$df
Sv4[idx + (1:sm_var3[[1]]$df), idx + (1:sm_var3[[1]]$df)] <- sm_var3[[1]]$S[[1]]; idx = idx + sm_var3[[1]]$df
Sv5[idx + (1:sm_var4[[1]]$df), idx + (1:sm_var4[[1]]$df)] <- sm_var4[[1]]$S[[1]]; idx = idx + sm_var4[[1]]$df
Sv6[idx + (1:sm_var5[[1]]$df), idx + (1:sm_var5[[1]]$df)] <- sm_var5[[1]]$S[[1]]; idx = idx + sm_var5[[1]]$df
Sv7[idx + (1:sm_var6[[1]]$df), idx + (1:sm_var6[[1]]$df)] <- sm_var6[[1]]$S[[1]]

## Splits the raster in chunks
bss <- bssize(nr, blks)
print(paste0('bss: ', bss$n))
d <- bss$n

## 'ag.m' is the function that calculates a chunk of the aggregated matrices based on two arguments (i = raster line; b = vector of coefficients)
### The initial vector of coefficients is a random guess, which is iteratively updated
ag.m <- function(i, b, v){
	print(paste0('received, i = ', i, ' : ', Sys.time()))

#	df.Y1 = df.Y[,1]

	## Then, a dataframe containing their values is created
	bl <- readValues(block, row = bss$row[i], nrows = bss$nrow[i])
	v.map <- readValues(map, row = bss$row[i], nrows = bss$nrow[i])
	v.cci <- readValues(cci, row = bss$row[i], nrows = bss$nrow[i])
	v.rbound <- readValues(rbound, row = bss$row[i], nrows = bss$nrow[i])
	v.cci = log(v.cci)
	v.cci[!is.finite(v.cci)] <- NA

	na = (v.map == 0 | v.rbound == 0)
	if(sum(na) == length(bl)) return(list(NULL, NULL, NULL))

	xy <- xyFromCell(block, cellFromRow(block, rownr = bss$row[i]:(bss$row[i] + bss$nrows[i] - 1)))

	tmp <- rep(list(NA), nrgs) ## Predictors
	lcc <- rep(list(NA), nrgs) ## Land cover, treated separately for being a time-series
	tmax <- rep(list(NA), nrgs) ## Maximum temperature, treated separately for the same reason
	mcwd <- rep(list(NA), nrgs) ## MCWD, treated separately for the same reason
	vpd <- rep(list(NA), nrgs) ## VPD, treated separately for the same reason

	for(k in seq_along(rgs)){
		# tmp
		j = rgs[k]
		tmp[[k]] <- cbind(xy, (j - 10) * 1e5 + bl, v.cci, v.map, v.rbound) # (j-10)*1e5 is a trick to have an unique id for each combination (coarse area, year)

		# lcc, tmax, mcwd, vpd
		lcc[[k]] = readValues(r.lcc[[k]], row = bss$row[i], nrows = bss$nrow[i])
		tmax[[k]] = readValues(r.tem0[[k]], row = bss$row[i], nrows = bss$nrow[i]) - readValues(r.tem1[[k]], row = bss$row[i], nrows = bss$nrow[i])
		mcwd[[k]] = readValues(r.mcw0[[k]], row = bss$row[i], nrows = bss$nrow[i]) - readValues(r.mcw1[[k]], row = bss$row[i], nrows = bss$nrow[i])
		vpd[[k]] = readValues(r.vpd0[[k]], row = bss$row[i], nrows = bss$nrow[i]) - readValues(r.vpd1[[k]], row = bss$row[i], nrows = bss$nrow[i])
	}	
	tmp <- do.call(rbind, tmp)
	lcc <- do.call(c, lcc)
	tmax <- do.call(c, tmax)
	mcwd <- do.call(c, mcwd)
	vpd <- do.call(c, vpd)
	colnames(tmp) <- c('x1', 'x2', 'block', 'cci', 'map', 'rbound')
	outval <- NA * tmp[,1]

	## The pixels with desert are treated as NA
	lcc[lcc == 10 & !is.na(lcc)] <- NA

	## The pixels that are not in the boundaries map are removed
#	nak <- (tmp[,'map'] == 0 & !is.na(tmp[,'map'])) | !(tmp[,'block'] %in% df.Y1)
	nak <- (tmp[,'map'] == 0 & !is.na(tmp[,'map'])) | (tmp[,'rbound'] == 0 & !is.na(tmp[,'rbound']))
	tmp[nak, 'block'] <- NA
	na <- (!complete.cases(tmp, tmax, mcwd, lcc))
	if(all(na)) return(list(NULL, NULL, NULL))

	## I convert the dataframe to a list (for practical reasons)
	if(sum(na) > 0){
		kp <- !na
		tmp <- list(x1 = tmp[,'x1'][kp], x2 = tmp[,'x2'][kp], block = tmp[,'block'][kp], cci = tmp[,'cci'][kp],
			lcc = lcc[kp], mcwd = mcwd[kp], tmax = tmax[kp], vpd = vpd[kp])
	} else {
		tmp <- list(x1 = tmp[,'x1'], x2 = tmp[,'x2'], block = tmp[,'block'], cci = tmp[,'cci'], 
			lcc = lcc, mcwd = mcwd, tmax = tmax, vpd = vpd)
	}
	
	## We save the ordering of groups in the aggregation function, which allows us identify which row belongs to each coarse pixel
	groups <- sort(unique(tmp$block))
	vk <- v[[1]]

	## Then, we evaluate the smoothers at the pixel values extracted, and create the model matrix at the fine scale ('tmpgam')
	### Based on this model matrix, we derive the matrices at the coarse scale, that allow us to recover the coefficients while enforcing the disaggregation constraints
	tmpgam_var1 <- mgcv:::Predict.matrix3(sm_var1[[1]], tmp)$X
	tmpgam_var1 <- tmpgam_var1 %*% qrQv1[,-1]
	i_st <- npar_v + 1
	i_en <- i_st + ncol(tmpgam_var1) - 1
	z. <- tpmv_vec_noblock(tmpgam_var1, vk[i_st:i_en])
	
	tmpgam_var2 <- mgcv:::Predict.matrix3(sm_var2[[1]], tmp)$X
	tmpgam_var2 <- tmpgam_var2 %*% qrQv2[,-1]
	i_st <- i_en + 1
	i_en <- i_st + ncol(tmpgam_var2) - 1
	z. <- tpmv_vec_noblock(tmpgam_var2, vk[i_st:i_en])

	tmpgam_var3 <- mgcv:::Predict.matrix3(sm_var3[[1]], tmp)$X
	tmpgam_var3 <- tmpgam_var3 %*% qrQv3[,-1]
	i_st <- i_en + 1
	i_en <- i_st + ncol(tmpgam_var3) - 1
	z. <- tpmv_vec_noblock(tmpgam_var3, vk[i_st:i_en])

	tmpgam_var4 <- mgcv:::Predict.matrix3(sm_var4[[1]], tmp)$X
	tmpgam_var4 <- tmpgam_var4 %*% qrQv4[,-1]
	i_st <- i_en + 1
	i_en <- i_st + ncol(tmpgam_var4) - 1
	z. <- tpmv_vec_noblock(tmpgam_var4, vk[i_st:i_en])

	tmpgam_var5 <- mgcv:::Predict.matrix3(sm_var5[[1]], tmp)$X
	tmpgam_var5 <- tmpgam_var5 %*% qrQv5[,-1]
	i_st <- i_en + 1
	i_en <- i_st + ncol(tmpgam_var5) - 1
	z. <- tpmv_vec_noblock(tmpgam_var5, vk[i_st:i_en])
	
	tmpgam_var6 <- mgcv:::Predict.matrix3(sm_var6[[1]], tmp)$X
	tmpgam_var6 <- tmpgam_var6 %*% qrQv6[,-1]
	i_st <- i_en + 1
	i_en <- i_st + ncol(tmpgam_var6) - 1
	z. <- tpmv_vec_noblock(tmpgam_var6, vk[i_st:i_en])
	
	z. <- z. + vk[1] # Adds the intercept

	vk. <- g(z.)

	if(sum(na) > 0){
		outval[!na] <- vk.
	} else {
		outval <- vk.
	}

	outval <- outval * rho[1]

	return(outval)
}

#####################################################
####### Step 2: handles the parallel cluster ########
#####################################################
save(list = c('bss', 'df.Y', 'rho', 'sm_var1', 'sm_var2', 'sm_var3', 'sm_var4', 'sm_var5', 'sm_var6', 'ag.m', 'cellFromRow', 'coordinates', 'lowres', 'id.min', 'id.max', 'g', 'dg', 'npar_v'),
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

		qrQv1 <- qr.Q(attr(sm_var1[[1]], 'qrc'), complete = TRUE)
		qrQv2 <- qr.Q(attr(sm_var2[[1]], 'qrc'), complete = TRUE)
		qrQv3 <- qr.Q(attr(sm_var3[[1]], 'qrc'), complete = TRUE)
		qrQv4 <- qr.Q(attr(sm_var4[[1]], 'qrc'), complete = TRUE)
		qrQv5 <- qr.Q(attr(sm_var5[[1]], 'qrc'), complete = TRUE)
		qrQv6 <- qr.Q(attr(sm_var6[[1]], 'qrc'), complete = TRUE)


		## This block loads all static spatial files
		if(lowres){
			block <- rast('../large_domain_allpixels/input_low/VOD/cellNumber-high.tif') #lowresolution
			map <- rast('../large_domain_allpixels/input_low/boundaries_high.tif') #lowresolution
			cci <- rast('../large_domain_allpixels/input_low/NCEO/AGB_map_2017v0-aligned.tif') #lowresolution
			rbound <- rast('../large_domain_allpixels/input_low/roughBoundary.tif') #lowresolution
		} else {
			block <- rast('../large_domain_allpixels/input/VOD/cellNumber-high.tif')
			map <- rast('../large_domain_allpixels/input/boundaries_high.tif')
			cci <- rast('../large_domain_allpixels/input/NCEO/AGB_map_2017v0-aligned.tif')
			rbound <- rast('../large_domain_allpixels/input/roughBoundary.tif')
		}
		readStart(block)
		readStart(map)
		readStart(cci)
		readStart(rbound)

		## This block loads all dynamic spatial files
		rgs <- id.min:id.max
		nrgs <- length(rgs)
		r.lcc <- rep(list(NA), nrgs)
		r.tem0 <- rep(list(NA), nrgs)
		r.tem1 <- rep(list(NA), nrgs)
		r.mcw0 <- rep(list(NA), nrgs)
		r.mcw1 <- rep(list(NA), nrgs)
		r.vpd0 <- rep(list(NA), nrgs)
		r.vpd1 <- rep(list(NA), nrgs)
		for(k in seq_along(rgs)){
			j = rgs[k]
			y = 2000 + j

			if(lowres){
				r.lcc[[k]] <- rast(paste0('../large_domain_allpixels/input_low/LandCover/LandCover', y-0, '.tif')) #lowresolution
				r.tem0[[k]] <- rast(paste0('../large_domain_allpixels/input_low/Temperature/MaxTemp_', y-0, '.tif'))[[1]] #lowresolution
				r.tem1[[k]] <- rast(paste0('../large_domain_allpixels/input_low/Temperature/MaxTemp_', y-1, '.tif'))[[1]] #lowresolution
				r.mcw0[[k]] <- rast(paste0('../large_domain_allpixels/input_low/TerraClimate/Max-TerraClimate_def_', y-0, '.tif')) #lowresolution
				r.mcw1[[k]] <- rast(paste0('../large_domain_allpixels/input_low/TerraClimate/Max-TerraClimate_def_', y-1, '.tif')) #lowresolution
				r.vpd0[[k]] <- rast(paste0('../large_domain_allpixels/input_low/TerraClimate/Max-TerraClimate_vpd_', y-0, '.tif')) #lowresolution
				r.vpd1[[k]] <- rast(paste0('../large_domain_allpixels/input_low/TerraClimate/Max-TerraClimate_vpd_', y-1, '.tif')) #lowresolution
			} else {
				r.lcc[[k]] <- rast(paste0('../large_domain_allpixels/input/LandCover/LandCover', y-0, '.tif'))
				r.tem0[[k]] <- rast(paste0('../large_domain_allpixels/input/Temperature/MaxTemp_', y-0, '.tif'))[[1]]
				r.tem1[[k]] <- rast(paste0('../large_domain_allpixels/input/Temperature/MaxTemp_', y-1, '.tif'))[[1]]
				r.mcw0[[k]] <- rast(paste0('../large_domain_allpixels/input/TerraClimate/Max-TerraClimate_def_', y-0, '.tif'))
				r.mcw1[[k]] <- rast(paste0('../large_domain_allpixels/input/TerraClimate/Max-TerraClimate_def_', y-1, '.tif'))
				r.vpd0[[k]] <- rast(paste0('../large_domain_allpixels/input/TerraClimate/Max-TerraClimate_vpd_', y-0, '.tif'))
				r.vpd1[[k]] <- rast(paste0('../large_domain_allpixels/input/TerraClimate/Max-TerraClimate_vpd_', y-1, '.tif'))
			}
			readStart(r.lcc[[k]])
			readStart(r.tem0[[k]])
			readStart(r.tem1[[k]])
			readStart(r.mcw0[[k]])
			readStart(r.mcw1[[k]])
			readStart(r.vpd0[[k]])
			readStart(r.vpd1[[k]])
		}

		return(TRUE)
	})
	return(cl)
}
stop_cl <- function(cl){
	clusterEvalQ(cl, {
		readStop(block)
		readStop(map)
		readStop(cci)
		readStop(rbound)
		for(k in seq_along(rgs)){
			readStop(r.lcc[[k]])
			readStop(r.tem0[[k]])
			readStop(r.tem1[[k]])
			readStop(r.mcw0[[k]])
			readStop(r.mcw1[[k]])
			readStop(r.vpd0[[k]])
			readStop(r.vpd1[[k]])
		}
		return(TRUE)
	})
	parallel::stopCluster(cl)
	return()
}

#####################################################
################# Step 3: exports ###################
#####################################################
exporta <- function(){
	for(i in 1:nodes){
		sendCall(cl[[i]], ag.m, list(i = i, b = NULL, v = list(v0)), tag = i)
	}

	for (i in 1:bss$n) {
		d <- recvOneData(cl)
		if (!d$value$success) {
			saveRDS(d, 'erro.rds')
			cat('erro no numero: ', d$value$tag, '\n'); flush.console();
			stop('cluster error')
		}

		ni <- nodes + i
		if (ni <= bss$n){
			sendCall(cl[[d$node]], ag.m, list(i = ni, b = NULL, v = list(v0)), tag = ni)
		}
		b <- d$value$tag
		print(paste0('recebido: ', b))
		saveRDS(d$value$value, paste0('output/tifs-', ifelse(lowres, 'low', 'high'), '/v_', b))
		rm(d)
	}
}
rt <- readRDS('output/1mtcs_sol-var.rds')
v0 <- as.numeric(rt[[6]])

cl <- start_cl(nodes)
exporta()
stop_cl(cl)
