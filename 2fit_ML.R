################### Script: 2fit_ML.R ##################################
## The following lines load the necessary libraries
library(disag.vAfrica)
library(numDeriv)
library(terra)
library(mgcv)
library(snow)
library(Rcpp)

lowres <- TRUE
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
	blks <- 1000 #lowresolution
} else {
	nr <- nrow(rast('../large_domain_allpixels/input/VOD/cellNumber-high.tif'))
	blks <- 12000
}

#auxiliary functions
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
source('utils/PredictMat_new.R')

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

## Loads the target coarse variables Y to a single dataframe 'df.Y'
for(j in id.min:id.max){
	k <- readRDS(paste0('../large_domain_allpixels/input/VOD/1df_20', j, '.rds'))
	k <- k[,c('blocks', 'rmean')]
	## Below is a trick to treat block data separately in time by assigning a unique ID. ID bl = 2010, ID bl+1e5 = 2011, ID bl + 2*1e5 = 2012 etc
	k[,1] <- k[,1] + (j-10) * 1e5
	assign(paste0('df.Y', j), k)
}
df.Y <- lapply(id.min:id.max, function(i){
	get(paste0('df.Y', i))
})
df.Y <- do.call(rbind, df.Y)
na <- which(df.Y$rmean <= 0) # Remove the VOD cells with negative AGC
if(length(na) > 0) df.Y <- df.Y[-na,]

if(lowres){
	if(file.exists('pk.rds')){
		df.Y <- df.Y[-readRDS('pk.rds'),] #lowresolution, pk is pre-computed
	}
}

# Reference smoothers for the log-likelihood calculation
sm = readRDS('output/1arable_sm1.rds')
sm1 = sm[[1]]; sm2 = sm[[2]]; sm3 = sm[[3]]; sm4 = sm[[4]]; sm5 = sm[[5]]
bd <- sapply(1:4, function(i) get(paste0('sm', i))[[1]]$df)
bd <- c(bd, rep(sm5[[1]]$df, 4))
dims <- 1 + cumsum(bd)
bd <- 1 + sum(bd) # 1 = intercept
S1 <- S2 <- S3 <- S4 <- S5 <- S6 <- S7 <- S8 <- S9 <- S10 <- matrix(0, bd, bd); idx <- 1
S1[idx + (1:sm1[[1]]$df), idx + (1:sm1[[1]]$df)] <- sm1[[1]]$S[[1]]
S2[idx + (1:sm1[[1]]$df), idx + (1:sm1[[1]]$df)] <- sm1[[1]]$S[[2]]; idx <- idx + sm1[[1]]$df
S3[idx + (1:sm2[[1]]$df), idx + (1:sm2[[1]]$df)] <- sm2[[1]]$S[[1]]; idx <- idx + sm2[[1]]$df
S4[idx + (1:sm3[[1]]$df), idx + (1:sm3[[1]]$df)] <- sm3[[1]]$S[[1]]
S5[idx + (1:sm3[[1]]$df), idx + (1:sm3[[1]]$df)] <- sm3[[1]]$S[[2]]; idx <- idx + sm3[[1]]$df
S6[idx + (1:sm4[[1]]$df), idx + (1:sm4[[1]]$df)] <- sm4[[1]]$S[[1]]; idx <- idx + sm4[[1]]$df
S7[idx + (1:sm5[[1]]$df), idx + (1:sm5[[1]]$df)] <- sm5[[1]]$S[[1]]; idx <- idx + sm5[[1]]$df
S8[idx + (1:sm5[[1]]$df), idx + (1:sm5[[1]]$df)] <- sm5[[1]]$S[[1]]; idx <- idx + sm5[[1]]$df
S9[idx + (1:sm5[[1]]$df), idx + (1:sm5[[1]]$df)] <- sm5[[1]]$S[[1]]; idx <- idx + sm5[[1]]$df
S10[idx + (1:sm5[[1]]$df), idx + (1:sm5[[1]]$df)] <- sm5[[1]]$S[[1]]

## Splits the raster in chunks and randomize to avoid overload
bss <- bssize(nr, blks)
print(paste0('bss: ', bss$n))
d <- bss$n

## 'ag.m' is the function that calculates a chunk of the aggregated matrices based on two arguments (i = raster line; b = vector of coefficients)
### The initial vector of coefficients is a random guess, which is iteratively updated
ag.m <- function(i, b){
	print(paste0('received, i = ', i, ' : ', Sys.time()))

	## This block loads all static spatial files
	if(lowres){
		block <- rast('../large_domain_allpixels/input_low/VOD/cellNumber-high.tif') #lowresolution
		cci <- rast('../large_domain_allpixels/input_low/CCIv6/CCI2010_AGB-t_ha.tif') #lowresolution
		map <- rast('../large_domain_allpixels/input_low/boundaries_high.tif') #lowresolution
	} else {
		block <- rast('../large_domain_allpixels/input/VOD/cellNumber-high.tif')
		cci <- rast('../large_domain_allpixels/input/CCIv6/CCI2010_AGB-t_ha.tif')
		map <- rast('../large_domain_allpixels/input/boundaries_high.tif')
	}
	readStart(block)
	readStart(cci)
	readStart(map)

	## This block loads all dynamic spatial files
	for(j in id.min:id.max){
		y = as.numeric(paste0('20', j))

		if(lowres){
			assign(paste0('lcc', j), rast(paste0('../large_domain_allpixels/input_low/LandCover/LandCover20', j, '.tif'))) #lowresolution
			assign(paste0('yfi', j), rast(paste0('../large_domain_allpixels/input_low/GABAM/YearsAfterFire_20', j, '.tif'))) #lowresolution
			assign(paste0('ydg', j), rast(paste0('../large_domain_allpixels/input_low/TMF/YearsAfterDegradation_20', j, '.tif'))) #lowresolution
			assign(paste0('lcc', j), rast(paste0('../large_domain_allpixels/input_low/LandCover/LandCover20', j, '.tif'))) #lowresolution
			assign(paste0('tem0', j), rast(paste0('../large_domain_allpixels/input_low/Temperature/MaxTemp_', y-0, '.tif'))[[1]]) #lowresolution
			assign(paste0('tem1', j), rast(paste0('../large_domain_allpixels/input_low/Temperature/MaxTemp_', y-1, '.tif'))[[1]]) #lowresolution
			assign(paste0('tem2', j), rast(paste0('../large_domain_allpixels/input_low/Temperature/MaxTemp_', y-2, '.tif'))[[1]]) #lowresolution
			assign(paste0('tem3', j), rast(paste0('../large_domain_allpixels/input_low/Temperature/MaxTemp_', y-3, '.tif'))[[1]]) #lowresolution
			assign(paste0('mcw0', j), rast(paste0('../large_domain_allpixels/input_low/MCWD/MCWD_', y-0, '.tif'))) #lowresolution
			assign(paste0('mcw1', j), rast(paste0('../large_domain_allpixels/input_low/MCWD/MCWD_', y-1, '.tif'))) #lowresolution
			assign(paste0('mcw2', j), rast(paste0('../large_domain_allpixels/input_low/MCWD/MCWD_', y-2, '.tif'))) #lowresolution
			assign(paste0('mcw3', j), rast(paste0('../large_domain_allpixels/input_low/MCWD/MCWD_', y-3, '.tif'))) #lowresolution
		} else {
			assign(paste0('lcc', j), rast(paste0('../large_domain_allpixels/input/LandCover/LandCover20', j, '.tif')))
			assign(paste0('yfi', j), rast(paste0('../large_domain_allpixels/input/GABAM/YearsAfterFire_20', j, '.tif')))
			assign(paste0('ydg', j), rast(paste0('../large_domain_allpixels/input/TMF/YearsAfterDegradation_20', j, '.tif')))
			assign(paste0('tem0', j), rast(paste0('../large_domain_allpixels/input/Temperature/MaxTemp_', y-0, '.tif'))[[1]])
			assign(paste0('tem1', j), rast(paste0('../large_domain_allpixels/input/Temperature/MaxTemp_', y-1, '.tif'))[[1]])
			assign(paste0('tem2', j), rast(paste0('../large_domain_allpixels/input/Temperature/MaxTemp_', y-2, '.tif'))[[1]])
			assign(paste0('tem3', j), rast(paste0('../large_domain_allpixels/input/Temperature/MaxTemp_', y-3, '.tif'))[[1]])
			assign(paste0('mcw0', j), rast(paste0('../large_domain_allpixels/input/MCWD/MCWD_', y-0, '.tif')))
			assign(paste0('mcw1', j), rast(paste0('../large_domain_allpixels/input/MCWD/MCWD_', y-1, '.tif')))
			assign(paste0('mcw2', j), rast(paste0('../large_domain_allpixels/input/MCWD/MCWD_', y-2, '.tif')))
			assign(paste0('mcw3', j), rast(paste0('../large_domain_allpixels/input/MCWD/MCWD_', y-3, '.tif')))
		}

		readStart(get(paste0('lcc', j)))
		readStart(get(paste0('yfi', j)))
		readStart(get(paste0('ydg', j)))
		readStart(get(paste0('tem0', j)))
		readStart(get(paste0('tem1', j)))
		readStart(get(paste0('tem2', j)))
		readStart(get(paste0('tem3', j)))
		readStart(get(paste0('mcw0', j)))
		readStart(get(paste0('mcw1', j)))
		readStart(get(paste0('mcw2', j)))
		readStart(get(paste0('mcw3', j)))
	}

	on.exit({
		readStop(cci)
		readStop(block)
		readStop(map)
		for(j in id.min:id.max){
			readStop(get(paste0('lcc', j)))
			readStop(get(paste0('yfi', j)))
			readStop(get(paste0('ydg', j)))
			readStop(get(paste0('tem0', j)))
			readStop(get(paste0('tem1', j)))
			readStop(get(paste0('tem2', j)))
			readStop(get(paste0('tem3', j)))
			readStop(get(paste0('mcw0', j)))
			readStop(get(paste0('mcw1', j)))
			readStop(get(paste0('mcw2', j)))
			readStop(get(paste0('mcw3', j)))
		}

		rm(list=ls())
		gc()
	})

	## Then, a dataframe containing their values is created
	bl <- readValues(block, row = bss$row[i], nrows = bss$nrow[i])
	xy <- xyFromCell(block, cellFromRow(block, rownr = bss$row[i]:(bss$row[i] + bss$nrows[i] - 1)))
	v.cci <- readValues(cci, row = bss$row[i], nrows = bss$nrow[i])
	v.map <- readValues(map, row = bss$row[i], nrows = bss$nrow[i])
	na = (v.map == 0)
	if(sum(na) == length(bl)) return(list(NULL, NULL, NULL, NULL, NULL))

	tmp <- lapply(id.min:id.max, function(j){
		v.yfi <- readValues(get(paste0('yfi', j)), row = bss$row[i], nrows = bss$nrow[i], mat = TRUE)
		v.ydg <- readValues(get(paste0('ydg', j)), row = bss$row[i], nrows = bss$nrow[i], mat = TRUE)
		v.lcc <- readValues(get(paste0('lcc', j)), row = bss$row[i], nrows = bss$nrow[i], mat = TRUE)

		yr <- rep(2000 + j, length(v.lcc))
		mat <- cbind(xy, (j-10)*1e5 + bl, v.map, v.lcc, v.cci, v.ydg, v.yfi, yr) # (j-10)*1e5 is a trick to have an unique id for each combination (coarse area, year)
		colnames(mat) <- c('x1', 'x2', 'block', 'map', 'lcc', 'cci', 'ydg', 'yfi', 'year')
		mat
	})
	tmp <- do.call(rbind, tmp)

	## Maximum temperature is treated separately for being a time-series
	tmax <- lapply(id.min:id.max, function(j){
		v.tem = data.frame(t0 = readValues(get(paste0('tem0', j)), row = bss$row[i], nrows = bss$nrow[i]),
		t1 = readValues(get(paste0('tem1', j)), row = bss$row[i], nrows = bss$nrow[i]),
		t2 = readValues(get(paste0('tem2', j)), row = bss$row[i], nrows = bss$nrow[i]),
		t3 = readValues(get(paste0('tem3', j)), row = bss$row[i], nrows = bss$nrow[i]))
		v.tem
	})
	tmax <- do.call(rbind, tmax)
	L <- 1 + 0 * tmax

	## MCWD is treated separately for the same reason
	mcwd <- lapply(id.min:id.max, function(j){
		v.mcw <- data.frame(t0 = readValues(get(paste0('mcw0', j)), row = bss$row[i], nrows = bss$nrow[i]),
		t1 = readValues(get(paste0('mcw1', j)), row = bss$row[i], nrows = bss$nrow[i]),
		t2 = readValues(get(paste0('mcw2', j)), row = bss$row[i], nrows = bss$nrow[i]),
		t3 = readValues(get(paste0('mcw3', j)), row = bss$row[i], nrows = bss$nrow[i]))
		v.mcw <- log(1 + abs(v.mcw))
	})
	mcwd <- do.call(rbind, mcwd)

	## The pixels without degradation are convert to NoData
	fx.ydg.1 <- which(tmp[,'ydg'] >= 24) # Some rare (around 2%) values affect the estimation, so we treat them separately.
	if(length(fx.ydg.1) > 0) tmp[fx.ydg.1,'ydg'] <- 24
	fx.ydg.2 <- which(is.na(tmp[,'ydg']))
	if(length(fx.ydg.2) > 0) tmp[fx.ydg.2, 'ydg'] <- 999 # 999 is used for degradation as NoData
	
	## The pixels without fire are convert to NoData
	fx.yfi.1 <- which(tmp[,'yfi'] >= 19) # Some rare (around 2%) values affect the estimation, so we treat them separately.
	if(length(fx.yfi.1) > 0) tmp[fx.yfi.1,'yfi'] <- 19
	fx.yfi.2 <- which(is.na(tmp[,'yfi']))
	if(length(fx.yfi.2) > 0) tmp[fx.yfi.2, 'yfi'] <- 999 # 999 is used for degradation as NoData

	## The pixels with desert are treated as NA
	fx.lcc.1 <- which(tmp[,'lcc'] == 8)
	if(length(fx.lcc.1) > 0) tmp[fx.lcc.1, 'lcc'] <- NA
	fx.lcc.2 <- which(tmp[,'lcc'] == 9)
	if(length(fx.lcc.2) > 0) tmp[fx.lcc.2, 'lcc'] <- 8 # Deciduous forest reclassified to 8

	## The pixels that are not in the boundaries map are removed
	nak <- which(tmp[,'map'] == 0)
	if(length(nak) > 0) tmp[nak, 'block'] <- NA
	na <- getRowIndex_NA(as.matrix(cbind(tmp, tmax, mcwd)))
	if(length(na) == nrow(tmp)) return(list(NULL, NULL, NULL, NULL, NULL))

	## This part is a little bit tricky. The coarse-level average is calculated based on the AGC of forest pixels + non-forest pixels
	### Since we assume that AGC(bare) = 0, we have to count how many non-forest pixels we have in order to add to the calculation
	### I call this "extra zeros", reliable pixels that are removed from the model, but not from the average calculation
	s0 <- data.frame(s0 = 0, block = tmp[,'block'])
	nak <- which(tmp[,'map'] != 0 & is.na(tmp[,'lcc']))
	if(length(nak) > 0) s0[nak, 1] <- 1
	s0 <- na.omit(s0)
	s0 <- cbind(sort(unique(s0[,2])), rowsum(s0[,1], group = s0[,2]))
	s0 <- s0[which(s0[,2] > 0),]

	## I convert the dataframe to a list (for practical reasons)
	if(length(na) > 0){
		tmp <- list(map = tmp[,'map'][-na], lcc = tmp[,'lcc'][-na], cci = tmp[,'cci'][-na],
			x1 = tmp[,'x1'][-na], x2 = tmp[,'x2'][-na], block = tmp[,'block'][-na],
			ydg = tmp[,'ydg'][-na], yfi = tmp[,'yfi'][-na],
			tmax = as.matrix(tmax[-na,]), mcwd = as.matrix(mcwd[-na,]), L = as.matrix(L[-na,]))
	} else {
		tmp <- list(map = tmp[,'map'], lcc = tmp[,'lcc'], cci = tmp[,'cci'], 
			x1 = tmp[,'x1'], x2 = tmp[,'x2'], block = tmp[,'block'],
			ydg = tmp[,'ydg'], yfi = tmp[,'yfi'],
			tmax = as.matrix(tmax), mcwd = as.matrix(mcwd), L = as.matrix(L))
	}

	## We save the ordering of groups in the aggregation function, which allows us identify which row belongs to each coarse pixel
	groups <- sort(unique(tmp$block))

	## Then, we evaluate the smoothers at the pixel values extracted, and create the model matrix at the fine scale ('tmpgam')
	### Based on this model matrix, we derive the matrices at the coarse scale, that allow us to recover the coefficients while enforcing the disaggregation constraints

	tmpgam.X1 <- mgcv:::Predict.matrix3(sm1[[1]], tmp)$X
	tmpgam.X1 <- tpmm_mat_noblock(tmpgam.X1, qrQ1[,2:ncol(tmpgam.X1)]) # tmpgam.X1 <- tmpgam.X1 %*% qrQ1[,2:ncol(tmpgam.X1)]
	i_st <- 2
	i_en <- i_st + ncol(tmpgam.X1) - 1
	z. <- tpmv_vec_noblock(tmpgam.X1, b[[1]][i_st:i_en]) # z. <- tmpgam.X1 %*% b[[1]][i_st:i_en]

	tmpgam.X2 <- mgcv:::Predict.matrix3(sm2[[1]], tmp)$X
	tmpgam.X2 <- tpmm_mat_noblock(tmpgam.X2, qrQ2[,2:ncol(tmpgam.X2)]) # tmpgam.X2 <- tmpgam.X2 %*% qrQ2[,2:ncol(tmpgam.X2)]
	i_st <- i_en + 1
	i_en <- i_st + ncol(tmpgam.X2) - 1
	z. <- z. + tpmv_vec_noblock(tmpgam.X2, b[[1]][i_st:i_en]) # z. <- z. + tmpgam.X2 %*% b[[1]][i_st:i_en]

	tmpgam.X3 <- FastPredictMat3(sm3[[1]], tmp, qrQ3)
	i_st <- i_en + 1
	i_en <- i_st + ncol(tmpgam.X3) - 1
	z. <- z. + tpmv_vec_noblock(tmpgam.X3, b[[1]][i_st:i_en]) # z. <- z. + tmpgam.X3 %*% b[[1]][i_st:i_en]

	tmpgam.X4 <- mgcv:::Predict.matrix3(sm4[[1]], tmp)$X
	tmpgam.X4 <- tpmm_mat_noblock(tmpgam.X4, qrQ4[,2:ncol(tmpgam.X4)]) # tmpgam.X4 <- tmpgam.X4 %*% qrQ4[,2:ncol(tmpgam.X4)]
	i_st <- i_en + 1
	i_en <- i_st + ncol(tmpgam.X4) - 1
	z. <- z. + tpmv_vec_noblock(tmpgam.X4, b[[1]][i_st:i_en]) # z. <- z. + tmpgam.X4 %*% b[[1]][i_st:i_en]

	# sm5 has the NoData encoded as 999, so it's handled separately
	tmpgam.sm5 <- mgcv:::Predict.matrix3(sm5[[1]], tmp)$X
#	tmpgam.sm5 <- tpmm_mat_noblock(tmpgam.sm5, qrQ5[,2:ncol(tmpgam.sm5)]) # tmpgam.sm5 <- tmpgam.sm5 %*% qrQ5[,2:ncol(tmpgam.sm5)]
	fx.yfi <- which(tmp$yfi == 999)
	if(length(fx.yfi) > 0) tmpgam.sm5[fx.yfi,] <- 0

	tmpgam.X5 <- tmpgam.sm5
	na <- which(tmp$lcc != 6) # Reclassified: 1=Undisturbed, 2=Disturbed, 3=Deforestation, 4=OtherCover, 5=Croplands, 6=Shrublands, 7=Grasslands, 8=Deciduous forest
	if(length(na) > 0) tmpgam.X5[na,] <- 0
	i_st <- i_en + 1
	i_en <- i_st + ncol(tmpgam.X5) - 1
	z. <- z. + tpmv_vec_noblock(tmpgam.X5, b[[1]][i_st:i_en]) # z. <- z. + tmpgam.X5 %*% b[[1]][i_st:i_en]

	tmpgam.X6 <- tmpgam.sm5
	na <- which(tmp$lcc != 7) # Reclassified: 1=Undisturbed, 2=Disturbed, 3=Deforestation, 4=OtherCover, 5=Croplands, 6=Shrublands, 7=Grasslands, 8=Deciduous forest
	if(length(na) > 0) tmpgam.X6[na,] <- 0
	i_st <- i_en + 1
	i_en <- i_st + ncol(tmpgam.X6) - 1
	z. <- z. + tpmv_vec_noblock(tmpgam.X6, b[[1]][i_st:i_en]) # z. <- z. + tmpgam.X6 %*% b[[1]][i_st:i_en]

	tmpgam.X7 <- tmpgam.sm5
	na <- which(!(tmp$lcc %in% c(1, 8))) # Reclassified: 1=Undisturbed, 2=Disturbed, 3=Deforestation, 4=OtherCover, 5=Croplands, 6=Shrublands, 7=Grasslands, 8=Deciduous forest
	if(length(na) > 0) tmpgam.X7[na,] <- 0
	i_st <- i_en + 1
	i_en <- i_st + ncol(tmpgam.X7) - 1
	z. <- z. + tpmv_vec_noblock(tmpgam.X7, b[[1]][i_st:i_en]) # z. <- z. + tmpgam.X7 %*% b[[1]][i_st:i_en]

	tmpgam.X8 <- tmpgam.sm5
	# We assume fire is already taken into account in the 'disturbed' class
	na <- which(!(tmp$lcc %in% c(3, 4, 5))) # Reclassified: 1=Undisturbed, 2=Disturbed, 3=Deforestation, 4=OtherCover, 5=Croplands, 6=Shrublands, 7=Grasslands, 8=Deciduous forest
	if(length(na) > 0) tmpgam.X8[na,] <- 0
	i_st <- i_en + 1
	i_en <- i_st + ncol(tmpgam.X8) - 1
	z. <- z. + tpmv_vec_noblock(tmpgam.X8, b[[1]][i_st:i_en]) # z. <- z. + tmpgam.X8 %*% b[[1]][i_st:i_en]

	z. <- z. + b[[1]][1] # Adds the intercept

	## From here on, we derive the matrices at the coarse scale. It is necessary to read the disaggregation paper to understand what they mean.
	### Please refer to its Section 2.2
	gz. <- g(z.)
	v.out <- rowsum(gz., group = tmp$block)
	tmp_Agz. <- v.out

	#MMt
	dg2 <- 0 * gz. + 1
	v.out <- rowsum(dg2, group = tmp$block)
	tmp_MMt <- v.out

	#MX
	tmpgam.X1 <- as.numeric(dg(z.)) * tmpgam.X1
	tmpgam.X2 <- as.numeric(dg(z.)) * tmpgam.X2
	tmpgam.X3 <- as.numeric(dg(z.)) * tmpgam.X3
	tmpgam.X4 <- as.numeric(dg(z.)) * tmpgam.X4
	tmpgam.X5 <- as.numeric(dg(z.)) * tmpgam.X5
	tmpgam.X6 <- as.numeric(dg(z.)) * tmpgam.X6
	tmpgam.X7 <- as.numeric(dg(z.)) * tmpgam.X7
	tmpgam.X8 <- as.numeric(dg(z.)) * tmpgam.X8
	v.out <- cbind(rowsum(dg(z.), group = tmp$block),
		rowsum(tmpgam.X1, group = tmp$block),
		rowsum(tmpgam.X2, group = tmp$block),
		rowsum(tmpgam.X3, group = tmp$block),
		rowsum(tmpgam.X4, group = tmp$block),
		rowsum(tmpgam.X5, group = tmp$block),
		rowsum(tmpgam.X6, group = tmp$block),
		rowsum(tmpgam.X7, group = tmp$block),
		rowsum(tmpgam.X8, group = tmp$block))
	tmp_AGX <- v.out

	return(list(tmp_Agz.[,1], tmp_AGX, tmp_MMt[,1], s0, groups))
}
environment(ag.m) <- .GlobalEnv # Prevents slow serialization: https://stackoverflow.com/a/59117823

## Here, the cluster is created and the variables, exported
clusterEvalQ(cl, rm(list=ls()))
clusterExport(cl, list('bss', 'sm1', 'sm2', 'sm3', 'sm4', 'sm5', 'df.Y', 'ag.m', 'g', 'dg', 'cellFromRow', 'coordinates', 'lowres', 'id.min', 'id.max'))
clusterApply(cl, seq_along(cl), function(i) sink(paste0('output/worker_', i)))
clusterEvalQ(cl, {
	library(disag.vAfrica)
	library(numDeriv)
	library(terra)
	library(mgcv)
	library(Rcpp)
	source('utils/PredictMat_new.R')

	qrQ1 <- qr.Q(attr(sm1[[1]], 'qrc'), complete = TRUE)
	qrQ2 <- qr.Q(attr(sm2[[1]], 'qrc'), complete = TRUE)
	qrQ3 <- qr.Q(attr(sm3[[1]], 'qrc'), complete = TRUE)
	qrQ4 <- qr.Q(attr(sm4[[1]], 'qrc'), complete = TRUE)
#	qrQ5 <- qr.Q(attr(sm5[[1]], 'qrc'), complete = TRUE)

	return(TRUE)
})

## While we don't calculate the analytical gradient, we use the numerical
grfitCpp <- function(x, ...) grad(fitCpp, x, ..., method = 'simple')
clusterExport(cl, list('grfitCpp'))

## This is the function to coordinate the generation of model matrices in parallel
out.fn <- function(b){
	on.exit({
		rm(list=ls())
		gc()
	})

	## All target variables are started with zero values, and then updated after the result of each raster block is received
	b <- list(matrix(b, ncol = 1))
	Agz. <- matrix(0, nrow = dim(df.Y)[1], ncol = 2)
	Agz.[,1] <- df.Y[,1]
	X_ <- matrix(0, nrow = dim(df.Y)[1], ncol = bd)
	diagMMt <- matrix(0, nrow = dim(df.Y)[1], ncol = 2)
	diagMMt[,1] <- df.Y[,1]

	## Sends the function to the workers, and updates the target variables
	for(i in 1:nodes){
		sendCall(cl[[i]], ag.m, list(i = i, b = b), tag = i)
	}
	for (i in 1:d){
		k <- recvOneData_fair(cl)

		ni <- nodes + i
		if (ni <= d) {
			sendCall(cl[[k$node]], ag.m, list(i = ni, b = b), tag = ni)
		}

		k <- k$value$value

		if(i == 1 | i %% round(d/5, 0) == 0) print(paste0('received: ', i, ' - ag.m - ', Sys.time()))

		if(!is.null(k[[5]])){
			groups <- k[[5]]
			gid <- match(df.Y[,1], groups)
			n.gid <- which(!is.na(gid))
			if(length(n.gid) > 0){
				v.out <- k[[1]]
				Agz.[n.gid,2] <- Agz.[n.gid,2] + v.out[gid[n.gid]]

				v.out <- k[[2]]
				X_[n.gid,] <- X_[n.gid,] + v.out[gid[n.gid],]

				v.out <- k[[3]]
				diagMMt[n.gid,2] <- diagMMt[n.gid,2] + v.out[gid[n.gid]]
			}

			# Again, for the extra zeros
			v.out <- k[[4]]
			groups <- v.out[,1]
			gid <- match(df.Y[,1], groups)
			n.gid <- which(!is.na(gid))
			if(length(n.gid) > 0){
				diagMMt[n.gid,2] <- diagMMt[n.gid,2] + v.out[gid[n.gid],2]
			}

			rm(k)
		}
	}
	print(paste0('received: all - ag.m - ', Sys.time()))

	b <- b[[1]]

	## After all coarse-scale model matrices are received, we have to derive yet another auxiliary matrices
	### Again, it is necessary to refer to the paper's Section 2.2.3
	## Any linear operation between the fine-scale and the coarse-scale data is supported (sum, mean, weighted sum etc.)
	### In the AGC case, the following block of lines is necessary because the coarse-scale values are the average of the fine-scale ones
	n <- diagMMt[,2]
	diagMMt[,2] <- diagMMt[,2] * (1/n)^2
	Agz.[,2] <- Agz.[,2] * (1/n)
	X_ <- X_ * (1/n)
	if(lowres){
		if(!file.exists('pk.rds')){
			if(sum(n == 0) > 0){ #lowresolution
				pk <- which(n == 0); saveRDS(pk, 'pk.rds') #lowresolution
			} #lowresolution
		}
	}

	## Calculates the auxiliary matrices
	Y_ <- df.Y[,2] - Agz.[,2]
	ql <- diagMMt[,2] > 0
	if(exists('ll.t2')){
		diagL <- sqrt(1/ll.t2)
	} else {
		diagL <- sqrt(1/diagMMt[,2])
		ll.t2 <- diagMMt[,2]
	}
	y. <- diagL * Y_
	X. <- diagL * X_
	X.qr <- qr(X., LAPACK = TRUE)
	Qty. <- qr.qty(X.qr, y.)
	R <- qr.R(X.qr, complete = TRUE)[,sort.list(X.qr$pivot)]

	r <- crossprod(y.) - crossprod(Qty.)
	N <- sum(ql)
	er <- crossprod(Y_)
	print(er) ## This will print the quadratic error between the true coarse values and the predicted ones. Ideally this decreases at each iteration

	saveRDS(list(N, n), 'output/1Nn.rds')
	saveRDS(list(Y_, X_, n), 'output/1gen_mtc-mu.rds')
	return(list(Qty., R, r, N, b, sum(log(ll.t2))))
}

## Here, we set some constants before effectivelly running the model:
### 'maxit': the maximum number of iterations
#### In practice I'm not using 'maxit' anymore since the model is converging with relatively few iterations, but it can be useful in more complex cases
### 'b0': the initial guess of model coefficients. I modify the first element so that feasiblity is more likely
### 'st': the initial guess of optimization parameters
maxit <- 50
ml <- aic <- bic <- rep(NA, maxit)
n.dif <- 1
res = FALSE
mtcs <- rep(list(NA), maxit)
it <- 0
set.seed(1)
if(file.exists('output/1mtcs_sol.rds')){
	b0 <- readRDS('output/1mtcs_sol.rds')[[5]]
	st <- readRDS('output/1mtcs_sol.rds')[[7]]
	np <- length(st)
} else {
	set.seed(1)
	b0 <- rnorm(bd, 0, 0.01); b0[1] = 4.5
	np <- 11
	st <- rep(0, np)
}
source('2fit_ML2.R')
