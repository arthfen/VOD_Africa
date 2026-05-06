################### Script: 1smoothers.R ##################################
## This script creates the reference smoothers for the mean parameter, based on a low-resolution sample of the data.

## The following lines load the necessary libraries
library(terra)
library(mgcv)

#####################################################
######## Step 1: setting up the environment #########
#####################################################
## In this step, we just define some necessary functions and loads the coarse variable
size <- 1e5 # Size of the spatial sample
yrmin <- 2010 # The first simulation year
yrmax <- 2022 # The last simulation year
conc.size <- 1e4 # Size of the sample to compute concurvity
blks <- nr <- nrow(rast('../large_domain_allpixels/input/VOD/cellNumber-high.tif'))

## The bssize and cellFromRow functions (see 0minmax.R)
### 'coordinates' returns the xy coordinates of grid cells
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

#####################################################
##### Step 2: creating the reference smoothers ######
#####################################################
## In this step, we define the basis functions of the multi and univariate smoothers that will be used
### During estimation step (Step 4, 2fit_ML2.R), we will find the coefficients that multiply these basis functions in order to generate the disaggregation curves

## Loads the rasters
cci <- rast('../large_domain_allpixels/input/NCEO/AGB_map_2017v0-aligned.tif')
lcc <- rast(paste0('../large_domain_allpixels/input/LandCover/LandCover', yrmax, '.tif'))
yfi <- rast(paste0('../large_domain_allpixels/input/GABAM/YearsAfterFire_', yrmax, '.tif'))
tem <- lapply(yrmin:yrmax, function(i) rast(paste0('../large_domain_allpixels/input/Temperature/MaxTemp_', i, '.tif')))
mcw <- lapply(yrmin:yrmax, function(i) rast(paste0('../large_domain_allpixels/input/TerraClimate/Max-TerraClimate_def_', i, '.tif')))
vpd <- lapply(yrmin:yrmax, function(i) rast(paste0('../large_domain_allpixels/input/TerraClimate/Max-TerraClimate_vpd_', i, '.tif')))

snd <- rast('../large_domain_allpixels/input/isdaAfrica/SoilSand.tif')
soc <- rast('../large_domain_allpixels/input/isdaAfrica/SOC.tif')
bdr <- rast('../large_domain_allpixels/input/boundaries_high.tif')
rbd <- rast('../large_domain_allpixels/input/roughBoundary.tif')
blk <- rast('../large_domain_allpixels/input/VOD/cellNumber-high.tif')

## Creates a dataframe 'tmp' containing the cell values of a lower resolution version of the rasters
### The smoother construction changes with the variable range and quantiles, so to have an approximation of the distribution of the original raster, we use a lower resolution
### The original raster (at the fine scale) can not be used because it contains too many pixels
set.seed(1)
tmpr <- spatSample(bdr, size, method = 'regular', as.raster = TRUE)
tmp <- data.frame(coordinates(tmpr[[1]]), readValues(tmpr, mat = TRUE))
colnames(tmp) <- c('x1', 'x2', 'bdr')

cel <- cellFromXY(cci, tmp[,1:2])

ndf <- c(cci, lcc, snd, soc, yfi, blk, rbd)[cel]
colnames(ndf) <- c('cci', 'lcc', 'snd', 'soc', 'yfi', 'blk', 'rbd') # format: lcc_<pixel value>

## Simulates the deltas
tem <- do.call(c, tem)
mcw <- do.call(c, mcw)
vpd <- do.call(c, vpd)
tss <- c(tem, mcw, vpd)[cel]
nc <- ncol(tss)/3
for(j in (0:2)*nc){
	for(i in rev(seq_along(yrmin:yrmax))){
		if(i == 1){
			tss[,j+i] = NA
		} else {
			tss[,j+i] = tss[,j+i] - tss[,j+i-1]
		}
	}
}
v.tem <- c(min(tss[,1:13], na.rm = TRUE), sample(na.omit(tss[,2]), nrow(tss)-2, replace = TRUE), max(tss[,1:13], na.rm = TRUE))
v.mcw <- c(min(tss[,14:26], na.rm = TRUE), sample(na.omit(tss[,16]), nrow(tss)-2, replace = TRUE), max(tss[,14:26], na.rm = TRUE))
v.vpd <- c(min(tss[,27:39], na.rm = TRUE), sample(na.omit(tss[,28]), nrow(tss)-2, replace = TRUE), max(tss[,27:39], na.rm = TRUE))


# Small reclassification to disregard deserts
na.lcc <- which(ndf$lcc == 10)
if(length(na.lcc) > 0){
	ndf$lcc[na.lcc] <- NA
}

## Log-transforms the CCI
ndf$cci <- log(ndf$cci) # warning: zeroes become NA, but we want to remove them
ndf$cci[!is.finite(ndf$cci)] <- NA

## Simulates a MCWD time-series
mcwd <- data.frame(t0 = v.mcw)
mcwd$t2 <- mcwd$t1 <- mcwd$t0

## Simulates a temperature time-series
tmax = data.frame(t0 = v.tem)
tmax$t2 <- tmax$t1 <- tmax$t0

## Simulates a vapor pressure deficit time-series
vpd = data.frame(t0 = v.vpd)
vpd$t2 <- vpd$t1 <- vpd$t0

## Simulates a land cover time-series
lac = data.frame(t0 = ndf$lcc)
lac$t2 <- lac$t1 <- lac$t0

## Simulates the years after fire and after degradation
set.seed(1)
tab.yfi <- table(ndf$yfi)
tmp.yfi <- sample(as.numeric(names(tab.yfi)), nrow(ndf), prob = tab.yfi/sum(tab.yfi), replace = TRUE)

## Simulates the year
set.seed(1)
tmp.year <- sample(2010:2022, nrow(ndf), replace = TRUE)
L <- matrix(rep(tmp.year, 3), ncol = 3, nrow = nrow(mcwd))
for(i in 2:3) L[,i] <- L[,(i-1)] - 1

tmp.sm <- data.frame(x1 = tmp$x1,
	x2 = tmp$x2,
	bdr = tmp$bdr,
	rbd = ndf$rbd,
	cci = ndf$cci,
	snd = ndf$snd,
	soc = ndf$soc,
	blk = ndf$blk,
	yfi = tmp.yfi,
	year = tmp.year)
tmp.sm$tmax <- as.matrix(tmax)
tmp.sm$mcwd <- as.matrix(mcwd)
tmp.sm$vpd <- as.matrix(vpd)
tmp.sm$lcc <- as.matrix(lac)
tmp.sm$L = as.matrix(L)

## Removes NAs
na <- (tmp.sm$bdr == 0 | tmp.sm$rbd == 0 | tmp.sm$lcc[,1] == 0 | tmp.sm$lcc[,2] == 0 | tmp.sm$lcc[,3] == 0 | tmp.sm$x2 > 17)
tmp.sm <- tmp.sm[!na,]
tmp.sm <- na.omit(tmp.sm)
saveRDS(tmp.sm, 'output/1tmp_sm.rds')

## The index for concurvity calculation
set.seed(1)
conc.idx <- sample(1:nrow(tmp.sm), conc.size, replace = FALSE)

## Here, we build the reference uni and multivariate smoothers to be used. We ask them to absorb the centering constraint because we add an intercept term later.
### The values set in file 1call.R are being used to define the dimensions.
sm1 <- smoothCon(te(x1, x2, year, bs = c(h.sp, h.sp, h.sp), m = list(h.m, h.m, h.m), k = c(h.b1, h.b2, h.b3)),
	data = tmp.sm, absorb.cons = TRUE) # s(long, lat)
sm1[[1]]$X <- sm1[[1]]$X[conc.idx,]
sm1[[1]]$margin[[1]]$X <- sm1[[1]]$margin[[2]]$X <- NULL

sm2 <- smoothCon(te(tmax, L, bs = c(h.sp, h.sp), m = h.m, k = c(h.b4, h.b5)),
	data = tmp.sm, absorb.cons = TRUE) # s(temperature)
sm2[[1]]$X <- sm2[[1]]$X[conc.idx,]
sm2[[1]]$margin[[1]]$X <- sm2[[1]]$margin[[2]]$X <- NULL
sm2[[1]]$ind <- NULL

sm3 <- smoothCon(te(mcwd, L, bs = c(h.sp, h.sp), m = h.m, k = c(h.b6, h.b7)),
	data = tmp.sm, absorb.cons = TRUE) # s(water deficit)
sm3[[1]]$X <- sm3[[1]]$X[conc.idx,]
sm3[[1]]$margin[[1]]$X <- sm3[[1]]$margin[[2]]$X <- NULL

sm4 <- smoothCon(te(vpd, L, bs = c(h.sp, h.sp), m = h.m, k = c(h.b8, h.b9)),
	data = tmp.sm, absorb.cons = TRUE) # s(vapor pressure deficit)
sm4[[1]]$X <- sm4[[1]]$X[conc.idx,]
sm4[[1]]$margin[[1]]$X <- sm4[[1]]$margin[[2]]$X <- NULL
sm4[[1]]$ind <- NULL

sm5 <- smoothCon(te(lcc, L, bs = c(h.sp, h.sp), m = h.m, k = c(h.b10, h.b11)),
	data = tmp.sm, absorb.cons = TRUE) # s(land cover)
sm5[[1]]$X <- sm5[[1]]$X[conc.idx,]
sm5[[1]]$margin[[1]]$X <- sm5[[1]]$margin[[2]]$X <- NULL
sm5[[1]]$ind <- NULL

set.seed(1)
tmp_lcc <- tmp.sm$lcc[,'t0']
lcc_yfi <- NA * tmp_lcc
lcc_yfi[tmp_lcc == 6] <- 1
lcc_yfi[tmp_lcc == 7] <- 2
lcc_yfi[tmp_lcc %in% c(1, 8, 9)] <- 3
lcc_yfi[tmp_lcc %in% c(2, 3, 4, 5)] <- 4
tmp.sm$lcfi <- lcc_yfi
sm6 <- smoothCon(te(lcfi, yfi, bs = c(h.sp, h.sp), m = c(h.m, h.m), k = c(h.b12, h.b13)),
	data = tmp.sm, absorb.cons = TRUE) # s(years after fire, by group of land cover class)
sm6[[1]]$X <- sm6[[1]]$X[conc.idx,]
sm6[[1]]$margin[[1]]$X <- sm6[[1]]$margin[[2]]$X <- NULL

sm7 <- smoothCon(s(cci, bs = h.sp, m = h.m, k = h.b14), data = tmp.sm, absorb.cons = TRUE)
sm7[[1]]$X <- sm7[[1]]$X[conc.idx,]
# reference AGC

sm8 <- smoothCon(s(snd, bs = h.sp, m = h.m, k = h.b15), data = tmp.sm, absorb.cons = TRUE)
sm8[[1]]$X <- sm8[[1]]$X[conc.idx,]
# soil sand

sm9 <- smoothCon(s(soc, bs = h.sp, m = h.m, k = h.b16), data = tmp.sm, absorb.cons = TRUE)
sm9[[1]]$X <- sm9[[1]]$X[conc.idx,]
# soil organic carbon

saveRDS(list(sm1, sm2, sm3, sm4, sm5, sm6, sm7, sm8, sm9), 'output/1arable_sm1.rds')
