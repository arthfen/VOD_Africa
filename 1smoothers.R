################### Script: 1smoothers.R ##################################
## This script creates the reference smoothers for the mean parameter, based on a low-resolution sample of the data.

## The following lines load the necessary libraries
library(terra)
library(mgcv)
library(Rcpp)
library(numDeriv)

#####################################################
######## Step 1: setting up the environment #########
#####################################################
## In this step, we just define some necessary functions and loads the coarse variable

set.seed(1)
nr <- nrow(rast('../large_domain_allpixels/input/VOD/cellNumber-high.tif'))
blks <- nr

size <- 1e5
yrmax <- 2022

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
cci <- rast('../large_domain_allpixels/input_low/CCIv6/CCI2010_AGB-t_ha.tif') # TODO: replace input_low
tem <- rast(paste0('../large_domain_allpixels/input_low/Temperature/MaxTemp_', yrmax, '.tif')) # TODO: replace input_low
lcc <- rast(paste0('../large_domain_allpixels/input_low/LandCover/LandCover', yrmax, '.tif')) # TODO: replace input_low
mcw <- rast(paste0('../large_domain_allpixels/input_low/MCWD/MCWD_', yrmax, '.tif')) # TODO: replace input_low
yfi <- rast(paste0('../large_domain_allpixels/input_low/GABAM/YearsAfterFire_', yrmax, '.tif')) # TODO: replace input_low
ydg <- rast(paste0('../large_domain_allpixels/input_low/TMF/YearsAfterDegradation_', yrmax, '.tif')) # TODO: replace input_low
block <- rast('../large_domain_allpixels/input_low/boundaries_high.tif') # TODO: replace input_low

## Creates a dataframe 'tmp' containing the cell values of a lower resolution version of the rasters
### The smoother construction changes with the variable range and quantiles, so to have an approximation of the distribution of the original raster, we use a lower resolution
### The original raster (at the fine scale) can not be used because it contains too many pixels
tmpr <- spatSample(block, size, method = 'regular', as.raster = TRUE)
tmp <- data.frame(coordinates(tmpr[[1]]), readValues(tmpr, mat = TRUE))
colnames(tmp) <- c('x1', 'x2', 'block')

sp.ndf <- spatSample(c(cci, tem, lcc, mcw, yfi, ydg), size, method = 'regular', as.raster = TRUE)
ndf <- data.frame(readValues(sp.ndf, mat = TRUE))
colnames(ndf) <- c('cci', 'tmax', 'lcc', 'mcwd', 'yfi', 'ydg') # format: lcc_<pixel value>

# Years after fire and after degradation have rare values that make estimation hard, we filter them out here
na.yfi <- which(ndf$yfi >= 19) # Around 2% of pixels in the land cover types we use
if(length(na.yfi) > 0){
	ndf$yfi[na.yfi] <- 19
}
na.ydg <- which(ndf$ydg >= 24) # Around 2% of pixels
if(length(na.ydg) > 0){
	ndf$ydg[na.ydg] <- 24
}

# Small reclassification to disregard deserts
na.lcc <- which(ndf$lcc == 8)
if(length(na.lcc) > 0){
	ndf$lcc[na.lcc] <- NA
}
na.lcc <- which(ndf$lcc == 9)
if(length(na.lcc) > 0){
	ndf$lcc[na.lcc] <- 8 # Deciduous forest reclassified to 8
}

na <- (tmp$block == 0 | ndf$lcc == 0)

## Simulates a MCWD time-series
mcwd = data.frame(t0 = log(1 + abs(ndf$mcwd)))
mcwd$t3 <- mcwd$t2 <- mcwd$t1 <- mcwd$t0
L <- 1 + mcwd * 0

## Simulates a temperature time-series
tmax = data.frame(t0 = ndf$tmax)
tmax$t3 <- tmax$t2 <- tmax$t1 <- tmax$t0

## Simulates the years after fire
set.seed(1)
tab.yfi <- table(ndf$yfi)
tmp.yfi <- sample(as.numeric(names(tab.yfi)), nrow(ndf), prob = tab.yfi/sum(tab.yfi), replace = TRUE)

tmp.sm <- data.frame(x1 = tmp$x1,
	x2 = tmp$x2,
	block = tmp$block,
	lcc = ndf$lcc,
	cci = ndf$cci,
	yfi = tmp.yfi)
tmp.sm$tmax <- as.matrix(tmax)
tmp.sm$mcwd <- as.matrix(mcwd)
tmp.sm$L = as.matrix(L)
saveRDS(tmp.sm, 'output/1tmp_sm.rds')

## Removes NAs
tmp.sm <- tmp.sm[!na,]
tmp.sm <- na.omit(tmp.sm)

## Here, we build the reference uni and multivariate smoothers to be used. We ask them to absorb the centering constraint because we add an intercept term later.
### The values set in file 1call.R are being used to define the dimensions.
sm1 <- smoothCon(te(x1, x2, bs = c(h.sp, h.sp), m = h.m, k = c(h.b0, h.b1)), data = tmp.sm, absorb.cons = TRUE) # s(long, lat)
sm1[[1]]$X <- sm1[[1]]$margin[[1]]$X <- sm1[[1]]$margin[[1]]$X <- NULL

sm2 <- smoothCon(s(cci, bs = h.sp, m = h.m, k = h.b2), data = tmp.sm, absorb.cons = TRUE) # s(initial AGC)
sm2[[1]]$X <- NULL

sm3 <- smoothCon(te(tmax, mcwd, by = L, bs = c(h.sp, h.sp), m = h.m, k = c(h.b3, h.b4)), data = tmp.sm, absorb.cons = TRUE) # s(climate)
sm3[[1]]$X <- sm3[[1]]$margin[[1]]$X <- sm3[[1]]$margin[[2]]$X <- NULL

sm4 <- smoothCon(s(lcc, bs = h.sp, m = h.m, k = h.b5), data = tmp.sm, absorb.cons = TRUE) # s(land cover)
sm4[[1]]$X <- NULL

sm5 <- smoothCon(s(yfi, bs = h.sp, m = h.m, k = h.b6), data = tmp.sm, absorb.cons = FALSE) # s(years after fire), different for each land cover - therefore, no centering constraint
sm5[[1]]$X <- NULL
## Note: we apply sm5 to Grasslands, Shrublands and Forests. Together, they account for ~ 90% of all fires.

saveRDS(list(sm1, sm2, sm3, sm4, sm5), 'output/1arable_sm1.rds')
sm1[[1]]$S <- sm2[[1]]$S <- sm3[[1]]$S <- sm4[[1]]$S <- sm5[[1]]$S <- NULL
for(j in 1:length(sm1[[1]]$margin)) sm1[[1]]$margin[[j]]$S <- sm1[[1]]$margin[[j]]$X <- NULL
rm(tmp)
