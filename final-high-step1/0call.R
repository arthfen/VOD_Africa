################### Script: 0call.R ##################################
## This script starts the processing cluster and define the basis dimensions
## It follows by calling the next scripts, '1smoothers.R' and '2fit_ML.R'

## Creates the dataframe with the basis dimensions for each smoother
### In order to understand this, go to file 2fit_ML.R
vf <- c(b1 = 5, b2 = 5, b3 = 7,
	b4 = 6, b5 = 7,
	b6 = 6, b7 = 7,
	b8 = 6, b9 = 7,
	b10 = 6, b11 = 7,
	b12 = 4, b13 = 10,
	b14 = 13, b15 = 13,
	b16 = 13, b17 = 13)
for(i in seq_along(vf)) assign(paste0('h.b', i), vf[paste0('b', i)])
h.sp <- 'cs'
h.m <- 2

#source('1smoothers.R')
source('2fit_ML.R')
