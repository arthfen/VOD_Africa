################### Script: 1call.R ##################################
## This script starts the processing cluster and define the basis dimensions
## It follows by calling the next scripts, '1smoothers.R' and '2fit_ML.R'

## Creates the parallel processing cluster of 'nodes' nodes
library(parallel)
nodes <- 19
cl <- makeCluster(nodes, type = 'PSOCK')

## Creates the dataframe with the basis dimensions for each smoother
### In order to understand this, go to file 2fit_ML.R
df <- data.frame(b0 = 12, b1 = 12, b2 = 12, b3 = 12, b4 = 12, b5 = 8, b6 = 12)
df$res <- 0

h.b0 = df[1, 'b0']
h.b1 = df[1, 'b1']
h.b2 = df[1, 'b2']
h.b3 = df[1, 'b3']
h.b4 = df[1, 'b4']
h.b5 = df[1, 'b5']
h.b6 = df[1, 'b6']
h.sp <- 'cs' # We are using cubic splines (see the documentation of the mgcv package)
h.m <- 3

source('1smoothers.R')
source('2fit_ML.R')
