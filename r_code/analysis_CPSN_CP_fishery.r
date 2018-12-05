# This is the R code to calculate the value of information for the fishery

# required libraries
library(R.matlab)
library(runjags)
library(coda)
library(snowfall)
library(parallel)
library(rjags)
library(vegan)
library(foreign)

#get functions source
source("functions.r")

# get results of evaluation - NW = no weights, W = weights (weights used in the manuscript)
MatFileExpW80 <- readMat("../problems/ConsMDP/results/ExpRc_FisheryWeights80Model3pr0.2pd0.2.mat")
MatFileExpNW80 <- readMat("../problems/ConsMDP/results/ExpRc_FisheryNoWeights80Model3pr0.2pd0.2.mat")

#calculate the % improvement in reward due to new information (i.e., the value of information)
BenefitW80 <- ((MatFileExpW80$ExpRc[,,,1]/MatFileExpW80$ExpRc[,,,2]) - 1) * 100
BenefitNW80 <- ((MatFileExpNW80$ExpRc[,,,1]/MatFileExpNW80$ExpRc[,,,2]) - 1) * 100
