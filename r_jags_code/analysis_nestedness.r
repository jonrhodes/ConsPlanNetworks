# libraries and functions
library(R.matlab)
library(runjags)
library(coda)
library(snowfall)
library(parallel)
library(rjags)
library(vegan)
library(foreign)

#set working directory and functions source
setwd("M:/Users/uqjrhode/work/networks_cons_plan/")
source("./code/functions.r")

#first calculate the nestedness of the fishery communities
Fishery95 <- read.csv("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/MDP_ConsPlanning/fishery_species95.csv",header=T,row.names=1)
Fishery80 <- read.csv("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/MDP_ConsPlanning/fishery_species80.csv",header=T,row.names=1)
Fishery50 <- read.csv("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/MDP_ConsPlanning/fishery_species50.csv",header=T,row.names=1)
Nested95 <- oecosimu(Fishery95,nesteddisc,"r0",nsimul=10000)
Nested80 <- oecosimu(Fishery80,nesteddisc,"r0",nsimul=10000)
Nested50 <- oecosimu(Fishery50,nesteddisc,"r0",nsimul=10000)

#set up matrix to hold nestedenss results
SimNest <- matrix(NA,nrow=10,ncol=2)
dimnames(SimNest) <- list(NULL,c("SimNest","Discr"))

# now calculate the nestedness of the simulated species
#loop through nestedness and replicates
for (i in 10:1)
{
	#get nestedness value
	SimNest[i,1] <- (i / 10)	
	
	Sum <- 0
	
	for (j in 1:20)
	{
		#get right .mat file
		File <- paste("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/ConsMDP/BR_library/BRSites6NumSp20Occ0.25Nested",1 - (i / 10),"Rep1.mat",sep="")
		MatFile <- readMat(File)
		#calculate discrepancy index
		Sum <- Sum + nesteddisc(MatFile$BR[[6]])$statistic
	}

	SimNest[i,2] <- Sum / 20
}

write.csv(SimNest,"M:/Users/uqjrhode/work/networks_cons_plan/figures/simnest.csv",row.names=F)
