# libraries and functions
library(R.matlab)
library(runjags)
library(coda)
library(snowfall)
library(parallel)
library(rjags)

#set working directory and functions source
setwd("M:/Users/uqjrhode/work/networks_cons_plan/")
source("./code/functions.r")

#GET THE DATA FOR EACH NETWORK SHAPE 

#STAR

#set up matrices to hold results
ResultsExp <- matrix(NA,ncol=9,nrow=1*1*10*5*10*10)
dimnames(ResultsExp) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN"))
ResultsSD <- matrix(NA,ncol=9,nrow=1*1*10*5*10*10)
dimnames(ResultsSD) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN"))

#create lists of variables
NumSpec <- c(20)
Occ <- c(0.25)
Nested <- seq(0,0.9,0.1)
PiR <- seq(0.1,1,0.1)
PiD <- seq(0.1,1,0.1)

RepCount <- 0

#loop through number of species
for (i in 1:length(NumSpec))
{
	ResultsExp[(((i - 1) * (1*10*5*10*10)) + 1):(i * (1*10*5*10*10)),"NumSp"] <- NumSpec[i]
	ResultsSD[(((i - 1) * (1*10*5*10*10)) + 1):(i * (1*10*5*10*10)),"NumSp"] <- NumSpec[i]
		
	#loop through occupancy
	for (j in 1:length(Occ))
	{
		ResultsExp[((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + 1):((((i - 1) * (1*10*5*10*10))) + (j * (10*5*10*10))),"Occ"] <- Occ[j]
		ResultsSD[((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + 1):((((i - 1) * (1*10*5*10*10))) + (j * (10*5*10*10))),"Occ"] <- Occ[j]
						
		#loop through nested
		for (k in 1:length(Nested))
		{
			StartVal <- ((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + (((k - 1) * (5*10*10))) + 1)
									
			#get the right .mat file
			# dimensions = [pir,pid,rep_id(BR replicate),method(1=CPSN,2=CP,3=SN,4=RND),number simulations] - no number simulations in the expectation and SD files 
			MatFileExp <- readMat(paste("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/ConsMDP/results/ExpRc_starModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
			MatFileSD <- readMat(paste("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/ConsMDP/results/SDRc_starModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
									
			for (a in 1:dim(MatFileExp$ExpRc)[3]) #BR rep
			{
				RepCount <- RepCount + 1
										
				for (b in 1:dim(MatFileExp$ExpRc)[1]) #pir
				{
					for (c in 1:dim(MatFileExp$ExpRc)[2]) #pid
					{
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"CPSN"] <- MatFileExp$ExpRc[b,c,a,1]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"CP"] <- MatFileExp$ExpRc[b,c,a,2]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"SN"] <- MatFileExp$ExpRc[b,c,a,3]
																		
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"CPSN"] <- MatFileSD$SDRc[b,c,a,1]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"CP"] <- MatFileSD$SDRc[b,c,a,2]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"SN"] <- MatFileSD$SDRc[b,c,a,3]						
					}
				}
			}
		}
	}	
}

#get data
#calculate log(E(CPSN)) - log(E(CP)) assuming the E(SN) and E(CP) are log-normally distributed
# Here, for the logged variables mu = log(m^2 / sqrt(v + m^2)) and sigma = sqrt(log((v/m^2) + 1))

SC <- log((ResultsExp[,"CPSN"]^2) / sqrt((ResultsSD[,"CPSN"]^2) + (ResultsExp[,"CPSN"]^2))) - log((ResultsExp[,"CP"]^2) / sqrt((ResultsSD[,"CP"]^2) + (ResultsExp[,"CP"]^2)))
TAUSC <- 1 / (log(((ResultsSD[,"CPSN"]^2) / (ResultsExp[,"CPSN"]^2)) + 1) + log(((ResultsSD[,"CP"]^2) / (ResultsExp[,"CP"]^2)) + 1))
N <- length(SC)
NUMBR <- length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3]

#get BR Rep level covariates 
NST <- as.vector(matrix(NA,nrow=length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpRc)[3],ncol=1))
for (i in 1:(nrow(ResultsExp) - 1))
{
	if ((ResultsExp[i,"BRRep"] != ResultsExp[i + 1,"BRRep"]) | (i == (nrow(ResultsExp) - 1)))
	{
		NST[ResultsExp[i,"BRRep"]] <- 1 - ResultsExp[i,"Nested"]
	}
}
#standardise
NST <- (NST - mean(NST)) / sd(NST)

#get within BRRep level covariates
BR <- ResultsExp[,"BRRep"]
PIR <- ResultsExp[,"PiR"]
PID <- ResultsExp[,"PiD"]
PIRD <- ResultsExp[,"PiR"] * ResultsExp[,"PiD"]
#standardise
PIR <- (PIR - mean(PIR)) / sd(PIR)
PID <- (PID - mean(PID)) / sd(PID)
PIRD <- (PIRD - mean(PIRD)) / sd(PIRD)

#complie data
DataStar <- list(NUMBR=NUMBR,N=N,SC=SC,TAUSC=TAUSC,NST=NST,PIR=PIR,PID=PID,BR=BR)

#LINE

#set up matrices to hold results
ResultsExp <- matrix(NA,ncol=9,nrow=1*1*10*5*10*10)
dimnames(ResultsExp) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN"))
ResultsSD <- matrix(NA,ncol=9,nrow=1*1*10*5*10*10)
dimnames(ResultsSD) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN"))

#create lists of variables
NumSpec <- c(20)
Occ <- c(0.25)
Nested <- seq(0,0.9,0.1)
PiR <- seq(0.1,1,0.1)
PiD <- seq(0.1,1,0.1)

RepCount <- 0

#loop through number of species
for (i in 1:length(NumSpec))
{
	ResultsExp[(((i - 1) * (1*10*5*10*10)) + 1):(i * (1*10*5*10*10)),"NumSp"] <- NumSpec[i]
	ResultsSD[(((i - 1) * (1*10*5*10*10)) + 1):(i * (1*10*5*10*10)),"NumSp"] <- NumSpec[i]
		
	#loop through occupancy
	for (j in 1:length(Occ))
	{
		ResultsExp[((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + 1):((((i - 1) * (1*10*5*10*10))) + (j * (10*5*10*10))),"Occ"] <- Occ[j]
		ResultsSD[((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + 1):((((i - 1) * (1*10*5*10*10))) + (j * (10*5*10*10))),"Occ"] <- Occ[j]
						
		#loop through nested
		for (k in 1:length(Nested))
		{
			StartVal <- ((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + (((k - 1) * (5*10*10))) + 1)
									
			#get the right .mat file
			# dimensions = [pir,pid,rep_id(BR replicate),method(1=CPSN,2=CP,3=SN,4=RND),number simulations] - no number simulations in the expectation and SD files 
			MatFileExp <- readMat(paste("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/ConsMDP/results/ExpRc_lineModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
			MatFileSD <- readMat(paste("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/ConsMDP/results/SDRc_lineModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
									
			for (a in 1:dim(MatFileExp$ExpRc)[3]) #BR rep
			{
				RepCount <- RepCount + 1
										
				for (b in 1:dim(MatFileExp$ExpRc)[1]) #pir
				{
					for (c in 1:dim(MatFileExp$ExpRc)[2]) #pid
					{
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"CPSN"] <- MatFileExp$ExpRc[b,c,a,1]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"CP"] <- MatFileExp$ExpRc[b,c,a,2]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"SN"] <- MatFileExp$ExpRc[b,c,a,3]
																		
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"CPSN"] <- MatFileSD$SDRc[b,c,a,1]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"CP"] <- MatFileSD$SDRc[b,c,a,2]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"SN"] <- MatFileSD$SDRc[b,c,a,3]						
					}
				}
			}
		}
	}	
}

#get data
#calculate log(E(CPSN)) - log(E(CP)) assuming the E(SN) and E(CP) are log-normally distributed
# Here, for the logged variables mu = log(m^2 / sqrt(v + m^2)) and sigma = sqrt(log((v/m^2) + 1))

SC <- log((ResultsExp[,"CPSN"]^2) / sqrt((ResultsSD[,"CPSN"]^2) + (ResultsExp[,"CPSN"]^2))) - log((ResultsExp[,"CP"]^2) / sqrt((ResultsSD[,"CP"]^2) + (ResultsExp[,"CP"]^2)))
TAUSC <- 1 / (log(((ResultsSD[,"CPSN"]^2) / (ResultsExp[,"CPSN"]^2)) + 1) + log(((ResultsSD[,"CP"]^2) / (ResultsExp[,"CP"]^2)) + 1))
N <- length(SC)
NUMBR <- length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3]

#get BR Rep level covariates 
NST <- as.vector(matrix(NA,nrow=length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpRc)[3],ncol=1))
for (i in 1:(nrow(ResultsExp) - 1))
{
	if ((ResultsExp[i,"BRRep"] != ResultsExp[i + 1,"BRRep"]) | (i == (nrow(ResultsExp) - 1)))
	{
		NST[ResultsExp[i,"BRRep"]] <- 1 - ResultsExp[i,"Nested"]
	}
}
#standardise
NST <- (NST - mean(NST)) / sd(NST)

#get within BRRep level covariates
BR <- ResultsExp[,"BRRep"]
PIR <- ResultsExp[,"PiR"]
PID <- ResultsExp[,"PiD"]
PIRD <- ResultsExp[,"PiR"] * ResultsExp[,"PiD"]
#standardise
PIR <- (PIR - mean(PIR)) / sd(PIR)
PID <- (PID - mean(PID)) / sd(PID)
PIRD <- (PIRD - mean(PIRD)) / sd(PIRD)

#complie data
DataLine <- list(NUMBR=NUMBR,N=N,SC=SC,TAUSC=TAUSC,NST=NST,PIR=PIR,PID=PID,BR=BR)

#WHEEL

#set up matrices to hold results
ResultsExp <- matrix(NA,ncol=9,nrow=1*1*10*5*10*10)
dimnames(ResultsExp) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN"))
ResultsSD <- matrix(NA,ncol=9,nrow=1*1*10*5*10*10)
dimnames(ResultsSD) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN"))

#create lists of variables
NumSpec <- c(20)
Occ <- c(0.25)
Nested <- seq(0,0.9,0.1)
PiR <- seq(0.1,1,0.1)
PiD <- seq(0.1,1,0.1)

RepCount <- 0

#loop through number of species
for (i in 1:length(NumSpec))
{
	ResultsExp[(((i - 1) * (1*10*5*10*10)) + 1):(i * (1*10*5*10*10)),"NumSp"] <- NumSpec[i]
	ResultsSD[(((i - 1) * (1*10*5*10*10)) + 1):(i * (1*10*5*10*10)),"NumSp"] <- NumSpec[i]
		
	#loop through occupancy
	for (j in 1:length(Occ))
	{
		ResultsExp[((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + 1):((((i - 1) * (1*10*5*10*10))) + (j * (10*5*10*10))),"Occ"] <- Occ[j]
		ResultsSD[((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + 1):((((i - 1) * (1*10*5*10*10))) + (j * (10*5*10*10))),"Occ"] <- Occ[j]
						
		#loop through nested
		for (k in 1:length(Nested))
		{
			StartVal <- ((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + (((k - 1) * (5*10*10))) + 1)
									
			#get the right .mat file
			# dimensions = [pir,pid,rep_id(BR replicate),method(1=CPSN,2=CP,3=SN,4=RND),number simulations] - no number simulations in the expectation and SD files 
			MatFileExp <- readMat(paste("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/ConsMDP/results/ExpRc_wheelModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
			MatFileSD <- readMat(paste("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/ConsMDP/results/SDRc_wheelModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
									
			for (a in 1:dim(MatFileExp$ExpRc)[3]) #BR rep
			{
				RepCount <- RepCount + 1
										
				for (b in 1:dim(MatFileExp$ExpRc)[1]) #pir
				{
					for (c in 1:dim(MatFileExp$ExpRc)[2]) #pid
					{
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"CPSN"] <- MatFileExp$ExpRc[b,c,a,1]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"CP"] <- MatFileExp$ExpRc[b,c,a,2]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"SN"] <- MatFileExp$ExpRc[b,c,a,3]
																		
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"CPSN"] <- MatFileSD$SDRc[b,c,a,1]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"CP"] <- MatFileSD$SDRc[b,c,a,2]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"SN"] <- MatFileSD$SDRc[b,c,a,3]						
					}
				}
			}
		}
	}	
}

#get data
#calculate log(E(CPSN)) - log(E(CP)) assuming the E(SN) and E(CP) are log-normally distributed
# Here, for the logged variables mu = log(m^2 / sqrt(v + m^2)) and sigma = sqrt(log((v/m^2) + 1))

SC <- log((ResultsExp[,"CPSN"]^2) / sqrt((ResultsSD[,"CPSN"]^2) + (ResultsExp[,"CPSN"]^2))) - log((ResultsExp[,"CP"]^2) / sqrt((ResultsSD[,"CP"]^2) + (ResultsExp[,"CP"]^2)))
TAUSC <- 1 / (log(((ResultsSD[,"CPSN"]^2) / (ResultsExp[,"CPSN"]^2)) + 1) + log(((ResultsSD[,"CP"]^2) / (ResultsExp[,"CP"]^2)) + 1))
N <- length(SC)
NUMBR <- length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3]

#get BR Rep level covariates 
NST <- as.vector(matrix(NA,nrow=length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpRc)[3],ncol=1))
for (i in 1:(nrow(ResultsExp) - 1))
{
	if ((ResultsExp[i,"BRRep"] != ResultsExp[i + 1,"BRRep"]) | (i == (nrow(ResultsExp) - 1)))
	{
		NST[ResultsExp[i,"BRRep"]] <- 1 - ResultsExp[i,"Nested"]
	}
}
#standardise
NST <- (NST - mean(NST)) / sd(NST)

#get within BRRep level covariates
BR <- ResultsExp[,"BRRep"]
PIR <- ResultsExp[,"PiR"]
PID <- ResultsExp[,"PiD"]
PIRD <- ResultsExp[,"PiR"] * ResultsExp[,"PiD"]
#standardise
PIR <- (PIR - mean(PIR)) / sd(PIR)
PID <- (PID - mean(PID)) / sd(PID)
PIRD <- (PIRD - mean(PIRD)) / sd(PIRD)

#complie data
DataWheel <- list(NUMBR=NUMBR,N=N,SC=SC,TAUSC=TAUSC,NST=NST,PIR=PIR,PID=PID,BR=BR)

#RING

#set up matrices to hold results
ResultsExp <- matrix(NA,ncol=9,nrow=1*1*10*5*10*10)
dimnames(ResultsExp) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN"))
ResultsSD <- matrix(NA,ncol=9,nrow=1*1*10*5*10*10)
dimnames(ResultsSD) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN"))

#create lists of variables
NumSpec <- c(20)
Occ <- c(0.25)
Nested <- seq(0,0.9,0.1)
PiR <- seq(0.1,1,0.1)
PiD <- seq(0.1,1,0.1)

RepCount <- 0

#loop through number of species
for (i in 1:length(NumSpec))
{
	ResultsExp[(((i - 1) * (1*10*5*10*10)) + 1):(i * (1*10*5*10*10)),"NumSp"] <- NumSpec[i]
	ResultsSD[(((i - 1) * (1*10*5*10*10)) + 1):(i * (1*10*5*10*10)),"NumSp"] <- NumSpec[i]
		
	#loop through occupancy
	for (j in 1:length(Occ))
	{
		ResultsExp[((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + 1):((((i - 1) * (1*10*5*10*10))) + (j * (10*5*10*10))),"Occ"] <- Occ[j]
		ResultsSD[((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + 1):((((i - 1) * (1*10*5*10*10))) + (j * (10*5*10*10))),"Occ"] <- Occ[j]
						
		#loop through nested
		for (k in 1:length(Nested))
		{
			StartVal <- ((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + (((k - 1) * (5*10*10))) + 1)
									
			#get the right .mat file
			# dimensions = [pir,pid,rep_id(BR replicate),method(1=CPSN,2=CP,3=SN,4=RND),number simulations] - no number simulations in the expectation and SD files 
			MatFileExp <- readMat(paste("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/ConsMDP/results/ExpRc_islandModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
			MatFileSD <- readMat(paste("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/ConsMDP/results/SDRc_islandModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
									
			for (a in 1:dim(MatFileExp$ExpRc)[3]) #BR rep
			{
				RepCount <- RepCount + 1
										
				for (b in 1:dim(MatFileExp$ExpRc)[1]) #pir
				{
					for (c in 1:dim(MatFileExp$ExpRc)[2]) #pid
					{
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"CPSN"] <- MatFileExp$ExpRc[b,c,a,1]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"CP"] <- MatFileExp$ExpRc[b,c,a,2]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"SN"] <- MatFileExp$ExpRc[b,c,a,3]
																		
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"CPSN"] <- MatFileSD$SDRc[b,c,a,1]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"CP"] <- MatFileSD$SDRc[b,c,a,2]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"SN"] <- MatFileSD$SDRc[b,c,a,3]						
					}
				}
			}
		}
	}	
}

#get data
#calculate log(E(CPSN)) - log(E(CP)) assuming the E(SN) and E(CP) are log-normally distributed
# Here, for the logged variables mu = log(m^2 / sqrt(v + m^2)) and sigma = sqrt(log((v/m^2) + 1))

SC <- log((ResultsExp[,"CPSN"]^2) / sqrt((ResultsSD[,"CPSN"]^2) + (ResultsExp[,"CPSN"]^2))) - log((ResultsExp[,"CP"]^2) / sqrt((ResultsSD[,"CP"]^2) + (ResultsExp[,"CP"]^2)))
TAUSC <- 1 / (log(((ResultsSD[,"CPSN"]^2) / (ResultsExp[,"CPSN"]^2)) + 1) + log(((ResultsSD[,"CP"]^2) / (ResultsExp[,"CP"]^2)) + 1))
N <- length(SC)
NUMBR <- length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3]

#get BR Rep level covariates 
NST <- as.vector(matrix(NA,nrow=length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpRc)[3],ncol=1))
for (i in 1:(nrow(ResultsExp) - 1))
{
	if ((ResultsExp[i,"BRRep"] != ResultsExp[i + 1,"BRRep"]) | (i == (nrow(ResultsExp) - 1)))
	{
		NST[ResultsExp[i,"BRRep"]] <- 1 - ResultsExp[i,"Nested"]
	}
}
#standardise
NST <- (NST - mean(NST)) / sd(NST)

#get within BRRep level covariates
BR <- ResultsExp[,"BRRep"]
PIR <- ResultsExp[,"PiR"]
PID <- ResultsExp[,"PiD"]
PIRD <- ResultsExp[,"PiR"] * ResultsExp[,"PiD"]
#standardise
PIR <- (PIR - mean(PIR)) / sd(PIR)
PID <- (PID - mean(PID)) / sd(PID)
PIRD <- (PIRD - mean(PIRD)) / sd(PIRD)

#complie data
DataRing <- list(NUMBR=NUMBR,N=N,SC=SC,TAUSC=TAUSC,NST=NST,PIR=PIR,PID=PID,BR=BR)

#RING-STAR

#set up matrices to hold results
ResultsExp <- matrix(NA,ncol=9,nrow=1*1*10*5*10*10)
dimnames(ResultsExp) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN"))
ResultsSD <- matrix(NA,ncol=9,nrow=1*1*10*5*10*10)
dimnames(ResultsSD) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN"))

#create lists of variables
NumSpec <- c(20)
Occ <- c(0.25)
Nested <- seq(0,0.9,0.1)
PiR <- seq(0.1,1,0.1)
PiD <- seq(0.1,1,0.1)

RepCount <- 0

#loop through number of species
for (i in 1:length(NumSpec))
{
	ResultsExp[(((i - 1) * (1*10*5*10*10)) + 1):(i * (1*10*5*10*10)),"NumSp"] <- NumSpec[i]
	ResultsSD[(((i - 1) * (1*10*5*10*10)) + 1):(i * (1*10*5*10*10)),"NumSp"] <- NumSpec[i]
		
	#loop through occupancy
	for (j in 1:length(Occ))
	{
		ResultsExp[((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + 1):((((i - 1) * (1*10*5*10*10))) + (j * (10*5*10*10))),"Occ"] <- Occ[j]
		ResultsSD[((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + 1):((((i - 1) * (1*10*5*10*10))) + (j * (10*5*10*10))),"Occ"] <- Occ[j]
						
		#loop through nested
		for (k in 1:length(Nested))
		{
			StartVal <- ((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + (((k - 1) * (5*10*10))) + 1)
									
			#get the right .mat file
			# dimensions = [pir,pid,rep_id(BR replicate),method(1=CPSN,2=CP,3=SN,4=RND),number simulations] - no number simulations in the expectation and SD files 
			MatFileExp <- readMat(paste("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/ConsMDP/results/ExpRc_islandstarModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
			MatFileSD <- readMat(paste("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/ConsMDP/results/SDRc_islandstarModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
									
			for (a in 1:dim(MatFileExp$ExpRc)[3]) #BR rep
			{
				RepCount <- RepCount + 1
										
				for (b in 1:dim(MatFileExp$ExpRc)[1]) #pir
				{
					for (c in 1:dim(MatFileExp$ExpRc)[2]) #pid
					{
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"CPSN"] <- MatFileExp$ExpRc[b,c,a,1]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"CP"] <- MatFileExp$ExpRc[b,c,a,2]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"SN"] <- MatFileExp$ExpRc[b,c,a,3]
																		
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"CPSN"] <- MatFileSD$SDRc[b,c,a,1]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"CP"] <- MatFileSD$SDRc[b,c,a,2]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"SN"] <- MatFileSD$SDRc[b,c,a,3]						
					}
				}
			}
		}
	}	
}

#get data
#calculate log(E(CPSN)) - log(E(CP)) assuming the E(SN) and E(CP) are log-normally distributed
# Here, for the logged variables mu = log(m^2 / sqrt(v + m^2)) and sigma = sqrt(log((v/m^2) + 1))

SC <- log((ResultsExp[,"CPSN"]^2) / sqrt((ResultsSD[,"CPSN"]^2) + (ResultsExp[,"CPSN"]^2))) - log((ResultsExp[,"CP"]^2) / sqrt((ResultsSD[,"CP"]^2) + (ResultsExp[,"CP"]^2)))
TAUSC <- 1 / (log(((ResultsSD[,"CPSN"]^2) / (ResultsExp[,"CPSN"]^2)) + 1) + log(((ResultsSD[,"CP"]^2) / (ResultsExp[,"CP"]^2)) + 1))
N <- length(SC)
NUMBR <- length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3]

#get BR Rep level covariates 
NST <- as.vector(matrix(NA,nrow=length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpRc)[3],ncol=1))
for (i in 1:(nrow(ResultsExp) - 1))
{
	if ((ResultsExp[i,"BRRep"] != ResultsExp[i + 1,"BRRep"]) | (i == (nrow(ResultsExp) - 1)))
	{
		NST[ResultsExp[i,"BRRep"]] <- 1 - ResultsExp[i,"Nested"]
	}
}
#standardise
NST <- (NST - mean(NST)) / sd(NST)

#get within BRRep level covariates
BR <- ResultsExp[,"BRRep"]
PIR <- ResultsExp[,"PiR"]
PID <- ResultsExp[,"PiD"]
PIRD <- ResultsExp[,"PiR"] * ResultsExp[,"PiD"]
#standardise
PIR <- (PIR - mean(PIR)) / sd(PIR)
PID <- (PID - mean(PID)) / sd(PID)
PIRD <- (PIRD - mean(PIRD)) / sd(PIRD)

#complie data
DataRingStar <- list(NUMBR=NUMBR,N=N,SC=SC,TAUSC=TAUSC,NST=NST,PIR=PIR,PID=PID,BR=BR)

#RING-LINE

#set up matrices to hold results
ResultsExp <- matrix(NA,ncol=9,nrow=1*1*10*5*10*10)
dimnames(ResultsExp) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN"))
ResultsSD <- matrix(NA,ncol=9,nrow=1*1*10*5*10*10)
dimnames(ResultsSD) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN"))

#create lists of variables
NumSpec <- c(20)
Occ <- c(0.25)
Nested <- seq(0,0.9,0.1)
PiR <- seq(0.1,1,0.1)
PiD <- seq(0.1,1,0.1)

RepCount <- 0

#loop through number of species
for (i in 1:length(NumSpec))
{
	ResultsExp[(((i - 1) * (1*10*5*10*10)) + 1):(i * (1*10*5*10*10)),"NumSp"] <- NumSpec[i]
	ResultsSD[(((i - 1) * (1*10*5*10*10)) + 1):(i * (1*10*5*10*10)),"NumSp"] <- NumSpec[i]
		
	#loop through occupancy
	for (j in 1:length(Occ))
	{
		ResultsExp[((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + 1):((((i - 1) * (1*10*5*10*10))) + (j * (10*5*10*10))),"Occ"] <- Occ[j]
		ResultsSD[((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + 1):((((i - 1) * (1*10*5*10*10))) + (j * (10*5*10*10))),"Occ"] <- Occ[j]
						
		#loop through nested
		for (k in 1:length(Nested))
		{
			StartVal <- ((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + (((k - 1) * (5*10*10))) + 1)
									
			#get the right .mat file
			# dimensions = [pir,pid,rep_id(BR replicate),method(1=CPSN,2=CP,3=SN,4=RND),number simulations] - no number simulations in the expectation and SD files 
			MatFileExp <- readMat(paste("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/ConsMDP/results/ExpRc_islandlineModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
			MatFileSD <- readMat(paste("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/ConsMDP/results/SDRc_islandlineModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
									
			for (a in 1:dim(MatFileExp$ExpRc)[3]) #BR rep
			{
				RepCount <- RepCount + 1
										
				for (b in 1:dim(MatFileExp$ExpRc)[1]) #pir
				{
					for (c in 1:dim(MatFileExp$ExpRc)[2]) #pid
					{
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"CPSN"] <- MatFileExp$ExpRc[b,c,a,1]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"CP"] <- MatFileExp$ExpRc[b,c,a,2]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"SN"] <- MatFileExp$ExpRc[b,c,a,3]
																		
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"CPSN"] <- MatFileSD$SDRc[b,c,a,1]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"CP"] <- MatFileSD$SDRc[b,c,a,2]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"SN"] <- MatFileSD$SDRc[b,c,a,3]						
					}
				}
			}
		}
	}	
}

#get data
#calculate log(E(CPSN)) - log(E(CP)) assuming the E(SN) and E(CP) are log-normally distributed
# Here, for the logged variables mu = log(m^2 / sqrt(v + m^2)) and sigma = sqrt(log((v/m^2) + 1))

SC <- log((ResultsExp[,"CPSN"]^2) / sqrt((ResultsSD[,"CPSN"]^2) + (ResultsExp[,"CPSN"]^2))) - log((ResultsExp[,"CP"]^2) / sqrt((ResultsSD[,"CP"]^2) + (ResultsExp[,"CP"]^2)))
TAUSC <- 1 / (log(((ResultsSD[,"CPSN"]^2) / (ResultsExp[,"CPSN"]^2)) + 1) + log(((ResultsSD[,"CP"]^2) / (ResultsExp[,"CP"]^2)) + 1))
N <- length(SC)
NUMBR <- length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3]

#get BR Rep level covariates 
NST <- as.vector(matrix(NA,nrow=length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpRc)[3],ncol=1))
for (i in 1:(nrow(ResultsExp) - 1))
{
	if ((ResultsExp[i,"BRRep"] != ResultsExp[i + 1,"BRRep"]) | (i == (nrow(ResultsExp) - 1)))
	{
		NST[ResultsExp[i,"BRRep"]] <- 1 - ResultsExp[i,"Nested"]
	}
}
#standardise
NST <- (NST - mean(NST)) / sd(NST)

#get within BRRep level covariates
BR <- ResultsExp[,"BRRep"]
PIR <- ResultsExp[,"PiR"]
PID <- ResultsExp[,"PiD"]
PIRD <- ResultsExp[,"PiR"] * ResultsExp[,"PiD"]
#standardise
PIR <- (PIR - mean(PIR)) / sd(PIR)
PID <- (PID - mean(PID)) / sd(PID)
PIRD <- (PIRD - mean(PIRD)) / sd(PIRD)

#complie data
DataRingLine <- list(NUMBR=NUMBR,N=N,SC=SC,TAUSC=TAUSC,NST=NST,PIR=PIR,PID=PID,BR=BR)

#STAR-LINE

#set up matrices to hold results
ResultsExp <- matrix(NA,ncol=9,nrow=1*1*10*5*10*10)
dimnames(ResultsExp) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN"))
ResultsSD <- matrix(NA,ncol=9,nrow=1*1*10*5*10*10)
dimnames(ResultsSD) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN"))

#create lists of variables
NumSpec <- c(20)
Occ <- c(0.25)
Nested <- seq(0,0.9,0.1)
PiR <- seq(0.1,1,0.1)
PiD <- seq(0.1,1,0.1)

RepCount <- 0

#loop through number of species
for (i in 1:length(NumSpec))
{
	ResultsExp[(((i - 1) * (1*10*5*10*10)) + 1):(i * (1*10*5*10*10)),"NumSp"] <- NumSpec[i]
	ResultsSD[(((i - 1) * (1*10*5*10*10)) + 1):(i * (1*10*5*10*10)),"NumSp"] <- NumSpec[i]
		
	#loop through occupancy
	for (j in 1:length(Occ))
	{
		ResultsExp[((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + 1):((((i - 1) * (1*10*5*10*10))) + (j * (10*5*10*10))),"Occ"] <- Occ[j]
		ResultsSD[((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + 1):((((i - 1) * (1*10*5*10*10))) + (j * (10*5*10*10))),"Occ"] <- Occ[j]
						
		#loop through nested
		for (k in 1:length(Nested))
		{
			StartVal <- ((((i - 1) * (1*10*5*10*10))) + (((j - 1) * (10*5*10*10))) + (((k - 1) * (5*10*10))) + 1)
									
			#get the right .mat file
			# dimensions = [pir,pid,rep_id(BR replicate),method(1=CPSN,2=CP,3=SN,4=RND),number simulations] - no number simulations in the expectation and SD files 
			MatFileExp <- readMat(paste("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/ConsMDP/results/ExpRc_starlineModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
			MatFileSD <- readMat(paste("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/ConsMDP/results/SDRc_starlineModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
									
			for (a in 1:dim(MatFileExp$ExpRc)[3]) #BR rep
			{
				RepCount <- RepCount + 1
										
				for (b in 1:dim(MatFileExp$ExpRc)[1]) #pir
				{
					for (c in 1:dim(MatFileExp$ExpRc)[2]) #pid
					{
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"CPSN"] <- MatFileExp$ExpRc[b,c,a,1]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"CP"] <- MatFileExp$ExpRc[b,c,a,2]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpRc)[1] * dim(MatFileExp$ExpRc)[2]) + ((b - 1) * dim(MatFileExp$ExpRc)[2]) + (c - 1)),"SN"] <- MatFileExp$ExpRc[b,c,a,3]
																		
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"CPSN"] <- MatFileSD$SDRc[b,c,a,1]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"CP"] <- MatFileSD$SDRc[b,c,a,2]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDRc)[1] * dim(MatFileSD$SDRc)[2]) + ((b - 1) * dim(MatFileSD$SDRc)[2]) + (c - 1)),"SN"] <- MatFileSD$SDRc[b,c,a,3]						
					}
				}
			}
		}
	}	
}

#get data
#calculate log(E(CPSN)) - log(E(CP)) assuming the E(SN) and E(CP) are log-normally distributed
# Here, for the logged variables mu = log(m^2 / sqrt(v + m^2)) and sigma = sqrt(log((v/m^2) + 1))

SC <- log((ResultsExp[,"CPSN"]^2) / sqrt((ResultsSD[,"CPSN"]^2) + (ResultsExp[,"CPSN"]^2))) - log((ResultsExp[,"CP"]^2) / sqrt((ResultsSD[,"CP"]^2) + (ResultsExp[,"CP"]^2)))
TAUSC <- 1 / (log(((ResultsSD[,"CPSN"]^2) / (ResultsExp[,"CPSN"]^2)) + 1) + log(((ResultsSD[,"CP"]^2) / (ResultsExp[,"CP"]^2)) + 1))
N <- length(SC)
NUMBR <- length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3]

#get BR Rep level covariates 
NST <- as.vector(matrix(NA,nrow=length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpRc)[3],ncol=1))
for (i in 1:(nrow(ResultsExp) - 1))
{
	if ((ResultsExp[i,"BRRep"] != ResultsExp[i + 1,"BRRep"]) | (i == (nrow(ResultsExp) - 1)))
	{
		NST[ResultsExp[i,"BRRep"]] <- 1 - ResultsExp[i,"Nested"]
	}
}
#standardise
NST <- (NST - mean(NST)) / sd(NST)

#get within BRRep level covariates
BR <- ResultsExp[,"BRRep"]
PIR <- ResultsExp[,"PiR"]
PID <- ResultsExp[,"PiD"]
PIRD <- ResultsExp[,"PiR"] * ResultsExp[,"PiD"]
#standardise
PIR <- (PIR - mean(PIR)) / sd(PIR)
PID <- (PID - mean(PID)) / sd(PID)
PIRD <- (PIRD - mean(PIRD)) / sd(PIRD)

#complie data
DataStarLine <- list(NUMBR=NUMBR,N=N,SC=SC,TAUSC=TAUSC,NST=NST,PIR=PIR,PID=PID,BR=BR)

#run jags model
#combine data
data <- list(DataStar,DataLine,DataRing,DataWheel,DataRingLine,DataRingStar,DataStarLine)

source("./jags_code/functions.r")

#run jags
sfInit(parallel=TRUE,cpus=7)

#export data, functions and libraries to workers
sfExportAll()
sfClusterEval(library(runjags))
sfClusterEval(library(coda))
sfClusterEval(library(rjags))
sfClusterEval(library(parallel))

Jags.Fits <- sfLapply(data,get.jags1)

#Jags.Fits <- get.jags1(data[[1]])

sfStop()

#save output
save(Jags.Fits,file="Jags_Fits.RData")
#save(Jags.Fits,file="Jags_Fits_GOF.RData")
