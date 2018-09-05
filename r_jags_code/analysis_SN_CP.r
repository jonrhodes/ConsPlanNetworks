# libraries and functions
library(R.matlab)
library(runjags)
library(coda)
library(snowfall)
library(parallel)
library(rjags)

#set working directory and functions source
setwd("D:/Users/jrhodes/networks_cons_plan_results")
source("./code/functions.r")

#GET THE DATA FOR EACH NETWORK SHAPE 

#STAR

#set up matrices to hold results
ResultsExp <- matrix(NA,ncol=10,nrow=2*2*11*5*11*11)
dimnames(ResultsExp) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN","RND"))
ResultsSD <- matrix(NA,ncol=10,nrow=2*2*11*5*11*11)
dimnames(ResultsSD) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN","RND"))

#create lists of variables
NumSpec <- c(20,60)
Occ <- c(0.1, 0.25)
Nested <- seq(0,1,0.1)
PiR <- seq(0,1,0.1)
PiD <- seq(0,1,0.1)

RepCount <- 0

#loop through number of species
for (i in 1:length(NumSpec))
{
	ResultsExp[(((i - 1) * (2*11*5*11*11)) + 1):(i * (2*11*5*11*11)),"NumSp"] <- NumSpec[i]
	ResultsSD[(((i - 1) * (2*11*5*11*11)) + 1):(i * (2*11*5*11*11)),"NumSp"] <- NumSpec[i]
		
	#loop through occupancy
	for (j in 1:length(Occ))
	{
		ResultsExp[((((i - 1) * (2*11*5*11*11))) + (((j - 1) * (11*5*11*11))) + 1):((((i - 1) * (2*11*5*11*11))) + (j * (11*5*11*11))),"Occ"] <- Occ[j]
		ResultsSD[((((i - 1) * (2*11*5*11*11))) + (((j - 1) * (11*5*11*11))) + 1):((((i - 1) * (2*11*5*11*11))) + (j * (11*5*11*11))),"Occ"] <- Occ[j]
						
		#loop through nested
		for (k in 1:length(Nested))
		{
			StartVal <- ((((i - 1) * (2*11*5*11*11))) + (((j - 1) * (11*5*11*11))) + (((k - 1) * (5*11*11))) + 1)
									
			#get the right .mat file
			# dimensions = [pir,pid,rep_id(BR replicate),method(1=CPSN,2=CP,3=SN,4=RND),number simulations] - no number simulations in the expectation and SD files 
			MatFileExp <- readMat(paste("./results/ExpR_starModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
			MatFileSD <- readMat(paste("./results/SDR_starModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
									
			for (a in 1:dim(MatFileExp$ExpR)[3]) #BR rep
			{
				RepCount <- RepCount + 1
										
				for (b in 1:dim(MatFileExp$ExpR)[1]) #pir
				{
					for (c in 1:dim(MatFileExp$ExpR)[2]) #pid
					{
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"CPSN"] <- MatFileExp$ExpR[b,c,a,1]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"CP"] <- MatFileExp$ExpR[b,c,a,2]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"SN"] <- MatFileExp$ExpR[b,c,a,3]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"RND"] <- MatFileExp$ExpR[b,c,a,4]
												
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"CPSN"] <- MatFileSD$SDR[b,c,a,1]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"CP"] <- MatFileSD$SDR[b,c,a,2]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"SN"] <- MatFileSD$SDR[b,c,a,3]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"RND"] <- MatFileSD$SDR[b,c,a,4]
					}
				}
			}
		}
	}	
}

#get data
#calculate log(E(SN)) - log(E(CP)) assuming the E(SN) and E(CP) are log-normally distributed
# Here, for the logged variables mu = log(m^2 / sqrt(v + m^2)) and sigma = sqrt(log((v/m^2) + 1))

SC <- log((ResultsExp[,"SN"]^2) / sqrt((ResultsSD[,"SN"]^2) + (ResultsExp[,"SN"]^2))) - log((ResultsExp[,"CP"]^2) / sqrt((ResultsSD[,"CP"]^2) + (ResultsExp[,"CP"]^2)))
TAUSC <- 1 / (log(((ResultsSD[,"SN"]^2) / (ResultsExp[,"SN"]^2)) + 1) + log(((ResultsSD[,"CP"]^2) / (ResultsExp[,"CP"]^2)) + 1))
N <- length(SC)
NUMBR <- length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3]

#get BR Rep level covariates 
SNUM <- as.vector(matrix(NA,nrow=length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3],ncol=1))
OCCP <- as.vector(matrix(NA,nrow=length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3],ncol=1))
NST <- as.vector(matrix(NA,nrow=length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3],ncol=1))
for (i in 1:(nrow(ResultsExp) - 1))
{
	if ((ResultsExp[i,"BRRep"] != ResultsExp[i + 1,"BRRep"]) | (i == (nrow(ResultsExp) - 1)))
	{
		SNUM[ResultsExp[i,"BRRep"]] <- ResultsExp[i,"NumSp"]
		OCCP[ResultsExp[i,"BRRep"]] <- ResultsExp[i,"Occ"]
		NST[ResultsExp[i,"BRRep"]] <- ResultsExp[i,"Nested"]
	}
}
#standardise
SNUM <- ifelse(SNUM==20,0,1)
OCCP <- ifelse(OCCP==0.1,0,1)
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
DataStar <- list(NUMBR=NUMBR,N=N,SC=SC,TAUSC=TAUSC,SNUM=SNUM,OCCP=OCCP,NST=NST,PIR=PIR,PID=PID,BR=BR)

#LINE

#set up matrices to hold results
ResultsExp <- matrix(NA,ncol=10,nrow=2*2*11*5*11*11)
dimnames(ResultsExp) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN","RND"))
ResultsSD <- matrix(NA,ncol=10,nrow=2*2*11*5*11*11)
dimnames(ResultsSD) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN","RND"))

#create lists of variables
NumSpec <- c(20,60)
Occ <- c(0.1, 0.25)
Nested <- seq(0,1,0.1)
PiR <- seq(0,1,0.1)
PiD <- seq(0,1,0.1)

RepCount <- 0

#loop through number of species
for (i in 1:length(NumSpec))
{
	ResultsExp[(((i - 1) * (2*11*5*11*11)) + 1):(i * (2*11*5*11*11)),"NumSp"] <- NumSpec[i]
	ResultsSD[(((i - 1) * (2*11*5*11*11)) + 1):(i * (2*11*5*11*11)),"NumSp"] <- NumSpec[i]
		
	#loop through occupancy
	for (j in 1:length(Occ))
	{
		ResultsExp[((((i - 1) * (2*11*5*11*11))) + (((j - 1) * (11*5*11*11))) + 1):((((i - 1) * (2*11*5*11*11))) + (j * (11*5*11*11))),"Occ"] <- Occ[j]
		ResultsSD[((((i - 1) * (2*11*5*11*11))) + (((j - 1) * (11*5*11*11))) + 1):((((i - 1) * (2*11*5*11*11))) + (j * (11*5*11*11))),"Occ"] <- Occ[j]
						
		#loop through nested
		for (k in 1:length(Nested))
		{
			StartVal <- ((((i - 1) * (2*11*5*11*11))) + (((j - 1) * (11*5*11*11))) + (((k - 1) * (5*11*11))) + 1)
									
			#get the right .mat file
			# dimensions = [pir,pid,rep_id(BR replicate),method(1=CPSN,2=CP,3=SN,4=RND),number simulations] - no number simulations in the expectation and SD files 
			MatFileExp <- readMat(paste("./results/ExpR_lineModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
			MatFileSD <- readMat(paste("./results/SDR_lineModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
									
			for (a in 1:dim(MatFileExp$ExpR)[3]) #BR rep
			{
				RepCount <- RepCount + 1
										
				for (b in 1:dim(MatFileExp$ExpR)[1]) #pir
				{
					for (c in 1:dim(MatFileExp$ExpR)[2]) #pid
					{
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"CPSN"] <- MatFileExp$ExpR[b,c,a,1]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"CP"] <- MatFileExp$ExpR[b,c,a,2]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"SN"] <- MatFileExp$ExpR[b,c,a,3]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"RND"] <- MatFileExp$ExpR[b,c,a,4]
												
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"CPSN"] <- MatFileSD$SDR[b,c,a,1]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"CP"] <- MatFileSD$SDR[b,c,a,2]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"SN"] <- MatFileSD$SDR[b,c,a,3]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"RND"] <- MatFileSD$SDR[b,c,a,4]
					}
				}
			}
		}
	}	
}

#get data
#calculate log(E(SN)) - log(E(CP)) assuming the E(SN) and E(CP) are log-normally distributed
# Here, for the logged variables mu = log(m^2 / sqrt(v + m^2)) and sigma = sqrt(log((v/m^2) + 1))

SC <- log((ResultsExp[,"SN"]^2) / sqrt((ResultsSD[,"SN"]^2) + (ResultsExp[,"SN"]^2))) - log((ResultsExp[,"CP"]^2) / sqrt((ResultsSD[,"CP"]^2) + (ResultsExp[,"CP"]^2)))
TAUSC <- 1 / (log(((ResultsSD[,"SN"]^2) / (ResultsExp[,"SN"]^2)) + 1) + log(((ResultsSD[,"CP"]^2) / (ResultsExp[,"CP"]^2)) + 1))
N <- length(SC)
NUMBR <- length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3]

#get BR Rep level covariates 
SNUM <- as.vector(matrix(NA,nrow=length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3],ncol=1))
OCCP <- as.vector(matrix(NA,nrow=length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3],ncol=1))
NST <- as.vector(matrix(NA,nrow=length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3],ncol=1))
for (i in 1:(nrow(ResultsExp) - 1))
{
	if ((ResultsExp[i,"BRRep"] != ResultsExp[i + 1,"BRRep"]) | (i == (nrow(ResultsExp) - 1)))
	{
		SNUM[ResultsExp[i,"BRRep"]] <- ResultsExp[i,"NumSp"]
		OCCP[ResultsExp[i,"BRRep"]] <- ResultsExp[i,"Occ"]
		NST[ResultsExp[i,"BRRep"]] <- ResultsExp[i,"Nested"]
	}
}
#standardise
SNUM <- ifelse(SNUM==20,0,1)
OCCP <- ifelse(OCCP==0.1,0,1)
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
DataLine <- list(NUMBR=NUMBR,N=N,SC=SC,TAUSC=TAUSC,SNUM=SNUM,OCCP=OCCP,NST=NST,PIR=PIR,PID=PID,BR=BR)

#RING

#set up matrices to hold results
ResultsExp <- matrix(NA,ncol=10,nrow=2*2*11*5*11*11)
dimnames(ResultsExp) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN","RND"))
ResultsSD <- matrix(NA,ncol=10,nrow=2*2*11*5*11*11)
dimnames(ResultsSD) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN","RND"))

#create lists of variables
NumSpec <- c(20,60)
Occ <- c(0.1, 0.25)
Nested <- seq(0,1,0.1)
PiR <- seq(0,1,0.1)
PiD <- seq(0,1,0.1)

RepCount <- 0

#loop through number of species
for (i in 1:length(NumSpec))
{
	ResultsExp[(((i - 1) * (2*11*5*11*11)) + 1):(i * (2*11*5*11*11)),"NumSp"] <- NumSpec[i]
	ResultsSD[(((i - 1) * (2*11*5*11*11)) + 1):(i * (2*11*5*11*11)),"NumSp"] <- NumSpec[i]
		
	#loop through occupancy
	for (j in 1:length(Occ))
	{
		ResultsExp[((((i - 1) * (2*11*5*11*11))) + (((j - 1) * (11*5*11*11))) + 1):((((i - 1) * (2*11*5*11*11))) + (j * (11*5*11*11))),"Occ"] <- Occ[j]
		ResultsSD[((((i - 1) * (2*11*5*11*11))) + (((j - 1) * (11*5*11*11))) + 1):((((i - 1) * (2*11*5*11*11))) + (j * (11*5*11*11))),"Occ"] <- Occ[j]
						
		#loop through nested
		for (k in 1:length(Nested))
		{
			StartVal <- ((((i - 1) * (2*11*5*11*11))) + (((j - 1) * (11*5*11*11))) + (((k - 1) * (5*11*11))) + 1)
									
			#get the right .mat file
			# dimensions = [pir,pid,rep_id(BR replicate),method(1=CPSN,2=CP,3=SN,4=RND),number simulations] - no number simulations in the expectation and SD files 
			MatFileExp <- readMat(paste("./results/ExpR_islandModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
			MatFileSD <- readMat(paste("./results/SDR_islandModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
									
			for (a in 1:dim(MatFileExp$ExpR)[3]) #BR rep
			{
				RepCount <- RepCount + 1
										
				for (b in 1:dim(MatFileExp$ExpR)[1]) #pir
				{
					for (c in 1:dim(MatFileExp$ExpR)[2]) #pid
					{
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"CPSN"] <- MatFileExp$ExpR[b,c,a,1]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"CP"] <- MatFileExp$ExpR[b,c,a,2]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"SN"] <- MatFileExp$ExpR[b,c,a,3]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"RND"] <- MatFileExp$ExpR[b,c,a,4]
												
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"CPSN"] <- MatFileSD$SDR[b,c,a,1]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"CP"] <- MatFileSD$SDR[b,c,a,2]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"SN"] <- MatFileSD$SDR[b,c,a,3]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"RND"] <- MatFileSD$SDR[b,c,a,4]
					}
				}
			}
		}
	}	
}

#get data
#calculate log(E(SN)) - log(E(CP)) assuming the E(SN) and E(CP) are log-normally distributed
# Here, for the logged variables mu = log(m^2 / sqrt(v + m^2)) and sigma = sqrt(log((v/m^2) + 1))

SC <- log((ResultsExp[,"SN"]^2) / sqrt((ResultsSD[,"SN"]^2) + (ResultsExp[,"SN"]^2))) - log((ResultsExp[,"CP"]^2) / sqrt((ResultsSD[,"CP"]^2) + (ResultsExp[,"CP"]^2)))
TAUSC <- 1 / (log(((ResultsSD[,"SN"]^2) / (ResultsExp[,"SN"]^2)) + 1) + log(((ResultsSD[,"CP"]^2) / (ResultsExp[,"CP"]^2)) + 1))
N <- length(SC)
NUMBR <- length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3]

#get BR Rep level covariates 
SNUM <- as.vector(matrix(NA,nrow=length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3],ncol=1))
OCCP <- as.vector(matrix(NA,nrow=length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3],ncol=1))
NST <- as.vector(matrix(NA,nrow=length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3],ncol=1))
for (i in 1:(nrow(ResultsExp) - 1))
{
	if ((ResultsExp[i,"BRRep"] != ResultsExp[i + 1,"BRRep"]) | (i == (nrow(ResultsExp) - 1)))
	{
		SNUM[ResultsExp[i,"BRRep"]] <- ResultsExp[i,"NumSp"]
		OCCP[ResultsExp[i,"BRRep"]] <- ResultsExp[i,"Occ"]
		NST[ResultsExp[i,"BRRep"]] <- ResultsExp[i,"Nested"]
	}
}
#standardise
SNUM <- ifelse(SNUM==20,0,1)
OCCP <- ifelse(OCCP==0.1,0,1)
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
DataRing <- list(NUMBR=NUMBR,N=N,SC=SC,TAUSC=TAUSC,SNUM=SNUM,OCCP=OCCP,NST=NST,PIR=PIR,PID=PID,BR=BR)

#WHEEL

#set up matrices to hold results
ResultsExp <- matrix(NA,ncol=10,nrow=2*2*11*5*11*11)
dimnames(ResultsExp) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN","RND"))
ResultsSD <- matrix(NA,ncol=10,nrow=2*2*11*5*11*11)
dimnames(ResultsSD) <- list(NULL,c("NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN","RND"))

#create lists of variables
NumSpec <- c(20,60)
Occ <- c(0.1, 0.25)
Nested <- seq(0,1,0.1)
PiR <- seq(0,1,0.1)
PiD <- seq(0,1,0.1)

RepCount <- 0

#loop through number of species
for (i in 1:length(NumSpec))
{
	ResultsExp[(((i - 1) * (2*11*5*11*11)) + 1):(i * (2*11*5*11*11)),"NumSp"] <- NumSpec[i]
	ResultsSD[(((i - 1) * (2*11*5*11*11)) + 1):(i * (2*11*5*11*11)),"NumSp"] <- NumSpec[i]
		
	#loop through occupancy
	for (j in 1:length(Occ))
	{
		ResultsExp[((((i - 1) * (2*11*5*11*11))) + (((j - 1) * (11*5*11*11))) + 1):((((i - 1) * (2*11*5*11*11))) + (j * (11*5*11*11))),"Occ"] <- Occ[j]
		ResultsSD[((((i - 1) * (2*11*5*11*11))) + (((j - 1) * (11*5*11*11))) + 1):((((i - 1) * (2*11*5*11*11))) + (j * (11*5*11*11))),"Occ"] <- Occ[j]
						
		#loop through nested
		for (k in 1:length(Nested))
		{
			StartVal <- ((((i - 1) * (2*11*5*11*11))) + (((j - 1) * (11*5*11*11))) + (((k - 1) * (5*11*11))) + 1)
									
			#get the right .mat file
			# dimensions = [pir,pid,rep_id(BR replicate),method(1=CPSN,2=CP,3=SN,4=RND),number simulations] - no number simulations in the expectation and SD files 
			MatFileExp <- readMat(paste("./results/ExpR_wheelModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
			MatFileSD <- readMat(paste("./results/SDR_wheelModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
									
			for (a in 1:dim(MatFileExp$ExpR)[3]) #BR rep
			{
				RepCount <- RepCount + 1
										
				for (b in 1:dim(MatFileExp$ExpR)[1]) #pir
				{
					for (c in 1:dim(MatFileExp$ExpR)[2]) #pid
					{
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"CPSN"] <- MatFileExp$ExpR[b,c,a,1]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"CP"] <- MatFileExp$ExpR[b,c,a,2]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"SN"] <- MatFileExp$ExpR[b,c,a,3]
						ResultsExp[(StartVal + ((a - 1) * dim(MatFileExp$ExpR)[1] * dim(MatFileExp$ExpR)[2]) + ((b - 1) * dim(MatFileExp$ExpR)[2]) + (c - 1)),"RND"] <- MatFileExp$ExpR[b,c,a,4]
												
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"Nested"] <- Nested[k]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"BRRep"] <- RepCount
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"PiR"] <- PiR[b]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"PiD"] <- PiD[c]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"CPSN"] <- MatFileSD$SDR[b,c,a,1]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"CP"] <- MatFileSD$SDR[b,c,a,2]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"SN"] <- MatFileSD$SDR[b,c,a,3]
						ResultsSD[(StartVal + ((a - 1) * dim(MatFileSD$SDR)[1] * dim(MatFileSD$SDR)[2]) + ((b - 1) * dim(MatFileSD$SDR)[2]) + (c - 1)),"RND"] <- MatFileSD$SDR[b,c,a,4]
					}
				}
			}
		}
	}	
}

#get data
#calculate log(E(SN)) - log(E(CP)) assuming the E(SN) and E(CP) are log-normally distributed
# Here, for the logged variables mu = log(m^2 / sqrt(v + m^2)) and sigma = sqrt(log((v/m^2) + 1))

SC <- log((ResultsExp[,"SN"]^2) / sqrt((ResultsSD[,"SN"]^2) + (ResultsExp[,"SN"]^2))) - log((ResultsExp[,"CP"]^2) / sqrt((ResultsSD[,"CP"]^2) + (ResultsExp[,"CP"]^2)))
TAUSC <- 1 / (log(((ResultsSD[,"SN"]^2) / (ResultsExp[,"SN"]^2)) + 1) + log(((ResultsSD[,"CP"]^2) / (ResultsExp[,"CP"]^2)) + 1))
N <- length(SC)
NUMBR <- length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3]

#get BR Rep level covariates 
SNUM <- as.vector(matrix(NA,nrow=length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3],ncol=1))
OCCP <- as.vector(matrix(NA,nrow=length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3],ncol=1))
NST <- as.vector(matrix(NA,nrow=length(NumSpec) * length(Occ) * length(Nested) * dim(MatFileExp$ExpR)[3],ncol=1))
for (i in 1:(nrow(ResultsExp) - 1))
{
	if ((ResultsExp[i,"BRRep"] != ResultsExp[i + 1,"BRRep"]) | (i == (nrow(ResultsExp) - 1)))
	{
		SNUM[ResultsExp[i,"BRRep"]] <- ResultsExp[i,"NumSp"]
		OCCP[ResultsExp[i,"BRRep"]] <- ResultsExp[i,"Occ"]
		NST[ResultsExp[i,"BRRep"]] <- ResultsExp[i,"Nested"]
	}
}
#standardise
SNUM <- ifelse(SNUM==20,0,1)
OCCP <- ifelse(OCCP==0.1,0,1)
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
DataWheel <- list(NUMBR=NUMBR,N=N,SC=SC,TAUSC=TAUSC,SNUM=SNUM,OCCP=OCCP,NST=NST,PIR=PIR,PID=PID,BR=BR)

#run jags model
#combine data
data <- list(DataStar,DataWheel,DataLine,DataRing)

source("./code/functions.r")

#run jags
sfInit(parallel=TRUE,cpus=4)

#export data, functions and libraries to workers
sfExportAll()
sfClusterEval(library(runjags))
sfClusterEval(library(coda))
sfClusterEval(library(rjags))
sfClusterEval(library(parallel))

Jags.Fits <- sfLapply(data,get.jags1)

#Jags.Fits <- get.jags2(data[[1]])

sfStop()

#save output
save(Jags.Fits,file="Jags_Fits_SC.RData")

