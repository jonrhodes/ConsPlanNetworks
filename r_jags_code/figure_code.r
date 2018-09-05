#get libraries
library(ggplot2)
library(rjags)
library(runjags)
library(coda)
library(R.matlab)
library(igraph)
library(easyGgplot2)
library(GGally)
library(network)
library(sna)
library(RColorBrewer)
library(intergraph) 
library(vegan)
library(foreign)

#set working directory and functions source
setwd("M:/Users/uqjrhode/work/networks_cons_plan/code")
source("./functions.r")

#generate shape figures
#create networks
Star <- graph_from_literal(1-2,1-3,1-4,1-5,1-6)
Wheel <- graph_from_literal(1-2,1-3,1-4,1-5,1-6,2-3,3-4,4-5,5-6,6-2)
Line <- graph_from_literal(1-2,2-3,3-4,4-5,5-6)
Ring <- graph_from_literal(1-2,2-3,3-4,4-5,5-6,6-1)
RingStar <- graph_from_literal(1-2,1-3,2-3,1-4,4-5,4-6)
StarLine <- graph_from_literal(1-2,1-3,1-4,4-5,5-6)
RingLine <- graph_from_literal(1-2,1-3,2-3,1-4,4-5,5-6)
CaseStudy <- graph_from_literal(1-3,1-2,2-3,3-4,5)

plot(Star,layout=layout_nicely,vertex.size=20,vertex.color="red",vertex.label=NA,edge.color="black",edge.size=10)
plot(Wheel,layout=layout_nicely,vertex.size=20,vertex.color="red",vertex.label=NA,edge.color="black",edge.size=10)
plot(Line,layout=layout_nicely,vertex.size=20,vertex.color="red",vertex.label=NA,edge.color="black",edge.size=10)
plot(Ring,layout=layout_nicely,vertex.size=20,vertex.color="red",vertex.label=NA,edge.color="black",edge.size=10)
plot(RingStar,layout=layout_nicely,vertex.size=20,vertex.color="red",vertex.label=NA,edge.color="black",edge.size=10)
plot(StarLine,layout=layout_nicely,vertex.size=20,vertex.color="red",vertex.label=NA,edge.color="black",edge.size=10)
plot(RingLine,layout=layout_nicely,vertex.size=20,vertex.color="red",vertex.label=NA,edge.color="black",edge.size=10)

#NETWORKS - FIGURE 1A

CentDens <- matrix(NA,nrow=10,ncol=3)
dimnames(CentDens) <- list(NULL,c("Shape","Metric","Value"))
CentDens <- as.data.frame(CentDens)
CentDens$Shape <- c("Ring","Ring","Line","Line","Ring-Star","Ring-Star","Wheel","Wheel","Star","Star")
CentDens$Shape <- factor(CentDens$Shape,levels=c(c("Ring","Line","Ring-Star","Wheel","Star")))
CentDens$Metric <- c("Closeness Centrality","Density","Closeness Centrality","Density","Closeness Centrality","Density","Closeness Centrality","Density","Closeness Centrality","Density")
CentDens$Value <- c(0.000001,0.4,0.29,0.33,0.43,0.4,0.64,0.67,1,0.33)
p <- ggplot(data = CentDens,aes(x = Shape, y = Value,fill = Metric))
p <- p + scale_fill_manual(values=c("gold","deepskyblue"),name="") + theme_bw()
p <- p + geom_bar(stat = "identity",position = position_dodge(0.9)) + labs(y = "Value",x="") + scale_y_continuous(breaks=seq(0,1.2,0.2),limits=c(0,1.2)) 
p <- p + theme(axis.text=element_text(size=14,face="bold",angle=45,vjust=0,hjust=1),legend.text=element_text(size=12,face="bold"),axis.title.y=element_text(size=14,face="bold"),legend.position="top") 
p

#BIODIVERSITY NESTEDNESS - FIGURE 1B
#get data
BRNested0.9 <- readMat("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/ConsMDP/BR_library/BRSites6NumSp20Occ0.25Nested0.9Rep1.mat")  
BRNested0.5 <- readMat("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/ConsMDP/BR_library/BRSites6NumSp20Occ0.25Nested0.5Rep1.mat")
BRNested0.1 <- readMat("M:/Users/uqjrhode/work/networks_cons_plan/spudd_final_run/ConsMDP/BR_library/BRSites6NumSp20Occ0.25Nested0.1Rep1.mat")
#nestedness plots
par(mfrow = c(3,1))
image(seq(1,20,1),seq(1,6,1),t(BRNested0.1$BR[[6]][sort(rowSums(BRNested0.1$BR[[6]]),index=T)$ix,]),col=c("white","red"),bty="l",xlab="",ylab="")
at <- seq(from = 0, to = 20, by = 1)
axis(side = 1, at = at)
title(ylab="Node Identity",cex.lab=1.5,line=2.5)
title(xlab="Species Identity",cex.lab=1.5,line=2.9)
title(main="Nested",cex.main=2)
box()
image(seq(1,20,1),seq(1,6,1),t(BRNested0.5$BR[[6]][sort(rowSums(BRNested0.5$BR[[6]]),index=T)$ix,]),col=c("white","red"),bty="l",xlab="",ylab="")
at <- seq(from = 0, to = 20, by = 1)
axis(side = 1, at = at)
title(ylab="Node Identity",cex.lab=1.5,line=2.5)
title(xlab="Species Identity",cex.lab=1.5,line=2.9)
title(main="Random",cex.main=2)
box()
image(seq(1,20,1),seq(1,6,1),t(BRNested0.9$BR[[6]][sort(rowSums(BRNested0.9$BR[[6]]),index=T)$ix,]),col=c("white","red"),bty="l",xlab="",ylab="")
at <- seq(from = 0, to = 20, by = 1)
axis(side = 1, at = at)
title(ylab="Node Identity",cex.lab=1.5,line=2.5)
title(xlab="Species Identity",cex.lab=1.5,line=2.9)
title(main="Un-nested",cex.main=2)
box()

#get MCMC output for the model fits
load("M:/Users/uqjrhode/work/networks_cons_plan/Jags_Fits.RData")
Star <- summary(Jags.Fits[[1]])
Line <- summary(Jags.Fits[[2]])
Ring <- summary(Jags.Fits[[3]])
Wheel <- summary(Jags.Fits[[4]])
RingLine <- summary(Jags.Fits[[5]])
RingStar <- summary(Jags.Fits[[6]])
StarLine <- summary(Jags.Fits[[7]])

#get MCMC output for predictionsPID
load("M:/Users/uqjrhode/work/networks_cons_plan/Jags_Fits_PredPID.RData")
StarPredPID <- summary(Jags.Fits[[1]])
LinePredPID <- summary(Jags.Fits[[2]])
RingPredPID <- summary(Jags.Fits[[3]])
WheelPredPID <- summary(Jags.Fits[[4]])
RingLinePredPID <- summary(Jags.Fits[[5]])
RingStarPredPID <- summary(Jags.Fits[[6]])
StarLinePredPID <- summary(Jags.Fits[[7]])

#get MCMC output for predictionsPIR
load("M:/Users/uqjrhode/work/networks_cons_plan/Jags_Fits_PredPIR.RData")
StarPredPIR <- summary(Jags.Fits[[1]])
LinePredPIR <- summary(Jags.Fits[[2]])
RingPredPIR <- summary(Jags.Fits[[3]])
WheelPredPIR <- summary(Jags.Fits[[4]])
RingLinePredPIR <- summary(Jags.Fits[[5]])
RingStarPredPIR <- summary(Jags.Fits[[6]])
StarLinePredPIR <- summary(Jags.Fits[[7]])

#VOI FOR SHAPE - FIGURE 2

DataAlpha <- as.data.frame(rbind(Ring["alpha",c("Mean","Lower95","Upper95")],Line["alpha",c("Mean","Lower95","Upper95")],RingStar["alpha",c("Mean","Lower95","Upper95")],Wheel["alpha",c("Mean","Lower95","Upper95")],Star["alpha",c("Mean","Lower95","Upper95")]))
#get percentage change
DataAlpha[,"Mean"] <- (exp(DataAlpha[,"Mean"]) - 1) * 100
DataAlpha[,"Lower95"] <- (exp(DataAlpha[,"Lower95"]) - 1) * 100
DataAlpha[,"Upper95"] <- (exp(DataAlpha[,"Upper95"]) - 1) * 100 

Names <- c("Ring","Line","Ring-Star","Wheel","Star")
Names <- factor(Names,levels = Names[order(c(1,2,3,4,5))])

p <- ggplot(data=DataAlpha,aes(x=Names,y=Mean))
p <- p + geom_bar(fill="red",stat="identity",position=position_dodge(0.9)) + xlab("") + ylab("Value of Network Information (%)") + scale_fill_discrete(name = "") + geom_errorbar(aes(ymin=DataAlpha[,"Lower95"],ymax=DataAlpha[,"Upper95"]),width=.25,position=position_dodge(.9)) + theme_bw() + theme(axis.text=element_text(size=12,face="bold"),legend.text=element_text(size=12),axis.title.y=element_text(size=12,face="bold"),axis.title.x=element_text(size=12,face="bold")) + scale_y_continuous(breaks=seq(0,16,2),limits=c(0,16))  
p

#COEFFICIENTS FOR NESTEDNESS AND TIE STRENGTH - FIGURE 3

DataCoeffs <- as.data.frame(rbind(Ring[c("nst","alpir","alpid","nspir","nspid"),c("Mean","Lower95","Upper95")],Line[c("nst","alpir","alpid","nspir","nspid"),c("Mean","Lower95","Upper95")],RingStar[c("nst","alpir","alpid","nspir","nspid"),c("Mean","Lower95","Upper95")],Wheel[c("nst","alpir","alpid","nspir","nspid"),c("Mean","Lower95","Upper95")],Star[c("nst","alpir","alpid","nspir","nspid"),c("Mean","Lower95","Upper95")]))
DataCoeffs$Shape <- c("Ring","Ring","Ring","Ring","Ring","Line","Line","Line","Line","Line","Ring-Star","Ring-Star","Ring-Star","Ring-Star","Ring-Star","Wheel","Wheel","Wheel","Wheel","Wheel","Star","Star","Star","Star","Star")
DataCoeffs$Shape <- factor(DataCoeffs$Shape,levels=c("Ring","Line","Ring-Star","Wheel","Star"))
DataCoeffs$Coeff <- c("Nestedness","Reserve Influence","Development Influence","Nest x Res Interaction","Nest x Dev Interaction","Nestedness","Reserve Influence","Development Influence","Nest x Res Interaction","Nest x Dev Interaction","Nestedness","Reserve Influence","Development Influence","Nest x Res Interaction","Nest x Dev Interaction","Nestedness","Reserve Influence","Development Influence","Nest x Res Interaction","Nest x Dev Interaction","Nestedness","Reserve Influence","Development Influence","Nest x Res Interaction","Nest x Dev Interaction")
DataCoeffs$Coeff <- factor(DataCoeffs$Coeff,levels=c("Nestedness","Reserve Influence","Development Influence","Nest x Res Interaction","Nest x Dev Interaction"))
p <- ggplot(data = DataCoeffs,aes(x = Coeff,y = Mean,fill = Shape))
p <- p + geom_bar(stat = "identity",position = position_dodge(0.9)) + geom_errorbar(aes(ymin=DataCoeffs[,"Lower95"],ymax=DataCoeffs[,"Upper95"]),width=.25,position=position_dodge(.9))
p <- p + ylab("Effect Size") + xlab("") + scale_fill_discrete(name = "") + theme_bw()
p <- p + theme(axis.text=element_text(size=12,face="bold",angle=45,vjust=0,hjust=1),legend.text=element_text(size=12),axis.title.y=element_text(size=12,face="bold"),axis.title.x=element_text(size=12,face="bold"),legend.position="top")
p

#VALUE OF INFORMATION FOR DIFFERENT SHAPES - PID - FIGURE 4
DataVOIRing <- (exp(as.data.frame(RingPredPID[,c("Mean","Lower95","Upper95")])) - 1) * 100
DataVOILine <- (exp(as.data.frame(LinePredPID[,c("Mean","Lower95","Upper95")])) - 1) * 100
DataVOIWheel <- (exp(as.data.frame(WheelPredPID[,c("Mean","Lower95","Upper95")])) - 1) * 100
DataVOIStar <- (exp(as.data.frame(StarPredPID[,c("Mean","Lower95","Upper95")])) - 1) * 100
DataVOIRing$Nestedness <- c("0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9")
DataVOIRing$Infl <- c(0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3,0.3,0.4,0.4,0.4,0.5,0.5,0.5,0.6,0.6,0.6,0.7,0.7,0.7,0.8,0.8,0.8,0.9,0.9,0.9,1.0,1.0,1.0)
DataVOILine$Nestedness <- c("0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9")
DataVOILine$Infl <- c(0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3,0.3,0.4,0.4,0.4,0.5,0.5,0.5,0.6,0.6,0.6,0.7,0.7,0.7,0.8,0.8,0.8,0.9,0.9,0.9,1.0,1.0,1.0)
DataVOIWheel$Nestedness <- c("0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9")
DataVOIWheel$Infl <- c(0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3,0.3,0.4,0.4,0.4,0.5,0.5,0.5,0.6,0.6,0.6,0.7,0.7,0.7,0.8,0.8,0.8,0.9,0.9,0.9,1.0,1.0,1.0)
DataVOIStar$Nestedness <- c("0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9")
DataVOIStar$Infl <- c(0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3,0.3,0.4,0.4,0.4,0.5,0.5,0.5,0.6,0.6,0.6,0.7,0.7,0.7,0.8,0.8,0.8,0.9,0.9,0.9,1.0,1.0,1.0)
p1 <- ggplot(data = DataVOIRing,aes(x=Infl,y=Mean,colour=Nestedness,group=Nestedness)) + geom_errorbar(aes(ymin=DataVOIRing[,"Lower95"],ymax=DataVOIRing[,"Upper95"]),colour="black",width=.1,position=position_dodge(0.1)) + geom_line(position=position_dodge(0.1)) + geom_point(position=position_dodge(0.1), size=1, shape=21, fill="white")+ theme_bw() + ggtitle("A.  Ring") + xlab("Development Influence Probability") + ylab("Value of Network Information (%)") + theme(plot.title = element_text(size=16, hjust=0),axis.text=element_text(size=8),legend.text=element_text(size=8),legend.title=element_text(size=8),axis.title.y=element_text(size=8),axis.title.x=element_text(size=8),legend.position="bottom") + scale_y_continuous(breaks=seq(-5,45,5),limits=c(-5,45))
p2 <- ggplot(data = DataVOILine,aes(x=Infl,y=Mean,colour=Nestedness,group=Nestedness)) + geom_errorbar(aes(ymin=DataVOILine[,"Lower95"],ymax=DataVOILine[,"Upper95"]),colour="black",width=.1,position=position_dodge(0.1)) + geom_line(position=position_dodge(0.1)) + geom_point(position=position_dodge(0.1), size=1, shape=21, fill="white")+ theme_bw() + ggtitle("B.  Line") + xlab("Development Influence Probability") + ylab("Value of Network Information (%)") + theme(plot.title = element_text(size=16, hjust=0),axis.text=element_text(size=8),legend.text=element_text(size=8),legend.title=element_text(size=8),axis.title.y=element_text(size=8),axis.title.x=element_text(size=8),legend.position="bottom") + scale_y_continuous(breaks=seq(-5,45,5),limits=c(-5,45))
p3 <- ggplot(data = DataVOIWheel,aes(x=Infl,y=Mean,colour=Nestedness,group=Nestedness)) + geom_errorbar(aes(ymin=DataVOIWheel[,"Lower95"],ymax=DataVOIWheel[,"Upper95"]),colour="black",width=.1,position=position_dodge(0.1)) + geom_line(position=position_dodge(0.1)) + geom_point(position=position_dodge(0.1), size=1, shape=21, fill="white")+ theme_bw() + ggtitle("C.  Wheel") + xlab("Development Influence Probability") + ylab("Value of Network Information (%)") + theme(plot.title = element_text(size=16, hjust=0),axis.text=element_text(size=8),legend.text=element_text(size=8),legend.title=element_text(size=8),axis.title.y=element_text(size=8),axis.title.x=element_text(size=8),legend.position="bottom") + scale_y_continuous(breaks=seq(-5,45,5),limits=c(-5,45)) 
p4 <- ggplot(data = DataVOIStar,aes(x=Infl,y=Mean,colour=Nestedness,group=Nestedness)) + geom_errorbar(aes(ymin=DataVOIStar[,"Lower95"],ymax=DataVOIStar[,"Upper95"]),colour="black",width=.1,position=position_dodge(0.1)) + geom_line(position=position_dodge(0.1)) + geom_point(position=position_dodge(0.1), size=1, shape=21, fill="white")+ theme_bw() + ggtitle("D.  Star") + xlab("Development Influence Probability") + ylab("Value of Network Information (%)") + theme(plot.title = element_text(size=16, hjust=0),axis.text=element_text(size=8),legend.text=element_text(size=8),legend.title=element_text(size=8),axis.title.y=element_text(size=8),axis.title.x=element_text(size=8),legend.position="bottom") + scale_y_continuous(breaks=seq(-5,45,5),limits=c(-5,45))
ggplot2.multiplot(p1,p2,p3,p4)

#VALUE OF INFORMATION FOR DIFFERENT SHAPES - PIR - FIGURE 5
DataVOIRing <- (exp(as.data.frame(RingPredPIR[,c("Mean","Lower95","Upper95")])) - 1) * 100
DataVOILine <- (exp(as.data.frame(LinePredPIR[,c("Mean","Lower95","Upper95")])) - 1) * 100
DataVOIWheel <- (exp(as.data.frame(WheelPredPIR[,c("Mean","Lower95","Upper95")])) - 1) * 100
DataVOIStar <- (exp(as.data.frame(StarPredPIR[,c("Mean","Lower95","Upper95")])) - 1) * 100
DataVOIRing$Nestedness <- c("0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9")
DataVOIRing$Infl <- c(0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3,0.3,0.4,0.4,0.4,0.5,0.5,0.5,0.6,0.6,0.6,0.7,0.7,0.7,0.8,0.8,0.8,0.9,0.9,0.9,1.0,1.0,1.0)
DataVOILine$Nestedness <- c("0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9")
DataVOILine$Infl <- c(0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3,0.3,0.4,0.4,0.4,0.5,0.5,0.5,0.6,0.6,0.6,0.7,0.7,0.7,0.8,0.8,0.8,0.9,0.9,0.9,1.0,1.0,1.0)
DataVOIWheel$Nestedness <- c("0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9")
DataVOIWheel$Infl <- c(0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3,0.3,0.4,0.4,0.4,0.5,0.5,0.5,0.6,0.6,0.6,0.7,0.7,0.7,0.8,0.8,0.8,0.9,0.9,0.9,1.0,1.0,1.0)
DataVOIStar$Nestedness <- c("0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9","0.1","0.5","0.9")
DataVOIStar$Infl <- c(0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3,0.3,0.4,0.4,0.4,0.5,0.5,0.5,0.6,0.6,0.6,0.7,0.7,0.7,0.8,0.8,0.8,0.9,0.9,0.9,1.0,1.0,1.0)
p1 <- ggplot(data = DataVOIRing,aes(x=Infl,y=Mean,colour=Nestedness,group=Nestedness)) + geom_errorbar(aes(ymin=DataVOIRing[,"Lower95"],ymax=DataVOIRing[,"Upper95"]),colour="black",width=.1,position=position_dodge(0.1)) + geom_line(position=position_dodge(0.1)) + geom_point(position=position_dodge(0.1), size=1, shape=21, fill="white")+ theme_bw() + ggtitle("A.  Ring") + xlab("Reserve Influence Probability") + ylab("Value of Network Information (%)") + theme(plot.title = element_text(size=16, hjust=0),axis.text=element_text(size=8),legend.text=element_text(size=8),legend.title=element_text(size=8),axis.title.y=element_text(size=8),axis.title.x=element_text(size=8),legend.position="bottom") + scale_y_continuous(breaks=seq(-5,35,5),limits=c(-5,35))
p2 <- ggplot(data = DataVOILine,aes(x=Infl,y=Mean,colour=Nestedness,group=Nestedness)) + geom_errorbar(aes(ymin=DataVOILine[,"Lower95"],ymax=DataVOILine[,"Upper95"]),colour="black",width=.1,position=position_dodge(0.1)) + geom_line(position=position_dodge(0.1)) + geom_point(position=position_dodge(0.1), size=1, shape=21, fill="white")+ theme_bw() + ggtitle("B.  Line") + xlab("Reserve Influence Probability") + ylab("Value of Network Information (%)") + theme(plot.title = element_text(size=16, hjust=0),axis.text=element_text(size=8),legend.text=element_text(size=8),legend.title=element_text(size=8),axis.title.y=element_text(size=8),axis.title.x=element_text(size=8),legend.position="bottom") + scale_y_continuous(breaks=seq(-5,35,5),limits=c(-5,35))
p3 <- ggplot(data = DataVOIWheel,aes(x=Infl,y=Mean,colour=Nestedness,group=Nestedness)) + geom_errorbar(aes(ymin=DataVOIWheel[,"Lower95"],ymax=DataVOIWheel[,"Upper95"]),colour="black",width=.1,position=position_dodge(0.1)) + geom_line(position=position_dodge(0.1)) + geom_point(position=position_dodge(0.1), size=1, shape=21, fill="white")+ theme_bw() + ggtitle("C.  Wheel") + xlab("Reserve Influence Probability") + ylab("Value of Network Information (%)") + theme(plot.title = element_text(size=16, hjust=0),axis.text=element_text(size=8),legend.text=element_text(size=8),legend.title=element_text(size=8),axis.title.y=element_text(size=8),axis.title.x=element_text(size=8),legend.position="bottom") + scale_y_continuous(breaks=seq(-5,35,5),limits=c(-5,35)) 
p4 <- ggplot(data = DataVOIStar,aes(x=Infl,y=Mean,colour=Nestedness,group=Nestedness)) + geom_errorbar(aes(ymin=DataVOIStar[,"Lower95"],ymax=DataVOIStar[,"Upper95"]),colour="black",width=.1,position=position_dodge(0.1)) + geom_line(position=position_dodge(0.1)) + geom_point(position=position_dodge(0.1), size=1, shape=21, fill="white")+ theme_bw() + ggtitle("D.  Star") + xlab("Reserve Influence Probability") + ylab("Value of Network Information (%)") + theme(plot.title = element_text(size=16, hjust=0),axis.text=element_text(size=8),legend.text=element_text(size=8),legend.title=element_text(size=8),axis.title.y=element_text(size=8),axis.title.x=element_text(size=8),legend.position="bottom") + scale_y_continuous(breaks=seq(-5,35,5),limits=c(-5,35))
ggplot2.multiplot(p1,p2,p3,p4)

#FISHERY NETWORK - FIGURE 6

#get fishery network and plot figure
Fishnet <- matrix(0,nrow=5,ncol=5)
#1 = gill net, 2 = spear gun, 3 = hand line, 4 = ring net, 5 = seine net
Fishnet[1,4]<-0.43
Fishnet[3,4]<-0.72
Fishnet[3,5]<-1.01
Fishnet[4,1]<-0.43
Fishnet[4,3]<-0.72
Fishnet[4,5]<-0.36
Fishnet[5,3]<-1.01
Fishnet[5,4]<-0.36
Fishnet <- Fishnet * 5
Fishnet <- network(Fishnet,directed = FALSE,ignore.eval=FALSE,names.eval="weight")
network.vertex.names(Fishnet) <- c("Gill net (9 species)","Spear gun (8 species)","Hand line (6 species)","Ring net (9 species)","Seine net (12 species)")
Fishnet %v% "species" <- c(9,8,6,9,12)
ggnet2(Fishnet,node.size = "species",node.color = "tomato", edge.color = "black",edge.size="weight",shape=15) + guides(color = FALSE, size = FALSE)

#NESTED FIGURE S1
par(mfrow = c(2,2))
#plot 1
plot(c(0,1),c(0,1),xlab="",ylab="",type="n",bty="l")
title(xlab=expression(phi),cex.lab=1.2,line=2.5)
title(ylab=expression(Un-normalized~italic(w[j])),cex.lab=1.2,line=2.5)
title("A. 2nd node",adj=0)
Phi <- seq(0,1,0.01)
lines(Phi,(0 - (1-Phi) * 0 + (1-Phi) * abs(0 - 1)),lwd=2,col=1) 
lines(Phi,(1 - (1-Phi) * 1 + (1-Phi) * abs(1 - 1)),lwd=2,col=2)
legend(0.3,1.1,legend=c("0","1"),lty=c(1,1),col=c(1,2),lwd=c(2,2),bty="n") 
#plot 2
plot(c(0,1),c(0,1),xlab="",ylab="",type="n",bty="l")
title(xlab=expression(phi),cex.lab=1.2,line=2.5)
title(ylab=expression(Un-normalized~italic(w[j])),cex.lab=1.2,line=2.5)
title("B. 3rd node",adj=0)
Phi <- seq(0,1,0.01)
lines(Phi,(0 - (1-Phi) * 0 + (1-Phi) * abs(0 - 1))^2,lwd=2,col=1) 
lines(Phi,(0 - (1-Phi) * 0 + (1-Phi) * abs(0 - 1)) * (1 - (1-Phi) * 1 + (1-Phi) * abs(1 - 1)),lwd=2,col=2)
lines(Phi,(1 - (1-Phi) * 1 + (1-Phi) * abs(1 - 1))^2,lwd=2,col=3)
legend(0.3,1.1,legend=c("0","1","2"),lty=c(1,1,1),col=c(1,2,3),lwd=c(2,2,2),bty="n") 
#plot 3
plot(c(0,1),c(0,1),xlab="",ylab="",type="n",bty="l")
title(xlab=expression(phi),cex.lab=1.2,line=2.5)
title(ylab=expression(Un-normalized~italic(w[j])),cex.lab=1.2,line=2.5)
title("C. 4th node",adj=0)
Phi <- seq(0,1,0.01)
lines(Phi,(0 - (1-Phi) * 0 + (1-Phi) * abs(0 - 1))^3,lwd=2,col=1) 
lines(Phi,(0 - (1-Phi) * 0 + (1-Phi) * abs(0 - 1))^2 * (1 - (1-Phi) * 1 + (1-Phi) * abs(1 - 1)),lwd=2,col=2)
lines(Phi,(0 - (1-Phi) * 0 + (1-Phi) * abs(0 - 1))^1 * (1 - (1-Phi) * 1 + (1-Phi) * abs(1 - 1))^2,lwd=2,col=3)
lines(Phi,(1 - (1-Phi) * 1 + (1-Phi) * abs(1 - 1))^3,lwd=2,col=4)
legend(0.3,1.1,legend=c("0","1","2","3"),lty=c(1,1,1,1),col=c(1,2,3,4),lwd=c(2,2,2,2),bty="n")
#plot 4
plot(c(0,1),c(0,1),xlab="",ylab="",type="n",bty="l")
title(xlab=expression(phi),cex.lab=1.2,line=2.5)
title(ylab=expression(Un-normalized~italic(w[j])),cex.lab=1.2,line=2.5)
title("D. 5th node",adj=0)
Phi <- seq(0,1,0.01)
lines(Phi,(0 - (1-Phi) * 0 + (1-Phi) * abs(0 - 1))^4,lwd=2,col=1) 
lines(Phi,(0 - (1-Phi) * 0 + (1-Phi) * abs(0 - 1))^3 * (1 - (1-Phi) * 1 + (1-Phi) * abs(1 - 1)),lwd=2,col=2)
lines(Phi,(0 - (1-Phi) * 0 + (1-Phi) * abs(0 - 1))^2 * (1 - (1-Phi) * 1 + (1-Phi) * abs(1 - 1))^2,lwd=2,col=3)
lines(Phi,(0 - (1-Phi) * 0 + (1-Phi) * abs(0 - 1))^1 * (1 - (1-Phi) * 1 + (1-Phi) * abs(1 - 1))^3,lwd=2,col=4)
lines(Phi,(1 - (1-Phi) * 1 + (1-Phi) * abs(1 - 1))^4,lwd=2,col=5)
legend(0.3,1.1,legend=c("0","1","2","3","4"),lty=c(1,1,1,1,1),col=c(1,2,3,4,5),lwd=c(2,2,2,2,2),bty="n")

#NESTED PERFORMANCE FIGURE S2

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

#write nested parameters and empirical estimates 
write.csv(SimNest,"M:/Users/uqjrhode/work/networks_cons_plan/figures/simnest.csv",row.names=F)

#reload simulated results
Nested <- read.csv("M:/Users/uqjrhode/work/networks_cons_plan/figures/simnest.csv")

#plot relationship
plot(Nested[,"SimNest"],Nested[,"Discr"],xlab="",ylab="",bty="l",xlim=c(0,1),ylim=c(0,18),yaxt="n",pch=22,col="red",bg="red")
title(xlab=expression(phi),cex.lab=1.2,line=2.5)
title(ylab="Discrepancy index",cex.lab=1.2,line=2.5)
axis(side = 1, at = c(seq(0,1,0.1)))
axis(side = 2, at = c(seq(0,18,3)))

#get correlation coefficient
cor.test(Nested[,"SimNest"],Nested[,"Discr"])






#VOI VARIATION WITH SHAPE AND NESTEDNESS - FIGURE 3

#GET DATA (Shapes: 1 = star, 2 = wheel, 3 = line, 4 = ring)

#STAR (1)

#set up matrices to hold results data
ResultsExp <- matrix(NA,ncol=11,nrow=2*2*11*5*11*11)
dimnames(ResultsExp) <- list(NULL,c("Shape","NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN","RND"))
ResultsSD <- matrix(NA,ncol=11,nrow=2*2*11*5*11*11)
dimnames(ResultsSD) <- list(NULL,c("Shape","NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN","RND"))

#create lists of variables
NumSpec <- c(20,60)
Occ <- c(0.1, 0.25)
Nested <- seq(0,1,0.1)
PiR <- seq(0,1,0.1)
PiD <- seq(0,1,0.1)

#set shape number
ResultsExp[,"Shape"] <- 1
ResultsSD[,"Shape"] <- 1

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
			# dimensions = [pir,pid,rep_id(BR replicate),method(1=CPSN,2=CP,3=SN,4=RND),number simulations] - note no "number simulations" in the expectation and SD files 
			MatFileExp <- readMat(paste("./results/ExpRc_starModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
			MatFileSD <- readMat(paste("./results/SDRc_starModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
									
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
ResultsExp1 <- ResultsExp
ResultsSD1 <- ResultsSD

#WHEEL (2)

#set up matrices to hold results
ResultsExp <- matrix(NA,ncol=11,nrow=2*2*11*5*11*11)
dimnames(ResultsExp) <- list(NULL,c("Shape","NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN","RND"))
ResultsSD <- matrix(NA,ncol=11,nrow=2*2*11*5*11*11)
dimnames(ResultsSD) <- list(NULL,c("Shape","NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN","RND"))

#create lists of variables
NumSpec <- c(20,60)
Occ <- c(0.1, 0.25)
Nested <- seq(0,1,0.1)
PiR <- seq(0,1,0.1)
PiD <- seq(0,1,0.1)

#set shape number
ResultsExp[,"Shape"] <- 2
ResultsSD[,"Shape"] <- 2

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
			MatFileExp <- readMat(paste("./results/ExpRc_wheelModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
			MatFileSD <- readMat(paste("./results/SDRc_wheelModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
									
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
ResultsExp2 <- ResultsExp
ResultsSD2 <- ResultsSD

#LINE (3)

#set up matrices to hold results
ResultsExp <- matrix(NA,ncol=11,nrow=2*2*11*5*11*11)
dimnames(ResultsExp) <- list(NULL,c("Shape","NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN","RND"))
ResultsSD <- matrix(NA,ncol=11,nrow=2*2*11*5*11*11)
dimnames(ResultsSD) <- list(NULL,c("Shape","NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN","RND"))

#create lists of variables
NumSpec <- c(20,60)
Occ <- c(0.1, 0.25)
Nested <- seq(0,1,0.1)
PiR <- seq(0,1,0.1)
PiD <- seq(0,1,0.1)

#set shape number
ResultsExp[,"Shape"] <- 3
ResultsSD[,"Shape"] <- 3

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
			MatFileExp <- readMat(paste("./results/ExpRc_lineModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
			MatFileSD <- readMat(paste("./results/SDRc_lineModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
									
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
ResultsExp3 <- ResultsExp
ResultsSD3 <- ResultsSD

#RING (4)

#set up matrices to hold results
ResultsExp <- matrix(NA,ncol=11,nrow=2*2*11*5*11*11)
dimnames(ResultsExp) <- list(NULL,c("Shape","NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN","RND"))
ResultsSD <- matrix(NA,ncol=11,nrow=2*2*11*5*11*11)
dimnames(ResultsSD) <- list(NULL,c("Shape","NumSp","Occ","Nested","BRRep","PiR","PiD","CPSN","CP","SN","RND"))

#create lists of variables
NumSpec <- c(20,60)
Occ <- c(0.1, 0.25)
Nested <- seq(0,1,0.1)
PiR <- seq(0,1,0.1)
PiD <- seq(0,1,0.1)

#set shape number
ResultsExp[,"Shape"] <- 4
ResultsSD[,"Shape"] <- 4

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
			MatFileExp <- readMat(paste("./results/ExpRc_islandModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
			MatFileSD <- readMat(paste("./results/SDRc_islandModel3nSites6NumSp",NumSpec[i],"Occ",Occ[j],"Nested",Nested[k],"pr0.2pd0.2.mat",sep=""))
									
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
ResultsExp4 <- ResultsExp
ResultsSD4 <- ResultsSD

#combine data
ResultsExp <- rbind(ResultsExp1,ResultsExp2,ResultsExp3,ResultsExp4)
ResultsSD <- rbind(ResultsSD1,ResultsSD2,ResultsSD3,ResultsSD4)


#aggregate data by shape and nestedness to find means across biodiversity replicates
ResultsExpAgg <- aggregate(ResultsExp[(ResultsExp[,"NumSp"] == 60) & (ResultsExp[,"Occ"] == 0.25),c("CPSN","CP")],by=list(Shape=ResultsExp[(ResultsExp[,"NumSp"] == 60) & (ResultsExp[,"Occ"] == 0.25),"Shape"],Nested=ResultsExp[(ResultsExp[,"NumSp"] == 60) & (ResultsExp[,"Occ"] == 0.25),"Nested"]),FUN=mean)

#few rare species
FRStar <- mean((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
FRStarMin <- min((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
FRStarMax <- max((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
FRWheel <- mean((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
FRWheelMin <- min((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
FRWheelMax <- max((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
FRLine <- mean((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
FRLineMin <- min((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
FRLineMax <- max((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
FRRing <- mean((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
FRRingMin <- min((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
FRRingMax <- max((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))

#few common species
FCStar <- mean((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
FCStarMin <- min((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
FCStarMax <- max((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
FCWheel <- mean((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
FCWheelMin <- min((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
FCWheelMax <- max((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
FCLine <- mean((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
FCLineMin <- min((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
FCLineMax <- max((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
FCRing <- mean((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))
FCRingMin <- min((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))
FCRingMax <- max((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))

#many rare species
MRStar <- mean((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
MRStarMin <- min((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
MRStarMax <- max((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
MRWheel <- mean((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
MRWheelMin <- min((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
MRWheelMax <- max((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
MRLine <- mean((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
MRLineMin <- min((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
MRLineMax <- max((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
MRRing <- mean((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
MRRingMin <- min((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
MRRingMax <- max((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))

#many common species
MCStar <- mean((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
MCStarMin <- min((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
MCStarMax <- max((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
MCWheel <- mean((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
MCWheelMin <- min((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
MCWheelMax <- max((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
MCLine <- mean((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))
MCLineMin <- min((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))
MCLineMax <- max((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==3) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))
MCRing <- mean((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))
MCRingMin <- min((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))
MCRingMax <- max((ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Shape"]==4) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))

#create figure
Star <- matrix(c(FRStar,FCStar,MRStar,MCStar) * 100,ncol=1,nrow=4,byrow=TRUE)
Wheel <- matrix(c(FRWheel,FCWheel,MRWheel,MCWheel) * 100,ncol=1,nrow=4,byrow=TRUE)
Line <- matrix(c(FRLine,FCLine,MRLine,MCLine) * 100,ncol=1,nrow=4,byrow=TRUE)
Ring <- matrix(c(FRRing,FCRing,MRRing,MCRing) * 100,ncol=1,nrow=4,byrow=TRUE)
PlotData <- rbind(Star,Wheel,Line,Ring)
SpeciesNames <- c("Few Rare Species","Few Common Species","Many Rare Species","Many Common Species","Few Rare Species","Few Common Species","Many Rare Species","Many Common Species","Few Rare Species","Few Common Species","Many Rare Species","Many Common Species","Few Rare Species","Few Common Species","Many Rare Species","Many Common Species")
dimnames(PlotData) <- list(NULL,c("Mean"))
PlotData <- as.data.frame(PlotData)
PlotData$Species <- SpeciesNames
PlotData$Species1 <- c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)
PlotData$Species2 <- with(PlotData, reorder(Species,Species1))
PlotData$Shape <- c("Star","Star","Star","Star","Wheel","Wheel","Wheel","Wheel","Line","Line","Line","Line","Ring","Ring","Ring","Ring")
PlotData$Shape1 <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)
PlotData$Shape2 <- with(PlotData, reorder(Shape,Shape1))
p <- ggplot(data = PlotData,aes(x = factor(Shape2), y = Mean,fill = factor(Species2)))
p <- p + geom_bar(stat = "identity",position = position_dodge(0.9)) + labs(x = "Shape", y = "Value of Social Network Information (% Improvement)") + scale_fill_discrete(name = "")
p <- p + theme(axis.text=element_text(size=12),legend.text=element_text(size=12),axis.title.y=element_text(size=14,margin=margin(0,10,0,0)),axis.title.x=element_text(size=14,margin=margin(10,0,0,0))) + ylim(0,20)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
p <- p + theme(legend.position = c(0.8, 0.9))
p

#VOI VARIATION WITH NESTEDNESS - FIGURE 4
#aggregare data by shape to find means across biodiversity replicates
ResultsExpAgg <- aggregate(ResultsExp[,c("CPSN","CP")],by=list(Nested=ResultsExp[,"Nested"],NumSp=ResultsExp[,"NumSp"],Occ=ResultsExp[,"Occ"]),FUN=mean)

#few rare species
FRN0_0 <- mean((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
FRN0_0Min <- min((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
FRN0_0Max <- max((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
FRN0_2 <- mean((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
FRN0_2Min <- min((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
FRN0_2Max <- max((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
FRN0_5 <- mean((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
FRN0_5Min <- min((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
FRN0_5Max <- max((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
FRN0_8 <- mean((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
FRN0_8Min <- min((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
FRN0_8Max <- max((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
FRN1_0 <- mean((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
FRN1_0Min <- min((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
FRN1_0Max <- max((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))

#few common species
FCN0_0 <- mean((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
FCN0_0Min <- min((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
FCN0_0Max <- max((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
FCN0_2 <- mean((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
FCN0_2Min <- min((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
FCN0_2Max <- max((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
FCN0_5 <- mean((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
FCN0_5Min <- min((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
FCN0_5Max <- max((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
FCN0_8 <- mean((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))
FCN0_8Min <- min((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))
FCN0_8Max <- max((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))
FCN1_0 <- mean((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))
FCN1_0Min <- min((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))
FCN1_0Max <- max((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==20) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))

#many rare species
MRN0_0 <- mean((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
MRN0_0Min <- min((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
MRN0_0Max <- max((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
MRN0_2 <- mean((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
MRN0_2Min <- min((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
MRN0_2Max <- max((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
MRN0_5 <- mean((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
MRN0_5Min <- min((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
MRN0_5Max <- max((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"])) 
MRN0_8 <- mean((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
MRN0_8Min <- min((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
MRN0_8Max <- max((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
MRN1_0 <- mean((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
MRN1_0Min <- min((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))
MRN1_0Max <- max((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.1),"CP"]))

#many common species
MCN0_0 <- mean((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
MCN0_0Min <- min((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
MCN0_0Max <- max((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
MCN0_2 <- mean((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
MCN0_2Min <- min((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
MCN0_2Max <- max((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.2) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"])) 
MCN0_5 <- mean((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))
MCN0_5Min <- min((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))
MCN0_5Max <- max((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.5) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))
MCN0_8 <- mean((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))
MCN0_8Min <- min((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))
MCN0_8Max <- max((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==0.8) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))
MCN1_0 <- mean((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))
MCN1_0Min <- min((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))
MCN1_0Max <- max((ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CPSN"] - ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]) /(ResultsExpAgg[(ResultsExpAgg[,"Nested"]==1) & (ResultsExpAgg[,"NumSp"]==60) & (ResultsExpAgg[,"Occ"]==0.25),"CP"]))

#create figure
N0_0 <- matrix(c(FRN0_0,FCN0_0,MRN0_0,MCN0_0) * 100,ncol=1,nrow=4,byrow=TRUE)
N0_2 <- matrix(c(FRN0_2,FCN0_2,MRN0_2,MCN0_2) * 100,ncol=1,nrow=4,byrow=TRUE)
N0_5 <- matrix(c(FRN0_5,FCN0_5,MRN0_5,MCN0_5) * 100,ncol=1,nrow=4,byrow=TRUE)
N0_8 <- matrix(c(FRN0_8,FCN0_8,MRN0_8,MCN0_8) * 100,ncol=1,nrow=4,byrow=TRUE)
N1_0 <- matrix(c(FRN1_0,FCN1_0,MRN1_0,MCN1_0) * 100,ncol=1,nrow=4,byrow=TRUE)
PlotData <- rbind(N0_0,N0_2,N0_5,N0_8,N1_0)
SpeciesNames <- c("Few Rare Species","Few Common Species","Many Rare Species","Many Common Species","Few Rare Species","Few Common Species","Many Rare Species","Many Common Species","Few Rare Species","Few Common Species","Many Rare Species","Many Common Species","Few Rare Species","Few Common Species","Many Rare Species","Many Common Species")
dimnames(PlotData) <- list(NULL,c("Mean"))
PlotData <- as.data.frame(PlotData)
PlotData$Species <- SpeciesNames
PlotData$Species1 <- c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)
PlotData$Species2 <- with(PlotData, reorder(Species,Species1))
PlotData$Shape <- c("Star","Star","Star","Star","Wheel","Wheel","Wheel","Wheel","Line","Line","Line","Line","Ring","Ring","Ring","Ring")
PlotData$Shape1 <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)
PlotData$Shape2 <- with(PlotData, reorder(Shape,Shape1))
p <- ggplot(data = PlotData,aes(x = factor(Shape2), y = Mean,fill = factor(Species2)))
p <- p + geom_bar(stat = "identity",position = position_dodge(0.9)) + labs(x = "Shape", y = "Value of Social Network Information (% Improvement)") + scale_fill_discrete(name = "")
p <- p + theme(axis.text=element_text(size=12),legend.text=element_text(size=12),axis.title.y=element_text(size=14,margin=margin(0,10,0,0)),axis.title.x=element_text(size=14,margin=margin(10,0,0,0))) + ylim(0,20)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
p <- p + theme(legend.position = c(0.8, 0.9))
p















#NETWORKS FIGURE S.2
barplot(c(0,0.29,0.64,1),names.arg=c("Ring","Line","Wheel","Star"),ylim=c(0,1.5),col="black")
title(xlab="Network Shape",cex.lab=1.2,line=2.5)
title(ylab="Closeness Centrality",cex.lab=1.2,line=2.5)

#GET DATA FOR REMAINING PLOTS 
load("M:/Users/uqjrhode/work/networks_cons_plan/Jags_Fits.RData")
Star <- summary(Jags.Fits[[1]])
Wheel <- summary(Jags.Fits[[2]])
Line <- summary(Jags.Fits[[3]])
Ring <- summary(Jags.Fits[[4]])


#EFFECT SIZE - PIR, PID - FIGURE 3
StarF1 <- Star[c("alpir","alpid"),c("Lower95","Upper95","Mean")]
WheelF1 <- Wheel[c("alpir","alpid"),c("Lower95","Upper95","Mean")]
LineF1 <- Line[c("alpir","alpid"),c("Lower95","Upper95","Mean")]
RingF1 <- Ring[c("alpir","alpid"),c("Lower95","Upper95","Mean")]
PlotDataF1 <- rbind(StarF1,WheelF1,LineF1,RingF1)
ParamNames <- c("Reserve Pr(influence)","Development Pr(influence)","Reserve Pr(influence)","Development Pr(influence)","Reserve Pr(influence)","Development Pr(influence)","Reserve Pr(influence)","Development Pr(influence)")
dimnames(PlotDataF1) <- list(NULL,c("Lower95","Upper95","Mean"))
PlotDataF1 <- as.data.frame(PlotDataF1)
PlotDataF1$Param <- ParamNames
PlotDataF1$Param1 <- c(1,2,1,2,1,2,1,2)
PlotDataF1$Param2 <- with(PlotDataF1, reorder(Param,Param1))
PlotDataF1$Shape <- c("Star","Star","Wheel","Wheel","Line","Line","Ring","Ring")
PlotDataF1$Shape1 <- c(4,4,3,3,2,2,1,1)
PlotDataF1$Shape2 <- with(PlotDataF1, reorder(Shape,Shape1))
limits <- aes(ymax = PlotDataF1$Upper95,ymin = PlotDataF1$Lower95)
p <- ggplot(data = PlotDataF1,aes(x = factor(Param), y = Mean,fill = factor(Shape2)))
p <- p + geom_bar(stat = "identity",position = position_dodge(0.9)) + geom_errorbar(limits, position = position_dodge(0.9),width = 0.25) + labs(x = "Parameter", y = "Effect on Benefit of Social Network Information") + scale_fill_discrete(name = "")
p <- p + theme(axis.text=element_text(size=12),legend.text=element_text(size=12),axis.title.y=element_text(size=15,face="bold",margin=margin(0,10,0,0)),axis.title.x=element_text(size=15,face="bold",margin=margin(10,0,0,0)))
p

#EFFECT SIZE - NUM SPECIES, RARITY, NESTEDNESS - FIGURE 4
StarF1 <- Star[c("sp","occ","nst"),c("Lower95","Upper95","Mean")]
WheelF1 <- Wheel[c("sp","occ","nst"),c("Lower95","Upper95","Mean")]
LineF1 <- Line[c("sp","occ","nst"),c("Lower95","Upper95","Mean")]
RingF1 <- Ring[c("sp","occ","nst"),c("Lower95","Upper95","Mean")]
PlotDataF1 <- rbind(StarF1,WheelF1,LineF1,RingF1)
ParamNames <- c("Number of Species","Species Rarity","Nestedness","Number of Species","Species Rarity","Nestedness","Number of Species","Species Rarity","Nestedness","Number of Species","Species Rarity","Nestedness")
dimnames(PlotDataF1) <- list(NULL,c("Lower95","Upper95","Mean"))
PlotDataF1 <- as.data.frame(PlotDataF1)
PlotDataF1$Param <- ParamNames
PlotDataF1$Param1 <- c(1,2,3,1,2,3,1,2,3,1,2,3)
PlotDataF1$Param2 <- with(PlotDataF1, reorder(Param,Param1))
PlotDataF1$Shape <- c("Star","Star","Star","Wheel","Wheel","Wheel","Line","Line","Line","Ring","Ring","Ring")
PlotDataF1$Shape1 <- c(4,4,4,3,3,3,2,2,2,1,1,1)
PlotDataF1$Shape2 <- with(PlotDataF1, reorder(Shape,Shape1))
limits <- aes(ymax = PlotDataF1$Upper95,ymin = PlotDataF1$Lower95)
p <- ggplot(data = PlotDataF1,aes(x = factor(Param2), y = Mean,fill = factor(Shape2)))
p <- p + geom_bar(stat = "identity",position = position_dodge(0.9)) + geom_errorbar(limits, position = position_dodge(0.9),width = 0.25) + labs(x = "Parameter", y = "Effect on Benefit of Social Network Information") + scale_fill_discrete(name = "")
p <- p + theme(axis.text=element_text(size=12),legend.text=element_text(size=12),axis.title.y=element_text(size=15,face="bold",margin=margin(0,10,0,0)),axis.title.x=element_text(size=15,face="bold",margin=margin(10,0,0,0)))
p

############################################



#EFFECT SIZE - INTERACTION WITH PIR - SPECIES, OCCUPANCY, NESTEDNESS - FIGURE 4
StarF1 <- Star[c("sppir","ocpir","nspir"),c("Lower95","Upper95","Mean")]
WheelF1 <- Wheel[c("sppir","ocpir","nspir"),c("Lower95","Upper95","Mean")]
LineF1 <- Line[c("sppir","ocpir","nspir"),c("Lower95","Upper95","Mean")]
RingF1 <- Ring[c("sppir","ocpir","nspir"),c("Lower95","Upper95","Mean")]
PlotDataF1 <- rbind(StarF1,WheelF1,LineF1,RingF1)
ParamNames <- c("Number of Species","Species Rarity","Nestedness","Number of Species","Species Rarity","Nestedness","Number of Species","Species Rarity","Nestedness","Number of Species","Species Rarity","Nestedness")
dimnames(PlotDataF1) <- list(NULL,c("Lower95","Upper95","Mean"))
PlotDataF1 <- as.data.frame(PlotDataF1)
PlotDataF1$Param <- ParamNames
PlotDataF1$Param1 <- c(1,2,3,1,2,3,1,2,3,1,2,3)
PlotDataF1$Param2 <- with(PlotDataF1, reorder(Param,Param1))
PlotDataF1$Shape <- c("Star","Star","Star","Wheel","Wheel","Wheel","Line","Line","Line","Ring","Ring","Ring")
PlotDataF1$Shape1 <- c(4,4,4,3,3,3,2,2,2,1,1,1)
PlotDataF1$Shape2 <- with(PlotDataF1, reorder(Shape,Shape1))
limits <- aes(ymax = PlotDataF1$Upper95,ymin = PlotDataF1$Lower95)
p <- ggplot(data = PlotDataF1,aes(x = factor(Param2), y = Mean,fill = factor(Shape2)))
p <- p + geom_bar(stat = "identity",position = position_dodge(0.9)) + geom_errorbar(limits, position = position_dodge(0.9),width = 0.25) + labs(x = "", y = "Interaction Effect with Pr(reserve influence)") + scale_fill_discrete(name = "")
p <- p + theme(axis.text=element_text(size=12,face="bold"),legend.text=element_text(size=12,face="bold"),axis.title.y=element_text(size=15,face="bold"))
p

#EFFECT SIZE - INTERACTION WITH PID - SPECIES, OCCUPANCY, NESTEDNESS - FIGURE 5
StarF1 <- Star[c("sppid","ocpid","nspid"),c("Lower95","Upper95","Mean")]
WheelF1 <- Wheel[c("sppid","ocpid","nspid"),c("Lower95","Upper95","Mean")]
LineF1 <- Line[c("sppid","ocpid","nspid"),c("Lower95","Upper95","Mean")]
RingF1 <- Ring[c("sppid","ocpid","nspid"),c("Lower95","Upper95","Mean")]
PlotDataF1 <- rbind(StarF1,WheelF1,LineF1,RingF1)
ParamNames <- c("Number of Species","Rarity","Nestedness","Number of Species","Rarity","Nestedness","Number of Species","Rarity","Nestedness","Number of Species","Rarity","Nestedness")
dimnames(PlotDataF1) <- list(NULL,c("Lower95","Upper95","Mean"))
PlotDataF1 <- as.data.frame(PlotDataF1)
PlotDataF1$Param <- ParamNames
PlotDataF1$Param1 <- c(1,2,3,1,2,3,1,2,3,1,2,3)
PlotDataF1$Param2 <- with(PlotDataF1, reorder(Param,Param1))
PlotDataF1$Shape <- c("Star","Star","Star","Wheel","Wheel","Wheel","Line","Line","Line","Ring","Ring","Ring")
PlotDataF1$Shape1 <- c(4,4,4,3,3,3,2,2,2,1,1,1)
PlotDataF1$Shape2 <- with(PlotDataF1, reorder(Shape,Shape1))
limits <- aes(ymax = PlotDataF1$Upper95,ymin = PlotDataF1$Lower95)
p <- ggplot(data = PlotDataF1,aes(x = factor(Param2), y = Mean,fill = factor(Shape2)))
p <- p + geom_bar(stat = "identity",position = position_dodge(0.9)) + geom_errorbar(limits, position = position_dodge(0.9),width = 0.25) + labs(x = "", y = "Interaction Effect with Pr(development influence)") + scale_fill_discrete(name = "")
p <- p + theme(axis.text=element_text(size=12,face="bold"),legend.text=element_text(size=12,face="bold"),axis.title.y=element_text(size=15,face="bold"))
p
