get.jags1 <- function(Data)
{
	#get initial values
	inits1 <- list(alpha=runif(1,-10,10),nst=runif(1,-10,10),alpir=runif(1,-10,10),nspir=runif(1,0,10),spir=runif(1,0,10),alpid=runif(1,-10,10),nspid=runif(1,0,10),spid=runif(1,0,10),sd=runif(1,0,10),TrueSC=runif(Data$N,-10,10))
	inits2 <- list(alpha=runif(1,-10,10),nst=runif(1,-10,10),alpir=runif(1,-10,10),nspir=runif(1,0,10),spir=runif(1,0,10),alpid=runif(1,-10,10),nspid=runif(1,0,10),spid=runif(1,0,10),sd=runif(1,0,10),TrueSC=runif(Data$N,-10,10))
	inits3 <- list(alpha=runif(1,-10,10),nst=runif(1,-10,10),alpir=runif(1,-10,10),nspir=runif(1,0,10),spir=runif(1,0,10),alpid=runif(1,-10,10),nspid=runif(1,0,10),spid=runif(1,0,10),sd=runif(1,0,10),TrueSC=runif(Data$N,-10,10))

	cl <- makeCluster(2)

	fit <- run.jags(model="jags_model1.txt",monitor=c("alpha","nst","sa","alpir","nspir","spir","alpid","nspid","spid","sd"),data=Data,n.chains=3,inits=list(inits1,inits2,inits3),burnin=20000,adapt=1000,sample=20000,jags="C:/Program Files/JAGS/JAGS-4.2.0/x64/bin/jags-terminal.exe",method="rjparallel",cl=cl)

	stopCluster(cl)

	return(fit)
}

# plotting JAGS output
plot.fit <- function(fit, cols=c("red", "blue", "green"))
{
  M <- as.mcmc(fit)
  ns <- length(M[,1])/3
  np <- length(M[1,])
  cn <- colnames(M)
  par(mfrow=c(np,2),mar=c(1,4,1,1))
  for (i in 1:np)
  {
    plot(c(1, ns), c(min(M[,i]), max(M[,i])), type="n", xlab="", ylab=cn[i])
    lines(c(1:ns), M[1:ns,i], col=cols[1])
    lines(c(1:ns), M[(ns+1):(2*ns),i], col=cols[2])
    lines(c(1:ns), M[(2*ns+1):(3*ns),i], col=cols[3])
    d1 <- density(M[1:ns,i])
    d2 <- density(M[(ns+1):(2*ns),i])
    d3 <- density(M[(2*ns+1):(3*ns),i])

    plot(c(min(c(d1$x, d2$x, d3$x)), max(c(d1$x, d2$x, d3$x))), c(0, max(c(d1$y, d2$y, d3$y))), type="n", main="", xlab="", ylab="")
    lines(d1, col=cols[1])
    lines(d2, col=cols[2])
    lines(d3, col=cols[3])
    rug(sample(M[1:ns,i], 100), col=cols[1])
    rug(sample(M[(ns+1):(2*ns),i], 100), col=cols[2])
    rug(sample(M[(2*ns+1):(3*ns),i], 100), col=cols[3])
  }
}
