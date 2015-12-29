rm(list=ls())

I=50    #Number of participants
J=100   #Number of items
K=4     #Number of response categories

iters=10000 #number of MCMC iterations
burnin=1000 #number of MCMC iterations to consider burnin


#Read in chains
len=5+2*(I+J+1)+I*(K-1)+1
zz <- file('an2.chn', "rb")
chains=readBin(zz, double(),n=len*iters)
close(zz)
dim(chains)=c(len,iters)
chains=t(chains)
#End readin


#define the columns of the chains
is.sig2s=1
is.sig2n=2
is.ds=2+1
is.alphas=3:(I+2)+1
is.betas=(I+3):(I+J+2)+1
is.dn=I+J+3+1
is.alphan=(I+J+4):(2*I+J+3)+1
is.betan=(2*I+J+4):(2*I+2*J+3)+1
is.sdeff=(2*I+2*J+4):(2*I+2*J+7)+1
is.crits=(2*I+2*J+9):(len)

cols=c("sig^2_s","sig^2_n","mu^s",paste("alpha^s_",1:I,sep=""),paste("beta^s_",1:J,sep=""),"mu^n",paste("alpha^n_",1:I,sep=""),paste("beta^n_",1:J,sep=""),
"sig2_alpha^n","sig2_alpha^s","sig2_beta^n","sig2_beta^s", paste("c_",as.vector(t(outer(1:I,1:(K-1),paste,sep="."))),sep=""))
colnames(chains)=cols


#define some new values
eta=sqrt(chains[,is.sig2s])/sqrt(chains[,is.sig2n])
dp=chains[,is.ds]/sqrt(chains[,is.sig2n])-chains[,is.dn]/sqrt(chains[,is.sig2n])
chains=cbind(chains,eta,dp)
is.eta=len+1
is.dp=len+2


#compute posterior means and credible intervals
estimates=colMeans(chains[(burnin+1):iters,])
dim(estimates)=c(1,length(estimates))

qs=apply(chains[(burnin+1):iters,],c(2),quantile,p=c(.025,.975),na.rm=T)
results=rbind(estimates,qs)
dimnames(results)[[1]]=c("Mean","2.5%","97.5%")
results=t(results)

#results contains posterior means and 95% central credible interval

#write chains to pdf, to check them

lo=c(1,2,1,3)
dim(lo)=c(2,2)

pdf("chains.pdf",width=10,height=8,one=T)
par(las=1,cex=1.4)
layout(lo)

plot(chains[,is.eta],t='l',xlab="Iteration",ylab="Eta",main="MCMC chain")
abline(h=results[is.eta,1],col="blue")
abline(h=results[is.eta,2:3],col="blue",lty=2)
abline(v=burnin,col="red")
plot(density(chains[(burnin+1):iters,is.eta]),xlab="Eta",main="Estimated Marginal Posterior")
abline(v=results[is.eta,1],col="blue")
abline(v=results[is.eta,2:3],col="blue",lty=2)
acf(chains[,is.eta])

plot(chains[,is.dp],t='l',xlab="Iteration",ylab="Overall d-prime",main="MCMC chain")
abline(h=results[is.dp,1],col="blue")
abline(h=results[is.dp,2:3],col="blue",lty=2)
abline(v=burnin,col="red")
plot(density(chains[(burnin+1):iters,is.dp]),xlab="Overall d-prime",main="Estimated Marginal Posterior")
abline(v=results[is.dp,1],col="blue")
abline(v=results[is.dp,2:3],col="blue",lty=2)
acf(chains[,is.dp])

plot(chains[,is.ds],t='l',xlab="Iteration",ylab="mu_s",main="MCMC chain")
abline(h=results[is.ds,1],col="blue")
abline(h=results[is.ds,2:3],col="blue",lty=2)
abline(v=burnin,col="red")
plot(density(chains[(burnin+1):iters,is.ds]),xlab="mu_s",main="Estimated Marginal Posterior")
abline(v=results[is.ds,1],col="blue")
abline(v=results[is.ds,2:3],col="blue",lty=2)
acf(chains[,is.ds])

plot(chains[,is.dn],t='l',xlab="Iteration",ylab="mu_n",main="MCMC chain")
abline(h=results[is.dn,1],col="blue")
abline(h=results[is.dn,2:3],col="blue",lty=2)
abline(v=burnin,col="red")
plot(density(chains[(burnin+1):iters,is.dn]),xlab="mu_n",main="Estimated Marginal Posterior")
abline(v=results[is.dn,1],col="blue")
abline(v=results[is.dn,2:3],col="blue",lty=2)
acf(chains[,is.dn])

dev.off()


#Examine item effects

pdf("items.pdf",width=10,height=8,one=T)
par(las=1,cex=1.4)

plot(results[is.betan,1],ylim=range(results[is.betan,]),pch=19,xlab="Item number",ylab="Item effect",main="Item effect beta_n")
segments(1:J,results[is.betan,2],1:J,results[is.betan,3])

plot(results[is.betas,1],ylim=range(results[is.betas,]),pch=19,xlab="Item number",ylab="Item effect",main="Item effect beta_s")
segments(1:J,results[is.betas,2],1:J,results[is.betas,3])

plot(results[is.betan,1],results[is.betas,1],xlab="beta_n",ylab="beta_s",pch=19)
abline(lm(results[is.betas,1]~results[is.betan,1]))

dev.off()

