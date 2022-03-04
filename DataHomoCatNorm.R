## data generation: true dist=norm, homogeneous ##
library(MASS)
lp=1000
ssize=200
NS=50
rho1=c(0,.1,.3,.5)
tau1=c(1,2) # tau-equivalency is met (1) or not (2);
dist1=1:4 # distribution of categorical items are all symmetric (1), all skewed in the same direction (2), mixed symmetric+skewednot (3), mixed left skewed and right skewed (4);
cat1=c(2,5) # number of categories for items;
I1=c(4,12)
cond=expand.grid(rho1,tau1,dist1,cat1,I1)
cond.idx=as.numeric(commandArgs(trailingOnly=TRUE))
rho=cond[cond.idx,1]
tau=cond[cond.idx,2]
dist=cond[cond.idx,3]
cat=cond[cond.idx,4]
I=cond[cond.idx,5]
## means and covariance matrix of the bi-variate variables in each study;
#mu=c(0,0) # true score means;
#varx=1 # true score variances;
#vary=1
#sigma=matrix(c(varx,rho*sqrt(varx*vary),rho*sqrt(varx*vary),vary),2,2) # fixed correlation in the population;
if (tau==1) {
  faloading=rep(.6,I)
  msmterrvarx=msmterrvary=rep(.55,I)
} 
if (tau==2) {
  faloading=rep(c(.3,.9),each=I/2)
  msmterrvarx=msmterrvary=rep(c(.91,.19),each=I/2)
}
# measurement error covariance matrix;
mumsmterrx <- mumsmterry <- rep(0,I)
sigmamsmterrx <- sigmamsmterry <- diag(msmterrvarx) # msmt error covariance matrix;
datadir=paste0('/gpfs/home/qzhang4/metacorr/Homogeneous/DataHomoCatNorm/cond',cond.idx+1000,'/')
#system(paste0('rm -r ',datadir))
dir.create(datadir)
# generate sample sizes;
nn=1
N=array(NA,NS) # N: sample sizes across studies;
while (nn <= NS) {
  X=rchisq(1,3)
  Ni=floor((X-3)/sqrt(2*3)*ssize/2+ssize)
  if (Ni<=3) {next}
  else {
    N[nn]=Ni
    nn=nn+1}
}
# replicate data sets, each with NS effect sizes;
flag=0 # flag problematic generated data sets for each condition;
for (i in 1:lp) {
  ro=ra=alphax=alphay=NULL
  nn=1
  while (nn <= NS) {
    # initialize matrices for measurement errors and items;
    msmterrx=msmterry=ix=iy=matrix(NA,N[nn],I)
    # generate the true scores:normal;
    xt=rnorm(N[nn])
    w2=rnorm(N[nn])
    yt=rho*xt+sqrt(1-rho^2)*w2
    msmterrx[,1:I] <- mvrnorm(N[nn],mumsmterrx,sigmamsmterrx) # generate msmt error for x;
    msmterry[,1:I] <- mvrnorm(N[nn],mumsmterry,sigmamsmterry) # generate msmt error for y;
    ix[,1:I]=xt%*%t(as.matrix(faloading))+msmterrx[,1:I] # item scores for x;
    iy[,1:I]=yt%*%t(as.matrix(faloading))+msmterry[,1:I] # item scores for y;
    if (cat==2) {
      if (dist==1) {
        for (j in 1:I) {
          # obtain the population var of observed variables
          popvar <- (faloading[j])^2+msmterrvarx[j] 
          # obtain quantiles at the population level given the normal distribution
          cutoff <- qnorm(.5,0,sqrt(popvar))
          xidx1 <- which(ix[,j]<cutoff)
          xidx2 <- which(ix[,j]>=cutoff)
          yidx1 <- which(iy[,j]<cutoff)
          yidx2 <- which(iy[,j]>=cutoff)
          ix[xidx1,j] <- 1
          ix[xidx2,j] <- 2
          iy[yidx1,j] <- 1
          iy[yidx2,j] <- 2
        }
      }
      if (dist==2) {
        for (j in 1:I) {
          # obtain the population var of observed variables
          popvar <- (faloading[j])^2+msmterrvarx[j] 
          # obtain quantiles at the population level given the normal distribution
          cutoff <- qnorm(.8,0,sqrt(popvar))
          xidx1 <- which(ix[,j]<cutoff)
          xidx2 <- which(ix[,j]>=cutoff)
          yidx1 <- which(iy[,j]<cutoff)
          yidx2 <- which(iy[,j]>=cutoff)
          ix[xidx1,j] <- 1
          ix[xidx2,j] <- 2
          iy[yidx1,j] <- 1
          iy[yidx2,j] <- 2
        }
      }
      if (dist==3) {
        for (j in 1:I/2) {
          # obtain the population var of observed variables
          popvar <- (faloading[j])^2+msmterrvarx[j] 
          # obtain quantiles at the population level given the normal distribution
          cutoff <- qnorm(.5,0,sqrt(popvar))
          xidx1 <- which(ix[,j]<cutoff)
          xidx2 <- which(ix[,j]>=cutoff)
          yidx1 <- which(iy[,j]<cutoff)
          yidx2 <- which(iy[,j]>=cutoff)
          ix[xidx1,j] <- 1
          ix[xidx2,j] <- 2
          iy[yidx1,j] <- 1
          iy[yidx2,j] <- 2
        }
        for (j in (I/2+1):I) {
          # obtain the population var of observed variables
          popvar <- (faloading[j])^2+msmterrvarx[j] 
          # obtain quantiles at the population level given the normal distribution
          cutoff <- qnorm(.8,0,sqrt(popvar))
          xidx1 <- which(ix[,j]<cutoff)
          xidx2 <- which(ix[,j]>=cutoff)
          yidx1 <- which(iy[,j]<cutoff)
          yidx2 <- which(iy[,j]>=cutoff)
          ix[xidx1,j] <- 1
          ix[xidx2,j] <- 2
          iy[yidx1,j] <- 1
          iy[yidx2,j] <- 2
        }
      }
      if (dist==4) {
        for (j in 1:I/2) {
          # obtain the population var of observed variables
          popvar <- (faloading[j])^2+msmterrvarx[j] 
          # obtain quantiles at the population level given the normal distribution
          cutoff <- qnorm(.8,0,sqrt(popvar))
          xidx1 <- which(ix[,j]<cutoff)
          xidx2 <- which(ix[,j]>=cutoff)
          yidx1 <- which(iy[,j]<cutoff)
          yidx2 <- which(iy[,j]>=cutoff)
          ix[xidx1,j] <- 1
          ix[xidx2,j] <- 2
          iy[yidx1,j] <- 1
          iy[yidx2,j] <- 2
        }
        for (j in (I/2+1):I) {
          # obtain the population var of observed variables
          popvar <- (faloading[j])^2+msmterrvarx[j] 
          # obtain quantiles at the population level given the normal distribution
          cutoff <- qnorm(.2,0,sqrt(popvar))
          xidx1 <- which(ix[,j]<cutoff)
          xidx2 <- which(ix[,j]>=cutoff)
          yidx1 <- which(iy[,j]<cutoff)
          yidx2 <- which(iy[,j]>=cutoff)
          ix[xidx1,j] <- 1
          ix[xidx2,j] <- 2
          iy[yidx1,j] <- 1
          iy[yidx2,j] <- 2
        }
      }
    }
    if (cat==5) {
      if (dist==1) {
        for (j in 1:I) {
          # obtain the population var of observed variables
          popvar <- (faloading[j])^2+msmterrvarx[j] 
          # obtain quantiles at the population level given the normal distribution
          cutoffs <- qnorm(c(.1,.3,.7,.9),0,sqrt(popvar))
          xidx1 <- which(ix[,j]<cutoffs[1])
          xidx2 <- which((ix[,j]>=cutoffs[1])&(ix[,j]<cutoffs[2]))
          xidx3 <- which((ix[,j]>=cutoffs[2])&(ix[,j]<cutoffs[3]))
          xidx4 <- which((ix[,j]>=cutoffs[3])&(ix[,j]<cutoffs[4]))
          xidx5 <- which(ix[,j]>=cutoffs[4])
          yidx1 <- which(iy[,j]<cutoffs[1])
          yidx2 <- which((iy[,j]>=cutoffs[1])&(iy[,j]<cutoffs[2]))
          yidx3 <- which((iy[,j]>=cutoffs[2])&(iy[,j]<cutoffs[3]))
          yidx4 <- which((iy[,j]>=cutoffs[3])&(iy[,j]<cutoffs[4]))
          yidx5 <- which(iy[,j]>=cutoffs[4])
          ix[xidx1,j] <- 1
          ix[xidx2,j] <- 2
          ix[xidx3,j] <- 3
          ix[xidx4,j] <- 4
          ix[xidx5,j] <- 5
          iy[yidx1,j] <- 1
          iy[yidx2,j] <- 2
          iy[yidx3,j] <- 3
          iy[yidx4,j] <- 4
          iy[yidx5,j] <- 5
        }
      }
      if (dist==2) {
        for (j in 1:I) {
          # obtain the population var of observed variables
          popvar <- (faloading[j])^2+msmterrvarx[j] 
          # obtain quantiles at the population level given the normal distribution
          cutoffs <- qnorm(c(.4,.65,.85,.95),0,sqrt(popvar))
          xidx1 <- which(ix[,j]<cutoffs[1])
          xidx2 <- which((ix[,j]>=cutoffs[1])&(ix[,j]<cutoffs[2]))
          xidx3 <- which((ix[,j]>=cutoffs[2])&(ix[,j]<cutoffs[3]))
          xidx4 <- which((ix[,j]>=cutoffs[3])&(ix[,j]<cutoffs[4]))
          xidx5 <- which(ix[,j]>=cutoffs[4])
          yidx1 <- which(iy[,j]<cutoffs[1])
          yidx2 <- which((iy[,j]>=cutoffs[1])&(iy[,j]<cutoffs[2]))
          yidx3 <- which((iy[,j]>=cutoffs[2])&(iy[,j]<cutoffs[3]))
          yidx4 <- which((iy[,j]>=cutoffs[3])&(iy[,j]<cutoffs[4]))
          yidx5 <- which(iy[,j]>=cutoffs[4])
          ix[xidx1,j] <- 1
          ix[xidx2,j] <- 2
          ix[xidx3,j] <- 3
          ix[xidx4,j] <- 4
          ix[xidx5,j] <- 5
          iy[yidx1,j] <- 1
          iy[yidx2,j] <- 2
          iy[yidx3,j] <- 3
          iy[yidx4,j] <- 4
          iy[yidx5,j] <- 5
        }
      }
      if (dist==3) {
        for (j in 1:I/2) {
          popvar <- (faloading[j])^2+msmterrvarx[j] 
          # obtain quantiles at the population level given the normal distribution
          cutoffs <- qnorm(c(.1,.3,.7,.9),0,sqrt(popvar))
          xidx1 <- which(ix[,j]<cutoffs[1])
          xidx2 <- which((ix[,j]>=cutoffs[1])&(ix[,j]<cutoffs[2]))
          xidx3 <- which((ix[,j]>=cutoffs[2])&(ix[,j]<cutoffs[3]))
          xidx4 <- which((ix[,j]>=cutoffs[3])&(ix[,j]<cutoffs[4]))
          xidx5 <- which(ix[,j]>=cutoffs[4])
          yidx1 <- which(iy[,j]<cutoffs[1])
          yidx2 <- which((iy[,j]>=cutoffs[1])&(iy[,j]<cutoffs[2]))
          yidx3 <- which((iy[,j]>=cutoffs[2])&(iy[,j]<cutoffs[3]))
          yidx4 <- which((iy[,j]>=cutoffs[3])&(iy[,j]<cutoffs[4]))
          yidx5 <- which(iy[,j]>=cutoffs[4])
          ix[xidx1,j] <- 1
          ix[xidx2,j] <- 2
          ix[xidx3,j] <- 3
          ix[xidx4,j] <- 4
          ix[xidx5,j] <- 5
          iy[yidx1,j] <- 1
          iy[yidx2,j] <- 2
          iy[yidx3,j] <- 3
          iy[yidx4,j] <- 4
          iy[yidx5,j] <- 5
        }
        for (j in (I/2+1):I) {
          popvar <- (faloading[j])^2+msmterrvarx[j] 
          # obtain quantiles at the population level given the normal distribution
          cutoffs <- qnorm(c(.4,.65,.85,.95),0,sqrt(popvar))
          xidx1 <- which(ix[,j]<cutoffs[1])
          xidx2 <- which((ix[,j]>=cutoffs[1])&(ix[,j]<cutoffs[2]))
          xidx3 <- which((ix[,j]>=cutoffs[2])&(ix[,j]<cutoffs[3]))
          xidx4 <- which((ix[,j]>=cutoffs[3])&(ix[,j]<cutoffs[4]))
          xidx5 <- which(ix[,j]>=cutoffs[4])
          yidx1 <- which(iy[,j]<cutoffs[1])
          yidx2 <- which((iy[,j]>=cutoffs[1])&(iy[,j]<cutoffs[2]))
          yidx3 <- which((iy[,j]>=cutoffs[2])&(iy[,j]<cutoffs[3]))
          yidx4 <- which((iy[,j]>=cutoffs[3])&(iy[,j]<cutoffs[4]))
          yidx5 <- which(iy[,j]>=cutoffs[4])          
          ix[xidx1,j] <- 1
          ix[xidx2,j] <- 2
          ix[xidx3,j] <- 3
          ix[xidx4,j] <- 4
          ix[xidx5,j] <- 5
          iy[yidx1,j] <- 1
          iy[yidx2,j] <- 2
          iy[yidx3,j] <- 3
          iy[yidx4,j] <- 4
          iy[yidx5,j] <- 5
        }
      }
      if (dist==4) {
        for (j in 1:I/2) {
          # obtain the population var of observed variables
          popvar <- (faloading[j])^2+msmterrvarx[j] 
          # obtain quantiles at the population level given the normal distribution
          cutoffs <- qnorm(c(.4,.65,.85,.95),0,sqrt(popvar))
          xidx1 <- which(ix[,j]<cutoffs[1])
          xidx2 <- which((ix[,j]>=cutoffs[1])&(ix[,j]<cutoffs[2]))
          xidx3 <- which((ix[,j]>=cutoffs[2])&(ix[,j]<cutoffs[3]))
          xidx4 <- which((ix[,j]>=cutoffs[3])&(ix[,j]<cutoffs[4]))
          xidx5 <- which(ix[,j]>=cutoffs[4])
          yidx1 <- which(iy[,j]<cutoffs[1])
          yidx2 <- which((iy[,j]>=cutoffs[1])&(iy[,j]<cutoffs[2]))
          yidx3 <- which((iy[,j]>=cutoffs[2])&(iy[,j]<cutoffs[3]))
          yidx4 <- which((iy[,j]>=cutoffs[3])&(iy[,j]<cutoffs[4]))
          yidx5 <- which(iy[,j]>=cutoffs[4])
          ix[xidx1,j] <- 1
          ix[xidx2,j] <- 2
          ix[xidx3,j] <- 3
          ix[xidx4,j] <- 4
          ix[xidx5,j] <- 5
          iy[yidx1,j] <- 1
          iy[yidx2,j] <- 2
          iy[yidx3,j] <- 3
          iy[yidx4,j] <- 4
          iy[yidx5,j] <- 5
        }
        for (j in (I/2+1):I) {
          # obtain the population var of observed variables
          popvar <- (faloading[j])^2+msmterrvarx[j] 
          # obtain quantiles at the population level given the normal distribution
          cutoffs <- qnorm(c(.05,.15,.35,.6),0,sqrt(popvar))
          xidx1 <- which(ix[,j]<cutoffs[1])
          xidx2 <- which((ix[,j]>=cutoffs[1])&(ix[,j]<cutoffs[2]))
          xidx3 <- which((ix[,j]>=cutoffs[2])&(ix[,j]<cutoffs[3]))
          xidx4 <- which((ix[,j]>=cutoffs[3])&(ix[,j]<cutoffs[4]))
          xidx5 <- which(ix[,j]>=cutoffs[4])
          yidx1 <- which(iy[,j]<cutoffs[1])
          yidx2 <- which((iy[,j]>=cutoffs[1])&(iy[,j]<cutoffs[2]))
          yidx3 <- which((iy[,j]>=cutoffs[2])&(iy[,j]<cutoffs[3]))
          yidx4 <- which((iy[,j]>=cutoffs[3])&(iy[,j]<cutoffs[4]))
          yidx5 <- which(iy[,j]>=cutoffs[4])
          ix[xidx1,j] <- 1
          ix[xidx2,j] <- 2
          ix[xidx3,j] <- 3
          ix[xidx4,j] <- 4
          ix[xidx5,j] <- 5
          iy[yidx1,j] <- 1
          iy[yidx2,j] <- 2
          iy[yidx3,j] <- 3
          iy[yidx4,j] <- 4
          iy[yidx5,j] <- 5
        }
      }
    }
    xo=apply(ix,1,sum)
    yo=apply(iy,1,sum)
    covx=cov(ix)
    covy=cov(iy)
    varsx=diag(covx)
    varsy=diag(covy)
    alphax1=I/(I-1)*(var(xo)-sum(varsx))/var(xo)
    alphay1=I/(I-1)*(var(yo)-sum(varsy))/var(yo)
    if ((sd(xo)!=0)&(sd(yo)!=0)&(alphax1<1)&(alphax1>0)&(alphay1<1)&(alphay1>0))
    {cor1=cor(xo,yo)
    ra1=cor1/sqrt(alphax1*alphay1)}
    else {
      flag=flag+1
      next}
    if (abs(ra1)<1)
    {ro=c(ro,cor1)
    ra=c(ra,ra1)
    alphax=c(alphax,alphax1)
    alphay=c(alphay,alphay1)
    nn=nn+1
    } else {
      flag=flag+1
      next}
  }
  dset=cbind(ro,ra,alphax,alphay,N)
  dataname=paste0(datadir,'data',i+1000,'.txt')
  write.table(dset,dataname,col.names=F,row.names=F)
}
flagname=paste0(datadir,'flag.txt')
write.table(flag,flagname,col.names=F,row.names=F)
