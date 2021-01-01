## data generation: true dist=skewed, heterogeneous ##
library(MASS)
lp=1000
ssize=200
NS=50
zrho1=c(0,.1028345,.3170305,.5618704)
tau1=c(1,2) # tau-equivalency is met (1) or not (2);
dist1=1:4 # distribution of categorical items are all symmetric (1), all skewed in the same direction (2), mixed symmetric+skewednot (3), mixed left skewed and right skewed (4);
cat1=c(2,5) # number of categories for items;
I1=c(4,12)
cond=expand.grid(zrho1,tau1,dist1,cat1,I1)
cond.idx=as.numeric(commandArgs(trailingOnly=TRUE))
zrho=cond[cond.idx,1]
tau=cond[cond.idx,2]
dist=cond[cond.idx,3]
cat=cond[cond.idx,4]
I=cond[cond.idx,5]
## means and covariance matrix of the bi-variate variables in each study;
mu=c(0,0) # true score means;
varx=1 # true score variances;
vary=1
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
datadir=paste0('/gpfs/home/qzhang4/metacorr/Heterogeneous/DataHeteCatSkew/cond',cond.idx+1000,'/')
dir.create(datadir)
### data generation ###
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
  zr=rnorm(NS,zrho,.16)
  rhos=(exp(2*zr)-1)/(exp(2*zr)+1)
  while (nn <= NS) {
    # initialize matrices for measurement errors and items;
    msmterrx=msmterry=ix=iy=matrix(NA,N[nn],I)
    z1u=rchisq(N[nn],1)
    z1=(z1u-mean(z1u))/sd(z1u)
    z2u=rchisq(N[nn],1)
    z2=(z2u-mean(z2u))/sd(z2u)
    xt=z1
    yt=rhos[nn]*z1+sqrt(1-rhos[nn]^2)*z2
    msmterrx[,1:I] <- mvrnorm(N[nn],mumsmterrx,sigmamsmterrx) # generate msmt error for x;
    msmterry[,1:I] <- mvrnorm(N[nn],mumsmterry,sigmamsmterry) # generate msmt error for y;
    ix[,1:I]=xt%*%t(as.matrix(faloading))+msmterrx[,1:I] # item scores for x;
    iy[,1:I]=yt%*%t(as.matrix(faloading))+msmterry[,1:I] # item scores for y;
    if (cat==2) {
      if (dist==1) {
        for (j in 1:I) {
          rankx <- rank(ix[,j])
          ranky <- rank(iy[,j])
          xidx1 <- which(rankx<(N[nn]*.50))
          xidx2 <- which(rankx>=(N[nn]*.50))
          yidx1 <- which(ranky<(N[nn]*.50))
          yidx2 <- which(ranky>=(N[nn]*.50))
          ix[xidx1,j] <- 1
          ix[xidx2,j] <- 2
          iy[yidx1,j] <- 1
          iy[yidx2,j] <- 2
        }
      }
      if (dist==2) {
        for (j in 1:I) {
          rankx <- rank(ix[,j])
          ranky <- rank(iy[,j])
          xidx1 <- which(rankx<(N[nn]*.80))
          xidx2 <- which(rankx>=(N[nn]*.80))
          yidx1 <- which(ranky<(N[nn]*.80))
          yidx2 <- which(ranky>=(N[nn]*.80))
          ix[xidx1,j] <- 1
          ix[xidx2,j] <- 2
          iy[yidx1,j] <- 1
          iy[yidx2,j] <- 2
        }
      }
      if (dist==3) {
        for (j in 1:I/2) {
          rankx <- rank(ix[,j])
          ranky <- rank(iy[,j])
          xidx1 <- which(rankx<(N[nn]*.50))
          xidx2 <- which(rankx>=(N[nn]*.50))
          yidx1 <- which(ranky<(N[nn]*.50))
          yidx2 <- which(ranky>=(N[nn]*.50))
          ix[xidx1,j] <- 1
          ix[xidx2,j] <- 2
          iy[yidx1,j] <- 1
          iy[yidx2,j] <- 2
        }
        for (j in (I/2+1):I) {
          rankx <- rank(ix[,j])
          ranky <- rank(iy[,j])
          xidx1 <- which(rankx<(N[nn]*.80))
          xidx2 <- which(rankx>=(N[nn]*.80))
          yidx1 <- which(ranky<(N[nn]*.80))
          yidx2 <- which(ranky>=(N[nn]*.80))
          ix[xidx1,j] <- 1
          ix[xidx2,j] <- 2
          iy[yidx1,j] <- 1
          iy[yidx2,j] <- 2
        }
      }
      if (dist==4) {
        for (j in 1:I/2) {
          rankx <- rank(ix[,j])
          ranky <- rank(iy[,j])
          xidx1 <- which(rankx<(N[nn]*.80))
          xidx2 <- which(rankx>=(N[nn]*.80))
          yidx1 <- which(ranky<(N[nn]*.80))
          yidx2 <- which(ranky>=(N[nn]*.80))
          ix[xidx1,j] <- 1
          ix[xidx2,j] <- 2
          iy[yidx1,j] <- 1
          iy[yidx2,j] <- 2
        }
        for (j in (I/2+1):I) {
          rankx <- rank(ix[,j])
          ranky <- rank(iy[,j])
          xidx1 <- which(rankx<(N[nn]*.20))
          xidx2 <- which(rankx>=(N[nn]*.20))
          yidx1 <- which(ranky<(N[nn]*.20))
          yidx2 <- which(ranky>=(N[nn]*.20))
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
          rankx <- rank(ix[,j])
          ranky <- rank(iy[,j])
          xidx1 <- which(rankx<(N[nn]*.10))
          xidx2 <- which((rankx>=N[nn]*.10)&(rankx<N[nn]*.3))
          xidx3 <- which((rankx>=N[nn]*.3)&(rankx<N[nn]*.7))
          xidx4 <- which((rankx>=N[nn]*.7)&(rankx<N[nn]*.9))
          xidx5 <- which(rankx>=N[nn]*.9)
          yidx1 <- which(ranky<(N[nn]*.10))
          yidx2 <- which((ranky>=N[nn]*.10)&(ranky<N[nn]*.3))
          yidx3 <- which((ranky>=N[nn]*.3)&(ranky<N[nn]*.7))
          yidx4 <- which((ranky>=N[nn]*.7)&(ranky<N[nn]*.9))
          yidx5 <- which(ranky>=N[nn]*.9)
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
          rankx <- rank(ix[,j])
          ranky <- rank(iy[,j])
          xidx1 <- which(rankx<(N[nn]*.4))
          xidx2 <- which((rankx>=N[nn]*.4)&(rankx<N[nn]*.65))
          xidx3 <- which((rankx>=N[nn]*.65)&(rankx<N[nn]*.85))
          xidx4 <- which((rankx>=N[nn]*.85)&(rankx<N[nn]*.95))
          xidx5 <- which(rankx>=N[nn]*.95)
          yidx1 <- which(ranky<(N[nn]*.4))
          yidx2 <- which((ranky>=N[nn]*.4)&(ranky<N[nn]*.65))
          yidx3 <- which((ranky>=N[nn]*.65)&(ranky<N[nn]*.85))
          yidx4 <- which((ranky>=N[nn]*.85)&(ranky<N[nn]*.95))
          yidx5 <- which(ranky>=N[nn]*.95)
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
          rankx <- rank(ix[,j])
          ranky <- rank(iy[,j])
          xidx1 <- which(rankx<(N[nn]*.4))
          xidx2 <- which((rankx>=N[nn]*.4)&(rankx<N[nn]*.65))
          xidx3 <- which((rankx>=N[nn]*.65)&(rankx<N[nn]*.85))
          xidx4 <- which((rankx>=N[nn]*.85)&(rankx<N[nn]*.95))
          xidx5 <- which(rankx>=N[nn]*.95)
          yidx1 <- which(ranky<(N[nn]*.4))
          yidx2 <- which((ranky>=N[nn]*.4)&(ranky<N[nn]*.65))
          yidx3 <- which((ranky>=N[nn]*.65)&(ranky<N[nn]*.85))
          yidx4 <- which((ranky>=N[nn]*.85)&(ranky<N[nn]*.95))
          yidx5 <- which(ranky>=N[nn]*.95)
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
          rankx <- rank(ix[,j])
          ranky <- rank(iy[,j])
          xidx1 <- which(rankx<(N[nn]*.10))
          xidx2 <- which((rankx>=N[nn]*.10)&(rankx<N[nn]*.3))
          xidx3 <- which((rankx>=N[nn]*.3)&(rankx<N[nn]*.7))
          xidx4 <- which((rankx>=N[nn]*.7)&(rankx<N[nn]*.9))
          xidx5 <- which(rankx>=N[nn]*.9)
          yidx1 <- which(ranky<(N[nn]*.10))
          yidx2 <- which((ranky>=N[nn]*.10)&(ranky<N[nn]*.3))
          yidx3 <- which((ranky>=N[nn]*.3)&(ranky<N[nn]*.7))
          yidx4 <- which((ranky>=N[nn]*.7)&(ranky<N[nn]*.9))
          yidx5 <- which(ranky>=N[nn]*.9)
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
          rankx <- rank(ix[,j])
          ranky <- rank(iy[,j])
          xidx1 <- which(rankx<(N[nn]*.4))
          xidx2 <- which((rankx>=N[nn]*.4)&(rankx<N[nn]*.65))
          xidx3 <- which((rankx>=N[nn]*.65)&(rankx<N[nn]*.85))
          xidx4 <- which((rankx>=N[nn]*.85)&(rankx<N[nn]*.95))
          xidx5 <- which(rankx>=N[nn]*.95)
          yidx1 <- which(ranky<(N[nn]*.4))
          yidx2 <- which((ranky>=N[nn]*.4)&(ranky<N[nn]*.65))
          yidx3 <- which((ranky>=N[nn]*.65)&(ranky<N[nn]*.85))
          yidx4 <- which((ranky>=N[nn]*.85)&(ranky<N[nn]*.95))
          yidx5 <- which(ranky>=N[nn]*.95)
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
          rankx <- rank(ix[,j])
          ranky <- rank(iy[,j])
          xidx1 <- which(rankx<(N[nn]*.05))
          xidx2 <- which((rankx>=N[nn]*.05)&(rankx<N[nn]*.15))
          xidx3 <- which((rankx>=N[nn]*.15)&(rankx<N[nn]*.35))
          xidx4 <- which((rankx>=N[nn]*.35)&(rankx<N[nn]*.6))
          xidx5 <- which(rankx>=N[nn]*.6)
          yidx1 <- which(ranky<(N[nn]*.05))
          yidx2 <- which((ranky>=N[nn]*.05)&(ranky<N[nn]*.15))
          yidx3 <- which((ranky>=N[nn]*.15)&(ranky<N[nn]*.35))
          yidx4 <- which((ranky>=N[nn]*.35)&(ranky<N[nn]*.6))
          yidx5 <- which(ranky>=N[nn]*.6)
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
    if ((sd(xo)!=0)&(sd(yo)!=0)&(alphax1<1)&(alphax1>0)&(alphay1<1)&(alphay1>0)) # problematic data sets are flagged;
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
