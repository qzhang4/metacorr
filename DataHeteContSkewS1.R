## data generation,continous items,true dist=normal, heterogeneous ##
library(MASS)
lp=1000
ssize1=c(30,50,100,200) # average sample size;
NS1=c(5,20,50) # number of studies;
Nn=expand.grid(ssize1,NS1)
zrho1=c(0,.1028345,.3170305,.5618704) # the average z_bar;
I1=c(4,12) # number of items per scale;
cond=expand.grid(zrho1,I1)
cond.idx=as.numeric(commandArgs(trailingOnly=TRUE))
zrho=cond[cond.idx,1]
I=cond[cond.idx,2]
## means and covariance matrix of the bi-variate variables in each study;
#mu=c(0,0) # true score means;
#varx=1 # true score variances;
#vary=1
# factor loading with equal or unequal loadings: reliability is the same (.72,.84,.89 respectively for I=4,8,12);
faloading=rep(.6,I)
msmterrvarx=msmterrvary=rep(.55,I)

# measurement error covariance matrix;
mumsmterrx <- mumsmterry <- rep(0,I)
sigmamsmterrx <- sigmamsmterry <- diag(msmterrvarx) # msmt error covariance matrix;
datadir=paste0('/gpfs/home/qzhang4/metacorr/Heterogeneous/DataHeteContSkewS1/cond',cond.idx+1000,'/')
dir.create(datadir)
for (sz in 1:nrow(Nn)) {
  ssize=Nn[sz,1]
  NS=Nn[sz,2]
  subdatadir=paste0(datadir,'n',ssize+1000,'NS',NS+1000,'/')
  dir.create(subdatadir)
  ### data generation ###
  # Generate sample sizes with NS studies; make this fixed across replications.
  nn=1
  N=array(NA,NS)
  while (nn <= NS) {
    X=rchisq(1,3)
    Ni=floor((X-3)/sqrt(2*3)*ssize/2+ssize)
    if (Ni<=3) {next}
    else {
      N[nn]=Ni
      nn=nn+1}
  }
  # replicate data sets, each with NS effect sizes;
  for (i in 1:lp) {
    ro=NULL
    nn=1
    zr=rnorm(NS,zrho,.16)
    rhos=(exp(2*zr)-1)/(exp(2*zr)+1)
    while (nn <= NS) {
      # initialize matrices for measurement errors and items;
      #sigma=matrix(c(varx,rhos[nn]*sqrt(varx*vary),rhos[nn]*sqrt(varx*vary),vary),2,2)
      msmterrx=msmterry=ix=iy=matrix(NA,N[nn],I)
      # generate the true scores:normal;
      z1u=rchisq(N[nn],1)
      xt=(z1u-1)/sqrt(2)
      z2u=rchisq(N[nn],1)
      z2=(z2u-1)/sqrt(2)
      yt=rhos[nn]*xt+sqrt(1-rhos[nn]^2)*z2
      msmterrx[,1:I] <- mvrnorm(N[nn],mumsmterrx,sigmamsmterrx) # generate msmt error for x;
      msmterry[,1:I] <- mvrnorm(N[nn],mumsmterry,sigmamsmterry) # generate msmt error for y;
      ix[,1:I]=xt%*%t(as.matrix(faloading))+msmterrx[,1:I] # item scores for x;
      iy[,1:I]=yt%*%t(as.matrix(faloading))+msmterry[,1:I] # item scores for y;
      xo=apply(ix,1,sum)
      yo=apply(iy,1,sum)
      cor1=cor(xo,yo)
      ro=c(ro,cor1)
      nn=nn+1
      } 
    dset=cbind(ro,N)
    dataname=paste0(subdatadir,'data',i+1000,'.txt')
    write.table(dset,dataname,col.names=F,row.names=F)
  }
}