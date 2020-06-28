## data generation,continous items,true dist=Skewed, homogeneous ##
library(MASS)
lp=1000
ssize1=c(30,50,100,200) # average sample size;
NS1=c(5,20,50) # number of studies;
Nn=expand.grid(ssize1,NS1)
rho1=c(0,.1,.3,.5)
tau1=c(1,2) # tau-equivalency is met (1) or not (2);
I1=c(4,8,12)
cond=expand.grid(rho1,tau1,I1)
cond.idx=as.numeric(commandArgs(trailingOnly=TRUE))
rho=cond[cond.idx,1]
tau=cond[cond.idx,2]
I=cond[cond.idx,3]
## means and covariance matrix of the bi-variate variables in each study;
mu=c(0,0) # true score means;
varx=1 # true score variances;
vary=1
sigma=matrix(c(varx,rho*sqrt(varx*vary),rho*sqrt(varx*vary),vary),2,2) # fixed correlation in the population;
# factor loading with equal or unequal loadings: reliability is the same (.72,.84,.89 respectively for I=4,8,12);
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
datadir=paste0('/gpfs/home/qzhang4/metacorr/Homogeneous/DataHomoContSkew/cond',cond.idx+1000,'/')
for (sz in 1:nrow(Nn)) {
  ssize=Nn[sz,1]
  NS=Nn[sz,2]
  subdatadir=paste0(datadir,'n',ssize+1000,'NS',NS+1000,'/')
  system(paste0('rm -r ',subdatadir))
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
  flag=0 # flag problematic generated data sets for each condition;
  for (i in 1:lp) {
    ro=ra=alphax=alphay=NULL
    nn=1
    while (nn <= NS) {
      # initialize matrices for measurement errors and items;
      msmterrx=msmterry=ix=iy=matrix(NA,N[nn],I)
      # generate the true scores:skewed;
      z1u=rchisq(N[nn],1)
      z1=(z1u-mean(z1u))/sd(z1u)
      z2u=rchisq(N[nn],1)
      z2=(z2u-mean(z2u))/sd(z2u)
      xt=z1
      yt=rho*z1+sqrt(1-rho^2)*z2
      msmterrx[,1:I] <- mvrnorm(N[nn],mumsmterrx,sigmamsmterrx) # generate msmt error for x;
      msmterry[,1:I] <- mvrnorm(N[nn],mumsmterry,sigmamsmterry) # generate msmt error for y;
      ix[,1:I]=xt%*%t(as.matrix(faloading))+msmterrx[,1:I] # item scores for x;
      iy[,1:I]=yt%*%t(as.matrix(faloading))+msmterry[,1:I] # item scores for y;
      xo=apply(ix,1,sum)
      yo=apply(iy,1,sum)
      covx=cov(ix)
      covy=cov(iy)
      varsx=diag(covx)
      varsy=diag(covy)
      alphax1=I/(I-1)*(var(xo)-sum(varsx))/var(xo) # cronbach's alpha for x;
      alphay1=I/(I-1)*(var(yo)-sum(varsy))/var(yo) # cronbach's alpha for y;
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
    dataname=paste0(subdatadir,'data',i+1000,'.txt')
    write.table(dset,dataname,col.names=F,row.names=F)
  }
  flagname=paste0(subdatadir,'flag.txt')
  write.table(flag,flagname,col.names=F,row.names=F)
}
