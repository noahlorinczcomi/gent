# matrix trace function
tr=function(x) sum(diag(x))
# effective number of independent SNPs
f_meff=function(ldlist) {
  if(!is.list(ldlist)) ldlist=list(ldlist)
  LD=matrix(0,nr=nrow(ldlist[[1]]),nc=ncol(ldlist[[1]]))
  for(i in 1:length(ldlist)) LD=LD+ldlist[[i]]
  LD=LD/length(ldlist)
  d=eigen(LD)$values
  sum(+(d>=1)+d*(d<=1))
}
# soft thresholing function
S=function(x,lambda) sapply(x,function(h) sign(h)*max(c(0,abs(h)-lambda)))
penmu=function(z,LD,upper=qnorm(0.975),reg_constant=1) {
  library(mvnfast)
  lams=seq(0,upper,0.05)
  pens=likes=regs=c()
  for(i in 1:length(lams)) {
    regz=S(z,lams[i])
    like=dmvn(z,regz,LD,log=TRUE)
    reg=reg_constant*sum(regz!=0)
    pens[i]=reg-like
    likes[i]=like
    regs[i]=reg
  }
  S(z,lams[which.min(pens)])
}
# calculate parameters of null or non-null gamma distribution
gampars=function(z,LD,null=TRUE) {
  mu=length(c(z))
  sigma=2*tr(LD%*%LD)
  if(!null) {
    muv=penmu(z,LD)
    mu=mu+sum(muv^2)
    sigma=sigma+4*t(muv)%*%LD%*%muv
  }
  rate=mu/sigma
  shape=mu^2/sigma
  list(shape=shape,rate=rate)  
}
# penalty function for mugent sel
sf=function(x) pos(2*(x-0.5))
# positive-adjust scalar
pos=function(x) ifelse(x<0,0,x)
# mixture deconvolution based on gamma distributions
mixp=function(z,LD) {
  p=ncol(z);m=nrow(z)
  if(!is.list(LD)) lapply(1:p,function(h) LD)
  alt_pars=list(); for(i in 1:p) alt_pars[[i]]=gampars(z[,i],LD[[i]],null=FALSE)
  null_pars=list(); for(i in 1:p) null_pars[[i]]=gampars(z[,i],LD[[i]],null=TRUE)
  alt_pars=do.call(rbind,alt_pars);null_pars=do.call(rbind,null_pars)
  alt_pars=matrix(c(unlist(alt_pars)),nr=p,nc=2);colnames(alt_pars)=c('shape','rate')
  null_pars=matrix(c(unlist(null_pars)),nr=p,nc=2);colnames(null_pars)=c('shape','rate')
  stats=colSums(z^2)
  p_alts=c()
  for(i in 1:nrow(alt_pars)) {
    p1=dgamma(stats[i],shape=alt_pars[i,1],rate=alt_pars[i,2])
    p0=dgamma(stats[i],shape=null_pars[i,1],rate=null_pars[i,2])
    p_alts[i]=p1/(p1+p0)
  }
  p_alts=sf(p_alts) # shrinkage
  pboth=prod(p_alts) # P(mu1!=0,mu2!=0)
  peither=sum(p_alts)-prod(p_alts) # P(mu1!=0 or mu2!=0)
  pone=peither-pboth # P(mu1!=0 or mu2!=0 and !(mu1!=0 and mu2!=0))
  pone # P(only one is associated)
}
# positive-adjust LD matrix
posadj=function(r) {
  lams=seq(1-1e-3,0.5-1e-3,length.out=20)
  I=diag(nrow(r))
  d=eigen(r)$values
  k=0
  while(any(d<0)) {
    k=k+1
    r=r*lams[k]+(1-lams[k])*I
    d=eigen(r)$values
  }
  r
}
# function to calculate the variance of vec[t(Z)%*%Z] under H0 of MuGenT
varmat=function(p,ldlist) {
  if(length(ldlist)!=p) stop(cat('The number of LD matrices in `ldlist` should be ',p,', not ',length(ldlist),'\n',sep=''))
  H=matrix(0,p^2,p^2)
  ix1=c(row(H[1:p,1:p]))
  ix2=c(col(H[1:p,1:p]))
  C=cbind(ix1,ix2)
  diags=seq(1,p^2,length.out=p) 
  for(i in 1:p^2) {
    for(j in 1:p^2) {
      left=C[i,]
      right=C[j,]
      both=c(left,right)
      boo1=length(unique(both))==2
      boo2=(length(unique(left))==2) & (length(unique(right))==2) 
      if(boo1 & boo2) H[i,j]=tr(ldlist[[left[1]]]%*%ldlist[[left[2]]])
      if(i==j & (i%in%diags)) H[i,j]=2*tr(ldlist[[left[1]]]%*%ldlist[[right[1]]])
    }
  }
  H
}