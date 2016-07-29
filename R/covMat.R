covMat<-function(pAll=1000, pie=0.05, lambda, corrDE, sigma){
  if(pie<=0|pie>=1|lambda<=0|corrDE<0|corrDE>1|sigma<0|sigma>1){
    stop("One of the variables pie, lambda, corrDE or sigma has a value out of the required range");
  }else{
    corr<-diag(1,pAll,pAll);
    eltsAll<-sum(upper.tri(corr)); #Number of cells in the upper triangular of the covariance matrix
    p2<-pie*pAll; p1<-floor(p2/2);
    eltsDE<-sum(upper.tri(corr[1:p2, 1:p2])); #cells of DE genes
    elts<-eltsAll-eltsDE  #Cells withwout DE genes
    varsmodule <- rexp(pAll, lambda); #Sample the vriances
    corrOther<-rnorm(elts, mean=0, sd=sigma); #Sample the non-DE genes correlations
    #Uniformly convert correlation values below -1 to  the interval  [-1, -0.5]
    dcoef<-runif(sum(corrOther<(-1)), min=-1, max=-0.5);
    corrOther[corrOther<(-1)]<-dcoef;
    #Uniformly convert correlation values above 1 to  the interval  [0.5, 1]
    ucoef<-runif(sum(corrOther>1), min=0.5, max=1)
    corrOther[corrOther>1]<-ucoef;
    
    #--Fill the upper-traingular of correlation matrix--#
    count=1;
    for(i in pAll:2){
      for(j in 1:(i-1)){
        if((i%in%1:p1 & j%in%1:p1)|(i%in%(p1+1):p2 & j%in%(p1+1):p2)){
          corr[i,j]<-corrDE; # Both up- and down-regulated genes take intra-class correlation values of vla1
        }else if((i%in%1:p1 & j%in%(p1+1):p2)|(i%in%(p1+1):p2 & j%in%1:p1)){
          corr[i,j]<-(-corrDE); # Up- and down-regulated genes take inter-class correlation values of  -vla1
        }else{
          corr[i,j]<-corrOther[count]; # Non-DE genes take both intra- and inter-class correlations from the sampled values
          count=count+1;
        }
      }
    }
    
    #--Fill the lower-traingular of correlation matrix--#
    all<-corr[lower.tri(corr)];
    corr<-t(corr)
    corr[lower.tri(corr)]<-all;
    
    #--Convert correlation matrix to covariance matrix--#
    diag(corr)<-sqrt(varsmodule);
    cov<-sdcor2cov(corr);
    restCov<-list(cov=cov, pie=pie);
    return(restCov);
  }
}