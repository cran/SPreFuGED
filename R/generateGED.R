generateGED<-function(covAll, nTrain, nTest, log2FC=1, niter=3, prob=0.5){
  allData<-list();
  samp<-rbinom(n=nTrain, size=1, prob=prob);
  for(i in 1:niter){
    pAll<-dim(covAll$cov)[1];
    p2<-covAll$pie*pAll; p1<-floor(p2/2);
    ##--Generate training data--##
    mean0<-meansmodule <- runif(pAll, min=6, max=10)
    deGenes<-1:p2; #--Define DE genes
    #--Simulated learning data for class 0
    n0<-length(samp[which(samp==0)]);
    data0<-suppressWarnings(rmvnorm(n0, mean=mean0, sigma=covAll$cov, method="chol"));
    n1<-length(samp[which(samp==1)]);
    #--Simulated learning data for class 1
    mean1<-mean0;
    mean1[deGenes[1:p1]]<-mean1[deGenes[1:p1]] + log2FC; #--Increment up-regulated expression by fc
    mean1[deGenes[(p1+1):p2]]<-mean1[deGenes[(p1+1):p2]] - log2FC #--Decrement down-regulated expression by fc
    data1<-suppressWarnings(rmvnorm(n1, mean=mean1, sigma=covAll$cov, method="chol"));
    #--Combine class 0 and 1 data to training data--#
    train.data<-t(rbind(data0, data1));
    rownames(train.data)<-paste("Gene", 1:pAll, sep="");
    colnames(train.data)<-paste("Sample", 1:nTrain, sep="");
    rm(list=c("data0", "data1"));
    
    #---Genenrate test set--#
    test.samp<-rbinom(n=nTest, size=1, prob=prob);
    test.samp1<-which(test.samp==1); test.samp0<-which(test.samp==0);
    test.n1<-length(test.samp1); test.n0<-length(test.samp0);
    test.data0<-suppressWarnings(rmvnorm(test.n0, mean=mean0, sigma=covAll$cov, method="chol"));
    test.data1<-suppressWarnings(rmvnorm(test.n1, mean=mean1, sigma=covAll$cov, method="chol"));
    test.data<-t(rbind(test.data0, test.data1));
    rownames(test.data)<-paste("Gene", 1:pAll, sep="");
    colnames(test.data)<-paste("Sample", 1:nTest, sep="");
    
    allData[[i]]<-list(trainData=train.data, trainLabels=samp, testData=test.data, testLabels=test.samp)
  }
  return(allData)
}