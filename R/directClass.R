directClass<-function(data, dataY=NULL, simulated=TRUE, fold=5){
  pam<-lda<-qda<-rf<-ridge<-lasso<-elas<-svm<-knn<-nnet<-numeric();
  topk<-cost<-alpha<-kn<-mtry<-nodesize<-size<-decay<-topP<-lambda<-norm.fraction<-delta<-numeric();
  if(simulated==TRUE){
    for(i in 1:length(data)){
      ##--Class data--### 
      trainData<-t(data[[i]]$trainData);
      trainDataY<-as.factor(data[[i]]$trainLabels);
      allData<-rbind(t(data[[i]]$trainData), t(data[[i]]$testData));
      allDataY<-as.factor(c(data[[i]]$trainLabels, data[[i]]$testLabels));
      
      #--Format training set---#
      train<-GenerateLearningsets(y=trainDataY, method="MCCV", niter=1, ntrain=length(trainDataY), strat=T);
      #---Rank genes---#
      method<-"limma";
      sel.train<-GeneSelection(X=trainData, y=trainDataY, learningsets=train, method=method);
      
      ##--Generate a  (5 fold) MCCV set to get optimal number of genes (k) for DAs & NNET---###
      train1<-GenerateLearningsets(y=trainDataY, method="MCCV", niter=5, ntrain=floor(2/3*length(trainDataY)), strat=T);
      sel.train1<-GeneSelection(X=trainData, y=trainDataY, learningsets=train1, method=method);
      
      #-----------------------------------#
      #--Build and evaluate classifiers---#
      #-----------------------------------#
      cat(sep="\n");
      cat("Training LDA");
      cat(sep="\n");
      ##--Linear Discriminant Ananlysis (LDA)--###
      len0<-sum(trainDataY==0);len1<-sum(trainDataY==1);
      if(len0>10 & len1>10){
        if(len0<=len1){topk<-seq(5, len0, 5)}else{topk<-seq(5, len1, 5)};
      }else{
        if(len0<=len1){topk<-seq(3, len0, 3)}else{topk<-seq(3, len1, 3)};
      }
      
      MER_lda<-MER_qda<-numeric();
      for(j in 1:length(topk)){
        cv.lda<-classification(X=trainData, y=trainDataY, learningsets=train1, genesel=sel.train1, nbgene=topk[j], classifier=ldaCMA, trace=F);
        MER_lda[j]<-summary(evaluation(cv.lda, measure="misclassification"))[4];
      }
      ##--With the opitmal k build LDA  and evaluate with the testset---###
      ngene_lda <- min(topk[which(MER_lda == min(MER_lda))]);
      class.lda<-classification(X=allData, y=allDataY, learningsets=train, genesel=sel.train, nbgene=ngene_lda, classifier=ldaCMA, trace=F);
      
      cat(sep="\n");
      cat("Training SVM");
      cat(sep="\n");
      #---SVM---#
      cost<-c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 25, 50, 100, 250, 500);
      tune.svm<-tune(X=trainData, y=trainDataY, learningsets=train, classifier=svmCMA, fold=5, strat=T, trace=F, grids=list(cost=cost));
      class.svm<-classification(X=allData, y=allDataY, learningsets=train, classifier=svmCMA, tuneres=tune.svm, trace=F);
      
      cat(sep="\n");
      cat("Training Elastic Net");
      cat(sep="\n");
      #---Elastic Net--#
      alpha = 2^{-4:4};
      norm.fraction = seq(from=0.1, to=0.9, length=9);
      tune.elas<-tune(X=trainData, y=trainDataY, learningsets=train, classifier=ElasticNetCMA, fold=5, strat=T, trace=F, grids=list(alpha = alpha, norm.fraction = norm.fraction));
      class.elas<-classification(X=allData, y=allDataY, learningsets=train, classifier=ElasticNetCMA, tuneres=tune.elas, trace=F);
      
      cat(sep="\n");
      cat("Training KNN");
      cat(sep="\n");
      #---KNN#--#
      if(2/3*length(trainDataY)>10){kn<-seq(3,10,2)}else{kn<-seq(3,5,2)};
      tune.knn<-tune(X=trainData, y=trainDataY, learningsets=train, classifier=knnCMA, fold=5, strat=T, trace=F, grids=list(k=kn));
      class.knn<-classification(X=allData, y=allDataY, learningsets=train, classifier=knnCMA, tuneres=tune.knn, trace=F);
      
      cat(sep="\n");
      cat("Training RF");
      cat(sep="\n");
      #--RF--#
      mtry<-ceiling(c(0.1, 0.25, 0.5, 1, 2)*sqrt(ncol(trainData))); nodesize <- c(1:3);
      tune.rf<-tune(X=trainData, y=trainDataY, learningsets=train, classifier=rfCMA, fold=5, strat=T, trace=F, grids=list(mtry=mtry, nodesize=nodesize));
      class.rf<-classification(X=allData, y=allDataY, learningsets=train, classifier=rfCMA, tuneres=tune.rf, trace=F);
      
      cat(sep="\n");
      cat("Training Ridge");
      cat(sep="\n");
      #--Ridge---#
      lambda = 2^{-4:4}
      tune.plr<-tune(X=trainData, y=trainDataY, learningsets=train, classifier=plrCMA, fold=5, strat=T, trace=F, grids=list(lambda=lambda));
      class.plr<-classification(X=allData, y=allDataY, learningsets=train, classifier=plrCMA, tuneres=tune.plr, trace=F);
      
      cat(sep="\n");
      cat("Training Lasso");
      cat(sep="\n");
      #--Lasso--#
      norm.fraction = seq(from=0.1, to=0.9, length=9);
      tune.lasso<-tune(X=trainData, y=trainDataY, learningsets=train, classifier=LassoCMA, fold=5, strat=T, trace=F, grids=list(norm.fraction=norm.fraction));
      class.lasso<-classification(X=allData, y=allDataY, learningsets=train, classifier=LassoCMA, tuneres=tune.lasso, trace=F);
      
      cat(sep="\n");
      cat("Training PAM");
      cat(sep="\n");
      #--PAM---#
      delta<-seq(0.1, 5, 0.1);
      tune.scda<-tune(X=trainData, y=trainDataY, learningsets=train, classifier=scdaCMA, fold=5, strat=T, trace=F, grids=list(delta=delta));
      class.scda<-classification(X=allData, y=allDataY,learningsets=train, classifier=scdaCMA, tuneres=tune.scda, trace=F);
      
      cat(sep="\n");
      cat("Training QDA");
      cat(sep="\n");
      #---QDA---#
      len0<-sum(trainDataY==0);len1<-sum(trainDataY==1);
      if(len0>10 & len1>10){
        if(len0<=len1){topk<-seq(5, ceiling(len0/2), 5)}else{topk<-seq(5, ceiling(len0/2), 5)};
      }else{
        if(len0<=len1){topk<-seq(2, ceiling(len0/2), 3)}else{topk<-seq(3, ceiling(len0/2), 3)};
      }
      MER_qda<-numeric();
      for(j in 1:length(topk)){
        cv.qda<-classification(X=trainData, y=trainDataY, learningsets=train1, genesel=sel.train1, nbgene=topk[j], classifier=qdaCMA, trace=F);
        MER_qda[j]<-summary(evaluation(cv.qda, measure="misclassification"))[4];
      }
      ##--With the opitmal k build the QDA classifiers and evalate with the testset---###
      ngene_qda<-min(topk[which(MER_qda == min(MER_qda))]);
      class.qda<-classification(X=allData, y=allDataY, learningsets=train, genesel=sel.train, nbgene=ngene_qda, classifier=qdaCMA, trace=F);
      
      cat(sep="\n");
      cat("Training NNET");
      cat(sep="\n");
      #---NNET----# 
      len0<-sum(trainDataY==0);len1<-sum(trainDataY==1);
      if(len0>10 & len1>10){
        if(len0<=len1){topP<-seq(5, len0, 5)}else{topP<-seq(5, len1, 5)};
      }else{
        if(len0<=len1){topP<-seq(3, len0, 3)}else{topP<-seq(3, len1, 3)};
      }
      size<-c(1:5); decay<-c(0, 2^{-(4:1)})
      nnetComp<-list(); t=1; MER_nnet<-numeric();
      for(j in 1:length(topP)){
        assign("top", topP[j], envir=.BaseNamespaceEnv)
        for(k in 1:length(size)){
          assign("siz", size[k], envir=.BaseNamespaceEnv)
          for(l in 1:length(decay)){
            assign("dec", decay[l], envir=.BaseNamespaceEnv)
            cv.nnet<-classification(X=trainData, y=trainDataY, learningsets=train1, genesel=sel.train1, nbgene=get("top", envir=.BaseNamespaceEnv), classifier=nnetCMA, size=get("siz", envir=.BaseNamespaceEnv), decay=get("dec", envir=.BaseNamespaceEnv));
            MER_nnet[t]<-summary(evaluation(cv.nnet, measure="misclassification"))[4];
            nnetComp[[t]]<-c(topP[j], size[k], decay[l], MER_nnet[t]);
            t<-t+1
          }
        }
      }
      ##--With the opitmal parameters build and evaluate NNET---###
      comb<-which(MER_nnet==min(MER_nnet))[1];
      assign("optComb", nnetComp[[comb]],  envir=.BaseNamespaceEnv);
      class.nnet<-classification(X=allData, y=allDataY, learningsets=train, genesel=sel.train, nbgene=get("optComb", envir=.BaseNamespaceEnv)[1], classifier=nnetCMA, size=get("optComb", envir=.BaseNamespaceEnv)[2], decay=get("optComb", envir=.BaseNamespaceEnv)[3], trace=F);
      
      ##--Compare and get the missclassification error rates---###
      classifierlist<- list(class.svm, class.elas, class.knn, class.rf, class.nnet, class.plr, class.scda, class.lda, class.lasso, class.qda);
      comparison<- compare(classifierlist, measure="misclassification");
      svm[i]<-comparison[1,];
      elas[i]<-comparison[2,];
      knn[i]<-comparison[3,];
      rf[i]<-comparison[4,];
      nnet[i]<-comparison[5,];
      ridge[i]<-comparison[6,];
      pam[i]<-comparison[7,];
      lda[i]<-comparison[8,];
      lasso[i]<-comparison[9,];
      qda[i]<-comparison[10,]; 
      
      cat(sep="\n");
      cat(paste("THIS IS SIMULATED DATA", i, sep=" "));
      cat(sep="\n"); cat(sep="\n");
    }
  }else if(!is.null(dataY) & dim(data)[1]>0){
    for(i in 1:fold){
      #--Generate learning/test sets and rank the genes using the learning set---###
      train<-GenerateLearningsets(y=dataY, method="MCCV", niter=1, ntrain=floor(2/3*length(dataY)), strat=TRUE);
      method<-"limma";
      if(dim(data)[1]!=length(dataY)){
        data<-t(data)
      }
      sel.train<-GeneSelection(X=data, y=dataY, learningsets=train, method=method);
      
      ##--Generate a CV (3 fold) set to get optimal number of genes (k) for DAs & NNET---###
      cv.data<-data[train@learnmatrix, ];
      cv.dataY<-dataY[train@learnmatrix];
      
      cv.train<-GenerateLearningsets(y=cv.dataY, method="MCCV", niter=3, ntrain=floor(3/4*length(cv.dataY)), strat=TRUE);    
      
      #--Replicate the ranked genes (as many as the fold) to use for cross validation--#
      rankings<-importance<-c(1:dim(data)[2])
      for(j in 1:cv.train@iter){
        rankings<-rbind(rankings, sel.train@rankings[[1]]);
        importance<-rbind(importance, sel.train@importance[[1]]);
      }
      #--Create a new genesel object with the replicated genes
      cv.sel.train<-new("genesel", rankings=list(rankings[-1, ]), importance=list(importance[-1,]), method=sel.train@method, scheme=sel.train@scheme);
      
      #-----------------------------------#
      #--Build and evaluate classifiers---#
      #-----------------------------------#
      cat(sep="\n");
      cat("Training LDA");
      cat(sep="\n");
      ##--Linear Discriminant Ananlysis (LDA)--###
      len0<-sum(cv.dataY==0);len1<-sum(cv.dataY==1);
      if(len0>10 & len1>10){
        if(len0<=len1){topk<-seq(5, len0, 5)}else{topk<-seq(5, len1, 5)};
      }else{
        if(len0<=len1){topk<-seq(2, len0, 2)}else{topk<-seq(2, len1, 2)};
      }
      
      MER_lda<-MER_qda<-numeric();
      for(j in 1:length(topk)){
        cv.lda<-classification(X=cv.data, y=cv.dataY, learningsets=cv.train, genesel=cv.sel.train, nbgene=topk[j], classifier=ldaCMA, trace=F);
        MER_lda[j]<-summary(evaluation(cv.lda, measure="misclassification"))[4];
      }
      ##--With the opitmal k build the DA classifiers and test with the testset---###
      ngene_lda <- min(topk[which(MER_lda == min(MER_lda))]);
      class.lda<-classification(X=data, y=dataY, learningsets=train, genesel=sel.train, nbgene=ngene_lda, classifier=ldaCMA, trace=F);
      
      cat(sep="\n");
      cat("Training SVM");
      cat(sep="\n");
      #---SVM---#
      cost<-c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 25, 50, 100, 250, 500);
      tune.svm<-suppressWarnings(tune(X=cv.data, y=cv.dataY, classifier=svmCMA, fold=3, strat=T, trace=F, grids=list(cost=cost)));
      class.svm<-classification(X=data, y=dataY, learningsets=train, classifier=svmCMA, tuneres=tune.svm, trace=F);
      
      cat(sep="\n");
      cat("Training Elastic Net");
      cat(sep="\n");
      #---Elastic Net--#
      alpha = 2^{-4:2};
      norm.fraction = seq(from=0.01, to=0.9, length=9);
      tune.elas<-suppressWarnings(tune(X=cv.data, y=cv.dataY, classifier=ElasticNetCMA, fold=3, strat=T, trace=F, grids=list(alpha = alpha, norm.fraction = norm.fraction)));
      class.elas<-classification(X=data, y=dataY, learningsets=train, classifier=ElasticNetCMA, tuneres=tune.elas, trace=F);
      
      cat(sep="\n");
      cat("Training NNET");
      cat(sep="\n");
      #---KNN#--#
      if(2/3*length(dataY)-1>10){kn<-seq(3,10,2)}else{kn<-seq(3,5,2)};
      tune.knn<-suppressWarnings(tune(X=cv.data, y=cv.dataY, classifier=knnCMA, fold=3, strat=T, trace=F, grids=list(k=kn)));
      class.knn<-classification(X=data, y=dataY, learningsets=train, classifier=knnCMA, tuneres=tune.knn, trace=F);
      
      cat(sep="\n");
      cat("Training RF");
      cat(sep="\n");
      #--RF--#
      mtry<-ceiling(c(0.1, 0.25, 0.5, 1, 2)*sqrt(ncol(data))); nodesize <- c(1:3);
      tune.rf<-suppressWarnings(tune(X=cv.data, y=cv.dataY, classifier=rfCMA, fold=3, strat=T, trace=F, grids=list(mtry=mtry, nodesize=nodesize)));
      class.rf<-classification(X=data, y=dataY, learningsets=train, classifier=rfCMA, tuneres=tune.rf, trace=F);
      
      cat(sep="\n");
      cat("Training Ridge");
      cat(sep="\n");
      #--Ridge---#
      lambda = 2^{-4:4}
      tune.plr<-suppressWarnings(tune(X=cv.data, y=cv.dataY, classifier=plrCMA, fold=3, strat=T, trace=F, grids=list(lambda=lambda)));
      class.plr<-classification(X=data, y=dataY, learningsets=train, classifier=plrCMA, tuneres=tune.plr, trace=F);
      
      cat(sep="\n");
      cat("Training Lasso");
      cat(sep="\n");
      #--Lasso--#
      norm.fraction = seq(from=0.01, to=0.9, length=9);
      tune.lasso<-suppressWarnings(tune(X=cv.data, y=cv.dataY, classifier=LassoCMA, fold=3, strat=T, trace=F, grids=list(norm.fraction=norm.fraction)));
      class.lasso<-classification(X=data, y=dataY, learningsets=train, classifier=LassoCMA, tuneres=tune.lasso, trace=F);
      
      cat(sep="\n");
      cat("Training PAM");
      cat(sep="\n");
      #--PAM---#
      delta<-seq(0.1, 5, 0.1);
      tune.scda<-suppressWarnings(tune(X=cv.data, y=cv.dataY, classifier=scdaCMA, fold=3, strat=T, trace=F, grids=list(delta=delta)));
      class.scda<-classification(X=data, y=dataY,learningsets=train, classifier=scdaCMA, tuneres=tune.scda, trace=F);
      
      cat(sep="\n");
      cat("Training QDA");
      cat(sep="\n");
      #---QDA---#
      len0<-sum(cv.dataY==0);len1<-sum(cv.dataY==1);
      if(len0>10 & len1>10){
        if(len0<=len1){topk<-seq(2, ceiling(len0/2), 3)}else{topk<-seq(2, ceiling(len1/2), 3)};
      }else{
        if(len0<=len1){topk<-seq(2, ceiling(len0/2), 1)}else{topk<-seq(2, ceiling(len1/2), 1)};
      }
      MER_qda<-numeric();
      for(j in 1:length(topk)){
        cv.qda<-classification(X=cv.data, y=cv.dataY, learningsets=cv.train, genesel=cv.sel.train, nbgene=topk[j], classifier=qdaCMA, trace=F);
        MER_qda[j]<-summary(evaluation(cv.qda, measure="misclassification"))[4];
      }
      ##--With the opitmal k build the QDA classifiers and test with the testset---###
      ngene_qda<-min(topk[which(MER_qda == min(MER_qda))]);
      class.qda<-classification(X=data, y=dataY, learningsets=train, genesel=sel.train, nbgene=ngene_qda, classifier=qdaCMA, trace=F);
      
      cat(sep="\n");
      cat("Training NNET");
      cat(sep="\n");
      #---NNET----# 
      len0<-sum(trainDataY==0);len1<-sum(trainDataY==1);
      if(len0>10 & len1>10){
        if(len0<=len1){topP<-seq(5, len0, 5)}else{topP<-seq(5, len1, 5)};
      }else{
        if(len0<=len1){topP<-seq(3, len0, 3)}else{topP<-seq(3, len1, 3)};
      }
      size<-c(1:5); decay<-c(0, 2^{-(4:1)})
      nnetComp<-list(); t=1; MER_nnet<-numeric();
      for(j in 1:length(topP)){
        assign("top", topP[j], envir=.BaseNamespaceEnv)
        for(k in 1:length(size)){
          assign("siz", size[k], envir=.BaseNamespaceEnv)
          for(l in 1:length(decay)){
            assign("dec", decay[l], envir=.BaseNamespaceEnv)
            cv.nnet<-classification(X=cv.data, y=cv.dataY, learningsets=cv.train, genesel=cv.sel.train, nbgene=get("top", envir=.BaseNamespaceEnv), classifier=nnetCMA, size=get("siz", envir=.BaseNamespaceEnv), decay=get("dec", envir=.BaseNamespaceEnv));
            MER_nnet[t]<-summary(evaluation(cv.nnet, measure="misclassification"))[4];
            nnetComp[[t]]<-c(topP[j], size[k], decay[l], MER_nnet[t]);
            t<-t+1
          }
        }
      }
      ##--With the opitmal parameters build and evaluate NNET---###
      comb<-which(MER_nnet==min(MER_nnet))[1];
      assign("optComb", nnetComp[[comb]],  envir=.BaseNamespaceEnv);
      class.nnet<-classification(X=data, y=dataY, learningsets=train, genesel=sel.train, nbgene=get("optComb", envir=.BaseNamespaceEnv)[1], classifier=nnetCMA, size=get("optComb", envir=.BaseNamespaceEnv)[2], decay=get("optComb", envir=.BaseNamespaceEnv)[3], trace=F);
      
      ##--Compare and get the missclassification rate---###
      classifierlist<- list(class.svm, class.elas, class.knn, class.rf, class.nnet, class.plr, class.scda, class.lda, class.lasso, class.qda);
      comparison<- compare(classifierlist, measure="misclassification");
      svm[i]<-comparison[1,];
      elas[i]<-comparison[2,];
      knn[i]<-comparison[3,];
      rf[i]<-comparison[4,];
      nnet[i]<-comparison[5,];
      ridge[i]<-comparison[6,];
      pam[i]<-comparison[7,];
      lda[i]<-comparison[8,];
      lasso[i]<-comparison[9,];
      qda[i]<-comparison[10,];
      
      cat(sep="\n");
      cat(paste("THIS IS FOLD NUMBER", i, sep=" "));
      cat(sep="\n"); cat(sep="\n");
    }
  }else{
    stop("Either data or dataY does not exist");
  }
  rest<-cbind(as.numeric(svm), as.numeric(elas), as.numeric(knn), as.numeric(rf), as.numeric(nnet),
              as.numeric(ridge), as.numeric(pam), as.numeric(lda), as.numeric(lasso), as.numeric(qda));
  colnames(rest)<-c("SVM", "PLR12", "KNN", "RF", "NNET", "PLR2", "PAM", "LDA", "PLR1", "QDA");
  return(rest);
}