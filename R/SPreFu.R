SPreFu<-function(dataCha, restModel){
  if(restModel$type=="Accuracy"){
    dataCha2<-dataCha
    for(i in 2:10){
      dataCha2<-rbind(dataCha2, dataCha);
    }
    dataCha2$class<-c("KNN", "LDA", "NNET", "PAM", "PLR1", "PLR12", "PLR2", "QDA", "RF", "SVM");
    avAccRest<-restModel$fitData
    #--Standardized data for predictions---#
    dataCha2$StdSampSize<-(dataCha2$sampSize - mean(avAccRest$sampSize))/sd(avAccRest$sampSize)
    dataCha2$StdPropDE<-(dataCha2$propDE - mean(avAccRest$propDE))/sd(avAccRest$propDE)
    dataCha2$StdVariance<-(dataCha2$varaince - mean(avAccRest$variance))/sd(avAccRest$variance)
    dataCha2$StdDECorr<-(dataCha2$deCorr - mean(avAccRest$deCorr))/sd(avAccRest$deCorr)
    dataCha2$StdOtherCorr<-(dataCha2$otherCorr - mean(avAccRest$otherCorr))/sd(avAccRest$otherCorr)
    dataCha2$StdLog2FC<-(dataCha2$log2FC - mean(avAccRest$log2FC))/sd(avAccRest$log2FC);
    
    #---Make predictions, convert to accuraces and save---#
    dataCha2$predVals<-predict(restModel$model, newdata=dataCha2, re.form=~(StdSampSize + StdPropDE + StdVariance + StdDECorr  + StdOtherCorr + StdLog2FC|class));
    dataCha2$Acc<-inv.logit(dataCha2$predVals);
  }else if(restModel$type=="Probability"){
    stop("At this moment, probabilistic classification has not been implemented yet");
  }else if(restModel$type=="Survival"){
    stop("At this moment, survival Prediction has not been implemented yet");
  }else{
    stop("restModel does not contain any of Accuracy, Probability or Survival");
  }
  restSPreFu<-list(dataCha=dataCha2, type=restModel$type)
  return(restSPreFu);
}