plotSPreFu<-function(restSPreFu){
  predData<-restSPreFu$dataCha;
  predData$classifier<-as.numeric(as.factor(predData$class));
  if(restSPreFu$type=="Accuracy"){
    predData$trsAcc<-as.numeric(predData$Acc)
    xyplot(trsAcc~classifier, data=predData,  type="p", grid=T, pch=20, cex=1.5, ylim=c(0,1),
           xlim=predData$class, xlab="Classifier", ylab="Predicted Accuracy",
           main=" Predicted Accuracy for each classification function");
  }else if(restSPreFu$type=="Probability"){
    stop("At this moment, probabilistic classification has not been implemented yet");
    #predData$trsBS<-as.numeric(predData$BS)
    #xyplot(trsBS~classifier, data=predData,  type="p", grid=T, pch=20, cex=1.5, ylim=c(0,1),
    #       xlim=predData$class, xlab="Classifier", ylab="Predicted transformed Brier score (trBS)",
    #       main=" Predicted transformed Brier score (trBS) for each classification function");
  }else if(restSPreFu$type=="Survival"){
    stop("At this moment, survival Prediction has not been implemented yet");
  }else{
    stop("restModel does not contain any of Accuracy, Probability or Survival");
  }
}