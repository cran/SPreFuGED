fitLMEModel<-function(type="Accuracy"){
  if(type=="Accuracy"){
    data(avAcc, envir=environment());
    #---Build the slected random effects model---#
    #--Standardized variables---#
    avAcc$StdVariance<-(avAcc$variance-mean(avAcc$variance))/sd(avAcc$variance);
    avAcc$StdDECorr<-(avAcc$deCorr-mean(avAcc$deCorr))/sd(avAcc$deCorr);
    avAcc$StdPropDE<-(avAcc$propDE-mean(avAcc$propDE))/sd(avAcc$propDE);
    avAcc$StdLog2FC<-(avAcc$log2FC-mean(avAcc$log2FC))/sd(avAcc$log2FC);
    avAcc$StdSampSize<-(avAcc$sampSize-mean(avAcc$sampSize))/sd(avAcc$sampSize);
    avAcc$StdOtherCorr<-(avAcc$otherCorr-mean(avAcc$otherCorr))/sd(avAcc$otherCorr);
    avAcc$logOddsAcc<-logit(1-(avAcc$acc+0.0001));
    #---Build model--#
    model<-lmer(logOddsAcc~StdSampSize + StdPropDE + StdVariance + StdDECorr + StdOtherCorr + StdLog2FC
                + StdSampSize*StdLog2FC + StdPropDE*StdLog2FC + StdVariance*StdLog2FC + StdDECorr*StdLog2FC + StdOtherCorr*StdLog2FC
                + StdSampSize*StdOtherCorr + StdPropDE*StdOtherCorr + StdVariance*StdOtherCorr +StdDECorr*StdOtherCorr
                + StdSampSize*StdDECorr + StdPropDE*StdDECorr + StdVariance*StdDECorr
                + StdSampSize*StdVariance + StdPropDE*StdVariance + StdSampSize*StdPropDE
                + (StdSampSize + StdPropDE + StdVariance + StdDECorr  + StdOtherCorr + StdLog2FC|class), data=avAcc, REML=T);
  }else if(type=="Probability"){
    stop("At this moment, survival Prediction has not been implemented yet");
  }else if(type=="Survival"){
    stop("At this moment, survival Prediction has not been implemented yet")
  }else{
    stop("Only one of: Accuracy, Probability or Survival is required")
  }
  restModel<-list(model=model, type=type, fitData=avAcc);
  return(restModel)
}