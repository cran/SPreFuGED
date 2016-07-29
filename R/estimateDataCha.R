estimateDataCha<-function(data, dataY, type="Binary"){
  if(type=="Binary"){
    #---Convert dataY to factor and check that the GED has samples in columns--#
    samples<-as.factor(dataY);
    if(dim(data)[2]!=length(samples)){
      data<-t(data)
    }
    
    #--Determine differentially expressed genes---#
    assaymat <- as.matrix(data);
    myexpset <- new("ExpressionSet", exprs = assaymat);
    eset <- exprs(myexpset);
    design <- model.matrix(~0 + samples);
    cont.mat <- makeContrasts(c(-1, 1), levels=design);
    fit <- lmFit(eset, design);
    fit.cont <- contrasts.fit(fit, cont.mat);
    fit2 <- eBayes(fit.cont);
    results <- topTable(fit2, number=dim(data)[1], lfc=1, p.value=1);
    
    if(dim(results)[1]<5){
      results <- topTable(fit2, number=dim(data)[1], lfc=0.5, p.value=1);
    }
    if(dim(results)[1]<5){
      results <- topTable(fit2, number=dim(data)[1], lfc=0.25, p.value=1);
    }
    if(dim(results)[1]<5){
      results <- topTable(fit2, number=dim(data)[1], lfc=0.125, p.value=1);
    }
    
    #---Variables---#
    p<-dim(data)[1]
    #--Informative Genes--#
    de<-dim(results)[1]*100/p;
    #--log fold change of informative genes--#
    logFC<-mean(abs(results$logFC));
    #--sample size--#
    n<-dim(data)[2]
    #--Variances--#
    variances<-apply(data, 1, var);
    lambda<-mean(variances);
    
    #--Correlation non-DE genes--#
    data2<-data[setdiff(rownames(data), rownames(results)), ]
    if(dim(data2)[1]>10000){
      warning("The number of non-DE genes is >10000. As such, a random sample of 10000 genes
              is being used to computed the correlation matrix");
      sel<-sample(1:dim(data2)[1], 10000, replace=F)
      corre<-cor(t(data2[sel, ]));
      otherCor<-corre[upper.tri(corre)]
      otherCorr<-sd(otherCor); #SD of the correlations
    }else{
      corre<-cor(t(data2));
      otherCor<-corre[upper.tri(corre)]
      otherCorr<-sd(otherCor); #SD of the correlations
    }
    rm(list=c("corre", "data2"));
    
    #--Correlation DE genes--#
    data2<-data[rownames(results), ];
    corre<-cor(t(data2));
    deCorr<-median(abs(corre[upper.tri(corre)])); #median of the absolute values of upper triangular
    rm(list=c("corre", "data2", "results"))
    
    dataCha<-data.frame(p, n, de, lambda, logFC, deCorr, otherCorr);
    names(dataCha)<-c("Genes", "sampSize", "propDE", "varaince", "log2FC", "deCorr", "otherCorr");
  }else if(type=="Survival"){
    stop("At this moment, survival Prediction has not been implemented yet")
  }else{
    stop("Only one of: Binary or Survival is required for the argument type")
  }
  return(dataCha)
}
