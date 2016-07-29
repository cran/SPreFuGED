plotDirectClass<-function(restDirectClass){
  classFxn<-rep(colnames(restDirectClass), each=dim(restDirectClass)[1])
  acc<-c(restDirectClass[,1], restDirectClass[,2], restDirectClass[,3], restDirectClass[,4], restDirectClass[,5],
         restDirectClass[,6], restDirectClass[,7], restDirectClass[,8], restDirectClass[,9], restDirectClass[,10]);
  obsAcc<-data.frame(classFxn, acc);
  obsAcc$trAcc<-1-obsAcc$acc;
  obsAcc$classifier<-as.numeric(as.factor(obsAcc$classFxn));
  class_name<-sort(unique(as.character(obsAcc$classFxn)))
  bwplot(trAcc~classifier, data=obsAcc,  grid=T, pch=20, cex=1.5, horizontal=F,
         xlim=class_name, xlab="Classifier", ylab="Accuracy",
         main="Accuracy of each classification function on test sets");
}
