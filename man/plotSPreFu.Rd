\name{plotSPreFu}
\alias{plotSPreFu}
\title{
A plotting function of the predicted performance of classification functions on a given dataset 
}
\description{
This function plots the predicted accuraces, Brier scores or Hazard ratios for each predictive function on a given dataset, as predicted by SPreFu.
}
\usage{
plotSPreFu(restSPreFu)
}
\arguments{
  \item{restSPreFu}{
an object returned by SPreFu
}
}
\value{
A plot of predictive funcions and their expected performance (Accuracy, Brier score or Integrated Brier score)
}
\references{
Jong VL, Novianti PW, Roes KCB & Eijkemans MJC. Selecting a classification function for class prediction with gene expression data. Bioinformatics (2016) 32(12): 1814-1822;
}
\author{
Victor Lih Jong
}
\seealso{
\code{\link{estimateDataCha}}, \code{\link{fitLMEModel}} and \code{\link{SPreFu}}
}
\examples{
#Let us consider a single simulated train data as our real-life dataset
myCov<-covMat(pAll=100, lambda=2, corrDE=0.75, sigma=0.25);
myData<-generateGED(covAll=myCov, nTrain=30, nTest=10);
data<-myData[[1]]$trainData;
dataY<-myData[[1]]$trainLabels;
myDataCha<-estimateDataCha(data, dataY);
\donttest{myFit<-fitLMEModel();}  #Takes roughly 250 Sec
\donttest{myPred<-SPreFu(myDataCha, myFit);}
\donttest{plotSPreFu(myPred);}
}
\keyword{Plot predicted performance}
