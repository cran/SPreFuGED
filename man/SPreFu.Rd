\name{SPreFu}
\alias{SPreFu}
\title{
A function for selecting an optimal predictive function for a given gene expression data.
}
\description{
This function uses the LME model and the estimated data characteristics to predict the accuracy (for binary direct classifiers) or transformed Brier score (for binary probabilistic classifiers) or transformed integrated Brier scores (for survival predictions).
}
\usage{
SPreFu(dataCha, restModel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataCha}{
an object returned by estimateDataCha and contains estimates of the data characteristics.}
  \item{restModel}{
an object returned by fitLMEModel and contains the fitted LME model to be used for predictions.}
}
\value{
A list containing:
\item{dataCha}{a data frame of the estimated data characteristics, the predictive functions and their predicted performance}
\item{type}{is the type of prediction (Accuracy, Probability or Survival) and determines what kind of plots to be produced by plotSPreFu}
}
\references{
Jong VL, Novianti PW, Roes KCB & Eijkemans MJC. Selecting a classification function for class prediction with gene expression data. Bioinformatics (2016) 32(12): 1814-1822;
}
\author{
Victor Lih Jong
}
\seealso{
\code{\link{estimateDataCha}}, \code{\link{fitLMEModel}} and \code{\link{plotSPreFu}}
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
}
\keyword{Select predictive function}
