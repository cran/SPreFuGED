\name{estimateDataCha}
\alias{estimateDataCha}
\title{
A function to estimate data characteristics
}
\description{
This function fits limma models (or univariate Cox's models t) to determine DE (informative) genes and then computes the proportion of DE (informative) genes, log2FC (coefficients or betas), pairwise correlation of DE (informative) and noisy genes, genes' variances, sample sizes and proportion of events (for survival data).
}
\usage{
estimateDataCha(data, dataY, type = "Binary")
}
\arguments{
  \item{data}{
a matrix of expression values with rows corresponding to genes and columns to samples
}
  \item{dataY}{
a binary vector of class labels or a survival outcome as produced by Surv. Its length must be equal to the number of columns of data.
}
  \item{type}{
takes Binary(Default) or Survival as values and correspond to binary classification or survival prediction
}
}
\details{
At the moment, only binary classification has been implemented.
}
\value{
A 1x7 (for Binary) or 1x8 (for Survival) matrix containing the estimates (row) of the data characteristics (columns)
}
\references{
Jong VL, Novianti PW, Roes KCB & Eijkemans MJC. Selecting a classification function for class prediction with gene expression data. Bioinformatics (2016) 32(12): 1814-1822
}
\author{
Victor Lih Jong
}
\seealso{
\code{\link{fitLMEModel}}, \code{\link{SPreFu}} and \code{\link{plotSPreFu}}
}
\examples{
#Let us consider a single simulated train data as our real-life dataset
myCov<-covMat(pAll=100, lambda=2, corrDE=0.75, sigma=0.25);
myData<-generateGED(covAll=myCov, nTrain=30, nTest=10);
data<-myData[[1]]$trainData;
dataY<-myData[[1]]$trainLabels;
myDataCha<-estimateDataCha(data, dataY);
}
\keyword{Estimate data characteristics}
