\name{generateGED}
\alias{generateGED}
\title{
A function to simulate two groups gene expression data (GED)
}
\description{
This function draws two sets of vectors (train and test samples' labels) from a binomial distribution and generates two gene expression datasets (train and test data) from a multivariate normal distribution with a mean vector U[6,10] and a given within-class covariance matrix, at each iteration. 
}
\usage{
generateGED(covAll, nTrain, nTest, log2FC = 1, niter = 3, prob = 0.5)
}
\arguments{
  \item{covAll}{
an object returned by covMax or a list containing a covariance matrix cov and the proportion of DE genes pie.
}
  \item{nTrain}{
the number of samples in the training set 
}
  \item{nTest}{
the number of samples in the test set
}
  \item{log2FC}{
the absolute Log2 fold changes (effect sizes) for DE genes (Default is 1)
}
  \item{niter}{
the number of iterations (train/test datasets to be generated). Default is 3
}
  \item{prob}{
the probability of success for the binomial sampling. Default is 0.5
}
}
\value{
A list of length niter. Each element of which is a list containing:
\item{trainData}{a matrix of the training data}
\item{trainLabels}{a binary vector of class labels of the training samples}
\item{testData}{a matrix of the test data}
\item{testLabels}{a binary vector of class labels of the test samples}
}
\references{
Jong VL, Novianti PW, Roes KCB & Eijkemans MJC. Selecting a classification function for class prediction with gene expression data. Bioinformatics (2016) 32(12): 1814-1822
}
\author{
Victor Lih Jong
}

\seealso{
\code{\link{covMat}}, \code{\link{directClass}} and \code{\link{plotDirectClass}}
}
\examples{
myCov<-covMat(pAll=100, lambda=2, corrDE=0.75, sigma=0.25);
myData<-generateGED(covAll=myCov, nTrain=30, nTest=10);
myFirstTrainData<-myData[[1]]$trainData; myFirstTrainLabels<-myData[[1]]$trainLabels;
myFirstTestData<-myData[[1]]$testData; myFirstTestLabels<-myData[[1]]$testLabels;
}
\keyword{Generate gene expression data }
