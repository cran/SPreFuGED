\name{directClass}
\alias{directClass}
\title{
A function to build and evaluate 10 different classification functions on a given gene expression data.
}
\description{
This function uses either simulated (train and test) or real-life gene expression data to build and evaluate binary direct classifiers with 10 different classification functions [LDA, KNN, NNET, PAM, QDA, RF, Ridge (PLR2), Lasso (PLR1), Elastic Net (PLR12) and SVM] by minimizing the misclassification error rates.
}
\usage{
directClass(data, dataY = NULL, simulated = TRUE, fold = 5)
}
\arguments{
  \item{data}{
an object returned by generateGED or a matrix of expression values containing genes in the rows and samples in the columns.
}
  \item{dataY}{
an optional vector of class lables that can be coerced to a factor. Must be supplied when the data argument is not simulated data. 
}
  \item{simulated}{
logical, indicating whether the data supplied is simulated data (Default is TRUE).
}
  \item{fold}{
the number of cross-validation times to divide the real-life data into 2/3 train and 1/3 test with stratisfication(Default is 5). For meaningful comparison we recommend fold=100.
}
}
\details{
See reference for detail on which classification functions and/or parameters are optimized in this function and how the classifiers are built and evaluated. Please note that for large datasets and large values of "niter" from generateGED or fold, this function might take quite sometime. Of course, if it takes time to train a single function, what more of training 10 functions at once?
}
\value{
An Lx10 matrix of misclassification error rates. Where L is the number of iterations (niter) when the data is simulated data from generateGED or the number of folds (fold=L) when the data is real-life data and 10 are the number of classification functions.
}
\references{
Jong VL, Novianti PW, Roes KCB & Eijkemans MJC. Selecting a classification function for class prediction with gene expression data. Bioinformatics (2016) 32(12): 1814-1822
}
\author{
Victor Lih Jong
}
\seealso{
\code{\link{covMat}}, \code{\link{generateGED}} and \code{\link{plotDirectClass}}
}
\examples{
#Let us use simulated data build the 10 classification functions
myCov<-covMat(pAll=100, lambda=2, corrDE=0.75, sigma=0.25);
myData<-generateGED(covAll=myCov, nTrain=30, nTest=10);
\donttest{myClassResultsSimulatedData<-directClass(data=myData);} #Takes roughly 60 Sec
}
\keyword{Construct Classifiers}
