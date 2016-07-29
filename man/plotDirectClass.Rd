\name{plotDirectClass}
\alias{plotDirectClass}
\title{
A plotting function for the performance (Accuracy) of direct classifiers
}
\description{
This function produces boxplots of the expected accuracies of each classification function from call to directClass.
}
\usage{
plotDirectClass(restDirectClass)
}
\arguments{
  \item{restDirectClass}{
an object returned by directClass.
}
}
\value{
A plot of classification funcions and their expected accuracies
}
\references{
Jong VL, Novianti PW, Roes KCB & Eijkemans MJC. Selecting a classification function for class prediction with gene expression data. Bioinformatics (2016) 32(12): 1814-1822;
}
\author{
Victor Lih Jong
}
\seealso{
\code{\link{covMat}}, \code{\link{generateGED}} and \code{\link{directClass}}
}
\examples{
#Let us use simulated data build the 10 classification functions
myCov<-covMat(pAll=100, lambda=2, corrDE=0.75, sigma=0.25);
myData<-generateGED(covAll=myCov, nTrain=30, nTest=10);
\donttest{myClassResultsSimulatedData<-directClass(data=myData);} #Takes roughly 60 Sec
\donttest{plotDirectClass(myClassResultsSimulatedData);}
}
\keyword{Plot expected accuracies}
