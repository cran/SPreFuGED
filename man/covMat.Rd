\name{covMat}
\alias{covMat}
\title{
A covariance matrix generating function.
}
\description{
For a given number of genes and a proportion of differentially expressed (informative) genes, the this function creates a covariance matrix by sampling variances from an exponential distribution with lambda and the correlation values corrDE and corrOther. Where corrOther is generated from a normal distribution with mean=0, and sigma as SD.
}
\usage{
covMat(pAll = 1000, pie = 0.05, lambda, corrDE, sigma)
}
\arguments{
  \item{pAll}{
the total number of genes (Default is 1000). For desktop users, we encourage pAll <=2500 for compuational reasons.
}
  \item{pie}{
a value in the interval (0, 1) and corresponds to the proportion of differentially expressed (informative) genes (Defualt is 0.05) 
}
  \item{lambda}{
a positive rate parameter for sampling variances from an exponential distribution. The smaller the value the larger the variances.
}
  \item{corrDE}{
a value in the interval [0, 1] specifying the correlation values of DE genes to each other. Half of which are up-regulated (positively associated to survival time) and the others are down-regulated (negatively associated to survival time). The inter-cluster (between up- and down-regulated) genes take negatively signed value of corrDE. The value 0 corresponds to complete independence of these DE genes.
}
  \item{sigma}{
a value in the interval [0, 1] specifying the distribution of correlations within noisy genes  and between noisy genes and informative genes. Where 0 means complete indipendence of noisy genes to each other and to informative genes.
}
}
\details{
This functions assumes three clusters of genes (up-regulated, down-regulated and noisy genes). While the pairwise correlations of the DE genes is a descrete value corrDE, the correlations of the non-DE genes are sampled from a normal distribution with mean zero and SD=sigma. Values beyond the interval [-1, 1] are unformly converted to that interval.
}
\value{
A list containing:
 \item{cov}{the covariance matrix generated}
 \item{pie}{the proportion of differentially expressed genes}
}
\references{
Jong VL, Novianti PW, Roes KCB & Eijkemans MJC. Selecting a classification function for class prediction with gene expression data. Bioinformatics (2016) 32(12): 1814-1822
}
\author{
Victor Lih Jong
}
\seealso{
\code{\link{generateGED}}
}
\examples{
myCov<-covMat(pAll=100, lambda=2, corrDE=0.75, sigma=0.25);
#Observe the covariance matrix of 6 genes, 2 each from up-regulated, down-regulated and non-DE
myCov$cov[c(1,2, 4,5, 30,31), c(1,2, 4,5, 30,31)];
myCov$pie;
}
\keyword{ covariance matrices}