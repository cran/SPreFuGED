\name{avAcc}
\alias{avAcc}
\docType{data}
\title{
Average accuracy data
}
\description{
A data frame of 11520 data points (rows) corresponding to the 1152 different simulation scenarios by 10 classification functions. The studied data variables and the average misclassification error rates (acc) of the classification functions on the columns.
}
\usage{data("avAcc")}
\format{
  A data frame with 11520 observations on the following 8 variables.
  \describe{
    \item{\code{class}}{a character vector of the predictive functions}
    \item{\code{variance}}{a numeric vector of the different values of the variance parameter}
    \item{\code{deCorr}}{a numeric vector of different pairwise correlation values within informative genes}
    \item{\code{propDE}}{a numeric vector different values of the proportion of informative genes}
    \item{\code{log2FC}}{a numeric vector of different values of absolute log2 fold changes}
    \item{\code{sampSize}}{a numeric vector of different values of absolute sample sizes}
    \item{\code{otherCorr}}{a numeric vector of the different values of the parameter for the pairwise correlations within non-informative and/or between non-informative and informative genes}
    \item{\code{acc}}{a numeric vector of misclassification error rates}
  }
}
\details{
This data contains the different values of the six studied factors and the error rates for the different fucntions in different simulation scenarios. And it is required for binary direct classification, to fit the LME model necessary to predict optimal function(s) for a given gene expression dataset.
}
\source{
Jong VL, Novianti PW, Roes KCB & Eijkemans MJC. Selecting a classification function for class prediction with gene expression data. Bioinformatics (2016) 32(12): 1814-1822
}
\examples{
data(avAcc)
## maybe str(avAcc) ; plot(avAcc) ...
}
\keyword{datasets}
