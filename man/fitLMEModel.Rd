\name{fitLMEModel}
\alias{fitLMEModel}
\title{
A function to fit a linear mixed effects (LME) model of performance measure on the studied variables
}
\description{
This function fits a LME model to the log-odds of accuracy (for binary direct classifiers), logit transformation of the transformed Brier Score (for binary probabilistic classifiers) or logit transformation of the transformed integrated Brier score (for survival data).
}
\usage{
fitLMEModel(type="Accuracy")
}
\arguments{
  \item{type}{
takes Accuracy (Default), for direct classifiers or Probability, for probabilistic classifiers or Survival, for survival predictions.
}
}
\details{
Depending of the value of type, this function uses either avAcc, avBS or avSurv data to build a LME model. Only the avAcc is available and hence LME model of log-odds accuracy is possible at the moment.
}
\value{
A list containing:
\item{model}{an object of class "lmer" for which several fucntions can be applied}
\item{type}{the type of predictions (Accuracy,  Probability or Survival)}
\item{fitData}{fitted data, contains the variables and their standardized values}
}
\references{
Jong VL, Novianti PW, Roes KCB & Eijkemans MJC. Selecting a classification function for class prediction with gene expression data. Bioinformatics (2016) 32(12): 1814-1822
}
\author{
Victor Lih Jong
}
\seealso{
\code{\link{estimateDataCha}}, \code{\link{SPreFu}} and \code{\link{plotSPreFu}}
}
\examples{
\donttest{myFit<-fitLMEModel();}  #Takes roughly 250 Sec
}
\keyword{ Fit LME model}
