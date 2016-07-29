\name{SPreFuGED-package}
\alias{SPreFuGED-package}
\alias{SPreFuGED}
\docType{package}
\title{
\packageTitle{SPreFuGED}
}
\description{
\packageDescription{SPreFuGED}
}
\details{
The package allows one to determine an optimal predictive function among several functions for either binary direct classification, binary probabilistic classification or survival prediction for a given gene expression data. It also presents an interface to simulate gene expression data and to compare classification survival prediction functions on a given data.
The most important workflow of this package is as follows:
\enumerate{
\item Estimate the gene expression data characteristics using estimateDataCha
\item Fit a specific linear mixed effects model using fitLMEModel
\item Predict the performance of the functions on the given data using SPreFu
\item Plot the results of step 3 using plotSPreFu.
\item Other functions are(were) used for simulation studies presented in most of the references.
}
}
\author{
\packageAuthor{SPreFuGED}

Maintainer: \packageMaintainer{SPreFuGED}
}
\references{
Jong VL, Novianti PW, Roes KCB & Eijkemans MJC. Selecting a classification function for class prediction with gene expression data. Bioinformatics (2016) 32(12): 1814-1822;
}