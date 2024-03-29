% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch_reduction.R
\name{fetch_reduction}
\alias{fetch_reduction}
\alias{fetch_reduction.Seurat}
\alias{fetch_reduction.SingleCellExperiment}
\alias{fetch_reduction.AnnDataR6}
\title{Get names of reduction keys}
\usage{
fetch_reduction(object, reduction, cells = NULL, dims = c(1, 2), ...)

\method{fetch_reduction}{Seurat}(object, reduction, cells = NULL, dims = c(1, 2))

\method{fetch_reduction}{SingleCellExperiment}(object, reduction, cells = NULL, dims = c(1, 2))

\method{fetch_reduction}{AnnDataR6}(object, reduction, cells = NULL, dims = c(1, 2))
}
\arguments{
\item{object}{a single-cell object. Currently, Seurat and
SingleCellExperiment objects are supported.}

\item{reduction}{the reduction to pull coordinates from}

\item{cells}{cells for which to pull reduction data}

\item{dims}{a two-element integer vector with the dimensions for which
data should be returned}

\item{...}{Currently unused.}
}
\description{
Given the name of a reduction and a set of dims, this function will return
the corresponding reduction data.
}
\section{Methods (by class)}{
\itemize{
\item \code{fetch_reduction(Seurat)}: Seurat objects

\item \code{fetch_reduction(SingleCellExperiment)}: SingleCellExperiment objects

\item \code{fetch_reduction(AnnDataR6)}: AnnDataR6 objects

}}
