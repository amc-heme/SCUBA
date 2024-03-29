% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reduction_dimnames.R
\name{reduction_dimnames}
\alias{reduction_dimnames}
\alias{reduction_dimnames.Seurat}
\alias{reduction_dimnames.SingleCellExperiment}
\alias{reduction_dimnames.AnnDataR6}
\title{Get names of reduction keys}
\usage{
reduction_dimnames(object, reduction, dims, ...)

\method{reduction_dimnames}{Seurat}(object, reduction, dims)

\method{reduction_dimnames}{SingleCellExperiment}(object, reduction, dims)

\method{reduction_dimnames}{AnnDataR6}(object, reduction, dims)
}
\arguments{
\item{object}{a single-cell object. Currently, Seurat and
SingleCellExperiment objects are supported.}

\item{reduction}{the reduction from which names should be formed}

\item{dims}{a two-element integer vector with the dimensions for which
names should be returned}

\item{...}{Currently unused.}
}
\description{
Given the name of a reduction and a set of dims, this function will return
the names of each dim as it appears in the reduction matrix. The output of
this function is passed to FetchData to pull information for reductions, and
it is also used to label reductions on DimPlots and feature plots.
}
\section{Methods (by class)}{
\itemize{
\item \code{reduction_dimnames(Seurat)}: Seurat objects

\item \code{reduction_dimnames(SingleCellExperiment)}: SingleCellExperiment objects

\item \code{reduction_dimnames(AnnDataR6)}: AnnDataR6 objects

}}
