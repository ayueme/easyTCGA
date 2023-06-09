% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_gene.R
\name{plot_KM}
\alias{plot_KM}
\title{K-M plot according to the expression of marker}
\usage{
plot_KM(exprset, marker, clin, optimal_cut = TRUE, return_data = TRUE)
}
\arguments{
\item{exprset}{expression matrix}

\item{marker}{marker you want to display}

\item{clin}{a data.frame with two columns: "time" and "event", and 1 for
live, 0 for dead.}

\item{optimal_cut}{use "optimal" cutpoint to do survival analysis. If FALSE,
median of expression will be used. Optimal cutpoint is calculated by
survminer::surv_cutpoint().}

\item{return_data}{return the data for plotting, default is TRUE.}
}
\value{
K-M plot and plotting data
}
\description{
K-M plot according to the expression of marker
}
