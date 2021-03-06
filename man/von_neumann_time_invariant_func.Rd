% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fisher_metrics.R
\name{von_neumann_time_invariant_func}
\alias{von_neumann_time_invariant_func}
\title{von_neumann_time_invariant_func}
\usage{
von_neumann_time_invariant_func(fisher_df)
}
\arguments{
\item{fisher_df}{A fisher output containing FI_means and time_windows}
}
\value{
a single row tibble with the time invariant von Neumann variance index
}
\description{
Calculates the von Neumann variance index accounting irregular intervals
}
