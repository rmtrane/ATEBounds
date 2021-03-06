% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/probability_manipulations.R
\name{gammas_from_tabp}
\alias{gammas_from_tabp}
\title{Calculate gammas from joint conditionals}
\usage{
gammas_from_tabp(tabp)
}
\arguments{
\item{tabp}{a \eqn{2x2xk} table given \eqn{P(Y = y, X = x | Z = z), y = 0,1, x = 0,1}, and \eqn{z = 0,1,..., k - 1}.}
}
\description{
Calculates \eqn{P(Y = 1 | Z = z)} from a table of \eqn{P(Y = y, X = x | Z = z)}
}
