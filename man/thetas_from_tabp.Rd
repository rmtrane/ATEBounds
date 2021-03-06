% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/probability_manipulations.R
\name{thetas_from_tabp}
\alias{thetas_from_tabp}
\title{Calculate thetas from joint conditionals}
\usage{
thetas_from_tabp(tabp)
}
\arguments{
\item{tabp}{a 2x2xn_z_levels table given P(Y = y, X = x | Z = z), y = 0,1, x = 0,1, and z = 0,1,..., n_z_levels - 1.}
}
\description{
Calculates P(X = 1 | Z = z) from a table of P(Y = y, X = x | Z = z)
}
