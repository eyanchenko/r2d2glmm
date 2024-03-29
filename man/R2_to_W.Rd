% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_functions.R
\name{R2_to_W}
\alias{R2_to_W}
\title{W for a given R2}
\usage{
R2_to_W(
  R2,
  b0 = 0,
  link = c("Gaussian", "Logistic", "Poisson", "Poisson_offsets", "ZIP", "NegBinom",
    "Arbitrary"),
  theta = 1,
  mu,
  sigma2,
  disp = FALSE
)
}
\arguments{
\item{R2}{value of R2}

\item{b0}{value of beta0}

\item{link}{link function used to compute R2. Default options are "Gaussian", "Logistic",
"Poisson", "Poisson_offsets", "ZIP" and "NegBinom". Also can have user-inputted link
function using "Arbitrary".}

\item{theta}{dispersion parameter (if necessary)}

\item{mu}{if link="Arbitrary", then this is the mean function of the linear predictor.
May depend on theta}

\item{sigma2}{if link="Arbitrary", then this is the variance function of the linear predictor.
May depend on theta}

\item{disp}{logical. TRUE if arbitrary link function depends on a dispersion parameter}
}
\value{
value of W
}
\description{
This function computes the value of W for
a given R2, b0 and link function
}
