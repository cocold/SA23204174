% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLRGibbs.R
\name{UpdateBetaHat}
\alias{UpdateBetaHat}
\title{Update the betahat}
\usage{
UpdateBetaHat(theta_array, matrix_z, covariates, y)
}
\arguments{
\item{theta_array}{the defined theta above}

\item{matrix_z}{the z matrix}

\item{covariates}{the new x subpopulation}

\item{y}{the dependent variable}
}
\value{
updated theta_array
}
\description{
Update the betahat
}
\examples{
\dontrun{
X<-matrix(rep(0,100),nrow=10,ncol=10)
y<-rep(1,100)
initialized<-init(3,X,y)
theta_array<-initialized$theta
new_Z<-UpdateArrayZ(theta_array,y,X,3)
covariates<-UpdateSubPop(theta_array,new_Z,X,3)
vector_beta_hat<-UpdateBetaHat(theta_array, new_Z, covariates)
}
}
