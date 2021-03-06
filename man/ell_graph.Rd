% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiagree.R
\name{ell_graph}
\alias{ell_graph}
\title{This function plots a confidence ellipse when 3 kappas are compared with the clustered bootstrap method}
\usage{
ell_graph(boot_val, a.level = 0.05, xlim = c(-1, 1), ylim = c(-1, 1),
  ylab, xlab)
}
\arguments{
\item{boot_val}{a ITN x 3 matrix with the bootstrapped values of the 3 kappa coefficients obtained with  \code{boot.pair}}

\item{a.level}{the significance level to draw the confidence ellipse}

\item{xlim}{the limits  for the x-axis}

\item{ylim}{the limits  for the y-axis}

\item{ylab}{the label for the y-axis}

\item{xlab}{the label for the x-axis}
}
\description{
This function plots a confidence ellipse when 3 kappas are compared with the clustered bootstrap method
}
\details{
This function plots differences between the bootstrapped kappa coefficients and a confidence ellipse. The triangle represents the coordinate (0,0) and the square the mean of the two differences
}
\examples{
#dataset (multilevel) (Vanbelle, xxx)

data(FEES)
attach(FEES)
dat<-cbind(val_CO,val_COR,val_MH,val_MHR,val_TB,val_TBR) #formating the data matrix
Kappa_val<-boot.pair(data=dat,cluster_id=subject,weight='equal',summary_k=FALSE)
ell_graph(boot_val=Kappa_val,xlim=c(-0.8,0.4),ylim=c(-0.2,0.4),xlab="K(obs1)-K(cons)",ylab="K(obs2)-K(cons)")
}
\author{
Sophie Vanbelle \email{sophie.vanbelle@maastrichtuniversity.nl}
}
