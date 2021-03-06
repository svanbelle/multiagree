% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiagree.R
\name{fleiss.pair}
\alias{fleiss.pair}
\title{Compare G independent pairwise kappa coefficients}
\usage{
fleiss.pair(rater1, rater2, cluster_id, weight = "equal", multilevel = T,
  a.level = 0.05, cov, ITN = 1000, meth)
}
\arguments{
\item{rater1}{a vector with the ratings of the first observer}

\item{rater2}{a vector with the ratings of the second observer}

\item{cluster_id}{a vector with the identification number of the clusters}

\item{weight}{the weighting scheme to be used for kappa coefficients. 'unweighted' for Cohen's kappa, 'equal' for linear weights and 'squared' for quadratic weights}

\item{multilevel}{a binary indicator equal to TRUE if the data are multilevel and FALSE otherwiwse.}

\item{a.level}{the significance level}

\item{cov}{the covariate determining the L groups to be compared}

\item{ITN}{number of bootstrap iterations needed if the bootstrap procedure is chosen}

\item{meth}{the method to be used to compute the standard error of the kappa coefficients: 'delta' for the delta method and 'boot' for the bootstrap method.}
}
\value{
$kappa the value of the L kappa coefficients and their standard error

$chi the value of the chi-squared statistic with L-1 degrees of freedom

$p the p-value
}
\description{
This functions compare L independent kappa coefficients using the method of Fleiss (1981)
}
\details{
This function compare L independent kappa coefficients using the method of Fleiss (1981). The data have to be entered in a vertical format.
}
\examples{
 
 
#dataset (multilevel) (Vanbelle, xxx)

data(FEES)
attach(FEES)

fleiss.pair(rater1=val_CO,rater2=val_COR,cluster_id=subject,
weight="unweighted",multilevel=TRUE,meth='delta',cov=group)
}
\references{
Fleiss, J. L. (1981). Statistical methods for rates and proportions (2nd ed.). New York: John Wiley

Vanbelle S. Comparing dependent agreement coefficients obtained on multilevel data. submitted

Vanbelle S. (2014) A New Interpretation of the Weighted Kappa Coefficients. Psychometrika. Advance online publication.  doi: 10.1007/s11336-014-9439-4
}
\author{
Sophie Vanbelle \email{sophie.vanbelle@maastrichtuniversity.nl}
}
