% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiagree.R
\name{boot.many1}
\alias{boot.many1}
\title{Compare dependent (multilevel) Fleiss kappas using the clustered bootstrap method}
\usage{
boot.many1(cluster_id, data, ncat = 2, a.level = 0.05, ITN = 1000,
  summary_k = T)
}
\arguments{
\item{cluster_id}{a vector of lenght N with the identification number of the clusters}

\item{data}{a N x sum(ncat_g) matrix representing the classification of the N items by the observers in group g in the ncat_g categories. For each group, the number of categories can vary}

\item{ncat}{a vector with G elements indicating how many categories are considered to compute each kappa coefficient. For example, c(3,5) means that the three first columns correspond to the classification of subjects on a 3 categorical scale by a group of observers and the five last columns correspond to the classification of subjects on a 5 categorical scale by a group of observers.}

\item{a.level}{the significance level}

\item{ITN}{the number of bootstrap iterations}

\item{summary_k, }{if true, Hotteling's T square test is performed, if false, only the bootstraped coefficients are returned}
}
\value{
$kappa a G x 2 matrix with the kappa coefficients in the first column and their corresponding standard error in the second column

$T_test a vector of length 2 with the value of Hotelling's T test as first element and the corresponding p-value as second element

$confidence confidence intervals for the pairwise comparisons of the measures

$cor the G x G correlation matrix for the kappa coefficients

$K when summary_k is false, the ITN x G matrix with the bootstrapped kappa coefficients is returned
}
\description{
This function performs Hotelling's T square test using a variance-covariance matrix based on the bootstrap method to compare dependent multi-observers kappa coefficients
}
\details{
This function compare several Fleiss kappa coefficients using Hotelling's T square with the variance-covariance matrix obtained by the clustered bootstrap method. If only one kappa is computed, it returns the estimate and confidence interval.
}
\examples{
 
#dataset (not multilevel) (Fleiss, 1971)
data(fleiss_psy)
attach(fleiss_psy)
set.seed(123) # to have the same results than in Vanbelle (submitted)
a<-boot.many1(data=fleiss_psy[,2:6],cluster_id=fleiss_psy[,1],ncat=c(5),a.level=0.05,ITN=1000)
}
\references{
Fleiss J.L. (1971). Measuring nominal scale agreement among many raters. Psychological Bulletin 76, 378-382.

Vanbelle S. (2017) Comparing dependent agreement coefficients obtained on multilevel data. Biometrical Journal, 59 (5):1016-1034

Vanbelle S. (submitted) On the asymptotic variability of (multilevel) multirater kappa coefficients
}
\author{
Sophie Vanbelle \email{sophie.vanbelle@maastrichtuniversity.nl}
}
