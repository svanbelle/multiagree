This R package permits the comparison of dependent kappa coefficients using a Hotelling's T square test.
Kappa coefficients could be obtained on multilevel data and computed between pairs of observers 
or between several observers.

To install the package, you have to type the following commands in Rstudio:

install.packages("devtools")
devtools::install_github("svanbelle/multiagree")

The function delta.pair compares dependent pairwise kappa coefficients (eventually on multilevel data) using the delta method to determine the variance-covariance matrix of the kappa coefficients

The function delta.many1 compares dependent Fleiss kappa coefficients obtained between several observers  (eventually on multilevel data) using the delta method to determine the variance-covariance matrix of the kappa coefficients

The function delta.many2 compares dependent Conger or Light kappa coefficients obtained between several observers  (eventually on multilevel data) using the delta method to determine the variance-covariance matrix of the kappa coefficients

The function boot.pair compares dependent pairwise kappa coefficients (eventually on multilevel data) using the bootstrap method to determine the variance-covariance matrix of the kappa coefficients

The function boot.many1 compares dependent Fleiss kappa coefficients obtained between several observers  (eventually on multilevel data) using the bootstrap method to determine the variance-covariance matrix of the kappa coefficients

The function boot.many2 compares dependent Conger kappa coefficients obtained between several observers  (eventually on multilevel data) using the bootstrap method to determine the variance-covariance matrix of the kappa coefficients


The function fleiss.pair compares independent pairwise kappa coefficients (eventually obtained on multilevel data) using the delta method to determine the variance-covariance matrix of the kappa coefficients


