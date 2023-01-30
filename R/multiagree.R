
#'  Compare dependent pairwise kappas (delta method)
#' 
#' This function performs Hotelling's T square test using a variance-covariance matrix based on the delta method to compare dependent pairwise kappa coefficients
#' 
#' @param data a N x R matrix representing the classification of the N items by the R observers. The kappa coefficients are computed pairwise between column (1,2), (3,4), etc....
#' @param cluster_id a vector of lenght N with the identification number of the K clusters
#' @param weight the weighting scheme to be used in the computation of the kappa coefficients: 'unweighted' for Cohen's kappa, 'equal' for linear weights and 'squared' for quadratic weights
#' @param multilevel a binary indicator equal to TRUE in the presence of multilevel data and FALSE otherwise
#' @param a.level the significance level
#' @return $kappa a G x 2 matrix with the G kappa coefficients to be compared in the first column and their corresponding standard error in the second column
#' @return $T_test a vector of length 2 with the value of Hotelling's T square test as first element and the p-value as second element
#' @return $confidence confidence intervals for the pairwise comparisons of kappa coefficients
#' @return $var the G x G correlation matrix of the kappa coefficients 
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function compare several dependent kappa coefficients obtained between pairs of observers. It uses Hotelling's T square test with the variance-covariance matrix obtained by the delta method. If only a single kappa coefficient is computed, the kappa coefficient and its standard error are returned.
#' @references Vanbelle S. and Albert A. (2008). A bootstrap method for comparing correlated kappa coefficients.  Journal of Statistical Computation and  Simulation, 1009-1015
#' @references Vanbelle S. (in press). Comparing dependent agreement coefficients obtained on multilevel data. Biometrical journal. doi: 10.1002/bimj.201600093
#' @references Vanbelle S. (2014). A New Interpretation of the Weighted Kappa Coefficients. Psychometrika. Advance online publication.  doi: 10.1007/s11336-014-9439-4
#' @export
#' @importFrom stats na.omit
#' @examples 
#'  
#' #dataset (not multilevel) (Vanbelle and Albert, 2008)
#' 
#' data(depression)
#' attach(depression)
#' delta.pair(data=cbind(diag,BDI,diag,GHQ),cluster_id=ID,weight='unweighted',multilevel=FALSE)
#' 
#' #dataset (multilevel) (Vanbelle, xxx)
#' 
#' data(FEES)
#' attach(FEES)
#' dat<-cbind(val_CO,val_COR,val_MH,val_MHR,val_TB,val_TBR) #formating the data matrix
#' delta.pair(data=dat,cluster_id=subject,weight='equal')

delta.pair<-function(data,cluster_id,weight,multilevel=T,a.level=0.05)
{
  if (a.level>1 | a.level<0) stop("the significance level should be between 0 and 1")
  if (ncol(data)%% 2!=0) stop("the data set should contain an even number of columns")
  
  #CREATE THE COMPLETE CASE DATASET
  data_<-na.omit(cbind(data,cluster_id))
  data2<-data_[,1:ncol(data)]
  c_id<-data_[,ncol(data)+1]
  
  
  #CREATE CONSECUTIVE ID FOR THE CLUSTERS
  id1<-unique(c_id);
  idd2<-rep(0,length(id1));
  for (i in 1:length(id1)){idd2<-ifelse(c_id==id1[i],i,idd2)}
  
  #CREATE FACTORS WITH THE CLASSIFICATION OF THE OBSERVERS (not usefull)
  all<-c(data2)
  all.f<-factor(all,levels=names(table(all)))
  ncat<-nlevels(all.f)  			#number of categories of the scale
  
  K<-max(idd2)						    #number of clusters
  nrater<-ncol(data)					#number of raters
  nkappa<-nrater/2  				  #number of kappa statistics to be computed
  
  
  
  #CREATE WEIGHT FOR THE ITEMS
  xi<-prop.table(table(idd2))			
  
  #CREATE THE WEIGHTED MATRX
  
if (weight=="equal"){
 w<-outer(1:ncat,1:ncat,function(i,j){1-abs(i-j)/(ncat-1)})
}

if (weight=="unweighted"){
    w<-diag(ncat)
  }
  
  
if (weight=="squared"){
    w<-outer(1:ncat,1:ncat,function(i,j){1-abs(i-j)^2/(ncat-1)^2})
  }
  
  #CREATE THE OBSERVED AGREEMENT
  Po<-replicate(nkappa,matrix(0,ncol=K,nrow=1), simplify=FALSE)
  ps<-replicate(nkappa,matrix(0,ncol=ncat,nrow=ncat), simplify=FALSE)
  
  for (t in 0:(nkappa-1)){
    for (s in 1:K){
      ps<-prop.table(table(factor(data2[,(2*t)+1][idd2==s],levels=names(table(all))),factor(data2[,(2*t)+2][idd2==s],levels=names(table(all)))))
      Po[[t+1]][s]<-sum(w*ps)
    }
  }
  
  
  
  #CREATE THE MARGINALS AT THE ITEM LEVEL
  p<-replicate(ncol(data),matrix(0,ncol=max(idd2),nrow=ncat), simplify=FALSE)
  
  
  for (i in 1:nrater){
    for (s in 1:K){
      p[[i]][,s]<-prop.table(table(factor(data2[,i][idd2==s],levels=names(table(all)))))
    }
  }
  
  #CREATE THE MATRIX OMEGA
  Omega<-(diag(K)-xi%*%matrix(1,ncol=K,nrow=1))%*%diag(xi^2)%*%(diag(K)-matrix(1,ncol=1,nrow=K)%*%t(xi))
  
  #CREATE THE MATRIX V
  V<-matrix(NA,ncol=(nkappa+nrater*ncat),nrow=(nkappa+nrater*ncat))
  
  #var-cov for observed agreement
  for (i in 1:nkappa){
    for (j in 1:nkappa){
      V[i,j]<-(K^2/(K-1))*Po[[i]]%*%Omega%*%t(Po[[j]])
    }
  }
  

  #var-cov for observed agreement versus marginals
  for (i in (1:nkappa)){
    for (j in (1:nrater)){
      V[i,(nkappa+1+(j-1)*ncat):(nkappa+j*ncat)]<-(K^2/(K-1))*Po[[i]]%*%Omega%*%t(p[[j]])
    }
  }
  
  #var-cov for marginals
  for (i in 1:nrater){
    for (j in (i:nrater)){
      V[(nkappa+1+(i-1)*ncat):(nkappa+i*ncat),(nkappa+1+(j-1)*ncat):(nkappa+j*ncat)]<-(K^2/(K-1))*p[[i]]%*%Omega%*%t(p[[j]])
    }
  }
  
  
  
  for (i in (2:(nkappa+nrater*ncat))){
    for (j in 1:(i-1)){
      V[i,j]<-V[j,i]
    }
  }
  
  
  
  #CREATE THE MARGINALS OVERALL
  M<-replicate(nrater,matrix(0,ncol=ncat,nrow=1), simplify=FALSE)
  
  for (i in 1:nrater){
    M[[i]]<-prop.table(table(factor(data2[,i],levels=names(table(all)))))
  }
  
  
  #CREATE THE MATRIX J
  J<-matrix(0,2*nkappa,(nkappa+nrater*ncat))
  
  if (nkappa>1){
    diag(J[1:nkappa,1:nkappa])<-1
  }
  if (nkappa==1){
    J[1:nkappa,1:nkappa]<-1
  }
  
  for (i in (0:(nkappa-1))){
    J[(nkappa+i+1),(nkappa+1+2*i*ncat):(nkappa+(1+i)*2*ncat)]<-c(t(M[[(2*i)+2]])%*%t(w),t(M[[(2*i)+1]])%*%w) 
  }
  
  #COMPUTE VARPSY
  varpsy<-J%*%V%*%t(J)
  
  #COMPUTE THE VECTOR OF EXPECTED AGREEMENT
  
  Pe<-matrix(NA,ncol=1,nrow=nkappa)
  P0<-matrix(NA,ncol=1,nrow=nkappa)
  
  for (i in 0:(nkappa-1)){
    Pe[i+1]<-t(M[[(2*i)+1]])%*%w%*%M[[(2*i)+2]]
    P0[i+1]<-sum(Po[[i+1]]*matrix(xi,ncol=K,nrow=1))               
    
  }
  
  s<-matrix(0,nkappa,2*nkappa)
  if (nkappa>1){
    diag(s[1:nkappa,1:nkappa])<-1/(1-Pe)
    diag(s[1:nkappa,(nkappa+1):(2*nkappa)])<-(P0-1)/((1-Pe)^2)
  }
  if (nkappa==1){
    s[1:nkappa,1:nkappa]<-1/(1-Pe)
    s[1:nkappa,(nkappa+1):(2*nkappa)]<-(P0-1)/((1-Pe)^2)
  }
  #COMPUTE THE VARIANCE-COVARIANCE MATRIX FOR THE KAPPAS
  
  if (multilevel==T){
    var_kappa<-s%*%varpsy%*%t(s)/K
    
  }
  
  if (multilevel==F){
    var_kappa<-s%*%varpsy%*%t(s)*(K-1)/K^2#correction factor for one unit per cluster
    
  }
  
  kappa<-matrix((P0-Pe)/(1-Pe),ncol=1)
  
  if (nkappa==1){
    bound<-c(NA,NA)
    bound[1]<-kappa-qnorm(1-a.level/2)*sqrt(var_kappa)
    bound[2]<-kappa+qnorm(1-a.level/2)*sqrt(var_kappa)
   
    result1<-matrix(c(kappa,sqrt(var_kappa)),nrow=1)
    colnames(result1) <-c("Kappa","SE(Kappa)")
    rownames(result1) <-""
    
    
    result2<-matrix(bound,nrow=1)
    colnames(result2) <-c("","")
    rownames(result2) <-""

    results<-list("kappa"=result1,"confidence"=result2,"level"=a.level)
    class(results)<-"kap1"
    return(results)
    
  }
  
  if (nkappa>1){hot.test(kappa=kappa,var_kappa=var_kappa,a.level=a.level,NN=K)}
}


#'   Compare dependent (multilevel) Fleiss kappas using the delta method
#' 
#'   This function performs Hotelling's T square test using a variance-covariance matrix based on the bootstrap method to compare dependent multi-observers kappa coefficients
#' 
#' @param cluster_id a vector of lenght N with the identification number of the clusters
#' @param data a N x sum(ncat_g) matrix representing the classification of the N items by the observers in group g in the ncat_g categories. For each group, the number of categories can vary 
#' @param ncat a vector with G elements indicating how many categories are considered to compute each kappa coefficient. For example, c(3,5) means that the three first columns correspond to the classification of subjects on a 3 categorical scale by a group of observers and the five last columns correspond to the classification of subjects on a 5 categorical scale by a group of observers.
#' @param a.level the significance level
#' @param multilevel a binary indicator equal to TRUE in the presence of multilevel data and FALSE otherwise
#' @return $kappa a G x 2 matrix with the kappa coefficients in the first column and their corresponding standard error in the second column
#' @return $T_test a vector of length 2 with the value of Hotelling's T test as first element and the corresponding p-value as second element
#' @return $confidence confidence intervals for the pairwise comparisons of the measures
#' @return $cor the G x G correlation matrix for the kappa coefficients
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function compare several Fleiss kappa coefficients using Hotelling's T square with the variance-covariance matrix obtained by the delta method. If only one kappa is computed, it returns the estimate and confidence interval. 
#' @export
#' @importFrom stats aggregate na.omit
#' @references Fleiss J.L. (1971). Measuring nominal scale agreement among many raters. Psychological Bulletin 76, 378-382.
#' @references Vanbelle S. (2017) Comparing dependent agreement coefficients obtained on multilevel data. Biometrical Journal, 59 (5):1016-1034
#' @references Vanbelle S. (submitted) On the asymptotic variability of (multilevel) multirater kappa coefficients
#' 

#' @examples  
#' #dataset (not multilevel) (Fleiss, 1971)

#' data(fleiss_psy)
#' attach(fleiss_psy)
#' delta.many1(data=fleiss_psy[,2:6],cluster_id=fleiss_psy[,1],ncat=c(5),a.level=0.05,multilevel=TRUE)
#' 

delta.many1<-function(data,cluster_id,ncat=c(2,3),multilevel=T,a.level=0.05)
{
  
  #CREATE THE COMPLETE CASE DATASET
  n_<-c(0,ncat)
  data_<-cbind(data,cluster_id)[rowSums(sapply(1:length(ncat),function(x){rowSums(data[,(1+sum(n_[1:x])):(sum(ncat[1:x]))])})>1)>=length(ncat),]
  data2<-data_[,1:ncol(data)]
  c_id<-data_[,ncol(data)+1]
  
  
  #CREATE CONSECUTIVE ID FOR THE CLUSTERS
  id1<-unique(c_id);
  idd2<-rep(0,length(id1));
  for (i in 1:length(id1)){idd2<-ifelse(c_id==id1[i],i,idd2)}
  
  #CREATE WEIGHT FOR THE ITEMS
  xi<-prop.table(table(idd2))			
  
  #CREATE FACTORS WITH THE CLASSIFICATION OF THE OBSERVERS (not usefull)
  K<-max(idd2)						    #number of clusters
  nrat<-sapply(1:length(ncat),function(x){rowSums(data2[,(1+sum(n_[1:x])):(sum(ncat[1:x]))])}) 				#number of raters per observation
  nkappa<-length(ncat)
  
  #CREATE THE OBSERVED AGREEMENT FOR EACH CLUSTER
  oprim_i<-sapply(simplify=FALSE,1:length(ncat),function(x){rowSums(data2[,(1+sum(n_[1:x])):(sum(ncat[1:x]))]*(data2[,(1+sum(n_[1:x])):(sum(ncat[1:x]))]-1)/(nrat[,x]*(nrat[,x]-1)))})
  Po<-t(data.matrix(unname(aggregate(oprim_i, by=list(idd2),FUN=mean, na.rm=TRUE)[,-1])))
  
  #CREATE THE MARGINALS AT THE ITEM LEVEL
  
  
  if (nkappa>1){p_ij<-sapply(1:length(ncat),function(x){sweep(data2[,(1+sum(n_[1:x])):(sum(ncat[1:x]))],1,nrat[,x],"/")},simplify=FALSE)}
  if (nkappa==1){p_ij<-sweep(data2,1,nrat,"/")}
  
  p_j<-data.matrix(aggregate(p_ij, by=list(idd2),FUN=mean, na.rm=TRUE)[,-1])
  p<-matrix(xi, ncol = K, nrow = 1)%*%p_j
  
  #CREATE THE MATRIX OMEGA
  Omega<-(diag(K)-xi%*%matrix(1,ncol=K,nrow=1))%*%diag(xi^2)%*%(diag(K)-matrix(1,ncol=1,nrow=K)%*%t(xi))
  
  
  #CREATE THE MATRIX V
  V<-matrix(NA,ncol=(nkappa+sum(ncat)),nrow=(nkappa+sum(ncat)))
  
  #var-cov for observed agreement
  V[1:nkappa,1:nkappa]<-K^2/(K-1)*Po%*%Omega%*%t(Po)
  
  
  #var-cov for observed agreement versus marginals
  V[1:nkappa,(nkappa+1):(nkappa+sum(ncat))]<-K^2/(K-1)*Po%*%Omega%*%p_j
  
  #var-cov for marginals
  V[(nkappa+1):(nkappa+sum(ncat)),(nkappa+1):(nkappa+sum(ncat))]<-(K^2/(K-1))*t(p_j)%*%Omega%*%p_j
  
  V[,1:nkappa]<-t(V[1:nkappa,])
  
  #compute the observed and expected agreement
  Po1<-Po%*%matrix(xi, ncol = 1, nrow = K)
  
  Pe1<-matrix(sapply(1:length(ncat),function(x){t(p[(1+sum(n_[1:x])):(sum(ncat[1:x]))]%*%p[(1+sum(n_[1:x])):(sum(ncat[1:x]))])}),ncol=1)
  
  #COMPUTE THE VARIANCE-COVARIANCE MATRIX FOR (Po,Pe)
  J <- matrix(0, 2 * nkappa, (nkappa + sum(ncat)))
  if (nkappa > 1) {
    diag(J[1:nkappa, 1:nkappa]) <- 1
  }
  if (nkappa == 1) {
    J[1:nkappa, 1:nkappa] <- 1
  }
  for (i in (1:nkappa)) {
    J[(nkappa+i), (nkappa+1+sum(n_[i])):(nkappa +sum(ncat[1:i]))] <-2*t(p[(1+sum(n_[i])):(sum(ncat[1:i]))]) 
  }
  #the formula of J was not correct prior to version 3.01
  var_PoPe<-J%*%V%*%t(J)
  
  #COMPUTE THE VARIANCE-COVARIANCE MATRIX FOR FLEISS KAPPA
  s <- matrix(0, nkappa, 2 * nkappa)
  if (nkappa > 1) {
    diag(s[1:nkappa, 1:nkappa]) <- 1/(1-Pe1)
    diag(s[1:nkappa, (nkappa+1):(2*nkappa)]) <- (Po1-1)/((1 - Pe1)^2)
  }
  if (nkappa == 1) {
    s[1:nkappa, 1:nkappa] <- 1/(1 - Pe1)
    s[1:nkappa, (nkappa + 1):(2 * nkappa)] <- (Po1 - 1)/((1-Pe1)^2)
  }
  
  
  
  if (multilevel==T){
    var_kappa<-s%*%var_PoPe%*%t(s)/K
    
  }
  
  if (multilevel==F){
    var_kappa<-s%*%var_PoPe%*%t(s)*(K-1)/K^2#correction factor for one unit per cluster
    
  }
  
  kappa<-(Po1-Pe1)/(1-Pe1)
  
  if (nkappa == 1) {
    bound <- c(NA, NA)
    bound[1] <- kappa - qnorm(1 - a.level/2) * sqrt(var_kappa)
    bound[2] <- kappa + qnorm(1 - a.level/2) * sqrt(var_kappa)
    result1 <- matrix(c(kappa, sqrt(var_kappa)), nrow = 1)
    colnames(result1) <- c("Kappa", "SE(Kappa)")
    rownames(result1) <- ""
    result2 <- matrix(bound, nrow = 1)
    colnames(result2) <- c("", "")
    rownames(result2) <- ""
    results <- list(kappa = result1, confidence = result2, 
                    level = a.level)
    class(results) <- "kap1"
    return(results)
  }
  if (nkappa > 1) {
    hot.test(kappa = kappa, var_kappa = var_kappa, a.level = a.level, 
             NN = K)
  }
  
}  


#'  Compare many (multilevel) Conger kappa coefficients using the delta method
#' 
#' This function performs Hotelling's T square test using a variance-covariance matrix based on the delta method to compare dependent Conger kappa coefficients
#' 
#' @param data a N x sum(Rg) matrix representing the classification of the N items by G groups of Rg observers (g=1,...,G). 
#' @param cluster_id a vector of lenght N with the identification number of the K clusters
#' @param nrat  a vector of lenght G indicating the number of observers in the G groups of observers
#' @param multilevel a binary indicator equal to TRUE in the presence of multilevel data and FALSE otherwise
#' @param a.level the significance level
#' @return $kappa a G x 2 matrix with the G kappa coefficients to be compared in the first column and their corresponding standard error in the second column
#' @return $T_test a vector of length 2 with the value of Hotelling's T square test as first element and the p-value as second element
#' @return $confidence confidence intervals for the pairwise comparisons of kappa coefficients
#' @return $var the G x G correlation matrix of the kappa coefficients 
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function compare several (multilevel) dependent Conger kappa coefficients. It uses Hotelling's T square test with the variance-covariance matrix obtained by the delta method. If only a single Conger kappa coefficient is computed, the kappa coefficient and its standard error are returned.
#' @references Vanbelle S. (2017) Comparing dependent agreement coefficients obtained on multilevel data. Biometrical Journal, 59 (5):1016-1034
#' @references Vanbelle S. (submitted) On the asymptotic variability of (multilevel) multirater kappa coefficients
#' @export
#' @importFrom stats aggregate na.omit
#' @examples 
#'  
#' #dataset (multilevel) (Vanbelle, 2008)
#' 
#' data(depression)
#' attach(depression)
#' delta.pair(data=cbind(diag,BDI,diag,GHQ),cluster_id=ID,weight='unweighted',multilevel=FALSE)
#' 
#' 
#' #dataset (multilevel) (Vanbelle, submitted)
#' data(CRACKLES)
#' attach(CRACKLES)
#' 
#'AGREEMENT<-matrix(NA,ncol=21,nrow=4)
#'
#'for (i in 1:7){
#'  AGREEMENT[1,((i-1)*3+1)]<-mean((rowSums(CRACKLES[UP==1,((i-1)*4+1):(i*4)])*(rowSums((CRACKLES[UP==1,((i-1)*4+1):(i*4)]))-1)+(4-rowSums(CRACKLES[UP==1,((i-1)*4+1):(i*4)]))*((4-rowSums((CRACKLES[UP==1,((i-1)*4+1):(i*4)])))-1))/12)
#'  AGREEMENT[1,((i-1)*3+2):(i*3)]<-delta.many2(cluster_id=patient[UP==1],data=CRACKLES[UP==1,((i-1)*4+1):(i*4)],nrat=c(4))$kappa
#'  AGREEMENT[2,((i-1)*3+1)]<-mean((rowSums(CRACKLES[LO==1,((i-1)*4+1):(i*4)])*(rowSums((CRACKLES[LO==1,((i-1)*4+1):(i*4)]))-1)+(4-rowSums(CRACKLES[LO==1,((i-1)*4+1):(i*4)]))*((4-rowSums((CRACKLES[LO==1,((i-1)*4+1):(i*4)])))-1))/12)
#'  AGREEMENT[2,((i-1)*3+2):(i*3)]<-delta.many2(cluster_id=patient[LO==1],data=CRACKLES[LO==1,((i-1)*4+1):(i*4)],nrat=c(4))$kappa
#'  AGREEMENT[3,((i-1)*3+1)]<-mean((rowSums(CRACKLES[UP!=1 & LO!=1,((i-1)*4+1):(i*4)])*(rowSums((CRACKLES[UP!=1 & LO!=1,((i-1)*4+1):(i*4)]))-1)+(4-rowSums(CRACKLES[UP!=1 & LO!=1,((i-1)*4+1):(i*4)]))*((4-rowSums((CRACKLES[UP!=1 & LO!=1,((i-1)*4+1):(i*4)])))-1))/12)
#'  AGREEMENT[3,((i-1)*3+2):(i*3)]<-delta.many2(cluster_id=patient[UP!=1 & LO!=1],data=CRACKLES[UP!=1 & LO!=1,((i-1)*4+1):(i*4)],nrat=c(4))$kappa
#'  AGREEMENT[4,((i-1)*3+1)]<-mean((rowSums(CRACKLES[,((i-1)*4+1):(i*4)])*(rowSums((CRACKLES[,((i-1)*4+1):(i*4)]))-1)+(4-rowSums(CRACKLES[,((i-1)*4+1):(i*4)]))*((4-rowSums((CRACKLES[,((i-1)*4+1):(i*4)])))-1))/12)
#'  AGREEMENT[4,((i-1)*3+2):(i*3)]<-delta.many2(cluster_id=patient,data=CRACKLES[,((i-1)*4+1):(i*4)],nrat=c(4))$kappa
#'}
#'
#'AGREEMENT2<-matrix(NA,ncol=14,nrow=4)
#'
#'for (i in 1:7)
#'{
#'  AGREEMENT2[,((i-1)*2+2)]<-paste(paste(paste(round(as.numeric(AGREEMENT[,((i-1)*3+2)]),2),' ('),round(as.numeric(AGREEMENT[,((i-1)*3+3)]),2),')'))
#'  AGREEMENT2[,((i-1)*2+1)]<-round(as.numeric(AGREEMENT[,((i-1)*3+1)]),2)
#'}
#' AGREEMENT2 #(table 3)



delta.many2<-function(data,cluster_id,nrat=c(6,6),multilevel=T,a.level=0.05)
{
  
  
  #CREATE THE COMPLETE CASE DATASET
  data_<-na.omit(cbind(data,cluster_id))
  data2<-data_[,1:ncol(data)]
  c_id<-data_[,ncol(data)+1]
  
  
  #CREATE CONSECUTIVE ID FOR THE CLUSTERS
  id1<-unique(c_id);
  idd2<-rep(0,length(id1));
  for (i in 1:length(id1)){idd2<-ifelse(c_id==id1[i],i,idd2)}
  
  #CREATE FACTORS WITH THE CLASSIFICATION OF THE OBSERVERS (not usefull)
  n_<-c(0,nrat)
  
  tab.l<-sapply(1:length(nrat),function(x){names(table(c(unlist(data2[,(1+sum(n_[1:x])):(sum(n_[1:(x+1)]))]))))},simplify=FALSE)
  ncat<-sapply(1:length(nrat),function(x){length(tab.l[[x]])},simplify=TRUE)
  
  K<-max(idd2)						    #number of clusters
  npair<-nrat*(nrat-1)/2  		#number of ditincts pairs within each group 
  nkappa<-length(nrat)        #number of kappas statistics
  npair_tot<-sum(npair)       #total number of pairs of observers
  
  #CREATE WEIGHT FOR THE ITEMS
  xi<-prop.table(table(idd2))
  
  
  #CREATE THE OBSERVED AGREEMENT VECTOR (each pair, each subject)
  tmp_<- sapply(1:nkappa,function(x){expand.grid((1+sum(n_[1:x])):(sum(n_[1:(x+1)])),(1+sum(n_[1:x])):(sum(n_[1:(x+1)])))},simplify=FALSE)
  tmp<-do.call(rbind,tmp_)
  tmp2<-tmp[tmp[,1]>tmp[,2],]
  RB<-tmp2[,1]#ID number of the rater 1 of pair i
  RA<-tmp2[,2]#ID number of the rater 2 of pair i
  
  a<-mapply(agree<-function(x,y){ifelse(data2[,x]==data2[,y],1,0)},tmp2[,2],tmp2[,1])
  
  Po<-unname(as.list(aggregate(a, by=list(idd2),FUN=mean, na.rm=TRUE)[,-1]))
  #if (nkappa>1){Po<-unname(as.list(aggregate(a, by=list(idd2),FUN=mean, na.rm=TRUE)[,-1]))}
  #if (nkappa==1){Po<-replicate(npair,matrix(0,ncol=K,nrow=1), simplify=FALSE);Po[[1]]<-unname(aggregate(a, by=list(idd2),FUN=mean, na.rm=TRUE)[,-1])}
  
  
  
  #CREATE THE MARGINALS AT THE ITEM LEVEL FOR EACH RATER
  nrat2<-rep(seq(1:length(nrat)),nrat)
  p<-sapply(1:sum(nrat),function(y){sapply(1:K,function(x){prop.table(table(factor(data2[, y][idd2==x], levels=tab.l[[nrat2[[y]]]])))})},simplify=F)
  
  
  #CREATE THE MATRIX OMEGA
  Omega<-(diag(K)-xi%*%matrix(1,ncol=K,nrow=1))%*%diag(xi^2)%*%(diag(K)-matrix(1,ncol=1,nrow=K)%*%t(xi))
  
  
  #CREATE THE MATRIX V
  V<-matrix(NA,ncol=(npair_tot+sum(nrat*ncat)),nrow=(npair_tot+sum(nrat*ncat)))
  
  #var-cov for observed agreement
  tmp = expand.grid(1:npair_tot,1:npair_tot)
  V[1:npair_tot,1:npair_tot]<-mapply(function(i, j) (K^2/(K-1))*matrix(Po[[i]],nrow=1)%*%Omega%*%t(matrix(Po[[j]],nrow=1)),tmp[,1],tmp[,2])
  
  #var-cov for observed agreement versus marginals
  tmp = expand.grid(1:sum(nrat),1:npair_tot)
  V[1:npair_tot,(npair_tot+1):(npair_tot+sum(nrat*ncat))]<-t(matrix(mapply(function(i, j)(K^2/(K-1))*matrix(Po[[i]],nrow=1)%*%Omega%*%t(p[[j]]),tmp[,2],tmp[,1]),ncol=npair_tot))
  
  #var-cov for marginals
  
  ncati<-c(0,rep(ncat,nrat))
  for (i in 1:sum(nrat)){
    for (j in (i:sum(nrat))){
      V[(npair_tot+1+sum(ncati[1:i])):(npair_tot+sum(ncati[1:(i+1)])),(npair_tot+1+sum(ncati[1:j])):(npair_tot+sum(ncati[1:(j+1)]))]<-(K^2/(K-1))*p[[i]]%*%Omega%*%t(p[[j]])
    }
  }
  
  for (i in ((npair_tot+1):(npair_tot+sum(nrat*ncat)))){
    for (j in 1:(i-1)){
      V[i,j]<-V[j,i]
    }
  }
  
  
  #CREATE THE MARGINALS OVERALL
  M<-sapply(1:sum(nrat),function(y){prop.table(table(factor(data2[, y],levels=tab.l[[nrat2[[y]]]])))},simplify=FALSE)
  
  #CREATE THE MATRIX J
  J<-matrix(0,2*npair_tot,(npair_tot+sum(nrat*ncat)))
  J[1:npair_tot, 1:npair_tot] <- diag(npair_tot)
  
  RAi<-c(0,RA)
  RBi<-c(0,RB)
  
  for (i in (1:npair_tot)){
    J[(npair_tot+i),(npair_tot+1+sum(ncati[1:RA[i]])):(npair_tot+sum(ncati[1:(RA[i]+1)]))]<-t(M[[RB[i]]])
    J[(npair_tot+i),(npair_tot+1+sum(ncati[1:RB[i]])):(npair_tot+sum(ncati[1:(RB[i]+1)]))]<-t(M[[RA[i]]])                                                  
  }
  
  #COMPUTE VARPSY
  varpsy<-J%*%V%*%t(J)
  
  #COMPUTE THE VECTOR OF EXPECTED AGREEMENT FOR EACH PAIR OF OBSERVER
  
  Pea<-matrix(NA,ncol=1,nrow=npair_tot)
  P0a<-matrix(NA,ncol=1,nrow=npair_tot)
  
  for (i in 0:(npair_tot-1)){
    
    P0a[i+1]<-sum(Po[[i+1]]*matrix(xi,ncol=K,nrow=1))               
    Pea[i+1]<-t(M[[RA[i+1]]])%*%M[[RB[i+1]]]
    
  }
  
  #COMPUTE THE VARIANCE-COVARIANCE MATRIX FOR THE PO and PE for each group
  npairi<-c(0,npair)
  
  L<-matrix(0,2*nkappa,2*npair_tot)
  for (i in 1:nkappa){
    L[i,(sum(npairi[1:i])+1):sum(npairi[1:(i+1)])]<-1/npair[i]
    L[i+nkappa,npair_tot+(sum(npairi[1:i])+1):sum(npairi[1:(i+1)])]<-1/npair[i]
    
  }
  
  var_bkappa<-L%*%varpsy%*%t(L)
  
  Pe<-matrix(NA,ncol=1,nrow=nkappa)
  P0<-matrix(NA,ncol=1,nrow=nkappa)
  
  for (i in 1:nkappa){
    
    P0[i]<-mean(P0a[nrat2==i])    
    Pe[i]<-mean(Pea[nrat2==i])
    
  }
  
  
  #COMPUTE THE VARIANCE-COVARIANCE MATRIX FOR THE KAPPAS
  
  s <- matrix(0, nkappa, 2 * nkappa)
  if (nkappa > 1) {
    diag(s[1:nkappa, 1:nkappa]) <- 1/(1 - Pe)
    diag(s[1:nkappa, (nkappa + 1):(2 * nkappa)]) <- (P0 - 1)/((1 - Pe)^2)
  }
  if (nkappa == 1) {
    s[1:nkappa, 1:nkappa] <- 1/(1 - Pe)
    s[1:nkappa, (nkappa + 1):(2 * nkappa)] <- (P0 - 1)/((1 - Pe)^2)}
  
  
  if (multilevel==T){
    var_kappa<-s%*%var_bkappa%*%t(s)/K
    
  }
  
  if (multilevel==F){
    var_kappa<-s%*%var_bkappa%*%t(s)*(K-1)/K^2#correction factor for one unit per cluster
    
  }
  
  
  kappa <- matrix((P0 - Pe)/(1 - Pe), ncol = 1)
  if (nkappa == 1) {
    bound <- c(NA, NA)
    bound[1] <- kappa - qnorm(1 - a.level/2) * sqrt(var_kappa)
    bound[2] <- kappa + qnorm(1 - a.level/2) * sqrt(var_kappa)
    result1 <- matrix(c(kappa, sqrt(var_kappa)), nrow = 1)
    colnames(result1) <- c("Kappa", "SE(Kappa)")
    rownames(result1) <- ""
    result2 <- matrix(bound, nrow = 1)
    colnames(result2) <- c("", "")
    rownames(result2) <- ""
    results <- list(kappa = result1, confidence = result2, 
                    level = a.level)
    class(results) <- "kap1"
    return(results)
  }
  if (nkappa > 1) {
    hot.test(kappa = kappa, var_kappa = round(var_kappa,10), a.level = a.level, 
             NN = K)
  }
  
  
}


#' Hotelling's T square test 
#' 
#' This function performs Hotelling's T square test when the vector of coefficients and the variance-covariance matrix are provided
#' 
#' @param kappa a vector of length G with the coefficients to be compared
#' @param var_kappa the G x G variance-covariance matrix of the coefficients. It can be provided in general or obtained with the \code{delta.pair}, the \code{boot.pair}, the \code{boot.many1}, the \code{boot.many2}, the \code{delta.many1} and the \code{delta.many2} functions in the kappa context. 
#' @param a.level the significance level
#' @param NN the number of clusters if the delta method is used to determine the variance-covariance matrix or the number of bootstrap iterations if the bootstrap method is used
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function performs Hotelling's T square test to compare dependent coefficients when the vector of coefficients and the variance-covariance matrix are provided. 
#' @return $kappa a G x 2 matrix with the coefficients in the first column and their corresponding standard error in the second column
#' @return $T_test a vector of length 2 with the value of Hotelling's T square test as first element and the corresponding p-value as second element
#' @return $confidence confidence intervals for the pairwise comparisons of the coefficients
#' @return $cor the G x G correlation matrix of the coefficients 
#' @export
#' @importFrom stats pchisq pf qf qnorm
#' @examples 
#'  
#' vect<-c(0.3,0.4,0.5)#vector of coefficients
#' v_c<-matrix(c(0.01,0.005,0.007,0.005,0.015,0.006,0.007,0.006,0.05),ncol=3)#variance-covariance matrix
#' hot.test(kappa=vect,var_kappa=v_c,a.level=0.05,NN=50)




hot.test<-function(kappa,var_kappa,a.level,NN){
  result4<-cov2cor(var_kappa)
  estimate<-cbind(kappa,diag(sqrt(var_kappa)))
   colnames(estimate)<-c("kappa","SE(kappa)")
  if (length(unique(kappa)) < length(kappa)){writeLines("Stop: The system is singular");return(list("kappa"=estimate,"cor"=result4));stop("The system is singular")}
  if (any(is.na(result4))|any(!is.na(diag(result4))==0) ){writeLines("Stop: There are missing values in the variance-covariance matrix of the kappa coefficients or null variances, Hotteling's test cannot be performed");
    return(list("kappa"=estimate,"cor"=result4));stop("Stop: There are missing values in the variance-covariance matrix of the kappa coefficients, Hotteling's test cannot be performed")}
  if (any(result4>(1+100 * .Machine$double.eps)) | any(result4<(-1-100 * .Machine$double.eps))){writeLines("Stop: The covariances provided are not compatible with the variances");return(list("kappa"=estimate,"cor"=result4));
    stop("The covariances provided are not compatible with the variances")}
  if (!isSymmetric(result4,tol = 100 * .Machine$double.eps)){writeLines("Stop: The variance-covariance matrix is not symmetric");return(list("kappa"=estimate,"cor"=result4)); stop("The variance-covariance matrix is not symmetric")}
  nkappa<-length(kappa)
  ncomb<-choose(nkappa,2) #number of pairwise comparisons
  labelk<-matrix(NA,nrow=ncomb,ncol=1)
  count<-0
 
  
  
  #MATRIX C
  C<-matrix(0,nkappa-1,nkappa)
  for (r in 1:(nkappa-1))
  {
    C[r,1]<-1
    C[r,r+1]<--1
  }
  
  CT<-t(C)
  f <- function(m) class(try(solve(m),silent=T))=="matrix"
  invert<-f(C%*%var_kappa%*%t(C))[1]
  if (invert==FALSE){stop("The variance-covariance matrix cannot be inversed")}
  if (invert==TRUE){
  T2<-t(C%*%kappa)%*%solve(C%*%var_kappa%*%t(C))%*%C%*%kappa
  quantile<-qf(1-a.level,nkappa-1,NN-nkappa+1)
  seuil<-quantile*(NN-1)*(nkappa-1)/(NN-nkappa+1)
  pvalue<-1-pf(T2*(NN-nkappa+1)/((NN-1)*(nkappa-1)), nkappa-1, NN-nkappa+1)
  
  
  Comb<-matrix(0,ncomb,nkappa)
  fin<-nkappa-1
  debut<-1
  r<-1
  p<-1
  for (p in 1:(nkappa-1))
  {
    for (i in debut:fin)
    {
      Comb[i,p]<-1
      Comb[i,r+p]<--1
      r<-r+1
    }
    debut<-fin+1
    fin<-fin+nkappa-(p+1)
    r<-1
  }
  
  limit<-matrix(0,ncomb,2)
  for (i in 1:(ncomb))
  {
    limit[i,1]<-Comb[i,]%*%kappa-sqrt(seuil)*sqrt(Comb[i,]%*%var_kappa%*%matrix(Comb[i,],nkappa,1))
    limit[i,2]<-Comb[i,]%*%kappa+sqrt(seuil)*sqrt(Comb[i,]%*%var_kappa%*%matrix(Comb[i,],nkappa,1))
  }
  

    for (m in 1:(nkappa-1)){
      for (p in (m+1):nkappa){
        count<-count+1
        labelk[count,1]<-c(paste("kappa",paste(m,paste("vs",paste("kappa",p)))))
      }
    }
    
    labelkappa<-matrix(0,nkappa,1)
    for (m in 1:nkappa){
      labelkappa[m,1]<-c(paste("kappa",m))
    }
    
    result1<-matrix(cbind(kappa,sqrt(diag(var_kappa))),ncol=2)
    colnames(result1) <- c("Kappa","SE")
    rownames(result1) <- labelkappa
    
    result2<-matrix(c(T2,pvalue),ncol=2)
    colnames(result2) <- c("T2","pvalue")
    rownames(result2) <- c("")
    
    result3<-matrix(limit,ncol=2)
    colnames(result3) <- c("Limit inf","Limit sup")
    rownames(result3) <- labelk
    
  colnames(result4) <- labelkappa
  rownames(result4) <- labelkappa
    
    results<-list("kappa"=result1,"T_test"=result2,"confidence"=result3,"cor"=result4)
    
    class(results) <-"kap"
    return(results)
  }
  
  
  
}




#' print objects of the class kap1
#' 
#' @param x Kappa coefficient and 95\% confidence interval
#' @param ... further arguments
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function prints kappa coefficients and 95\% confidence intervals
#' @export
#' @examples  
#' #dataset (Vanbelle and Albert, 2008)
#' 
#' data(depression)
#' attach(depression)
#' a<-boot.pair(data=cbind(diag,BDI),cluster_id=ID,weight="unweighted",ITN=2000)
#' print(a)
#' 
#
print.kap1 <- function(x,...) {
  cat("Value of the kappa coefficients\n\n")
  print(round(x$kappa,3))
  cat("\n")
  cat(100*(1-x$level),"% Confidence interval\n")
  print(round(x$confidence,3))
  writeLines("\n \n Warning: on binary scales, only the intraclass kappa coefficient measures reliability.\n On nominal scales, the intraclass kappa coefficient has to be determined for each category separately.\n For ordinal scales, only the quadratic weighted kappa coefficient is a reliability measure.\n Other types of kappa coefficients are relative agreement measures.")
}


#' print objects of the class kap
#' 
#' @param x outputs from the Hotelling's T^2 test
#' @param ... further arguments
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function prints the outputs from Hotelling's T square test
#' @export
#' @examples  
#' #dataset (not multilevel) (Vanbelle and Albert,2008)
#' 
#' set.seed(103) #to get the same results as in the paper
#' data(depression)
#' attach(depression)
#' a<-boot.pair(data=cbind(diag,BDI,diag,GHQ),cluster_id=ID,weight='unweighted',ITN=1000)
#' print(a)
#' 
#



print.kap <- function(x,...) {
  cat("Value of the kappa coefficients\n\n")
  print(round(x$kappa,3))
  cat("\n")
  cat("Hotellings T^2\n\n")
  print(x$T_test)
  cat("\n")
  cat("Pairwise comparisons, 95% CI\n\n")
  print(round(x$confidence,3))
  cat("\n")
  cat("Correlation matrix between the kappa statistics\n\n")
  print(round(x$cor,3))
  writeLines("\n \n Warning: on binary scales, only the intraclass kappa coefficient measures reliability.\n On nominal scales, the intraclass kappa coefficient has to be determined for each category separately.\n For ordinal scales, only the quadratic weighted kappa coefficient is a reliability measure.\n Other types of kappa coefficients are relative agreement measures.")
}

#' 
#' Compare dependent pairwise kappas (bootstrap method)
#' 
#' This function performs Hotelling's T square test using a variance-covariance matrix based on the bootstrap method to compare dependent pairwise kappa coefficients
#' 
#' @param cluster_id a vector of lenght N with the identification number of the clusters
#' @param data a N x R matrix representing the classification of the N items by the R observers. The kappa coefficients are computed between column (1,2), (3,4), etc....
#' @param weight the weighting scheme to be used for kappa coefficients. 'unweighted' for Cohen's kappa, 'equal' for linear weights and 'squared' for quadratic weights
#' @param a.level significance level
#' @param ITN the number of bootstrap iterations
#' @param summary_k, if true, Hotteling's T square test is performed, if false, only the bootstraped kappa coefficients are returned
#' @return $kappa a G x 2 matrix with the G kappa coefficients to be compared in the first column and their corresponding standard error in the second column
#' @return $T_test a vector of length 2 with the value of Hotelling's T square test as first element and the corresponding p-value as second element
#' @return $confidence confidence intervals for the pairwise comparisons of the kappa coefficients
#' @return $cor the G x G correlation matrix of the kappa coefficients
#' @return $K when summary_k is false, the ITN x G matrix with the bootstrapped kappa coefficients 
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function compares several dependent pairwise kappa coefficients using Hotelling's T square with the variance-covariance matrix obtained by the bootstrap method. If only one kappa is computed, it returns the estimate and confidence interval.
#' @export
#' @importFrom irr kappa2
#' @importFrom stats cor qnorm var
#' @references Vanbelle S. and Albert A. (2008). A bootstrap method for comparing correlated kappa coefficients.  Journal of Statistical Computation and  Simulation, 1009-1015
#' @references Vanbelle S. Comparing dependent agreement coefficients obtained on multilevel data. submitted
#' @references Vanbelle S. (2014) A New Interpretation of the Weighted Kappa Coefficients. Psychometrika. Advance online publication.  doi: 10.1007/s11336-014-9439-4
#' 
#' @examples  
#' #dataset (not multilevel) (Vanbelle and Albert, 2008)
#' 
#' set.seed(103) #to get the same results as in the paper
#' data(depression)
#' attach(depression)
#' a<-boot.pair(data=cbind(diag,BDI,diag,GHQ),cluster_id=ID,weight='unweighted')
#' 
#' 
#' #dataset (multilevel) (Vanbelle, xxx)
#' 
#' data(FEES)
#' attach(FEES)
#' dat<-cbind(val_CO,val_COR,val_MH,val_MHR,val_TB,val_TBR) #formating the data matrix
#' boot.pair(data=dat,cluster_id=subject,weight='equal',summary_k=FALSE)
#' 


boot.pair<-function(cluster_id,data,weight='equal',a.level=0.05,ITN=1000,summary_k=T){
  
  id1 <- unique(cluster_id)
  idd2 <- rep(0, length(id1))
  
  for (i in 1:length(id1)) {
    idd2 <- ifelse(cluster_id == id1[i], i, idd2)
  }
  
 
  ID_ind<-tapply(1:nrow(data),idd2,function(x){x})
  
  nk<-sapply(ID_ind,length)
  
  ncluster <- max(idd2)
  nrater <- ncol(data)
  
  nobs <- length(idd2)
  nkappa<-nrater/2
  nk<-as.vector(table(idd2))

  K<-matrix(NA,ncol=nkappa,nrow=ITN)
  
  labelk<-matrix(NA,ncol=1,nrow=choose(nkappa,2))
  
  for (i in 1:ITN)
  {
    
    num<-sample(x=1:ncluster,size=ncluster,replace=T)
    
    nknew <- nk[num]
    Nnew <- sum(nk[num])
    rex <- matrix(NA, ncol = nrater, nrow = Nnew)
    all_boot_IND<-rep(NA,Nnew)
    s<-1
    endIND=0
    for (s in 1:length(num)){
      
      startIND=endIND+1
      endIND=startIND+nk[num[s]]-1
      all_boot_IND[startIND:endIND]=ID_ind[[num[s]]]
    }
    rex<-data[all_boot_IND,]
    deg<-sapply(1:nrater,f1<-function(j){length(table(rex[,j]))})
    if (nkappa==1){
      K[i,]<-ifelse(deg[1]==1|deg[2]==1,NA,kappa2(rex, weight = weight)$value)
    }
    if (nkappa>1){
    K[i,]<-sapply(seq(0:(nkappa-1))-1,k_kap<-function(idd){ifelse((deg[((2 *idd) + 1)]==1|deg[((2 *idd) + 2)]==1),NA,kappa2(rex[, ((2 *idd) + 1):((2 *idd) + 2)], weight = weight)$value)})}
  }
  
  na.kappa<-sum(is.na(K))
 if (na.kappa>0){writeLines(paste("\n \n Warning: the number of missing kappa estimates is ",na.kappa))}
 
var_kappa<-var(K,na.rm=T)
kappa<-matrix(colMeans(K,na.rm=T),nkappa,1)

if (summary_k==F){return(K)}
if (nkappa==1 & summary_k==T){
  
  bound<-c(NA,NA)
  bound[1]<-kappa-qnorm(1-a.level/2)*sqrt(var_kappa)
  bound[2]<-kappa+qnorm(1-a.level/2)*sqrt(var_kappa)
  
  result1<-matrix(c(kappa,sqrt(var_kappa)),nrow=1)
  result2<-matrix(bound,nrow=1)
  
  colnames(result2) <-c('','')
  rownames(result2) <-''
  
 
   
    colnames(result1) <-c('Kappa','SE')
    rownames(result1) <-""
  
    results<-list('kappa'=result1,'confidence'=result2,'level'=a.level)
    class(results)<-'kap1'
    return(results)  
    
  
}
  if (nkappa>1 & summary_k==T)  {hot.test(kappa=kappa,var_kappa=var_kappa,a.level=a.level,NN=ncluster)}


}

#'   Compare dependent multi-observers kappas (Conger, light, fleiss) with the bootstrap method (bootstrap method)
#' 
#'   This function performs Hotelling's T square test using a variance-covariance matrix based on the bootstrap method to compare dependent multi-observers kappa coefficients
#' 
#' @param cluster_id a vector of lenght N with the identification number of the clusters
#' @param data a N x R matrix representing the classification of the N items by the R observers. 
#' @param group a vector with G elements indicating how many observers are considered to compute each kappa coefficient. For example, c(3,5) means that a kappa coefficient between the first 3 observersand a kappa coefficient between the 5 last observers will be computed 
#' @param method The type of kappa to be computed. 'light' for Light's kappa coefficient, 'conger' for Conger's kappa coefficient and 'fleiss' for Fleiss' kappa coefficient. In case of Fleiss kappa, it works only if each item is rated by the same number of observers.
#' @param a.level the significance level
#' @param ITN the number of bootstrap iterations
#' @param summary_k, if true, Hotteling's T square test is performed, if false, only the bootstraped coefficients are returned
#' @return $kappa a G x 2 matrix with the kappa coefficients in the first column and their corresponding standard error in the second column
#' @return $T_test a vector of length 2 with the value of Hotelling's T test as first element and the corresponding p-value as second element
#' @return $confidence confidence intervals for the pairwise comparisons of the measures
#' @return $cor the G x G correlation matrix for the kappa coefficients
#' @return $K when summary_k is false, the ITN x G matrix with the bootstrapped kappa coefficients is returned
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function compare several kappa coefficients for several observers using Hotelling's T square with the variance-covariance matrix obtained by the bootstrap method. If only one kappa is computed, it returns the estimate and confidence interval. 
#' @export
#' @importFrom irr kappam.light kappam.fleiss
#' @importFrom stats cor qnorm var
#' @references Conger A.J. (1980). Integration and generalization of kappas for multiple raters. Psychological Bulletin 88, 322-328.
#' @references Fleiss J.L. (1971). Measuring nominal scale agreement among many raters. Psychological Bulletin 76, 378-382.
#' @references Light R.J. (1971). Measures of response agreement for qualitative data: some generalizations and alternatives. Psychological Bulletin 76, 365-377.
#' @references Vanbelle S. and Albert A. (2008). A bootstrap method for comparing correlated kappa coefficients.  Journal of Statistical Computation and  Simulation, 1009-1015
#' @references Vanbelle S. Comparing dependent agreement coefficients obtained on multilevel data. submitted

#' @examples  
#' #dataset (not multilevel) (Vanbelle and Albert, 2008)

#' data(depression)
#' attach(depression)
#' a<-boot.many2(data=cbind(diag,BDI,GHQ,BDI,GHQ),cluster_id=ID,method='light',group=c(3,2),summary_k=TRUE)

 

boot.many2<-function(cluster_id,data,group,method='fleiss',a.level=0.05,ITN=1000,summary_k=T){
  id1<-unique(cluster_id);
  idd2<-rep(0,length(id1));
  for (i in 1:length(id1)){idd2<-ifelse(cluster_id==id1[i],i,idd2)}
  
  nrater<-ncol(data)
  ncluster<-max(idd2)
  nobs<-length(idd2)
  nkappa<-length(group)
  
  nk<-as.vector(table(idd2))
  
  K<-matrix(NA,ncol=nkappa,nrow=ITN)
  
  labelk<-matrix(NA,ncol=1,nrow=choose(nkappa,2))
  
  
  for (i in 1:ITN)
  {
    
    
    num<-sample(unique(as.vector(idd2)),replace=T)
    nknew<-nk[num]
    Nnew<-sum(nknew)
   
    rex<-matrix(NA,ncol=nrater,nrow=Nnew)
    
    count<-1
    for (j in 1:length(num)) 
    {
      
      rex[(count:(count+nknew[j]-1)),]<-data[idd2==num[j],]
      count<-count+nk[num[j]]  			 
    }
    
    count<-0
    for (m in 1:(nkappa))
    {
      count<-count+1
      start<-1
      if (m>1){start<-sum(group[1:(m-1)])+1}
      
        
        if (method=="light"){
        K[i,count]<-kappam.light(rex[,start:sum(group[1:m])])$value
        }
        if (method=="fleiss"){
          K[i,count]<-kappam.fleiss(rex[,start:sum(group[1:m])],exact=FALSE)$value
        }
        if (method=="conger"){
          K[i,count]<-kappam.fleiss(rex[,start:sum(group[1:m])],exact=TRUE)$value
        }
        

    }
    
  }
  
  
  var_kappa<-var(K,na.rm=T)
  kappa<-matrix(colMeans(K,na.rm=T),nkappa,1)
  if (summary_k==F){return(K)}
  if (nkappa==1 & summary_k==T){
    
    bound<-c(NA,NA)
    bound[1]<-kappa-qnorm(1-a.level/2)*sqrt(var_kappa)
    bound[2]<-kappa+qnorm(1-a.level/2)*sqrt(var_kappa)
    
    result1<-matrix(c(kappa,sqrt(var_kappa)),nrow=1)
    result2<-matrix(bound,nrow=1)
    
    colnames(result2) <-c("","")
    rownames(result2) <-""
    
 
      
      colnames(result1) <-c("Kappa","SE")
      rownames(result1) <-""
      
      results<-list("kappa"=result1,"confidence"=result2,"level"=a.level)
      class(results)<-"kap1"
      return(results)   
    
  }
  if (nkappa>1 & summary_k==T){hot.test(kappa=kappa,var_kappa=var_kappa,a.level=a.level,NN=ncluster)}
}



#'   Compare dependent (multilevel) Fleiss kappas using the clustered bootstrap method
#' 
#'   This function performs Hotelling's T square test using a variance-covariance matrix based on the bootstrap method to compare dependent multi-observers kappa coefficients
#' 
#' @param cluster_id a vector of lenght N with the identification number of the clusters
#' @param data a N x sum(ncat_g) matrix representing the classification of the N items by the observers in group g in the ncat_g categories. For each group, the number of categories can vary 
#' @param ncat a vector with G elements indicating how many categories are considered to compute each kappa coefficient. For example, c(3,5) means that the three first columns correspond to the classification of subjects on a 3 categorical scale by a group of observers and the five last columns correspond to the classification of subjects on a 5 categorical scale by a group of observers.
#' @param a.level the significance level
#' @param ITN the number of bootstrap iterations
#' @param summary_k, if true, Hotteling's T square test is performed, if false, only the bootstraped coefficients are returned
#' @return $kappa a G x 2 matrix with the kappa coefficients in the first column and their corresponding standard error in the second column
#' @return $T_test a vector of length 2 with the value of Hotelling's T test as first element and the corresponding p-value as second element
#' @return $confidence confidence intervals for the pairwise comparisons of the measures
#' @return $cor the G x G correlation matrix for the kappa coefficients
#' @return $K when summary_k is false, the ITN x G matrix with the bootstrapped kappa coefficients is returned
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function compare several Fleiss kappa coefficients using Hotelling's T square with the variance-covariance matrix obtained by the clustered bootstrap method. If only one kappa is computed, it returns the estimate and confidence interval. 
#' @export
#' @importFrom stats cor qnorm var
#' @references Fleiss J.L. (1971). Measuring nominal scale agreement among many raters. Psychological Bulletin 76, 378-382.
#' @references Vanbelle S. (2017) Comparing dependent agreement coefficients obtained on multilevel data. Biometrical Journal, 59 (5):1016-1034
#' @references Vanbelle S. (submitted) On the asymptotic variability of (multilevel) multirater kappa coefficients
#' 

#' @examples  
#' #dataset (not multilevel) (Fleiss, 1971)

#' data(fleiss_psy)
#' attach(fleiss_psy)
#' set.seed(123) # to have the same results than in Vanbelle (submitted)
#' a<-boot.many1(data=fleiss_psy[,2:6],cluster_id=fleiss_psy[,1],ncat=c(5),a.level=0.05,ITN=1000)

boot.many1<-function(cluster_id,data,ncat=2,a.level=0.05, ITN = 1000,summary_k=T) 
{
  ncluster <- length(unique(cluster_id))
  nkappa <- ncol(data)/ncat
  
  K <- matrix(NA, ncol = nkappa, nrow = ITN)
  labelk <- matrix(NA, ncol = 1, nrow = choose(nkappa, 2))
  
  ID_ind<-sapply(1:ncluster, function(x) which(as.character(cluster_id) %in% x ) )
  
  for (i in 1:ITN) {
    
    num<-sample(x=1:ncluster,size=ncluster,replace=T)
    
    if (length(cluster_id)>ncluster){ind2<-sapply(1:ncluster, function(x) ID_ind[[num[x]]] )}
    if (length(cluster_id)==ncluster){ind2<-ID_ind[num]}
    
    rex<- data[unlist(ind2),]
    
    if (nkappa == 1) {
      K[i,] <- fleiss.sum(rex)
    }
    if (nkappa > 1) {
      K[i, ] <- sapply(seq(0:(nkappa - 1))-1,k_kap <- function(idd)fleiss.sum(rex[, ((ncat*idd) +1):(ncat*(idd+1))]))
    }
  }
  
  var_kappa <- var(K, na.rm = T)
  kappa <- matrix(colMeans(K, na.rm = T), nkappa, 1)
  if (summary_k == F) {
    return(K)
  }
  if (nkappa == 1 & summary_k == T) {
    bound <- c(NA, NA)
    bound[1] <- kappa - qnorm(1 - a.level/2) * sqrt(var_kappa)
    bound[2] <- kappa + qnorm(1 - a.level/2) * sqrt(var_kappa)
    result1 <- matrix(c(kappa, sqrt(var_kappa)), nrow = 1)
    result2 <- matrix(bound, nrow = 1)
    colnames(result2) <- c("", "")
    rownames(result2) <- ""
    colnames(result1) <- c("Kappa", "SE")
    rownames(result1) <- ""
    results <- list(kappa = result1, confidence = result2, 
                    level = a.level)
    class(results) <- "kap1"
    return(results)
  }
  if (nkappa > 1 & summary_k == T) {
    hot.test(kappa = kappa, var_kappa = var_kappa, a.level = a.level,NN = ncluster)
  }
}








#' Compare G independent pairwise kappa coefficients
#' 
#' This functions compare L independent kappa coefficients using the method of Fleiss (1981)
#' 
#' @param cluster_id a vector with the identification number of the clusters
#' @param rater1 a vector with the ratings of the first observer
#' @param rater2 a vector with the ratings of the second observer
#' @param meth the method to be used to compute the standard error of the kappa coefficients: 'delta' for the delta method and 'boot' for the bootstrap method. 
#' @param cov the covariate determining the L groups to be compared 
#' @param weight the weighting scheme to be used for kappa coefficients. 'unweighted' for Cohen's kappa, 'equal' for linear weights and 'squared' for quadratic weights
#' @param a.level the significance level
#' @param multilevel a binary indicator equal to TRUE if the data are multilevel and FALSE otherwiwse. 
#' @param ITN number of bootstrap iterations needed if the bootstrap procedure is chosen
#' @return $kappa the value of the L kappa coefficients and their standard error
#' @return $chi the value of the chi-squared statistic with L-1 degrees of freedom
#' @return $p the p-value 
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function compare L independent kappa coefficients using the method of Fleiss (1981). The data have to be entered in a vertical format.
#' @export
#' @importFrom utils write.table
#' @references Fleiss, J. L. (1981). Statistical methods for rates and proportions (2nd ed.). New York: John Wiley
#' @references Vanbelle S. Comparing dependent agreement coefficients obtained on multilevel data. submitted
#' @references Vanbelle S. (2014) A New Interpretation of the Weighted Kappa Coefficients. Psychometrika. Advance online publication.  doi: 10.1007/s11336-014-9439-4
#' @examples  
#'  
#' #dataset (multilevel) (Vanbelle, xxx)
#' 
#' data(FEES)
#' attach(FEES)
#' 
#' fleiss.pair(rater1=val_CO,rater2=val_COR,cluster_id=subject,
#' weight="unweighted",multilevel=TRUE,meth='delta',cov=group)






fleiss.pair<-function(rater1,rater2,cluster_id,weight="equal",multilevel=T,a.level=0.05,cov,ITN=1000,meth)
{
  if (a.level>1 | a.level<0) stop("the significance level should be between 0 and 1")
 
 
  cov.f <- factor(cov)
  cov.f2<-cov.f
  ncat<-nlevels(cov.f)
  levels(cov.f) <- seq(from=1,to=ncat)
  kap<-matrix(NA,ncol=2,nrow=ncat)
 for (i in 1:ncat){
if (meth=='delta'){kap[i,]<-delta.pair(data=cbind(rater1[cov.f==i],rater2[cov.f==i]),cluster_id=cluster_id[cov.f==i],weight=weight,multilevel=multilevel,a.level=a.level)$kappa}
if (meth=='boot'){kap[i,]<-boot.pair(data=cbind(rater1[cov.f==i],rater2[cov.f==i]),cluster_id=cluster_id[cov.f==i],weight=weight,ITN=ITN,a.level=a.level)$kappa}
}
  
  wg<-1/(kap[,2]*kap[,2])
  k_ass<-sum(wg*kap[,1])/sum(wg)
  
  chi_tot<-sum(wg*kap[,1]*kap[,1])
  chi_ass<-sum(wg*kap[,1])^2/sum(wg)
  chi_hom<-chi_tot-chi_ass
  p.val<-1-pchisq(chi_hom,df=length(kap[,1])-1)
  

resultp1<-as.table(rbind(c("Category","Kappa","SE"),cbind(levels(cov.f2),round(kap[,1],3),round(kap[,2],3))))
cat("\n")
cat("Value of the kappa coefficients\n")
cat("\n")
write.table(format(resultp1, justify="right"),row.names=F,col.names=F, quote=F)

cat("\n")
cat("Chi-square test\n")
cat("\n")
resultp2<-as.table(rbind(c("chi^2","p-value"),cbind(round(chi_hom,3),round(p.val,5))))
write.table(format(resultp2, justify="right"),row.names=F,col.names=F, quote=F)

return(list("kappa"=kap,"Test"=c(chi_hom,p.val)))

writeLines("\n \n Warning: on binary scales, only the intraclass kappa coefficient measures reliability.\n On nominal scales, the intraclass kappa coefficient has to be determined for each category separately.\n For ordinal scales, only the quadratic weighted kappa coefficient is a reliability measure.\n Other types of kappa coefficients are relative agreement measures.")

  
}
  
#' This function plots a confidence ellipse when 3 kappas are compared with the clustered bootstrap method
#'
#' @param boot_val a ITN x 3 matrix with the bootstrapped values of the 3 kappa coefficients obtained with  \code{boot.pair}
#' @param a.level the significance level to draw the confidence ellipse
#' @param xlim the limits  for the x-axis
#' @param ylim the limits  for the y-axis
#' @param xlab the label for the x-axis
#' @param ylab the label for the y-axis
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function plots differences between the bootstrapped kappa coefficients and a confidence ellipse. The triangle represents the coordinate (0,0) and the square the mean of the two differences
#' @export
#' @importFrom ellipse ellipse 
#' @importFrom graphics lines plot points
#' @examples
#' #dataset (multilevel) (Vanbelle, xxx)
#'
#' data(FEES)
#' attach(FEES)
#' dat<-cbind(val_CO,val_COR,val_MH,val_MHR,val_TB,val_TBR) #formating the data matrix
#' Kappa_val<-boot.pair(data=dat,cluster_id=subject,weight='equal',summary_k=FALSE)
#' ell_graph(boot_val=Kappa_val,xlim=c(-0.8,0.4),ylim=c(-0.2,0.4),xlab="K(obs1)-K(cons)",ylab="K(obs2)-K(cons)")


ell_graph<-function(boot_val,a.level=0.05,xlim=c(-1,1),ylim=c(-1,1),ylab,xlab){
  ITN<-nrow(boot_val)
  quantile<-qf(1-a.level,3-1,ITN-3+1)
  seuil<-quantile*(ITN-1)*(3-1)/(ITN-3+1)
  d1<-boot_val[,1]-boot_val[,3];
  d2<-boot_val[,2]-boot_val[,3];
  correlation<-cor(d1,d2);
  plot(d1,d2,pch=20,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,col=8);
  lines(ellipse(correlation,c(sqrt(var(d1)),sqrt(var(d2))),c(mean(d1),mean(d2)),0.95,sqrt(seuil)),type="l",col=1,lwd=3);
  points(mean(d1), mean(d2), pch=15,col=1,cex=1.2);
  points(0,0, pch=17,col=1,cex=1.2);
}






#' return the value of Fleiss kappa coefficient
#' 
#' @param data NxK dataset to compute Fleiss kappa coefficient. For each of the N subjects (row), the number of observers classifying the subject in the K-categories of the scale (columns)
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function compute Fleiss kappa coefficient
#' @return the value of Fleiss kappa coefficent
#' @export

fleiss.sum<-function(data){
  
  N<-nrow(data)            #number of subjects
  ncat<-ncol(data)         #number of categories
  nrat<-rowSums(data)      #number of raters per subject
  idd2<-seq(1,nrow(data))  #subject ID
  
  #OBSERVED AGREEMENT
  oprim_i<-rowSums(data*(data-1))/(nrat*(nrat-1))
  P_o<-mean(oprim_i)
  
  #MARGINAL PROBABILITY DISTRIBUTION
  p_j<-colMeans(sweep(data,1,nrat,"/"))
  eprim_i<-rowSums(sweep(data,2,p_j,"*"))/nrat
  
  #EXPECTED AGREEMENT
  P_e1<-mean(eprim_i)		

  #FLEISS KAPPA COEFFICIENT
  kappa_1<-(P_o-P_e1)/(1-P_e1)
  return(c(kappa_1))
  
}



#' convert a variance-covariance matrix in a correlation matrix
#' 
#' @param V a G x G matrix
#' @author Sophie Vanbelle \email{sophie.vanbelle@@maastrichtuniversity.nl}
#' @details This function covert a variance-covariance matrix in a correlation matrix
#' @return a G x G correlation matrix
#' @export


cov2cor<-function (V) 
{
  D<-diag(1/(sqrt(diag(V))))
  CO<-D %*% V %*%D
  
  if (any(is.na(CO))| any(is.na(CO))) {return(CO);stop("There are some missing value in the vaiance-covariance matrix")}
  if (any(!is.na(CO)>(1+10*.Machine$double.eps)) | any(!is.na(CO)<(-1-10*.Machine$double.eps))){return(CO); stop("The covariances you provided are not compatible with the variances")}
  if (!isSymmetric(CO,tol = 10 * .Machine$double.eps)) stop("The variance-covariance matrix is not symmetric")
  return(CO)
}


#' Diagnosis of depression
#'
#' A dataset containing the diagnosis of depression of 50 subjects on three binary scales: diagnosis, BDI and GHQ.  The variables are as follows:
#' \itemize{
#'   \item ID. subject ID (1--50)
#'   \item diag. depression diagnosis (0: no, 1: yes)
#'   \item BDI. depression diagnosis following BDI scale (0: no, 1: yes)
#'   \item GHQ. depression diagnosis following PHQ scale (0: no, 1: yes)
#' }
#'
#' @format A data frame with 50 rows and 4 variables
#' @name depression
#' @docType data
NULL

#' FEES
#'
#' A dataset containing the valleculae score of 20 subjects obtained with thin and thick liquid.  The variables are as follows:
#' \itemize{
#'   \item subject. subject ID (1--37)
#'   \item swallow. type of swallowing (1=thin liquid, 2=thick liquid)
#'   \item val_TB. valleculae score given by rater TB on occasion  1 (0: no, 1: <50\%, 2: >=50\%)
#'   \item val_TBR. valleculae score given by rater TB on occasion  2 (0: no, 1: <50\%, 2: >=50\%)
#'   \item val_MH. valleculae score given by rater MH on occasion  1 (0: no, 1: <50\%, 2: >=50\%)
#'   \item val_MHR. valleculae score given by rater MH on occasion  2 (0: no, 1: <50\%, 2: >=50\%)
#'   \item val_CO. valleculae score given in consensus on occasion  1 (0: no, 1: <50\%, 2: >=50\%)   
#'   \item val_COR. valleculae score given in consensus on occasion  2 (0: no, 1: <50\%, 2: >=50\%)
#' }
#'
#' @format A data frame with 40 rows and 8 variables
#' @name FEES
#' @docType data
NULL

#Load the data
# fleiss_psy<-matrix(scan('C:/Users/sophie.vanbelle/surfdrive/fleiss/fleiss.dat'),ncol=6,byrow=T)
# colnames(fleiss_psy)<-c("subject","cat1","cat2","cat3","cat4","cat5")
# 
# fleiss_psy<-data.frame(fleiss_psy)
# 
# save(fleiss_psy,file='fleiss_psy.rda')

#load("fleiss_psy.rda")


#' Psychiatric diagnosis (Fleiss)
#'
#' A dataset containing the diagnosis of depression of 50 subjects on three binary scales: diagnosis, BDI and GHQ.  The variables are as follows:
#' \itemize{
#'   \item subject. subject ID (1--30)
#'   \item cat1. number of observers classifying the subject in category 1 (depression)
#'   \item cat2. number of observers classifying the subject in category 2 (personality disorder)
#'   \item cat3. number of observers classifying the subject in category 3 (schizophrenia)
#'   \item cat4. number of observers classifying the subject in category 4 (neurosis)
#'   \item cat5. number of observers classifying the subject in category 5 (other psychiatric disorder)
#' }
#'
#' @format A data frame with 30 rows and 6 variables
#' @name fleiss_psy
#' @docType data
NULL


#' CRACKLES (Tromso study)
#'
#' A dataset containing the assessment of the presence of crackles in 20 subjects by 28 observers. The variables are as follows:
#' \itemize{
#'   \item EXP1...4. Classification made by experts 1...4
#'   \item NOR1...4. Classification made by general practitioners 1...4 from Norway 
#'   \item RUS1...4. Classification made by general practitioners 1...4 from Russia 
#'   \item WAL1...4. Classification made by general practitioners 1...4 from Wales 
#'   \item NLD1...4. Classification made by general practitioners 1...4 from The Netherlands 
#'   \item PUL1...4. Classification made by pulmonologists 1...4
#'   \item STU1...4. Classification made by students 1...4
#'   \item patient. subject ID (1--20)
#'   \item UP. indicator variable. 1: examination of the upper posterior thorax; 0: otherwise
#'   \item LO. indicator variable. 1: examination of the lower posterior thorax; 0: otherwise  
#' }
#'
#' @format A data frame with 120 rows and 31 variables
#' @name CRACKLES
#' @docType data
NULL