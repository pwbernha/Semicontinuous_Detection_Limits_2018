##############GenIMS Application Code#####################

###################Necessary Packages########################
library(evd)
library(flexsurv)
library(mcmc)
library(brglm)
library(flexmix)
library(numDeriv)


####################Data Read and Clean-Up####################

#Reading in  GenIMS longitudinal data from laptop
GenIMS <- read.csv("C://Users/pbernh01/Documents/Data.csv")

#Reading in GenIMS single time indepdent data from laptop
GenIMSInd <- read.csv("C://Users/pbernh01/Documents/Data2.csv")

#Attaching sex, race, apache, and age covariates from time independent data and restricting 
#to first day measurements (note apache not used in final analysis)
GenIMS1 <- cbind(GenIMS[GenIMS[,2]==1,],GenIMSInd[,c(4,5,12,69)])

#Removing patients who did not actually have CAP 
GenIMS1 <- GenIMS1[GenIMSInd[,68]==1,]

#Removing patients who were immediately discharged
GenIMS1 <- GenIMS1[GenIMS1[,66]=="INPATIENT",]

#Defining Response (survival time and survival indicator)
T <- as.numeric(GenIMS1[,60])
T[which(T>90)]<- 90
cutoff <- 90
Y <- T

#Defining Covariate matrix: covariate 1 is TNF, covariate 2 is IL-6, covariate 3 is IL-10, 
#covariate 4 is sex (1=male), covariate 5 is race, covariate 6 is apache score, covariate 7 is age
X <- cbind(1, as.numeric(GenIMS1[,22]),as.numeric(GenIMS1[,23]),as.numeric(GenIMS1[,24]),as.numeric(GenIMS1[,82]-1),as.numeric(GenIMS1[,83]),as.numeric(GenIMS1[,84]), as.numeric(GenIMS1[,85]))

#Eliminating those without first day observations (should not affect analysis greatly since 
#lack of first day observation occured only when patient arrived on a weekend or holiday)
ResCov <- cbind(Y,X)
ResCov <- na.omit(ResCov)
Y <- ResCov[,1]
X <- ResCov[,2:9]

#Defining n=number of individuals (1418)
n <- length(ResCov[,1])

#creating an n x 3 matrix of the detection limits for convenience below (all rows are the same except for those individuals with censoring at 2 for IL-6
d <- c(4, 5, 5)
dmatrix <- matrix(log(d),n,3,byrow=TRUE)
dmatrix[which(X[,3]==-2) ,2] <- log(2)

#Taking log tranformation of covariates to make them approximately normal; defining race as a binary variable
for(i in 1:n){
  if(X[i,2]>0) X[i,2] <- log(X[i,2]) else X[i,2] <- log(-X[i,2])
  if(X[i,3]>0) X[i,3] <- log(X[i,3]) else X[i,3] <- log(-X[i,3])
  if(X[i,4]>0) X[i,4] <- log(X[i,4]) else X[i,4] <- log(-X[i,4])
  if(X[i,6]>1) X[i,6] <- 0   #1 if white and 0 otherwise
}


###############User Defined Functions for Analysis Below############

#unnormalized density function used to generate imputations
XgenNA <- function(x,Ymat,Xmat,BetaA,GammaA,PiA,AlphaAA,AlphaBA,SigmaA,Sigma2A,ShapeA){
  log((dgev(log(Ymat), loc=c(Xmat[1],x,Xmat[3:5])%*%BetaA, scale=Sigma2A, shape=ShapeA)/pgev(log(cutoff), loc=c(Xmat[1],x,Xmat[3:5])%*%BetaA, scale=Sigma2A, shape=ShapeA)*(1+exp(c(Xmat[1],x,Xmat[3:5])%*%GammaA))^(-1)*(Ymat<cutoff)+(1+exp(-c(Xmat[1],x,Xmat[3:5])%*%GammaA))^(-1)*(Ymat==cutoff))*(PiA[1]*dnorm(x,c(Xmat[1],Xmat[3:5])%*%AlphaAA,SigmaA[1])+PiA[2]*dnorm(x,c(Xmat[1],Xmat[3:5])%*%AlphaBA,SigmaA[2])))
  
}

#full-data likelihood that is maximized with appropriate weights (need to approximate conditional expected log-likelihood
#in E-step of EM algorithm;  note that logistic regression and AFT parameters can be estimated separately if
#numerical issues arise
MaxNA <- function(Y,X){
  weight <- X[,6]
  X <- X[,1:5]
  like <- function(par){
    llike <- rep(0,length(Y))
    llike <- -weight*log((dgev(log(Y), loc=X%*%par[1:5], scale=abs(par[11]), shape=par[12])/pgev(log(cutoff), loc=X%*%par[1:5], scale=abs(par[11]), shape=par[12])*(1+exp(X%*%par[6:10]))^(-1)*(Y<cutoff)+(1+exp(-X%*%par[6:10]))^(-1)*(Y==cutoff)))
    return(sum(llike))
  }
  like
}	

#Rather than used Louis method for estimating the variance of the MLEs, this observed data likelihood can be used
#to obtain variance estimates if it is numerically stable; note, some issues did exist if large number of imputations
#weren't used in the MCEM algorithm, usually originating in the scale parameter of the generalized gamma distribution
Var <- function(par,Y,X,dvals){
  llike <- rep(0,length(Y))
  for(i in 1:length(Y)){if(X[i,2]!=dvals[i]) llike[i] <- -log(((dgev(log(Y[i]), loc=X[i,]%*%par[1:5], scale=abs(par[11]), shape=par[12])/pgev(log(cutoff), loc=X[i,]%*%par[1:5], scale=abs(par[11]), shape=par[12])*(1+exp(X[i,]%*%par[6:10]))^(-1)*(Y[i]<cutoff)+(1+exp(-X[i,]%*%par[6:10]))^(-1)*(Y[i]==cutoff)))*(par[13]*dnorm(X[i,2],X[i,-2]%*%par[14:17],par[22])+(1-par[13])*dnorm(X[i,2],X[i,-2]%*%par[18:21],par[23])))
  if(X[i,2]==dvals[i]) llike[i] <- -log(integrate(integrand,lower=-Inf,upper=dvals[i], par=par,Y=Y[i],X=X[i,])$value)}
  return(sum(llike))
}

#function that requires integration as part of the observed data likelihood
integrand <- function(z,par,Y,X){
  (dgev(log(Y),(par[1] + z*par[2] + X[3:5]%*%par[3:5]),scale=abs(par[11]),shape=par[12])/pgev(log(cutoff),(par[1] + z*par[2] + X[3:5]%*%par[3:5]),scale=abs(par[11]),shape=par[12])*(1+exp((par[6] + z*par[7] + X[3:5]%*%par[8:10])))^(-1)*(Y<cutoff)+(1+exp(-(par[6] + z*par[7] + X[3:5]%*%par[8:10])))^(-1)*(Y==cutoff))*(par[13]*dnorm(z,X[-2]%*%par[14:17],par[22])+(1-par[13])*dnorm(z,X[-2]%*%par[18:21],par[23]))
}


###############Semicontinuous Survival Analysis of GenIMS Data#################

#Note:
#Three separate codes (for each analysis are provided)...could use a function to make code much shorter, though
#since codes changes were minimum, was faster time-wise to simply copy the code and make appropriate changes



#####Analysis for First Cytokine Biomarker (TNF)#####
X1 <- X[,c(1,2,5,6,8)]	#Defining observed covariate matrix

#Finding starting values for MC-EM algorithm
X1star <- X1		#Creating second copy of covariate matrix

#Lazy way of replacing censored values for the purpose of getting starting values (DL/2 + standard normal error)
X1star[which(X1star[,2]==dmatrix[,1]),2] <- dmatrix[which(X1star[,2]==dmatrix[,1]),1]/2+rnorm(length(dmatrix[which(X1star[,2]==dmatrix[,1]),1]))

#Obtaining mixture of normals regression parameter estimates for covariate distribution
F <- flexmix(X1star[,2] ~ X1star[,3] + X1star[,4] + X1star[,5], k = 2)
if(which(parameters(F)[1,]==max(parameters(F)[1,]))==1){ PiI <- prior(F) 
SigmaI <- as.matrix(parameters(F))[5,]}
if(which(parameters(F)[1,]==max(parameters(F)[1,]))!=1) { PiI <- rev(prior(F))
SigmaI <- rev(as.matrix(parameters(F))[5,])}
AlphaAI <- as.matrix(parameters(F))[1:4,which(parameters(F)[1,]==max(parameters(F)[1,]))]
AlphaBI <- as.matrix(parameters(F))[1:4,which(parameters(F)[1,]==min(parameters(F)[1,]))]

#Obtaining generalized gamma AFT parameter estimates for a model that treats all "90" responses as censored at 90
G <- flexsurvreg(Surv(Y,(Y<90)) ~X1star[,2] + X1star[,3] + X1star[,4] + X1star[,5], dist="gengamma")
BetaI <- c(G$coefficients[1],G$coefficients[4:7])
Sigma2I <- 10 #G$coefficients[2]
ShapeI <- 0.1  #G$coefficients[3]

#Obtaining logistic regression parameter estimates for point mass probility (brglm reduces bias in estimates)
H <- brglm((Y==90) ~ X1star[,2] + X1star[,3] + X1star[,4] + X1star[,5] , dist="binomial")
GammaI <- H$coefficients

#Defining matrix to store estimates at each EM iteration, with the first row being 0's (purely for programming 
#purposes regarding defining convergence in the while loop below) and the second row the starting values 
ParEstsA <- rbind(rep(0,12),c(BetaI,GammaI,Sigma2I, ShapeI))	#parameters of main interest
ParEstsB <- rbind(rep(0,12),c(PiI,AlphaAI,AlphaBI,SigmaI))		#nuisance parameters


###EM algorithm update loop###

j=2 #index defining current row of matrix

#Continues updating estimates until consecutive sets of parameter estimates are close; a stricter requirement
#is placed on the parameters of main interest in ParEsts A
while(sum(abs(ParEstsA[j,]-ParEstsA[j-1,]))>0.005 | sum(abs(ParEstsB[j,]-ParEstsB[j-1,]))>0.1 ){
  j <- j+1
  n.MC <- 30 		   #number of imputations for approximately E-step expecations; starting low
  if(j>7) n.MC <- 200  #increasing number of imputations
  if(j>19) n.MC <- 1000	#further increasing imputations (in practice, after convergence, I ran
  #another iteration with 10,000 imputations
  
  #alternatively, if estimates are close, rampe up the imputations earlier
  if(sum(abs(ParEstsA[j-1,]-ParEstsA[j-2,]))<0.01 & sum(abs(ParEstsB[j-1,]-ParEstsB[j-2,]))<0.2 & j<20) n.MC <- 500
  
  Xnew <- NULL	#X matrix with imputed biomarker values (when censoring occured); last column is weight 1 or 1/n.MC
  Ynew <- NULL	#Y matrix with n.MC repeated Y[i] values for each individual with censored biomarker
  
  #For each individual, if biomark is censored, replaces observation with n.MC imputations
  for(k in 1:n){
    if(X1[k,2] != dmatrix[k,1]){ Xnew <- rbind(Xnew, c(X1[k,], 1))
    Ynew <- c(Ynew, Y[k])}
    
    Xvals <-NULL	#vector to temporarily store imputed values for a single individual
    if(X1[k,2]==dmatrix[k,1]){ 
      
      #loop that iterates until enough imputations are generated (this is used since the imputations
      #must be less than the detection limit, i.e. an easy way to generate data from a truncated dist.)
      while(length(Xvals)<n.MC){
        (Xvals <- c(Xvals, metrop(XgenNA, initial=0, nspac=5, nbatch = 2*n.MC, Ymat=Y[k], Xmat=X1[k,], BetaA=BetaI, GammaA=GammaI, PiA=PiI, AlphaAA=AlphaAI, AlphaBA=AlphaBI, SigmaA=SigmaI,Sigma2A=Sigma2I, ShapeA=ShapeI)$batch))
        Xvals <- Xvals[which(Xvals<dmatrix[k,1])]	
      }
      Xnew <- rbind(Xnew, cbind(rep(1,n.MC),Xvals[1:n.MC],matrix(X1[k,3:5],n.MC,3,byrow=TRUE), rep(1/n.MC,n.MC)))
      Ynew <- c(Ynew, rep(Y[k],n.MC))
    }
  }
  
  
  #Obtaining estimates in mixture of normals regression based on imputed data with appropriate weighting
  F <- flexmix(Xnew[,2] ~ Xnew[,3] + Xnew[,4] + Xnew[,5], weights=as.integer(n.MC*Xnew[,6]), k = 2)
  if(which(parameters(F)[1,]==max(parameters(F)[1,]))==1){ PiI <- prior(F) 
  SigmaI <- as.matrix(parameters(F))[5,]}
  if(which(parameters(F)[1,]==max(parameters(F)[1,]))!=1) { PiI <- rev(prior(F))
  SigmaI <- rev(as.matrix(parameters(F))[5,])}
  AlphaAI <- as.matrix(parameters(F))[1:4,which(parameters(F)[1,]==max(parameters(F)[1,]))]
  AlphaBI <- as.matrix(parameters(F))[1:4,which(parameters(F)[1,]==min(parameters(F)[1,]))]
  
  #Obtaining estimates for logistic and AFT regression (could be found separately if numerical issues arose)
  vals <- tryCatch(optim(c(BetaI,GammaI,Sigma2I,ShapeI),MaxNA(Ynew,Xnew),method="BFGS")$par, error=function(...) nlm(MaxNA(Ynew,Xnew), c(BetaI,GammaI,Sigma2I,ShapeI))$estimate)
  BetaI <- vals[1:5]
  GammaI <- vals[6:10]
  Sigma2I <- vals[11]
  ShapeI <- vals[12]
  ParEstsA <- rbind(ParEstsA,c(BetaI,GammaI,Sigma2I,ShapeI))
  ParEstsB <- rbind(ParEstsB,c(PiI,AlphaAI,AlphaBI,SigmaI))
}

#Final MLEs and variance matrix
ParEsts1 <- c(ParEstsA[j,],ParEstsB[j,])
hessb <- hessian(Var,ParEsts1[-14],Y=Y,X=X1,dvals=dmatrix[,1])
VarEsts1 <- solve(hessb)






#####Analysis for Second Cytokine Biomarker (IL-6)#####
X1 <- X[,c(1,3,5,6,8)]	#Defining observed covariate matrix

#Finding starting values for MC-EM algorithm
X1star <- X1		#Creating second copy of covariate matrix

#Lazy way of replacing censored values for the purpose of getting starting values (DL/2 + standard normal error)
X1star[which(X1star[,2]==dmatrix[,2]),2] <- dmatrix[which(X1star[,2]==dmatrix[,2]),1]/2+rnorm(length(dmatrix[which(X1star[,2]==dmatrix[,2]),1]))

#Obtaining mixture of normals regression parameter estimates for covariate distribution
F <- flexmix(X1star[,2] ~ X1star[,3] + X1star[,4] + X1star[,5], k = 2)
if(which(parameters(F)[1,]==max(parameters(F)[1,]))==1){ PiI <- prior(F) 
SigmaI <- as.matrix(parameters(F))[5,]}
if(which(parameters(F)[1,]==max(parameters(F)[1,]))!=1) { PiI <- rev(prior(F))
SigmaI <- rev(as.matrix(parameters(F))[5,])}
AlphaAI <- as.matrix(parameters(F))[1:4,which(parameters(F)[1,]==max(parameters(F)[1,]))]
AlphaBI <- as.matrix(parameters(F))[1:4,which(parameters(F)[1,]==min(parameters(F)[1,]))]

#Obtaining generalized gamma AFT parameter estimates for a model that treats all "90" responses as censored at 90
G <- flexsurvreg(Surv(Y,(Y<90)) ~X1star[,2] + X1star[,3] + X1star[,4] + X1star[,5], dist="gengamma")
BetaI <- c(G$coefficients[1],G$coefficients[4:7])
Sigma2I <- 10 #G$coefficients[2]
ShapeI <- 0.1  #G$coefficients[3]

#Obtaining logistic regression parameter estimates for point mass probility (brglm reduces bias in estimates)
H <- brglm((Y==90) ~ X1star[,2] + X1star[,3] + X1star[,4] + X1star[,5] , dist="binomial")
GammaI <- H$coefficients

#Defining matrix to store estimates at each EM iteration, with the first row being 0's (purely for programming 
#purposes regarding defining convergence in the while loop below) and the second row the starting values 
ParEstsA <- rbind(rep(0,12),c(BetaI,GammaI,Sigma2I, ShapeI))	#parameters of main interest
ParEstsB <- rbind(rep(0,12),c(PiI,AlphaAI,AlphaBI,SigmaI))		#nuisance parameters


###EM algorithm update loop###

j=2 #index defining current row of matrix

#Continues updating estimates until consecutive sets of parameter estimates are close; a stricter requirement
#is placed on the parameters of main interest in ParEsts A
while(sum(abs(ParEstsA[j,]-ParEstsA[j-1,]))>0.005 | sum(abs(ParEstsB[j,]-ParEstsB[j-1,]))>0.1 ){
  j <- j+1
  n.MC <- 30 		   #number of imputations for approximately E-step expecations; starting low
  if(j>7) n.MC <- 200  #increasing number of imputations
  if(j>19) n.MC <- 1000	#further increasing imputations (in practice, after convergence, I ran
  #another iteration with 10,000 imputations
  
  #alternatively, if estimates are close, rampe up the imputations earlier
  if(sum(abs(ParEstsA[j-1,]-ParEstsA[j-2,]))<0.01 & sum(abs(ParEstsB[j-1,]-ParEstsB[j-2,]))<0.2 & j<20) n.MC <- 500
  
  Xnew <- NULL	#X matrix with imputed biomarker values (when censoring occured); last column is weight 1 or 1/n.MC
  Ynew <- NULL	#Y matrix with n.MC repeated Y[i] values for each individual with censored biomarker
  
  #For each individual, if biomark is censored, replaces observation with n.MC imputations
  for(k in 1:n){
    if(X1[k,2] != dmatrix[k,2]){ Xnew <- rbind(Xnew, c(X1[k,], 1))
    Ynew <- c(Ynew, Y[k])}
    
    Xvals <-NULL	#vector to temporarily store imputed values for a single individual
    if(X1[k,2]==dmatrix[k,2]){ 
      
      #loop that iterates until enough imputations are generated (this is used since the imputations
      #must be less than the detection limit, i.e. an easy way to generate data from a truncated dist.)
      while(length(Xvals)<n.MC){
        (Xvals <- c(Xvals, metrop(XgenNA, initial=0, nspac=5, nbatch = 2*n.MC, Ymat=Y[k], Xmat=X1[k,], BetaA=BetaI, GammaA=GammaI, PiA=PiI, AlphaAA=AlphaAI, AlphaBA=AlphaBI, SigmaA=SigmaI,Sigma2A=Sigma2I, ShapeA=ShapeI)$batch))
        Xvals <- Xvals[which(Xvals<dmatrix[k,2])]	
      }
      Xnew <- rbind(Xnew, cbind(rep(1,n.MC),Xvals[1:n.MC],matrix(X1[k,3:5],n.MC,3,byrow=TRUE), rep(1/n.MC,n.MC)))
      Ynew <- c(Ynew, rep(Y[k],n.MC))
    }
  }
  
  
  #Obtaining estimates in mixture of normals regression based on imputed data with appropriate weighting
  F <- flexmix(Xnew[,2] ~ Xnew[,3] + Xnew[,4] + Xnew[,5], weights=as.integer(n.MC*Xnew[,6]), k = 2)
  if(which(parameters(F)[1,]==max(parameters(F)[1,]))==1){ PiI <- prior(F) 
  SigmaI <- as.matrix(parameters(F))[5,]}
  if(which(parameters(F)[1,]==max(parameters(F)[1,]))!=1) { PiI <- rev(prior(F))
  SigmaI <- rev(as.matrix(parameters(F))[5,])}
  AlphaAI <- as.matrix(parameters(F))[1:4,which(parameters(F)[1,]==max(parameters(F)[1,]))]
  AlphaBI <- as.matrix(parameters(F))[1:4,which(parameters(F)[1,]==min(parameters(F)[1,]))]
  
  #Obtaining estimates for logistic and AFT regression (could be found separately if numerical issues arose)
  vals <- tryCatch(optim(c(BetaI,GammaI,Sigma2I,ShapeI),MaxNA(Ynew,Xnew),method="BFGS")$par, error=function(...) nlm(MaxNA(Ynew,Xnew), c(BetaI,GammaI,Sigma2I,ShapeI))$estimate)
  BetaI <- vals[1:5]
  GammaI <- vals[6:10]
  Sigma2I <- vals[11]
  ShapeI <- vals[12]
  ParEstsA <- rbind(ParEstsA,c(BetaI,GammaI,Sigma2I,ShapeI))
  ParEstsB <- rbind(ParEstsB,c(PiI,AlphaAI,AlphaBI,SigmaI))
}

#Final MLEs and variance matrix
ParEsts2 <- c(ParEstsA[j,],ParEstsB[j,])
hessb <- hessian(Var,ParEsts2[-14],Y=Y,X=X1,dvals=dmatrix[,2])
VarEsts2 <- solve(hessb)








#####Analysis for Third Cytokine Biomarker (IL-10)#####
X1 <- X[,c(1,4,5,6,8)]	#Defining observed covariate matrix

#Finding starting values for MC-EM algorithm
X1star <- X1		#Creating second copy of covariate matrix

#Lazy way of replacing censored values for the purpose of getting starting values (DL/2 + standard normal error)
X1star[which(X1star[,2]==dmatrix[,3]),2] <- dmatrix[which(X1star[,2]==dmatrix[,3]),1]/2+rnorm(length(dmatrix[which(X1star[,2]==dmatrix[,3]),1]))

#Obtaining mixture of normals regression parameter estimates for covariate distribution
F <- flexmix(X1star[,2] ~ X1star[,3] + X1star[,4] + X1star[,5], k = 2)
if(which(parameters(F)[1,]==max(parameters(F)[1,]))==1){ PiI <- prior(F) 
SigmaI <- as.matrix(parameters(F))[5,]}
if(which(parameters(F)[1,]==max(parameters(F)[1,]))!=1) { PiI <- rev(prior(F))
SigmaI <- rev(as.matrix(parameters(F))[5,])}
AlphaAI <- as.matrix(parameters(F))[1:4,which(parameters(F)[1,]==max(parameters(F)[1,]))]
AlphaBI <- as.matrix(parameters(F))[1:4,which(parameters(F)[1,]==min(parameters(F)[1,]))]

#Obtaining generalized gamma AFT parameter estimates for a model that treats all "90" responses as censored at 90
G <- flexsurvreg(Surv(Y,(Y<90)) ~X1star[,2] + X1star[,3] + X1star[,4] + X1star[,5], dist="gengamma")
BetaI <- c(G$coefficients[1],G$coefficients[4:7])
Sigma2I <- 10 #G$coefficients[2]
ShapeI <- 0.1  #G$coefficients[3]

#Obtaining logistic regression parameter estimates for point mass probility (brglm reduces bias in estimates)
H <- brglm((Y==90) ~ X1star[,2] + X1star[,3] + X1star[,4] + X1star[,5] , dist="binomial")
GammaI <- H$coefficients

#Defining matrix to store estimates at each EM iteration, with the first row being 0's (purely for programming 
#purposes regarding defining convergence in the while loop below) and the second row the starting values 
ParEstsA <- rbind(rep(0,12),c(BetaI,GammaI,Sigma2I, ShapeI))	#parameters of main interest
ParEstsB <- rbind(rep(0,12),c(PiI,AlphaAI,AlphaBI,SigmaI))		#nuisance parameters


###EM algorithm update loop###

j=2 #index defining current row of matrix

#Continues updating estimates until consecutive sets of parameter estimates are close; a stricter requirement
#is placed on the parameters of main interest in ParEsts A
while(sum(abs(ParEstsA[j,]-ParEstsA[j-1,]))>0.005 | sum(abs(ParEstsB[j,]-ParEstsB[j-1,]))>0.1 ){
  j <- j+1
  n.MC <- 30 		   #number of imputations for approximately E-step expecations; starting low
  if(j>7) n.MC <- 200  #increasing number of imputations
  if(j>19) n.MC <- 1000	#further increasing imputations (in practice, after convergence, I ran
  #another iteration with 10,000 imputations
  
  #alternatively, if estimates are close, rampe up the imputations earlier
  if(sum(abs(ParEstsA[j-1,]-ParEstsA[j-2,]))<0.01 & sum(abs(ParEstsB[j-1,]-ParEstsB[j-2,]))<0.2 & j<20) n.MC <- 500
  
  Xnew <- NULL	#X matrix with imputed biomarker values (when censoring occured); last column is weight 1 or 1/n.MC
  Ynew <- NULL	#Y matrix with n.MC repeated Y[i] values for each individual with censored biomarker
  
  #For each individual, if biomark is censored, replaces observation with n.MC imputations
  for(k in 1:n){
    if(X1[k,2] != dmatrix[k,3]){ Xnew <- rbind(Xnew, c(X1[k,], 1))
    Ynew <- c(Ynew, Y[k])}
    
    Xvals <-NULL	#vector to temporarily store imputed values for a single individual
    if(X1[k,2]==dmatrix[k,3]){ 
      
      #loop that iterates until enough imputations are generated (this is used since the imputations
      #must be less than the detection limit, i.e. an easy way to generate data from a truncated dist.)
      while(length(Xvals)<n.MC){
        (Xvals <- c(Xvals, metrop(XgenNA, initial=0, nspac=5, nbatch = 2*n.MC, Ymat=Y[k], Xmat=X1[k,], BetaA=BetaI, GammaA=GammaI, PiA=PiI, AlphaAA=AlphaAI, AlphaBA=AlphaBI, SigmaA=SigmaI,Sigma2A=Sigma2I, ShapeA=ShapeI)$batch))
        Xvals <- Xvals[which(Xvals<dmatrix[k,3])]	
      }
      Xnew <- rbind(Xnew, cbind(rep(1,n.MC),Xvals[1:n.MC],matrix(X1[k,3:5],n.MC,3,byrow=TRUE), rep(1/n.MC,n.MC)))
      Ynew <- c(Ynew, rep(Y[k],n.MC))
    }
  }
  
  
  #Obtaining estimates in mixture of normals regression based on imputed data with appropriate weighting
  F <- flexmix(Xnew[,2] ~ Xnew[,3] + Xnew[,4] + Xnew[,5], weights=as.integer(n.MC*Xnew[,6]), k = 2)
  if(which(parameters(F)[1,]==max(parameters(F)[1,]))==1){ PiI <- prior(F) 
  SigmaI <- as.matrix(parameters(F))[5,]}
  if(which(parameters(F)[1,]==max(parameters(F)[1,]))!=1) { PiI <- rev(prior(F))
  SigmaI <- rev(as.matrix(parameters(F))[5,])}
  AlphaAI <- as.matrix(parameters(F))[1:4,which(parameters(F)[1,]==max(parameters(F)[1,]))]
  AlphaBI <- as.matrix(parameters(F))[1:4,which(parameters(F)[1,]==min(parameters(F)[1,]))]
  
  #Obtaining estimates for logistic and AFT regression (could be found separately if numerical issues arose)
  vals <- tryCatch(optim(c(BetaI,GammaI,Sigma2I,ShapeI),MaxNA(Ynew,Xnew),method="BFGS")$par, error=function(...) nlm(MaxNA(Ynew,Xnew), c(BetaI,GammaI,Sigma2I,ShapeI))$estimate)
  BetaI <- vals[1:5]
  GammaI <- vals[6:10]
  Sigma2I <- vals[11]
  ShapeI <- vals[12]
  ParEstsA <- rbind(ParEstsA,c(BetaI,GammaI,Sigma2I,ShapeI))
  ParEstsB <- rbind(ParEstsB,c(PiI,AlphaAI,AlphaBI,SigmaI))
}

#Final MLEs and variance matrix
ParEsts3 <- c(ParEstsA[j,],ParEstsB[j,])
hessb <- hessian(Var,ParEsts1[-14],Y=Y,X=X1,dvals=dmatrix[,3])
VarEsts3 <- solve(hessb)




#######Getting 2 df test (for IL-10 in example here)######
params3 <- ParEsts3[-14]
variance3 <- VarEsts3
r <- matrix(c(0,1,rep(0,21),0,0,0,0,0,0,1,rep(0,16)),2,23,byrow=TRUE)
t(r%*%params3)%*%solve(r%*%(variance3)%*%t(r))%*%(r%*%params3)
