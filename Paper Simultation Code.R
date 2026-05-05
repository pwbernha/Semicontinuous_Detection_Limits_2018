#####################################################################################
#NOTE: PLEASE DO NOT USE CODE FOR ANY PUBLICATION PURPOSES WITHOUT FIRST CONTACTING #
#      PAUL BERNHARDT AT PAUL.BERNHARDT@VILLANOVA.EDU                               #
#####################################################################################


####################Necessary Packages for Simulation######################
library(MCMCpack)
library(mvtnorm)
library(msm)
library(mcmc)
library(flexmix)


##################Simulation Parameters######################
set.seed(8089205)					#seed number
n <- 200						#number of individuals per dataset
N <- 500						#number of simulated datasets
cutoff <- 90					#cure time
d <- 1.24 #20% censoring			#detection limit
#d <- 2.23 #40% censoring
Beta <- c(8,-0.5,-0.03,-0.01,-1)		#AFT regression parameters
Gamma <- c(1.35,-0.05,-0.02,-0.005,-0.5)	#logistic regression parameters
AlphaA <- c(3,0.05,-0.03,1)			#Mixture of normals mean1
AlphaB <- c(2,-0.1,0.08,1)			#Mixture of normals mean2
Pi <- c(0.3,0.7)					#Mixture of normals probabilities
Sigmamix <- c(1,1)



###################Necessary Program Functions####################

#unnormalized log likelihood function used to generate imputations for censored values using mixture of normals f(x|z)
Xgen <- function(x,Ymat,Xmat,BetaA,GammaA,PiA,AlphaAA,AlphaBA,SigmaA,Sigma2A){
  log((dnorm(log(Ymat),c(Xmat[1],x,Xmat[3:5])%*%BetaA,Sigma2A)/pnorm(log(cutoff),c(Xmat[1],x,Xmat[3:5])%*%BetaA,Sigma2A)*(1+exp(c(Xmat[1],x,Xmat[3:5])%*%GammaA))^(-1)*(Ymat<cutoff)+(1+exp(-c(Xmat[1],x,Xmat[3:5])%*%GammaA))^(-1)*(Ymat==cutoff))*(PiA[1]*dnorm(x,c(Xmat[1],Xmat[3:5])%*%AlphaAA,SigmaA[1])+PiA[2]*dnorm(x,c(Xmat[1],Xmat[3:5])%*%AlphaBA,SigmaA[2])))
}

#unnormalized log likelihood function used to generate imputations for censored values using normal f(x|z)
Xgen2 <- function(x,Ymat,Xmat,BetaA,GammaA,BetaNormA,SigmaNormA, Sigma2A){
  log((dnorm(log(Ymat),c(Xmat[1],x,Xmat[3:5])%*%BetaA,Sigma2A)/pnorm(log(cutoff),c(Xmat[1],x,Xmat[3:5])%*%BetaA,Sigma2A)*(1+exp(c(Xmat[1],x,Xmat[3:5])%*%GammaA))^(-1)*(Ymat<cutoff)+(1+exp(-c(Xmat[1],x,Xmat[3:5])%*%GammaA))^(-1)*(Ymat==cutoff))*(dnorm(x,c(Xmat[1],Xmat[3:5])%*%BetaNormA,SigmaNormA)))
}

#log-likelihood function (with weighting to approximate expectation in E-step of EM) for proposed model,
# to be maximized after imputations
Max <- function(Y,X){
  weight <- X[,6]
  X <- X[,1:5]
  like <- function(par){
    llike <- rep(0,length(Y))
    llike <- -weight*log((dnorm(log(Y),X%*%par[1:5],par[11])/pnorm(log(cutoff),X%*%par[1:5],par[11])*(1+exp(X%*%par[6:10]))^(-1)*(Y<cutoff)+(1+exp(-X%*%par[6:10]))^(-1)*(Y==cutoff)))
    return(sum(llike))
  }
  like
}	

#log-likelihood function (with weighting to approximate expectation in E-step of EM) for proposed model except 
#using normal f(x|z), to be maximized after imputations
Max2 <- function(Y,X){
  weight <- X[,6]
  X <- X[,1:5]
  like <- function(par){
    llike <- rep(0,length(Y))
    llike <- -weight*log(dnorm(X[,2],X[,c(1,3:5)]%*%par[1:4],par[5]))
    return(sum(llike))
  }
  like
}	

#Funcion to generate complete set of simulation data (for N data sets) based on simulation parameters
GenData <- function(N,n,Pi,AlphaA,AlphaB,Sigmamix,Gamma,Beta,cutoff,d){
  Xgenerate <- NULL
  Ygenerate <- NULL
  for(i in 1:N){
    ###Generate X matrix###
    X <- matrix(c(1,rep(0,4)),n,5,byrow=TRUE)
    age <- rbeta(n,3,2)*84+18					#generate age variable
    gender <- rbinom(n,1,0.49)					#generate sex variable
    X[,3] <- meanX[2] + rnorm(n,30+5*gender+0.4*age,30)	#generate continuous covariate
    X[,4:5] <- cbind(age,gender)
    prop <- rbinom(n,1,Pi[1])
    
    #generate cytokine subject to DL
    X[,2] <- prop*rnorm(n, X[,c(1,3:5)]%*%AlphaA, Sigmamix[1])+(1-prop)*rnorm(n, X[,c(1,3:5)]%*%AlphaB, Sigmamix[2])
    
    #Add censoring below LOD to X1#
    X[which(X[,2]<d),2] <- d
    
    ###Generate Y vector###
    p <- rbinom(n, 1, (1+exp(-X%*%Gamma))^(-1))
    logY <- rtnorm(n,X%*%Beta,2,upper=log(cutoff))
    Y <- exp(logY)
    Y[which(p==1)]<- cutoff
    
    ###Adding X,Y to simulated data matrix###
    Ygenerate <- c(Ygenerate,Y)
    Xgenerate <- rbind(Xgenerate,X)
  }
  
  return(list(Ygenerate,Xgenerate))
}




##############Simulation##################

#####Generate Data for Simulations ####
Data <- GenData(N,n,Pi,AlphaA,AlphaB,Sigmamix,Gamma,Beta,cutoff,d)
Yvector <- Data[[1]]
Xmatrix <- Data[[2]]

####Matrices to Store Final Parameter Estimates####
Ests <- matrix(0,N,23)	#Proposed Method
Ests2 <- matrix(0,N,16)	#Proposed Method using Normal instead of Mixture of Noramls to Model f(x|z)
Ests3 <- matrix(0,N,11)	#Full data analysis with DL/sqrt(2) replacing censored x's
Ests4 <- matrix(0,N,11)	#Complete data analysis

#Looping through each data set to get parameter estimates for each method
for(i in 1:N){
  
  #Defining response vector and covariance matrix for the i^th dataset
  Y <- Yvector[((i-1)*n+1):(i*n)]
  X <- Xmatrix[((i-1)*n+1):(i*n),]
  
  ############Sim #1: Proposed Method############
  BetaI <- Beta
  GammaI <- Gamma
  PiI <- Pi
  AlphaAI <- AlphaA
  AlphaBI <- AlphaB
  SigmaI <- Sigmamix
  Sigma2I <- 2
  
  #matrices to store EM estimate updates (for programming convenience using while statement below, first row is 0s 
  #and second row is simply the true values, which may or may not be better than finding starting values)
  ParEstsA <- rbind(rep(0,11),c(BetaI,GammaI,Sigma2I))
  ParEstsB <- rbind(rep(0,12),c(PiI,AlphaAI,AlphaBI,SigmaI))
  
  ####EM Algorithm Estimate Updates#####
  j=2 #index for current EM iteration
  
  #Continue Iterating until consecutive parameter estimates are close
  while(sum(abs(ParEstsA[j,]-ParEstsA[j-1,]))>0.01 | sum(abs(ParEstsB[j,]-ParEstsB[j-1,]))>0.1 ){
    j <- j+1 #increasing iteration index by 1
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
      if(X[k,2] != d[1]){ Xnew <- rbind(Xnew, c(X[k,], 1))
      Ynew <- c(Ynew, Y[k])}
      
      Xvals <-NULL	#vector to temporarily store imputed values for a single individual
      if(X[k,2]==d[1]){ 
        #loop that iterates until enough imputations are generated (this is used since the imputations
        #must be less than the detection limit, i.e. an easy way to generate data from a truncated dist.)
        while(length(Xvals)<n.MC){
          (Xvals <- c(Xvals, metrop(Xgen, initial=0, nspac=5, nbatch = 2*n.MC, Ymat=Y[k], Xmat=X[k,], BetaA=BetaI, GammaA=GammaI, PiA=PiI, AlphaAA=AlphaAI, AlphaBA=AlphaBI, SigmaA=SigmaI,Sigma2A=Sigma2I)$batch))
          Xvals <- c(Xvals,Xvals[-c(1:500)]) #allowing burn-in
          Xvals <- Xvals[which(Xvals<d[1])]	
        }
        
        Xnew <- rbind(Xnew, cbind(rep(1,n.MC),Xvals[1:n.MC],matrix(X[k,3:5],n.MC,3,byrow=TRUE), rep(1/n.MC,n.MC)))
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
    
    vals <- tryCatch(nlm(Max(Ynew,Xnew), c(BetaI,GammaI,Sigma2I))$estimate, error=function(...) optim(c(BetaI,GammaI,Sigma2I),Max(Ynew,Xnew),method="BFGS")$par)
    BetaI <- vals[1:5]
    GammaI <- vals[6:10]
    Sigma2I <- vals[11]
    ParEstsA <- rbind(ParEstsA,c(BetaI,GammaI,Sigma2I))
    ParEstsB <- rbind(ParEstsB,c(PiI,AlphaAI,AlphaBI,SigmaI))
  }
  
  #Final Estimate for ith generated data set, proposed method
  Ests[i,] <- c(ParEstsA[j,],ParEstsB[j,])
  
  
  
  
  ##############Sim #2: Proposed Method using normal for f(x|z)#############
  BetaI <- Beta
  GammaI <- Gamma
  Sigma2I <- 2
  BetaNormI <- (AlphaA+AlphaB)/2
  SigmaNormI <- 3
  vals <- tryCatch(nlm(Max2(Ynew,Xnew), c(BetaNormI, SigmaNormI))$estimate, error=function(...) optim( c(BetaNormI, SigmaNormI),Max2(Ynew,Xnew),method="BFGS")$par)
  BetaNormI <- vals[1:4]
  SigmaNormI <- vals[5]
  
  #matrices to store EM estimate updates (for programming convenience using while statement below, first row is 0s 
  #and second row is simply the true values, which may or may not be better than finding starting values)
  ParEstsA <- rbind(rep(0,11),c(BetaI,GammaI,Sigma2I))
  ParEstsB <- rbind(rep(0,5),c(BetaNormI, SigmaNormI))
  
  ####EM Algorithm Estimate Updates#####
  j=2 #index for current EM iteration
  
  #Continue Iterating until consecutive parameter estimates are close
  while(sum(abs(ParEstsA[j,]-ParEstsA[j-1,]))>0.02 | sum(abs(ParEstsB[j,]-ParEstsB[j-1,]))>0.1 ){
    
    j <- j+1 #increasing iteration index by 1
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
      if(X[k,2] != d[1]){ Xnew <- rbind(Xnew, c(X[k,], 1))
      Ynew <- c(Ynew, Y[k])}
      
      Xvals <-NULL	#vector to temporarily store imputed values for a single individual
      
      if(X[k,2]==d[1]){ 
        
        #loop that iterates until enough imputations are generated (this is used since the imputations
        #must be less than the detection limit, i.e. an easy way to generate data from a truncated dist.)
        while(length(Xvals)<n.MC){
          (Xvals <- c(Xvals, metrop(Xgen2, initial=0, nspac=5, nbatch = 2*n.MC, Ymat=Y[k], Xmat=X[k,], BetaA=BetaI, GammaA=GammaI, BetaNormA=BetaNormI, SigmaNormA=SigmaNormI, Sigma2A=Sigma2I)$batch))
          Xvals <- c(Xvals,Xvals[-c(1:500)])   #allowing for burn-in
          Xvals <- Xvals[which(Xvals<d[1])]	
        }
        Xnew <- rbind(Xnew, cbind(rep(1,n.MC),Xvals[1:n.MC],matrix(X[k,3:5],n.MC,3,byrow=TRUE), rep(1/n.MC,n.MC)))
        Ynew <- c(Ynew, rep(Y[k],n.MC))
      }
    }
    
    #Obtaining estimates for logistic and AFT regression (could be found separately if numerical issues arose)
    vals <- tryCatch(nlm(Max(Ynew,Xnew), c(BetaI,GammaI,Sigma2I))$estimate, error=function(...) optim(c(BetaI,GammaI,Sigma2I),Max(Ynew,Xnew),method="BFGS")$par)
    BetaI <- vals[1:5]
    GammaI <- vals[6:10]
    Sigma2I <- vals[11]
    
    vals <- tryCatch(nlm(Max2(Ynew,Xnew), c(BetaNormI, SigmaNormI))$estimate, error=function(...) optim( c(BetaNormI, SigmaNormI),Max2(Ynew,Xnew),method="BFGS")$par)
    BetaNormI <- vals[1:4]
    SigmaNormI <- vals[5]
    ParEstsA <- rbind(ParEstsA,c(BetaI,GammaI,Sigma2I))
    ParEstsB <- rbind(ParEstsB,c(BetaNormI,SigmaNormI))
    
  }
  
  #Final Estimate for ith generated data set, proposed method with normal covariate distribution
  Ests2[i,] <- c(ParEstsA[j,], ParEstsB[j,])
  
  
  
  
  ######Sim #3: Single Imputation by Detection Limit/sqrt(2) #######
  
  Ynew <- Y
  Xnew <- X
  Xnew[which(Xnew[,2]==d),2] <- d/sqrt(2)   #replacing censored values by DL/sqrt(2)
  Xnew <- cbind(Xnew,1)
  
  #Obtaining estimates for logistic and AFT regression (could be found separately)
  vals <- tryCatch(nlm(Max(Ynew,Xnew), c(BetaI,GammaI,Sigma2I))$estimate, error=function(...) optim(c(BetaI,GammaI,Sigma2I),Max(Ynew,Xnew),method="BFGS")$par)
  
  #Final Estimate for ith generated data set, DL method
  Ests3[i,] <- vals
  
  
  
  
  ######Sim #4: Complete Case Method#######
  
  Ynew <- Y
  Xnew <- X
  
  #throwing out individuals with censored covariate value
  Xnew <- Xnew[-which(Xnew[,2]==d),]    
  Ynew <- Ynew[-which(X[,2]==d)] 
  Xnew <- cbind(Xnew,1)
  
  #Obtaining estimates for logistic and AFT regression (could be found separately)	
  vals <- tryCatch(nlm(Max(Ynew,Xnew), c(BetaI,GammaI,Sigma2I))$estimate, error=function(...) optim(c(BetaI,GammaI,Sigma2I),Max(Ynew,Xnew),method="BFGS")$par)
  
  #Final Estimate for ith generated data set, Complete Case Method
  Ests4[i,] <- vals
  
}


