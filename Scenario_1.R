
rm(list = ls())

working_directory = "/Users/oscar/Documents/Network_code"
setwd(working_directory)


source("auxiliary_functions.R")
# 
source("functions.R")
source("BinarySegmentation.R")
source("localrefine.R")

# load packages
# ==========
#library('hdbinseg')
library('InspectChangepoint')
library('R.utils')
source("Util.R")
library(foreach)
library(doParallel)
library(matrixcalc)
#library(Rmisc)
library(Matrix)
#library(gSeg)
#library(ade4)

n_grid  =  c(300,200,100)
T_grid =   c(150)

NMC = 100
ind_n  = 1
ind_t = 1

hausdorff_ks = array(0,c(length(n_grid),NMC,  length(T_grid)))  #matrix(0,NMC,  length(T_grid))
hausdorff2_ks =  array(0,c(length(n_grid),NMC,  length(T_grid)))
error_ks =  array(0,c(length(n_grid),NMC,  length(T_grid)))

hausdorff_op = array(0,c(length(n_grid),NMC,  length(T_grid)))  #matrix(0,NMC,  length(T_grid))
hausdorff2_op =  array(0,c(length(n_grid),NMC,  length(T_grid)))
error_op =  array(0,c(length(n_grid),NMC,  length(T_grid)))

hausdorff_MNBS = array(0,c(length(n_grid),NMC,  length(T_grid)))  #matrix(0,NMC,  length(T_grid))
hausdorff2_MNBS =  array(0,c(length(n_grid),NMC,  length(T_grid)))
error_MNBS =  array(0,c(length(n_grid),NMC,  length(T_grid)))

rho_a = 0.9

for(ind_n in 1:length(n_grid))
{
  n = n_grid[ind_n]
  
  print("n")
  print(n)
  
  for(ind_t in 1: length(T_grid) )
  {
    T = T_grid[ind_t]
    
    print("T")
    print(T)
    #####
    v =  c( floor(T/3)+1,2*floor(T/3)+1 )
    
    
    for(iter  in 1:NMC)
    {
      ### generate data
      d= 10
      xhat =  array(0,c(T,n,d))
      
      data  = matrix(0,n^2, T)
      
      E =  list()
      for(t in 1:T)
      {
        ###
        if(t==1 ||  t== v[2]+1 )
        {
          P =  matrix(0.3,n,n)
          P[1:floor(n/4), 1:floor(n/4)] = 0.5
          P[(1+floor(n/4)):(2*floor(n/4)),(1+floor(n/4)):(2*floor(n/4)) ] = 0.5
          P[(1+2*floor(n/4)):(3*floor(n/4)),(1+2*floor(n/4)):(3*floor(n/4)) ] = 0.5
          P[(1+3*floor(n/4)):n,(1+3*floor(n/4)):n ] = 0.5
          diag(P) = 0
          
          
          A = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),P),n,n)
          aux = drop(A)
          aux[lower.tri(aux)] = aux[lower.tri(aux)]  
          diag(aux) = 0
          
          data[,t] = drop(matrix(A,n^2,1 ))
          
          temp = svd(A)
          xhat[t,,] =  temp$u[,1:d] %*%  diag( sqrt(temp$d[1:d]) ) 
        }
        if( (t > 1      && t <=  v[1])  ||  ( t >v[2]+1 ) )
        {
           aux1 = P +  (1-P)*rho_a
           aux2 = P*(1-rho_a)
           
           aux1 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux1),n,n)
           aux2 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux2),n,n)
           A =  aux1*A + aux2*(1-A)
           
           aux = drop(A)
           aux[lower.tri(aux)] = aux[lower.tri(aux)]  
           diag(aux) = 0
           
           data[,t] = drop(matrix(A,n^2,1 ))
           
           temp = svd(A)
           xhat[t,,] =  temp$u[,1:d] %*%  diag( sqrt(temp$d[1:d]) )
        }
        if(t ==  v[1]+1)
        {
          Q =  matrix(0.2,n,n)
          Q[1:floor(n/4), 1:floor(n/4)] = 0.45
          Q[(1+floor(n/4)):(2*floor(n/4)),(1+floor(n/4)):(2*floor(n/4)) ] = 0.45
          Q[(1+2*floor(n/4)):(3*floor(n/4)),(1+2*floor(n/4)):(3*floor(n/4)) ] = 0.45
          Q[(1+3*floor(n/4)):n,(1+3*floor(n/4)):n ] = 0.45
          diag(Q) = 0
        
          A = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),Q),n,n)
          aux = drop(A)
          aux[lower.tri(aux)] = aux[lower.tri(aux)]  
          diag(aux) = 0
          data[,t] = drop(matrix(A,n^2,1 ))
          
          temp = svd(A)
          xhat[t,,] =  temp$u[,1:d] %*%  diag( sqrt(temp$d[1:d]) ) 
        }
        if(t > v[1]+1 && t <= v[2] )
        {
          aux1 = Q +  (1-Q)*rho_a
          aux2 = Q*(1-rho_a)
          
          aux1 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux1),n,n)
          aux2 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux2),n,n)
          A =  aux1*A + aux2*(1-A)
          
          aux = drop(A)
          aux[lower.tri(aux)] = aux[lower.tri(aux)]  
          diag(aux) = 0
          
          data[,t] = drop(matrix(A,n^2,1 ))
          
          temp = svd(A)
          xhat[t,,] =  temp$u[,1:d] %*%  diag( sqrt(temp$d[1:d]) )
        }
        ##################
       
        E[[t]] = A
      }####for t
      ########
      ################ WBS
      
      yhat = matrix(0, T, floor(n/2)   )
      
      
      for(t in 1:T)
      {
        phat =    drop(xhat[t,,] %*% t(xhat[t,,]))
        ind =  sample(1:n,n, replace = FALSE)
        #aux = phat[ind[2*(1:floor(n/2))], ind[2*(1:floor(n/2)) -1]  ]
        
        for(i in 1:floor(n/2) )
        {
          yhat[t,i ] = phat[ind[2*i], ind[2*i-1] ]
        }
      }
      
      
      M =   120
      alpha =  sample.int(size =M  , n = T,replace = TRUE)
      beta =   sample.int(size =M  , n = T,replace = TRUE)#alpha + floor((T- alpha)*runif(M))
      #
      for(j in 1:M)
      {
        aux =  alpha[j]
        aux2 =  beta[j]
        #
        alpha[j] = min(aux,aux2)
        beta[j] = max(aux,aux2)
      }
      grid =   seq(min(yhat),max(yhat), length = M)
      N =  rep(floor(n/2),T)
      S = NWBS_estimator(yhat[,1:floor(n/2)],grid,Ntau =20,gam=10,N,alpha,beta)
      #Model_selection_V4(yhat[,1:floor(n/2)],grid,Ntau =20,gam=10,N,alpha,beta)
      S =  sort(unique(S))
      hausdorff_ks[ind_n,iter,ind_t] =  dist_change_points(S,v)
      hausdorff2_ks[ind_n,iter,ind_t] =  dist_change_points(v,S)
      error_ks[ind_n,iter,ind_t] =  length(v) - length(S) 
      
      #################
      data = data[gen.lower.coordinate(n), ]
      data.1 = data[, seq(2, T, 2)]
      data.2 = data[, seq(1, T - 1, 2)]
      
      rho.hat = quantile(rowSums(data) / T, 0.95)
      population.change = c(1, 1 : T)
      nbs = Binary.Segmentation(data.1, data.2, tau = n * rho.hat * log(T)^2/20)
      
      hausdorff_op[ind_n,iter,ind_t] =  dist_change_points(nbs,v)
      hausdorff2_op[ind_n,iter,ind_t] =  dist_change_points(v,nbs)
      error_op[ind_n,iter,ind_t] =  length(v) - length(nbs) 
      
      #################################3
      ###MNBS
      ##  The following  needs  code  from Zhang et. al 2019 
      #D0 <- 0.25 # scaling parameter for threshold
      #beta <- 0.5 # MNBS omega D0*(log(n))^(1+delta)/n^0.5/h^beta
      #C0 <- 3 # neighborhood size for MNBS C*log(n)/n^0.5/h^beta
      #delta0 <- 0.1 # MNBS thresholding
      #h <- max(10, floor(sqrt(T))) # window size for SaRa
      #threshold_MNBS <- D0*(log(n))^(1/2+delta0)/n^0.5/h^beta # threshold
      
      #T0 <- h
      #T1 <- T-h
      #d2_MNBS <- rep(0, T)
      #for (t in T0:T1){
       # E1 <- E[(t-h+1):t]
       # E2 <- E[(t+1):(t+h)]
       # p_MNBS_1 <- mnbs(beta,E1,C=C0)$P 
       # p_MNBS_2 <- mnbs(beta,E2,C=C0)$P
       # d2_MNBS[t] <- (d2infty_norm(p_MNBS_1-p_MNBS_2))^2 # d2infty norm
      #}
      #MNBS_cp <- SaRa(d2_MNBS, h, threshold_MNBS)
      
      
      #hausdorff_MNBS[ind_n,iter,ind_t] =  dist_change_points(MNBS_cp,v)
      #hausdorff2_MNBS[ind_n,iter,ind_t] =  dist_change_points(v,MNBS_cp)
      #error_MNBS[ind_n,iter,ind_t] =  length(v) - length(MNBS_cp) 
      
      ## print
      print("iter")
      
      print("ks")
      print(median(hausdorff_ks[ind_n,iter,ind_t]))
      print(median(hausdorff2_ks[ind_n,iter,ind_t]))
      print(mean(abs(error_ks[ind_n,iter,ind_t])))
      
      print("op")
      print(median(hausdorff_op[ind_n,iter,ind_t]))
      print(median(hausdorff2_op[ind_n,iter,ind_t]))
      print(mean(abs(error_op[ind_n,iter,ind_t])))
      
      print("MNBS")
      print(median(hausdorff_MNBS[ind_n,iter,ind_t]))
      print(median(hausdorff2_MNBS[ind_n,iter,ind_t]))
      print(mean(abs(error_MNBS[ind_n,iter,ind_t])))
    }##### close for iter
    
    print("Averaging")
    print("ks")
    print(median(hausdorff_ks[ind_n,,ind_t]))
    print(median(hausdorff2_ks[ind_n,,ind_t]))
    print(mean(abs(error_ks[ind_n,,ind_t])))
    
    print("op")
    print(median(hausdorff_op[ind_n,,ind_t]))
    print(median(hausdorff2_op[ind_n,,ind_t]))
    print(mean(abs(error_op[ind_n,,ind_t])))
    
    print("MNBS")
    print(median(hausdorff_MNBS[ind_n,,ind_t]))
    print(median(hausdorff2_MNBS[ind_n,,ind_t]))
    print(mean(abs(error_MNBS[ind_n,,ind_t])))
  }###########
}
