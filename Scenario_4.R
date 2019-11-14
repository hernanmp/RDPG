
rm(list = ls())

working_directory = "/Users/oscar/Documents/Network_code"
setwd(working_directory)



source("auxiliary_functions.R")
# 
source("functions.R")
source("BinarySegmentation.R")
source("localrefine.R")
library(MCMCpack)
# load packages
# ==========
#library('hdbinseg')
library('InspectChangepoint')
library('R.utils')
source("Util.R")
library(foreach)
library(doParallel)
library(matrixcalc)
library(Rmisc)
library(Matrix)
library(gSeg)
library(ade4)

n_grid  =  c(300,200,100)
T_grid =   c(150)

NMC =  25
ind_n  = 1
ind_t = 1
delt =  0.3
rho_a = 0.9


hausdorff_ks = array(0,c(length(n_grid),NMC,  length(T_grid)))  #matrix(0,NMC,  length(T_grid))
hausdorff2_ks =  array(0,c(length(n_grid),NMC,  length(T_grid)))
error_ks =  array(0,c(length(n_grid),NMC,  length(T_grid)))

hausdorff_op = array(0,c(length(n_grid),NMC,  length(T_grid)))  #matrix(0,NMC,  length(T_grid))
hausdorff2_op =  array(0,c(length(n_grid),NMC,  length(T_grid)))
error_op =  array(0,c(length(n_grid),NMC,  length(T_grid)))

hausdorff_MNBS = array(0,c(length(n_grid),NMC,  length(T_grid)))  #matrix(0,NMC,  length(T_grid))
hausdorff2_MNBS =  array(0,c(length(n_grid),NMC,  length(T_grid)))
error_MNBS =  array(0,c(length(n_grid),NMC,  length(T_grid)))

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
      x = array(0,c(T,n,5))
      data  = matrix(0,n^2, T)
      E =  list()
      for(t in 1:T)
      {
        ###
        if(t ==1   || t == v[2]+1 )
        {
          x[t,,] = rdirichlet(n,rep(1,5))
          P =   drop( x[t,,]%*% t(x[t,,] ))
          
          A = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),P),n,n)
          aux = drop(A)
          aux[lower.tri(aux)] = aux[lower.tri(aux)]  
          diag(aux) = 0
          
          data[,t] = drop(matrix(A,n^2,1 ))
          
          temp = svd(A)
          xhat[t,,] =  temp$u[,1:d] %*%  diag( sqrt(temp$d[1:d]) )
          phat =   xhat[t,,] %*% t(xhat[t,,])
        }
        
        if( (t > 1      && t <=  v[1])  ||  ( t >v[2]+1 ) )
        {
          for(i in 1:n)
          {
            coin = rbinom(1,1,rho_a)
            
            if(coin ==1 )
            {
              x[t,i,] = x[t-1,i,] 
            }
            if(coin ==0 )
            {
              x[t,i,] =  rdirichlet(1,rep(1,5))
            }
          }
          ######
          P=   drop(x[t,,]%*% t(x[t,,]  ))
          aux = P
          diag(aux) = 0
          ind =  which(aux>1)
          aux[ind] = 1
          P = aux
          #P=   #0.6*exp(temp)/(1+exp(temp)    )
          diag(P) = 0
          
          A = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),P),n,n)
          aux = drop(A)
          aux[lower.tri(aux)] = aux[lower.tri(aux)]  
          diag(aux) = 0
          
          data[,t] = drop(matrix(A,n^2,1 ))
          
          temp = svd(A)
          xhat[t,,] =  temp$u[,1:d] %*%  diag( sqrt(temp$d[1:d]) )
        }
        
        if(t==v[1] +1)
        {
          x[t,,] =  rdirichlet(n,rep(1,5))
          x[t,1:floor(n*delt),] =  rdirichlet(floor(n*delt),rep(500,5))
          P =   x[t,,]%*%t(x[t,,])
          #0.6*matrix(rbeta(n^2,100,100),n,n)
          diag(P) = 0
          #mean(P)
          A = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),P),n,n)
          aux = drop(A)
          aux[lower.tri(aux)] = aux[lower.tri(aux)]  
          diag(aux) = 0
          
          data[,t] = drop(matrix(A,n^2,1 ))
          
          temp = svd(A)
          xhat[t,,] =  temp$u[,1:d] %*%  diag( sqrt(temp$d[1:d]) )
        }
        if(t >v[1]+1  &&  t <=v[2])
        {
          for(i in 1:n)
          {
            coin = rbinom(1,1,rho_a)   
            
            if(coin ==1)
            {
              x[t,i,] = x[t-1,i,]
            }
            if(coin ==0)
            {
              x[t,i,] = rdirichlet(1,rep(1,5))
              if(i<= floor(n*delt))
              {
                x[t,i,] = rdirichlet(1,rep(500,5))
              }
              #        x[t,,] =  rdirichlet(n,rep(1,5))
              #x[t,1:floor(n*delt),] =  rdirichlet(n,rep(50,5))
            }
          }
          
          P =   x[t,,]%*%t(x[t,,])
          diag(P) = 0
          
          A = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),P),n,n)
          aux = drop(A)
          aux[lower.tri(aux)] = aux[lower.tri(aux)]  
          diag(aux) = 0
          
          data[,t] = drop(matrix(A,n^2,1 ))
          
          temp = svd(A)
          xhat[t,,] =  temp$u[,1:d] %*%  diag( sqrt(temp$d[1:d]) )
          
        }
        E[[t]] =  A
      }
      ####for t
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
      
      #############################################
      ##MNBS
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
        #E1 <- E[(t-h+1):t]
        #E2 <- E[(t+1):(t+h)]
        #p_MNBS_1 <- mnbs(beta,E1,C=C0)$P 
        #p_MNBS_2 <- mnbs(beta,E2,C=C0)$P
        #d2_MNBS[t] <- (d2infty_norm(p_MNBS_1-p_MNBS_2))^2 # d2infty norm
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



# pdf("example4_before.pdf")
# image(A,axes =FALSE)
# #axis(2)
# axis(1, at=seq(0,1.0,0.25), labels=c(seq(0,1.0,0.25)*300))
# axis(2, at=seq(0,1.0,0.25), labels=c(seq(0,1.0,0.25)*300))
# dev.off()
# # #
# pdf("example4_after.pdf")
# image(A,axes =FALSE)
# #axis(2)
# axis(1, at=seq(0,1.0,0.25), labels=c(seq(0,1.0,0.25)*300))
# axis(2, at=seq(0,1.0,0.25), labels=c(seq(0,1.0,0.25)*300))
# dev.off()