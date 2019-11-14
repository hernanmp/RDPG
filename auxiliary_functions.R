
BS = function(y,tau, gam,s,e,flag,S, N)
{
  if(e-s <   1+  3*gam  ||  flag ==1)
  {
    S =  c(S,NULL)
    return(S);
  }
  else{
    
    ###  calculate  statistics
    a =  rep(0,e-s+1  )
    
    for(t  in    (s+1):(e-1)  )
    {
      
      a[t-s ] =  Delta_se_t(y,s,e,t,N)
    }
    
    
    best_t  =  which.max( a)
    
    if(a[best_t]  <tau  )
    {
      return(S)
    }
    print(s+best_t)
    best_t =   s+ best_t 
    S1 = BS(y,tau, gam,s,best_t-1,flag,S,N)
    S2 = BS(y,tau, gam,best_t+1,e,flag, S,N)
    
    S =  c(S,best_t)
    
    
    S  =   c(S,S1,S2)
    return(S)
  }
}



Delta_se_t = function(y,s,e,t,N)
{
  #T =   dim(y)[2]
  n =  dim(y)[2]
  
  n_st = sum(N[s:t])  #n*(t-s+1)
  n_se = sum(N[s:e])  #n*(e-s+1)
  n_te =sum(N[(t+1):e]) #n*(e-(t+1) +1)
  
  aux =  as.vector(y[s:t,])
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  vec_y =  as.vector(y[s:e,])
  vec_y = vec_y[which(is.na(vec_y)==FALSE)]
  Fhat_st =  temp(vec_y)# temp(grid)
  
  aux = y[(t+1):e,]
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  Fhat_te =  temp(vec_y)# temp(grid)
  
  temp =  sqrt( n_st*n_te / n_se   ) *max(abs(Fhat_te - Fhat_st  ))
  
  return(temp )
}

compute_F0 =   function(grid, p0 )
{
  cdf  = cumsum(p0)  
  loc  = (1:length(p0))/length(p0)
  F0 =  rep(0, length(grid))
  
  for(j in 1:length(grid))
  {
    ind =  which.min(abs(grid[j] - loc   )) 
    F0[j] =   cdf[ind[1]]
  }
  
  return(F0)
}

######################################


compute_AIC =  function(y,m=50,tau_grid,s,e,gam = 20,N)
{
 T =  dim(y)[1] 
 n =  dim(y)[2]
 min_y = min(y)
 max_y = max(y) 
 y = (y -min_y)/(max_y - min_y)
 
 AIC_score =  rep(0, length(tau_grid))
#  N = 100
  break_p =  (1:m)/m
  break_p = c(0,break_p)
  #mid =  (mid[2:N]+ mid[1:(N-1)])/2
  
  y_counts  = matrix(0, T, m)
  
  for(t in 1:T)
  {
    y_counts[t,] = hist(y[t,1:N[t]],breaks = break_p)$counts
  }
  
  
  for(j in 1:length(tau_grid))
  { 
    S = BS(y,tau_grid[j],gam,s,e,0,NULL,N)
    S =  sort(S)
   # print(S)
    
    par  =  matrix(0,length(S)+1,m)
    
    if(length(S)>0)
    {
      
      
      temp = hist(y[1:S[1],],breaks = break_p)$counts
      par[1,] =   temp/sum(temp)
      i = 1
      
      if(length(S)>1)
      {
        for(i in 1:(length(S)-1) )
        {
          temp = hist(y[(S[i]+1):S[i+1],],breaks = break_p)$counts
          par[i+1,] = temp/sum(temp)
        }
      }
      temp = hist(y[(S[length(S)]+1):T,],breaks = break_p)$counts
      par[length(S)+1,] =   temp/sum(temp)
      
    }
    ###################################
    minus_log_lik =   0
    if(length(S)==0)
    {
      temp = hist(y,breaks = break_p)$counts
      par[1,] =   temp/sum(temp)
      
      ind = which(par[1,] >0)
      temp  =  apply(y_counts,1,  function(x){ -sum(x[ind] *log(par[1,ind]) )   }    )
      minus_log_lik =   minus_log_lik + sum(temp)
    }
    if(length(S)>0)
    {
      ind = which(par[1,] >0)
      temp  =  apply(y_counts[1:S[1],],1,  function(x){ -sum(x[ind] *log(par[1,ind]) )   }    )
      minus_log_lik =   minus_log_lik + sum(temp)
      
      if(length(S)>1)
      {
        for(i in 1: (length(S)-1) )
        {
          ind = which(par[i+1,] >0)
          temp  =  apply(y_counts[(S[i]+1):S[i+1],  ],1,  function(x){ -sum(x[ind] *log(par[i+1,ind]) )   }    )
          minus_log_lik =   minus_log_lik + sum(temp)
        }
      }  
      
      ind = which(par[length(S)+1,] >0)
      temp  =  apply(y_counts[(S[length(S)]+1):T,],1,  function(x){ -sum(x[ind] *log(par[length(S)+1,ind]) )   }    )
      minus_log_lik =   minus_log_lik + sum(temp)
      
      
    }
    
    #########
    AIC_score[j] =  2*minus_log_lik  +   2*(  length(S)+1 )*(m-1)
  }
  return(AIC_score)
}

#######################################3

WBS =  function(y,tau, gam,s,e,flag,S, alpha,beta,N)
{
  
  
  alpha_new =  pmax(alpha,s)
  beta_new = pmin(beta,e)
  # 
  #  ind = which( beta_new- alpha_new > 1+  3*gam)
  #  alpha_new =  alpha_new[ind]
  #  beta_new=   beta_new[ind]
  # # 
  ind = which( beta_new- alpha_new >1)
  alpha_new =  alpha_new[ind]
  beta_new=   beta_new[ind]
  M =  length(alpha_new)

  xi = 1/8
  alpha_new2 =alpha_new
  beta_new2  =  beta_new
  alpha_new =   ceiling((1-xi)*alpha_new2+  xi*beta_new2)
  beta_new =    ceiling((1-xi)*beta_new2 +  xi*alpha_new2)
  ind = which( beta_new- alpha_new >1)
  alpha_new =  alpha_new[ind]
   beta_new=   beta_new[ind]
  M =  length(alpha_new)

 # print(S)
  
  if(M ==  0)
  {
    return(NULL) 
  }
  
  b  =  rep(0,M)
  a =  rep(0,M)
  
  #print(beta_new)
  #print(alpha_new)
#  print()
  for( m in 1:M  )
  {
     temp  =  rep(0,beta_new[m]-alpha_new[m]+1)
     for(t  in    (alpha_new[m]+1):(beta_new[m]-1 )  )
     {
       temp[t-(alpha_new[m]) ] =  Delta_se_t(y,alpha_new[m],beta_new[m],t,N)
     }
     best_ind  =  which.max(temp)
     a[m] =  alpha_new[m] +  best_ind
     b[m] =   temp[best_ind]
  }
  best_ind =  which.max(b)
  
  
#  print(b)
  if(b[best_ind]  <tau  )
  {
    return(S)
  }
 # print(S)
  S1 = WBS(y,tau, gam,s,a[best_ind]-1,flag, S,alpha, beta,N)
  S2 = WBS(y,tau, gam,a[best_ind]+1,e,flag, S,alpha, beta,N)
  
  S =  c(S,a[best_ind])
  
  
  S  =   c(S,S1,S2)
    
  
    return(S)
  
}



######################################


compute_AIC_wbs =  function(y,m=50,tau_grid,s,e,gam = 20,alpha, beta,N)
{
  T =  dim(y)[1] 
  n =  dim(y)[2]
  min_y = min(y)
  max_y = max(y) 
  y = (y -min_y)/(max_y - min_y)
  
  AIC_score =  rep(0, length(tau_grid))
  #  N = 100
  break_p =  (1:m)/m
  break_p = c(0,break_p)
  #mid =  (mid[2:N]+ mid[1:(N-1)])/2
  
  y_counts  = matrix(0, T, m)
  
  for(t in 1:T)
  {
    y_counts[t,] = hist(y[t,1:N[t]],breaks = break_p)$counts
  }
  
  
  for(j in 1:length(tau_grid))
  { 
    S = WBS(y,tau_grid[j],gam,s,e,0,NULL,alpha,beta, N)
      #WBS(y,tau_grid[j],gam,s,e,0,grid,NULL)
    S =  sort(S)
  #  print(S)
    
    par  =  matrix(0,length(S)+1,m)
    
    if(length(S)>0)
    {
      
      
      temp = hist(y[1:S[1],],breaks = break_p)$counts
      par[1,] =   temp/sum(temp)
      i = 1
      
      
      if(length(S)>1)
      {
        for(i in 1:(length(S)-1) )
        {
          temp = hist(y[(S[i]+1):S[i+1],],breaks = break_p)$counts
          par[i+1,] = temp/sum(temp)
        }
      }
      temp = hist(y[(S[length(S)]+1):T,],breaks = break_p)$counts
      par[length(S)+1,] =   temp/sum(temp)
      
    }
    ###################################
    minus_log_lik =   0
    if(length(S)==0)
    {
      temp = hist(y,breaks = break_p)$counts
      par[1,] =   temp/sum(temp)
      
      ind = which(par[1,] >0)
      temp  =  apply(y_counts,1,  function(x){ -sum(x[ind] *log(par[1,ind]) )   }    )
      minus_log_lik =   minus_log_lik + sum(temp)
    }
    if(length(S)>0)
    {
      ind = which(par[1,] >0)
      temp  =  apply(y_counts[1:S[1],],1,  function(x){ -sum(x[ind] *log(par[1,ind]) )   }    )
      minus_log_lik =   minus_log_lik + sum(temp)
      
      if(length(S)>1)
      {
        for(i in 1: (length(S)-1) )
        {
          ind = which(par[i+1,] >0)
          temp  =  apply(y_counts[(S[i]+1):S[i+1],  ],1,  function(x){ -sum(x[ind] *log(par[i+1,ind]) )   }    )
          minus_log_lik =   minus_log_lik + sum(temp)
        }
      }  
      
      ind = which(par[length(S)+1,] >0)
      temp  =  apply(y_counts[(S[length(S)]+1):T,],1,  function(x){ -sum(x[ind] *log(par[length(S)+1,ind]) )   }    )
      minus_log_lik =   minus_log_lik + sum(temp)
      
      
    }
    
    #########
    AIC_score[j] =  2*minus_log_lik  +   2*(  length(S)+1 )*(m-1)
  }
  return(AIC_score)
}

#############################################3


mean_WBS  = function(y,tau, gam,s,e,alpha,beta,flag =0,S= NULL)
{
  
  alpha_new =  pmax(alpha,s)
  beta_new = pmin(beta,e)
  
  ind = which( beta_new- alpha_new > 1+  3*gam)
  alpha_new =  alpha_new[ind]
  beta_new=   beta_new[ind]
  M =  length(alpha_new)
  
  #print(S)
  
  if(M ==  0)
  {
    return(NULL) 
  }
  
  b  =  rep(0,M)
  a =  rep(0,M)
  
  
  for( m in 1:M  )
  {
    temp  =  rep(0,beta_new[m]-alpha_new[m]+1)
    for(t  in    (alpha_new[m]+gam+1):(beta_new[m]-  gam-1 )  )
    {
      temp[t] =  Delta_mean_se_t(y,alpha_new[m],beta_new[m],t)
    }
    best_ind  =  which.max(temp)
    a[m] =  best_ind
    b[m] =   temp[best_ind]
  }
  best_ind =  which.max(b)
  
  
  #  print(b)
  if(b[best_ind]  <tau  )
  {
    return(S)
  }
  
  S1 = mean_WBS(y,tau, gam,s,a[best_ind]-1,flag,S)
  S2 = mean_WBS(y,tau, gam,a[best_ind]+1,e,flag, S)
  
  S =  c(S,a[best_ind])
  
  
  S  =   c(S,S1,S2)
  
  
  return(S)
}  

########################################################



Delta_mean_se_t = function(y,s,e,t)
{
  T =   dim(y)[2]
  
  
  n_st =  (t-s+1)
  n_se =  (e-s+1)
  n_te = (e-(t+1) +1)
  
 
  temp  = sqrt(n_te/(  n_st*  n_se )  )*sum(y[s:t,1])     -  sqrt(n_st/(  n_te*  n_se )  )*sum(y[(t+1):e,1])
  temp =  abs(temp)  
  return(temp )
}



mean_BS  = function(y,tau, gam,s,e,flag =0,S = NULL)
{
  if(e-s <   1+  3*gam  ||  flag ==1)
  {
    S =  c(S,NULL)
    return(S);
  }
  else{
    
    ###  calculate  statistics
    a =  rep(0,e-s+1  )
    
    for(t  in    (s+gam+1):(e-  gam-1 )  )
    {
      
      a[t] =  Delta_mean_se_t(y,s,e,t)
    }
    
    
    best_t  =  which.max( a)
    
    if(a[best_t]  <tau  )
    {
      return(S)
    }
    
    
    S1 = mean_BS(y,tau, gam,s,best_t-1,flag,S)
    S2 = mean_BS(y,tau, gam,best_t+1,e,flag,S)
    
    S =  c(S,best_t)
    
    
    S  =   c(S,S1,S2)
    return(S)
  }
}
############################33

L_obj =  function(y,S,grid,N)
{
  T = dim(y)[1]
   y_vec =  as.vector(y)
   y_vec = y_vec[which(is.na(y_vec)==FALSE)]
  # 
   K =  length(S)
  # 
   val = rep(0,length(grid)) 
    
      for(k in 0:(K))
      {
        if(k== 0)
        {
          s =  1
          if(K>0)
          {
            e =  S[1]
          }
          if(K==0)
          {
            e = T
          }
        }
        ########
        ##3
        if(k   == K && k > 0)
        {
          s =  S[K]+1
          e =  T
        }
        ##############3
        if( 0 < k &&  k < K)
        {
          s= S[k]+1
          e= S[k+1]
        }
        ######################33
        aux =  as.vector(y[s:e,])
        aux =  aux[which(is.na(aux)==FALSE)]
        Fhat = ecdf(  aux )
        F_hat_z =  Fhat(grid)
        
        for(j in 1:length(grid))
        {
          aux = y[s:e,] - grid[j]
          aux = as.vector(aux)
          aux =  aux[which(is.na(aux)==FALSE)]
          aux2 =  which(aux > 0)
          aux3 = which(aux <= 0)
          aux[aux2] = 0
          aux[aux3] = 1
          aux  =  aux - F_hat_z[j]
          aux = sum(aux^2)
          val[j] =  val[j] + aux
        }
      }
     return(max(val) + median(grid)*15*log(sum(N))*length(S) )
}
########################################


L_obj_V2 =  function(y,z,S,grid,N)
{
  T = dim(y)[1]
  y_vec =  as.vector(y)
  y_vec = y_vec[which(is.na(y_vec)==FALSE)]
  
  z_vec =  as.vector(z)
  z_vec = z_vec[which(is.na(z_vec)==FALSE)]
  # 
  K =  length(S)
  # 
  val = rep(0,length(grid)) 
  cost = 0
  
  for(k in 0:(K))
  {
    if(k== 0)
    {
      s =  1
      if(K>0)
      {
        e =  S[1]
      }
      if(K==0)
      {
        e = T
      }
    }
    ########
    ##3
    if(k   == K && k > 0)
    {
      s =  S[K]+1
      e =  T
    }
    ##############3
    if( 0 < k &&  k < K)
    {
      s= S[k]+1
      e= S[k+1]
    }
    ######################33
    aux =  as.vector(z[s:e,])
    aux =  aux[which(is.na(aux)==FALSE)]
    # aux = mean(aux)
    # aux = sum((aux - y[s:e,])^2)
    # cost =  cost +  aux
    Fhat = ecdf(  aux )
    F_hat_z =  Fhat(grid)

    for(j in 1:length(grid))
    {
      aux = y[s:e,] - grid[j]
      aux = as.vector(aux)
      aux =  aux[which(is.na(aux)==FALSE)]
      aux2 =  which(aux > 0)
      aux3 = which(aux <= 0)
      aux[aux2] = 0
      aux[aux3] = 1
      aux  =  aux-F_hat_z[j]
      aux = sum(aux^2)
      val[j] =  val[j] + aux
    }
    
  }
  return(max(val))#+ .015*log(sum(N))*length(S))
  #return(max(cost))
  # + .15*log(sum(N))*length(S) )
}



#######################################


L_obj2 =  function(y,S,grid,N)
{
  T = dim(y)[1]
  y_vec =  as.vector(y)
  y_vec = y_vec[which(is.na(y_vec)==FALSE)]
  # 
  K =  length(S)
  # 
  val = rep(0,length(grid)) 
  
  for(k in 0:(K))
  {
    if(k== 0)
    {
      s =  1
      if(K>0)
      {
        e =  S[1]
      }
      if(K==0)
      {
        e = T
      }
    }
    ########
    ##3
    if(k   == K && k > 0)
    {
      s =  S[K]+1
      e =  T
    }
    ##############3
    if( 0 < k &&  k < K)
    {
      s= S[k]+1
      e= S[k+1]
    }
    ######################33
    aux =  as.vector(y[s:e,])
    aux =  aux[which(is.na(aux)==FALSE)]
    Fhat = ecdf(  aux )
    F_hat_z =  Fhat(grid)
    
    for(j in 1:length(grid))
    {
      aux = y[s:e,] - grid[j]
      aux = as.vector(aux)
      aux =  aux[which(is.na(aux)==FALSE)]
      aux2 =  which(aux > 0)
      aux3 = which(aux <= 0)
      aux[aux2] = 0
      aux[aux3] = 1
      aux  =  aux - F_hat_z[j]
      aux = sum(aux^2)
      val[j] =  val[j] + aux
    }
  }
  return(max(val) + .08*log(sum(N))*length(S) )
}
############################################

Model_selection=  function(y,tau_grid,gam,alpha, beta,N)
{
  n = max(N)
  T = dim(y)[1]
  y1 =  matrix(NA, T,floor(n/2) )
  y2 = matrix(NA, T, floor(n/2))
  
  for(t in 1:T)
  {
    y1[t,1:floor(N[t]/2)] =  y[t,1:floor(N[t]/2)]
    
    y2[t,1:floor(N[t]/2)] =  y[t,(1+floor(N[t]/2)):(2* floor(N[t]/2) ) ]
  }
  S =  list()
  for(j in 1:length(tau_grid))
  {
    temp = WBS(y1,tau = tau_grid[j],gam,1,T,0,NULL,alpha, beta,floor(N/2))
  
    
    if(length(temp)==0)
    {
      break;
    }
    S[[j]] = sort(temp)
   # print(S[[j]])
    #print(L_obj(y2,S[[j]],grid,floor(N/2)))
  }
  
  score =  rep(0,length(S)+1)
  if(length(S)>0)
  {
    for(j in 1:length(S))
    {
      score[j] = L_obj(y2,S[[j]],grid,floor(N/2))
    }
  }
  score[j+1] = L_obj(y2,NULL,grid,floor(N/2))
  best_ind = which.min(score)
  
  if(length(S) < best_ind)
  {
    #print(L_obj(y2,NULL,grid,floor(N/2)))
    return(NULL)
  }
  return(S[[best_ind]])
#  return()
}
######################################3

Model_selection_V2=  function(y,tau_grid,gam,alpha, beta,N)
{
   
  n = max(N)
  T = dim(y)[1]
  y1 =  matrix(NA, T,floor(n/2) )
  y2 = matrix(NA, T, floor(n/2))
 
  for(t in 1:T)
  {
    y1[t,1:floor(N[t]/2)] =  y[t,1:floor(N[t]/2)]
    
    y2[t,1:floor(N[t]/2)] =  y[t,(1+floor(N[t]/2)):(2* floor(N[t]/2) ) ]
  }
 #grid =  grid[which(grid < min(max(y1),max(y2)))]
  #grid =  grid[which(grid > max(min(y1),min(y2)))] 
  S =  list()
  for(j in 1:length(tau_grid))
  {
    temp = WBS(y1,tau = tau_grid[j],gam,1,T,0,NULL,alpha, beta,floor(N/2))
    
    
    if(length(temp)==0)
    {
      break;
    }
    S[[j]] = sort(temp)
    # print(S[[j]])
    #print(L_obj(y2,S[[j]],grid,floor(N/2)))
  }
  
  score =  rep(0,length(S)+1)
  if(length(S)>0)
  {
    for(j in 1:length(S))
    {
      score[j] = L_obj_V2(y2,y1,S[[j]],grid,floor(N/2))# +  2*length(S[[j]])
    }
  }
  score[j+1] = L_obj_V2(y2,y1,NULL,grid,floor(N/2))
  #best_ind = which.min(score)
  
  ####################################
  S =  list()
  for(j in 1:length(tau_grid))
  {
    temp = WBS(y2,tau = tau_grid[j],gam,1,T,0,NULL,alpha, beta,floor(N/2))
    
    
    if(length(temp)==0)
    {
      break;
    }
    S[[j]] = sort(temp)
    # print(S[[j]])
    #print(L_obj(y2,S[[j]],grid,floor(N/2)))
  }
  
  score2 =  rep(0,length(S)+1)
  if(length(S)>0)
  {
    for(j in 1:length(S))
    {
      score2[j] = L_obj_V2(y1,y2,S[[j]],grid,floor(N/2))# +  2*length(S[[j]])
    }
  }
  score2[j+1] = L_obj_V2(y1,y2,NULL,grid,floor(N/2))
  
  len =  min(length(score),length(score2) )
  score3 =  score[1:len] + score2[1:len]
  best_ind = which.min(score3)
  
  # T = 2*T
  # y = rep(0,T)
  # y[2*(1:(T/2))-1] = y1
  # y[2*(1:(T/2))] = y2
  # N_old = rep(0,T)
  # N_old[2*(1:(T/2))-1] = floor(N/2)
  # N_old[2*(1:(T/2))] = floor(N/2)
  S =  WBS(y1,tau = tau_grid[best_ind],gam,1,T,0,NULL,alpha, beta,floor(N/2))
  S2 =  WBS(y2,tau = tau_grid[best_ind],gam,1,T,0,NULL,alpha, beta,floor(N/2))
  if(length(S2)< length(S))
  {
    return(S2)
  }
  return(S)
  ######################################
  
  if(length(S) < best_ind)
  {
    #print(L_obj(y2,NULL,grid,floor(N/2)))
    return(NULL)
  }
  return(S[[best_ind]])
  #  return()
}




###################################################
Model_selection_bs =  function(y,tau_grid,gam,N)
{
  n = max(N)
  T = dim(y)[1]
  y1 =  matrix(NA, T,floor(n/2) )
  y2 = matrix(NA, T, floor(n/2))
  
  for(t in 1:T)
  {
    y1[t,1:floor(N[t]/2)] =  y[t,1:floor(N[t]/2)]
    
    y2[t,1:floor(N[t]/2)] =  y[t,(1+floor(N[t]/2)):(2* floor(N[t]/2) ) ]
  }
  S =  list()
  for(j in 1:length(tau_grid))
  {
    temp = BS(y1,tau = tau_grid[j],gam,1,T,0,NULL,floor(N/2))
    
    
    if(length(temp)==0)
    {
      break;
    }
    S[[j]] = sort(temp)
    #print(S[[j]])
    #print(L_obj(y2,S[[j]],grid,floor(N/2)))
  }
  
  score =  rep(0,length(S)+1)
  if(length(S)>0)
  {
    for(j in 1:length(S))
    {
      score[j] = L_obj2(y2,S[[j]]-1,grid,floor(N/2))
    }
  }
  score[j+1] = L_obj2(y2,NULL,grid,floor(N/2))
  best_ind = which.min(score)
  
  if(length(S) < best_ind)
  {
   # print(L_obj(y2,NULL,grid,floor(N/2)))
    return(NULL)
  }
  return(S[[best_ind]])
  #  return()
}
#####################################################################3



Model_selection_bs_V2 =  function(y,tau_grid,gam,N)
{
  n = max(N)
  T = dim(y)[1]
  y1 =  matrix(NA, T,floor(n/2) )
  y2 = matrix(NA, T, floor(n/2))
  
  for(t in 1:T)
  {
    y1[t,1:floor(N[t]/2)] =  y[t,1:floor(N[t]/2)]
    
    y2[t,1:floor(N[t]/2)] =  y[t,(1+floor(N[t]/2)):(2* floor(N[t]/2) ) ]
  }
  #grid =  grid[which(grid < min(max(y1),max(y2)))]
  #grid =  grid[which(grid > max(min(y1),min(y2)))] 
  S =  list()  
  for(j in 1:length(tau_grid))
  {
    temp = BS(y1,tau = tau_grid[j],gam,1,T,0,NULL,floor(N/2))
    
    
    if(length(temp)==0)
    {
      break;
    }
    S[[j]] = sort(temp)
    #print(S[[j]])
    #print(L_obj(y2,S[[j]],grid,floor(N/2)))
  }
  
  score =  rep(0,length(S)+1)
  if(length(S)>0)
  {
    for(j in 1:length(S))
    {
      score[j] = L_obj_V2(y2,y1,S[[j]],grid,floor(N/2))
    }
  }
  score[j+1] = L_obj_V2(y2,y1,NULL,grid,floor(N/2))
  #best_ind = which.min(score)
  
  ####################################
  S =  list()
  for(j in 1:length(tau_grid))
  {
    temp = BS(y2,tau = tau_grid[j],gam,1,T,0,NULL,floor(N/2))
    
    
    if(length(temp)==0)
    {
      break;
    }
    S[[j]] = sort(temp)
    # print(S[[j]])
    #print(L_obj(y2,S[[j]],grid,floor(N/2)))
  }
  
  score2 =  rep(0,length(S)+1)
  if(length(S)>0)
  {
    for(j in 1:length(S))
    {
      score2[j] = L_obj_V2(y1,y2,S[[j]],grid,floor(N/2))# +  2*length(S[[j]])
    }
  }
  score2[j+1] = L_obj_V2(y1,y2,NULL,grid,floor(N/2))
  
  len =  min(length(score),length(score2) )
  score3 =  score[1:len] + score2[1:len]
  best_ind = which.min(score3)
  
  #y =  rep(0,2*T)
  # T = 2*T
  # y = rep(0,T)
  # y[2*(1:(T/2))-1] = y1
  # y[2*(1:(T/2))] = y2
  # N_old = rep(0,T)
  # N_old[2*(1:(T/2))-1] = floor(N/2)
  # N_old[2*(1:(T/2))] = floor(N/2)
  S =  BS(y1,tau = tau_grid[best_ind],gam,1,T,0,NULL,floor(N/2))
  S2 =  BS(y2,tau = tau_grid[best_ind],gam,1,T,0,NULL,floor(N/2))
  
  if(length(S2)< length(S))
  {
    return(S2)
  }
  
  return(S)
  
  #################################
  if(length(S) < best_ind)
  {
    # print(L_obj(y2,NULL,grid,floor(N/2)))
    return(NULL)
  }
  return(S[[best_ind]])
  #  return() 
}


#############################################################3
radiological_example =  function(K, celsium  = 1, T)
{
  working_directory = "C:/Users/oscar/Desktop/Anomaly_detection/supplementary_materials/supplementary_materials/Experiments"
  
  
  setwd(paste(working_directory,"/Code_J",sep=""))
  source('utils.R')
  setwd(paste(working_directory,"/Data_J",sep=""))
  
  # How strong is the background rate, and how large is the source in milliCuries?
  background_cps = 39  
  source_size = 100#650  #mCi
  # Distance to the source
  dist_grid = 20*(1:K)
  i = 1## index of distance used, if for example i = 1, then distance to source is 50m
  
  ##################################3
  
  
 # K =  length(dist_grid)
  m = 2048
  
  f0 =    matrix(0, K+1,m)
  
  
  #f0[2,] =  sim_spectrum
  #f0[3,] =  sim_spectrum2
  
  
  ###########################
  
  for(i in 1:K)
  {
    
    # Read in background
    background_train = as.numeric(read.csv('background_train_mean.csv', header=FALSE))
    plot(background_train, type='l')
    n_bins = length(background_train)
    
    # choose anomaly
    #celsium  = 1 # if celsium = 1 then the anomaly is celsium is celsium-137. If set to zero, we use cobalt
    
    # Read in anomaly
    
    if(celsium==1)
    {
      cs137_spectrum = as.numeric(read.csv('2013-cs137-5cm.csv', header=FALSE))
    }
    if(celsium ==0)
    {
      cs137_spectrum = as.numeric(read.csv('co60-4min-320cps-c7-20cm.csv', header=FALSE))
    }
    
    
    # Winsorize: assign all counts after bin 2048 to bin 2048, and then truncate
    cs137_spectrum[2048] = sum(cs137_spectrum[2048:4096])
    cs137_spectrum = head(cs137_spectrum, 2048)
    
    
    # Normalize and plot
    cs137_spectrum = cs137_spectrum/sum(cs137_spectrum)
    plot(cs137_spectrum, type='l')
    
    
    # Create composite spectrum from background + anomaly
    this_dist =  dist_grid[i]
    sim_spectrum = inject_source(background_train, background_cps, cs137_spectrum, source_size, this_dist)
    sim_spectrum = sim_spectrum/sum(sim_spectrum)
    
    f0[i+1,] =  sim_spectrum
  }
  f0[1,] =  background_train
  return(list(f0 = f0, background_train = background_train))
}

dist_change_points =  function(Shat,S0)
{
  
  if(length(Shat)==0)
  {return(Inf)}
  
  if(length(S0)==0)
  {return(-Inf)}
  
  temp =rep(0,length(S0))
  for(j in 1:length(S0))
  {
    temp[j] = min(abs(S0[j] - Shat))
  }
  return( max(temp) )
}

correct = function(S,gam)
{
      if(length(S)==0)
      {return(NULL)}
  
       S_c = c()
       for(j in 1:length(S))
       {
           if(j ==1)
           {
             S_c = c(S[1])
           }
           if( j >1)
           {
             
              if(abs(S[j] -  S[j-1]) > gam )
              {
                S_c = c(S_c ,S[j])
              }
           }
       }
       return(S_c)
         
}

########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################

new_BS = function(y, gam,s,e,flag,S,Dval,pos, N)
{
  if(e-s <   1+  3*gam  ||  flag ==1)
  {
    S =  c(S,NULL)
    Dval =  c(Dval, NULL)
    pos  = c(pos,NULL)
    return(S);
  }
  else{
    e = e- gam
    s = s + gam 
    ###  calculate  statistics
    a =  rep(0,e-s+1  )
    
    for(t  in    (s+1):(e-1)  )
    {
      
      a[t-s ] =  Delta_se_t(y,s,e,t,N)
    }
    
    
    best_t  =  which.max( a)
    
    #  if(a[best_t]  <tau  )
    #  {
    #3   return(S)
    #  }
    #   print(s+best_t)
    best_t =   s+ best_t 
    pos1 =  pos
    pos2 = pos
    pos1[length(pos)] = pos[length(pos)]+1
    pos2[length(pos)] = pos[length(pos)]+1
    temp1 = new_BS(y,gam,s,best_t-1,flag,S,Dval,pos1,N)
    temp2 = new_BS(y, gam,best_t+1,e,flag, S,Dval,pos2,N)
    S1 = temp1$S
    Dval1 = temp1$Dval    
    pos1 = temp1$pos 
    S2 = temp2$S
    Dval2 = temp2$Dval 
    pos2 = temp2$pos 
    
    #pos = c(pos,pos[length(pos)]+1) 
    S =  c(S,best_t)
    Dval = c(Dval,a[best_t-s])
    
    
    S  =   c(S,S1,S2)
    Dval =  c(Dval,Dval1,Dval2)
    pos = c(pos,pos1,pos2)
    return(list(S=S,Dval = Dval,pos=pos))
  }
}
###########################################

BS_path = function(y, gam,s,e,flag,S,Dval,pos, N)
{
  if(e-s <   1+  2*gam  ||  flag ==1)#3*gam 
  {
    S =  c(S,NULL)
    Dval =  c(Dval, NULL)
    pos  = c(pos,NULL)
    return(S);
  }
  else{
    #e = e- gam
    #s = s + gam 
    ###  calculate  statistics
    a =  rep(0,e-s+1  )
    
    for(t  in   (s+gam):(e-gam) )#(s+gam):(e-gam)  )
    {
      
      a[t-s ] =  Delta_se_t(y,s,e,t,N)# Delta_se_t(y,s+gam,e-gam,t,N)
        #Delta_se_t(y,s+gam,e-gam,t,N)
    }
    
    
    best_t  =  which.max( a)
    
    #  if(a[best_t]  <tau  )
    #  {
    #3   return(S)
    #  }
    #   print(s+best_t)
    best_t =   s+ best_t 
    pos1 =  pos
    pos2 = pos
    pos1[length(pos)] = pos[length(pos)]+1
    pos2[length(pos)] = pos[length(pos)]+1
    temp1 = BS_path(y,gam,s,best_t,flag,S,Dval,pos1,N)#0-1
    temp2 = BS_path(y, gam,best_t,e,flag, S,Dval,pos2,N)#+1
    S1 = temp1$S
    Dval1 = temp1$Dval    
    pos1 = temp1$pos 
    S2 = temp2$S
    Dval2 = temp2$Dval 
    pos2 = temp2$pos 
    
    #pos = c(pos,pos[length(pos)]+1) 
    S =  c(S,best_t)
    Dval = c(Dval,a[best_t-s])
    
    
    S  =   c(S,S1,S2)
    Dval =  c(Dval,Dval1,Dval2)
    pos = c(pos,pos1,pos2)
    return(list(S=S,Dval = Dval,pos=pos))
  }
}


##########################################
new_BS_threshold =  function(temp,tau,p)
{
  ind = which(temp$Dval >  tau)
  
  Shat =  c()
  
  if(length(ind)==0)
  {
    return(NULL)
  }
  
  for( j in 1:length(ind))
  {
    if(p[ind[j]]==0)
    {
      Shat = c(Shat,temp$S[ind[j]])
    }
    if(p[ind[j]] > 0  && min(abs(Shat - temp$S[p[ind[j]]] ))==0   )
    {
      Shat = c(Shat,temp$S[ind[j]])
    }
  }
  
  return(Shat)
}
parent =  function(temp)
{
  p= rep(1,length(temp$pos))
  p[1] = 0
  for(i  in 2:length(temp$pos))
  {
    ind = which(temp$pos[1:(i-1)] == (temp$pos[i]-1))
    p[i] =  ind[length(ind)]
  }
  return(p)
}

#S = Model_selection_bs_V3      (y_new,grid,Ntau =15,gam=1,rep(2,T/2))
Model_selection_bs_V3 =  function(y,grid,Ntau =15,gam,N)
{
  n = max(N)
  T = dim(y)[1]
  y1 =  matrix(NA, T,floor(n/2) )
  y2 = matrix(NA, T, floor(n/2))
  
  for(t in 1:T)
  {
    y1[t,1:floor(N[t]/2)] =  y[t,1:floor(N[t]/2)]
    
    y2[t,1:floor(N[t]/2)] =  y[t,(1+floor(N[t]/2)):(2* floor(N[t]/2) ) ]
  }
  #grid = 
  temp1 = new_BS(y1, gam,1,T,0,NULL,NULL,1, floor(N/2))  
  Dval = temp1$Dval
  p1 =  parent(temp1)
  aux = sort(Dval,decreasing = TRUE)
  tau_grid = rev(aux[1:100]-10^{-4})
  tau_grid = c(tau_grid,10)
  
  S =  c()
  for( j in 1:length(tau_grid))
  {
    aux = new_BS_threshold(temp1,tau_grid[j],p1)
    #if(length(aux)==0)
    #{
    #  break;
    #}
    S[[j]] = sort(aux)
  }
  
   T= dim(y1)[1]
   S = unique(S)
   score  = rep(0,length(S)+1)
   
   for(j in 1:length(S))
   {
    #  print(S[[j]])
      # print(y2)
    # print(dim(y2))
      score[j] = L_obj(y2,S[[j]],sort(y2),floor(N/2))
        #L_obj_V2(y2,y1,S[[j]],sort(y1),floor(N/2))
        #L_obj(y2,S[[j]],sort(y2),floor(N/2))
        #L_obj_V2(y2,y1,S[[j]],grid,N)
        #L_obj_v3(y1,y2,S[[j]])
   }
   score[j+1] = L_obj(y2,S[[j]],sort(y2),floor(N/2))
     #L_obj_V2(y2,y1,NULL,sort(y1),floor(N/2))
     #L_obj_V2(y2,y1,NULL,grid,N)
     #L_obj_v3(y1,y2,NULL)
   
   best_ind = which.min(score)
   
   if(best_ind==length(score))
   {
     return(NULL)
   }
   return( S[[best_ind]])
  # 
  # T= dim(y1)[1]
  # S = unique(S)
  # for(j in rev(1:length(S)))
  # {
  #   Shat = S[[j]]
  #   
  #   if(j==length(S))
  #   {
  #     dif =  S[[j]]
  #   }
  #   if(j < length(S))
  #   {  dif =  setdiff(S[[j]],S[[j+1]])}
  #   
  #   if(j == 1)
  #   { break;}
  #   
  #   ns =  floor(min(diff(c(1,S[[j]],T)))/2)
  #   ns = min(20,ns)
  #   st =rep(0,length(dif))
  #   
  #   for(s in 1:length(dif))
  #   {
  #     low =  dif[s]-ns
  #     up = dif[s] +ns
  #     t =  dif[s]
  #     
  #     aux =  as.vector(y2[low:(t-floor(ns/2)),])
  #     aux =  aux[which(is.na(aux)==FALSE)]
  #     Fhat = ecdf(aux)
  #     Fhat = Fhat(grid)
  #     
  #     aux =  as.vector(y2[(t+floor(ns/2)):up,])
  #     aux =  aux[which(is.na(aux)==FALSE)]
  #     Fhat2 = ecdf(aux)
  #     Fhat2 = Fhat2(grid)
  #     
  #     st[s] =  max(abs(Fhat-Fhat2)) * sqrt(ns/2)
  #   }
  #   
  #   if(max(st)<1.51)
  #     break;
  # }
  # Shat = NULL
  # if( j < length(S))
  # {
  #   Shat = S[[j+1]]
  # }
  # return(Shat)
}

#######################################3

new_WBS =  function(y, gam,s,e,flag,S,Dval,pos,alpha,beta,N)
{
  
  
  alpha_new =  pmax(alpha,s)
  beta_new = pmin(beta,e)
  # 
  ind = which( beta_new- alpha_new > 1+ 3 )#3*gam)
  alpha_new =  alpha_new[ind] + gam
  beta_new=   beta_new[ind] - gam
  M =  length(alpha_new)
  
  xi = 0#1/8
  alpha_new2 =alpha_new
  beta_new2  =  beta_new
  alpha_new =   ceiling((1-xi)*alpha_new2+  xi*beta_new2)
  beta_new =    ceiling((1-xi)*beta_new2 +  xi*alpha_new2)
  ind = which( beta_new- alpha_new >1)
  alpha_new =  alpha_new[ind]
  beta_new=   beta_new[ind]
  M =  length(alpha_new)
  
  # print(S)
  temp2 = cbind(alpha_new,beta_new)
  temp2 =unique(temp2)
  M  =  dim(temp2)[1]
  if(M>0)
  {
    alpha_new =  temp2[,1]
    beta_new =  temp2[,2]
  }
  
  #  beta_new=   beta_new[ind]
  # # 
#  ind = whic
  if(M ==  0 ||  pos[length(pos)]>12)
  {
    S =  c(S,NULL)
    Dval =  c(Dval, NULL)
    pos  = NULL#c(pos,NULL)
    return(list(S=S,pos = pos,Dval = Dval))
  }
  
  b  =  rep(0,M)
  a =  rep(0,M)
  
  #print(beta_new)
  #print(alpha_new)
  #  print()
  for( m in 1:M  )
  {
    temp  =  rep(0,beta_new[m]-alpha_new[m]+1)
    for(t  in    (alpha_new[m]+1):(beta_new[m]-1 )  )
    {
      temp[t-(alpha_new[m]) ] =  Delta_se_t(y,alpha_new[m],beta_new[m],t,N)
    }
    best_ind  =  which.max(temp)
    a[m] =  alpha_new[m] +  best_ind
    b[m] =   temp[best_ind]
  }
  best_ind =  which.max(b)
  
  
  #  print(b)
#  if(b[best_ind]  <tau  )
 # {
  #  return(S)
  #}
  # print(S)
  #best_t =   s+ best_t 
  pos1 =  pos
  pos2 = pos
  pos1[length(pos)] = pos[length(pos)]+1
  pos2[length(pos)] = pos[length(pos)]+1
  
  #print(c(a[best_ind],b[best_ind]))
  temp2 = new_WBS(y,gam,a[best_ind]+1,e,flag, S,Dval,pos2,alpha, beta,N)
  temp1 = new_WBS(y,gam,s,a[best_ind]-1,flag, S,Dval,pos1,alpha, beta,N)
  S1 = temp1$S 
  Dval1 = temp1$Dval     
  pos1 = temp1$pos  
  S2 = temp2$S 
  Dval2 = temp2$Dval 
  pos2 = temp2$pos 
  
  S =  c(S,a[best_ind])
  Dval = c(Dval,b[best_ind])
  
  
  S  =   c(S,S1,S2)
  Dval =  c(Dval,Dval1,Dval2)
  pos = c(pos,pos1,pos2)
  
  return(list(S=S,Dval = Dval,pos=pos))
  # S  =   c(S,S1,S2)
  # 
  # 
  # return(S)
  
}

#######################


Model_selection_V3 =  function(y,grid,Ntau =15,gam,N,alpha,beta)
{
  n = max(N)
  T = dim(y)[1]
  y1 =  matrix(NA, T,floor(n/2) )
  y2 = matrix(NA, T, floor(n/2))
  
  for(t in 1:T)
  {
    y1[t,1:floor(N[t]/2)] =  y[t,1:floor(N[t]/2)]
    
    y2[t,1:floor(N[t]/2)] =  y[t,(1+floor(N[t]/2)):(2* floor(N[t]/2) ) ]
  }
  #grid = 
  temp1 = new_NWBS(y1, gam,1,T,0,NULL,NULL,1,alpha,beta,rep(1,T))
    #new_WBS(y1, gam=1,1,T,0,NULL,NULL,1, N)  
  Dval = temp1$Dval
  p1 =  parent(temp1)
  aux = sort(Dval,decreasing = TRUE)
  tau_grid = rev(aux[1:50]-10^{-4})
  tau_grid = c(tau_grid,10)
  
  S =  c()
  for( j in 1:length(tau_grid))
  {
    aux = new_BS_threshold(temp1,tau_grid[j],p1)
    #if(length(aux)==0)
    #{
    #  break;
    #}
    S[[j]] = sort(aux)
  }
  
  T= dim(y1)[1]
  S = unique(S)
  score  = rep(0,length(S)+1)
  
  for(j in 1:length(S))
  {
    #  print(S[[j]])
    # print(y2)
    # print(dim(y2))
    score[j] = L_obj_v3(y1,y2,S[[j]])
  }
  score[j+1] = L_obj_v3(y1,y2,NULL)
  
  best_ind = which.min(score)
  
  if(best_ind==length(score))
  {
    return(NULL)
  }
  return( S[[best_ind]])
  
  # T= dim(y1)[1]
  # S = unique(S)
  # for(j in rev(1:length(S)))
  # {
  #   Shat = S[[j]]
  #   
  #   if(j==length(S))
  #   {
  #     dif =  S[[j]]
  #   }
  #   if(j < length(S))
  #   {  dif =  setdiff(S[[j]],S[[j+1]])}
  #   
  #   if(j == 1)
  #   { break;}
  #   
  #   ns =  floor(min(diff(c(1,S[[j]],T)))/2)
  #   st =rep(0,length(dif))
  #   ns  = min(20,ns)
  #   for(s in 1:length(dif))
  #   {
  #     low =  dif[s]-ns
  #     up = dif[s] +ns
  #     t =  dif[s]
  #     
  #     aux =  as.vector(y2[low:(t-1),])#floor(ns/2)
  #     aux =  aux[which(is.na(aux)==FALSE)]
  #     Fhat = ecdf(aux)
  #     Fhat = Fhat(grid)
  #     
  #     aux =  as.vector(y2[(t+1):up,])
  #     aux =  aux[which(is.na(aux)==FALSE)]
  #     Fhat2 = ecdf(aux)
  #     Fhat2 = Fhat2(grid)
  #     
  #     st[s] =  max(abs(Fhat-Fhat2)) * sqrt(ns/2)
  #   }
  #   
  #   if(max(st)<1.22)
  #     break;
  # }
  # Shat = NULL
  # if( j < length(S))
  # {
  #   Shat = S[[j+1]]
  # }
  # return(Shat)
}
L_obj_v3 =  function(y1,y2,S)
{
   K =  length(S)
  # # 
   T =  dim(y1)[1]
   val = rep(0,length(grid)) 
   cost = 0
   sorted_y2 = sort(y2)
  # 
   for(k in 0:(K))
   {
     if(k== 0)
     {
       s =  1
       if(K>0)
       {
         e =  S[1]
       }
       if(K==0)
       {
         e = T
       }
     }
     ########
     ##3
     if(k   == K && k > 0)
     {
       s =  S[K]+1
       e =  T
     }
     ##############3
     if( 0 < k &&  k < K)
     {
       s= S[k]+1
       e= S[k+1]
     }
     ######################33
     
     aux =  as.vector(y2[s:e,])
     aux =  aux[which(is.na(aux)==FALSE)]
     
     Fhat = ecdf(  aux )
     F_hat_z = Fhat(sorted_y2)
     
     ny =  length(sorted_y2)
     a = 1:ny
     a = a*(ny-a)
     
     ind =  which(F_hat_z==0)
     if(length(ind)>0)
     {
       F_hat_z = F_hat_z[-ind]
       a = a[-ind] 
     }
     ind =  which(F_hat_z==1)
     if(length(ind)>0)
     {
       F_hat_z = F_hat_z[-ind]
       a = a[-ind] 
     }
     cost = cost +   ny*(e-s+1)*sum( (F_hat_z*log(F_hat_z)+(1-F_hat_z)*log(1-F_hat_z))/a     ) 
   
     
   }
   #0.5*(1/log(10))
  cost = -cost +  0.25*length(S)*(log(ny))^{2.1}
  return(cost   )
  
  
  
  # K =  length(S)
  # # 
  # T =  dim(y1)[1]
  # val = rep(0,length(grid)) 
  # cost = 0
  # sorted_y2 = sort(y2)
  # 
  # u =  c()
  # for(k in 0:(K))
  # {
  #   if(k== 0)
  #   {
  #     s =  1
  #     if(K>0)
  #     {
  #       e =  S[1]
  #     }
  #     if(K==0)
  #     {
  #       e = T
  #     }
  #   }
  #   ########
  #   ##3
  #   if(k   == K && k > 0)
  #   {
  #     s =  S[K]+1
  #     e =  T
  #   }
  #   ##############3
  #   if( 0 < k &&  k < K)
  #   {
  #     s= S[k]+1
  #     e= S[k+1]
  #   }
  #   ######################33
  #   
  #   aux =  as.vector(y2[s:e,])
  #   aux =  aux[which(is.na(aux)==FALSE)]
  #   Fhat = ecdf(  aux )
  #   F_hat_z =  Fhat(grid)
  # 
  #   
  # }
}



Model_selection_bs_V4 =  function(y,grid,Ntau =15,gam,N)
{
  n = max(N)
  T = dim(y)[1]

  temp1 = BS_path(y, gam,1,T,0,NULL,NULL,1, N)  
    #new_BS(y, gam,1,T,0,NULL,NULL,1, N)  
  Dval = temp1$Dval
  p1 =  parent(temp1)
  aux = sort(Dval,decreasing = TRUE)
  tau_grid = rev(aux[1:100]-10^{-4})
  tau_grid = c(tau_grid,10)
  
  S =  c()
  for( j in 1:length(tau_grid))
  {
    aux = new_BS_threshold(temp1,tau_grid[j],p1)
    #if(length(aux)==0)
    #{
    #  break;
    #}
    S[[j]] = sort(aux)
  }
  
  T= dim(y)[1]
  S = unique(S)
  score  = rep(0,length(S)+1)
  
  for(j in 1:length(S))
  {
    #  print(S[[j]])
    # print(y2)
    # print(dim(y2))
    score[j] = L_obj_v3(as.matrix(y[,1]),as.matrix(y[,1]),S[[j]])
    #L_obj_V2(y2,y1,S[[j]],grid,N)
    #L_obj_v3(y1,y2,S[[j]])
  }
  score[j+1] = L_obj_v3(as.matrix(y[,1]),as.matrix(y[,1]),NULL)
  #L_obj_V2(y2,y1,NULL,grid,N)
  #L_obj_v3(y1,y2,NULL)
  
  best_ind = which.min(score)
  
  if(best_ind==length(score))
  {
    return(NULL)
  }
  return( S[[best_ind]])
  # 

}
#############################################3


Model_selection_V4 =  function(y,grid,Ntau =15,gam,N,alpha,beta)
{
  n = max(N)
  T = dim(y)[1]
  #grid = 
  temp1 = new_WBS(y, gam,1,T,0,NULL,NULL,1,alpha,beta,rep(1,T))
  #new_WBS(y1, gam=1,1,T,0,NULL,NULL,1, N)  
  Dval = temp1$Dval
  p1 =  parent(temp1)
  aux = sort(Dval,decreasing = TRUE)
  tau_grid = rev(aux[1:50]-10^{-4})
  tau_grid =  tau_grid[which(is.na(tau_grid)==FALSE)] ### *
  tau_grid = c(tau_grid,10) 
  
  S =  c() 
  for( j in 1:length(tau_grid))
  {
    aux = new_BS_threshold(temp1,tau_grid[j],p1)

    if(length(aux)==0)  ##*
      break;###*
    
    S[[j]] = sort(aux)
  }
  
  T= dim(y)[1]
  S = unique(S)
  score  = rep(0,length(S)+1)
  
  for(j in 1:length(S))
  {

    score[j] = L_obj_v3(as.matrix(y[,1]),as.matrix(y[,1]),S[[j]])
  }
  score[j+1] = L_obj_v3(as.matrix(y[,1]),as.matrix(y[,1]),NULL)
  
  best_ind = which.min(score)
  
  if(best_ind==length(score))
  {
    return(NULL)
  }
  return( S[[best_ind]])
}


#############################################################3
################################################################


Model_selection_bs_V5 =  function(y,z,gam,N)
{
  n = max(N)
  T = dim(y)[1]
  # 
  # if(length(z[,1])<1000 && dim(z)[2] == 1  )
  # {
  #   ys =matrix(0,c(2*T,1 ))
  #   ys[2*(1:T)-1,1] = y[,1]
  #   ys[2*(1:T),1] = z[,1]   
  #   temp1 = new_BS(ys,gam= 10,1,2*T,0,NULL,NULL,1, rep(1,2*T))  
  #   Dval = temp1$Dval
  #   p1 =  parent(temp1)
  #   aux = sort(Dval,decreasing = TRUE)
  #   tau_grid = rev(aux[1:min(60,length(Dval))]-10^{-5})
  #   tau_grid = c(tau_grid,10)
  #   
  #   S =  c()
  #   for( j in 1:length(tau_grid))
  #   {
  #     aux = new_BS_threshold(temp1,tau_grid[j],p1)
  #     #if(length(aux)==0)
  #     #{
  #     #  break;
  #     #}
  #     if(length(aux)==0)
  #       break;
  #     
  #     S[[j]] = floor(sort(aux)/2)
  #   }
  # }
  #   
 # if(length(z[,1])>=1000 || dim(z)[2] > 1  )
  #{
  temp1 = new_BS(z, gam,1,T,0,NULL,NULL,1, N)  
  Dval = temp1$Dval
  p1 =  parent(temp1)
  aux = sort(Dval,decreasing = TRUE)
  tau_grid = rev(aux[1:min(20,length(Dval))]-10^{-5})
  tau_grid = c(tau_grid,10)
  
  S =  c()
  for( j in 1:length(tau_grid))
  {
    aux = new_BS_threshold(temp1,tau_grid[j],p1)
    if(length(aux)==0)
    {
      break;
    }
    S[[j]] = sort(aux)
  }
  #}
  T= dim(y)[1]
  S = unique(S)
  if(length(S)==0)
  {
    return(NULL)
  }  
  lamb =log(sum(N))/1.5#2.5#1.5#2#2.555#
  for(j in 1:length(S))#)
  {
    if(length(S[[j]])==0)
    {
      j = j+1;
      
      if(j>length(S))
        break;
    }
    
      B2  =  S[[j]]
      if(j==length(S))
      {
        B1 = NULL
      }
       if(j< length(S))
       {
         B1 = S[[j+1]]
       }
      temp = setdiff(B2,B1)

      st =  -10^15
      #Delta_se_t(z,eta1,eta2,eta,N)^2 
      for(l in 1:length(temp))
      {
         eta =  temp[l]
         
         if(j == length(S))
         {
           eta1 = 1
           eta2 = T
         }
         if(j < length(S))
         {
           for(k in 1:length(S[[j+1]]))
           {
             if(S[[j+1]][k]> eta  )
               break;
           }
           if(S[[j+1]][k]> eta )
           {
             eta2 = S[[j+1]][k]
             
             if(k ==1)
               eta1 = 1
             
             if(k > 1)
               eta1 = S[[j+1]][k-1]+1
           }
           if(S[[j+1]][k]< eta )
           {
             eta1 = S[[j+1]][k]+1
             eta2 = T
           }
         }
         st_aux = Delta_se_t(y,eta1,eta2,eta,N)^2
          # print(st_au)
         if(st_aux> st)
         {
           st = st_aux
         }
      }###  close for defining  eta1 and eta2

    
      # print(c1 - c2 - Delta_se_t(z,eta1+1,eta2,eta,N)^2 + lamb)
      if(st >   lamb)
      {
        #B1 = B2
        return(B2)
      }
     # print(st)
  }
  #c1 - c2 - Delta_se_t(z,eta1+1,eta2,eta,N)^2 + lamb
  
  return(B1)
  # 
  
}
#############################################3

arg_max_Delta_se_t = function(y,s,e,t,N)
{
  #T =   dim(y)[2]
  n =  dim(y)[2]
  
  n_st = sum(N[s:t])  #n*(t-s+1)
  n_se = sum(N[s:e])  #n*(e-s+1)
  n_te =sum(N[(t+1):e]) #n*(e-(t+1) +1)
  
  aux =  as.vector(y[s:t,])
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  vec_y =  as.vector(y[s:e,])
  vec_y = vec_y[which(is.na(vec_y)==FALSE)]
  Fhat_st =  temp(vec_y)# temp(grid)
  
  aux = y[(t+1):e,]
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  Fhat_te =  temp(vec_y)# temp(grid)
  
  #temp =  sqrt( n_st*n_te / n_se   ) *max(abs(Fhat_te - Fhat_st  ))
  ind =  which.max(abs(Fhat_te - Fhat_st  ))
  z_hat = min(vec_y[ ind])
  #print(min(abs(Fhat_te - Fhat_st  )))
  return(z_hat)
}
################################################3



Model_selection_V5 =  function(y,z,gam,N,alpha,beta)
{
  
  n = max(N)
  T = dim(y)[1]
  #grid = 
  temp1 = new_WBS(z, gam,1,T,0,NULL,NULL,1,alpha,beta,N)
  #new_WBS(y1, gam=1,1,T,0,NULL,NULL,1, N)  
  Dval = temp1$Dval
  p1 =  parent(temp1)
  aux = sort(Dval,decreasing = TRUE)
  tau_grid = rev(aux[1:min(20,length(Dval))]-10^{-4})
  tau_grid = c(tau_grid,10)
  
  S =  c()
  for( j in 1:length(tau_grid))
  {
    aux = new_BS_threshold(temp1,tau_grid[j],p1)
    
    S[[j]] = sort(aux)
  }
  
  T= dim(y)[1]
  S = unique(S)
  if(length(S)==0)
  {
    return(NULL)
  }  

  lamb = log(sum(N))/1.5#2.5#1.25#1.5#1.45#1.4#*0.2
  for(j in 1:(length(S)))
  {
    if(length(S[[j]]) == 0)
    {
      j = 2
    }
    B2  =  S[[j]]
    if(j < length(S))
    {
      B1 = S[[j+1]]
    }
     if(j  == length(S))
     {
       B1 = NULL
     }
    temp = setdiff(B2,B1)
    
    st =  -10^15
    #Delta_se_t(z,eta1,eta2,eta,N)^2 
    count = 0
    for(l in 1:length(temp))
    {
      eta =  temp[l]
      
      if( length(B1)==0)
      {
        eta1 = 1
        eta2 = T
      }
      if( length(B1)>0)
      {
        for(k in 1:length(B1))
        {
          if(B1[k]> eta  )
            break;
        }
        if(B1[k]> eta )
        {
          eta2 = B1[k]
          
          if(k ==1)
            eta1 = 1
          
          if(k > 1)
            eta1 = B1[k-1]+1
        }
        if(B1[k]< eta )
        {
          eta1 = B1[k]+1
          eta2 = T
        }
      }
      if( length(temp) > 1     )
      {
        st_aux = max(wbs_Delta_se_t(y,eta1,eta2,eta,N,alpha,beta)^2)  
        if(st_aux == 0)
        {
          st_aux = Delta_se_t(y,eta1,eta2,eta,N)^2
        }
      }
      if(length(temp)==1){st_aux = max(Delta_se_t(y,eta1,eta2,eta,N)^2)  }
                   #,Delta_se_t(z,eta1,eta2,eta,N)^2,Delta_se_t(u,eta1,eta2,eta,N)^2)
    #  print(st_aux)
      if(st_aux> st)
      {
        st = st_aux
    #    B1 = B2
      }
    }  

    if(st >  lamb)#  st > lamb
    {
    # B1 = B2
      return(B2)
    }
#     print(st)
  }#
  #c1 - c2 - Delta_se_t(z,eta1+1,eta2,eta,N)^2 + lamb
  return(B1)
  # 
  
}

wbs_Delta_se_t  =  function(y,s,e,t,N,alpha, beta)
{
  alpha_new =  pmax(alpha,s)
  beta_new = pmin(beta,e)
  # 
  #  ind = which( beta_new- alpha_new > 1+  3*gam)
  #  alpha_new =  alpha_new[ind]
  #  beta_new=   beta_new[ind]
  # # 
  ind = which( beta_new- alpha_new >20)
  alpha_new =  alpha_new[ind]
  beta_new=   beta_new[ind]
  M =  length(alpha_new)
  
  xi = 1/8
  alpha_new2 =alpha_new
  beta_new2  =  beta_new
  alpha_new =   ceiling((1-xi)*alpha_new2+  xi*beta_new2)
  beta_new =    ceiling((1-xi)*beta_new2 +  xi*alpha_new2)
  ind = which( beta_new- alpha_new >1)
  alpha_new =  alpha_new[ind]
  beta_new=   beta_new[ind]
  M =  length(alpha_new)
  
  # print(S)
  
  if(M ==  0)
  {
    return(0) 
  }
  
  b  =  rep(0,M)
  a =  rep(0,M)
  
  #print(beta_new)
  #print(alpha_new)
  #  print()
  
  for( m in 1:M  )
  {
  #  temp  =  rep(0,beta_new[m]-alpha_new[m]+1)
  #  for(t  in    (alpha_new[m]+1):(beta_new[m]-1 )  )
  #  {
      if(alpha_new[m]<t &&  t  < beta_new[m]  )
      {
        b[m] =  Delta_se_t(y,alpha_new[m],beta_new[m],t,N)
      }
  #  }
   # best_ind  =  which.max(temp)
  #  a[m] =  alpha_new[m] +  best_ind
  #  b[m] =   temp[best_ind]
  }
  #best_ind =  which.max(b)
  return(max(b))
  #return(b[best_ind])
}

####################################3
##############################################3
#########################



NWBS_estimator =  function(y,grid,Ntau =15,gam,N,alpha,beta)
{
  
  n = max(N)
  T = dim(y)[1]
  #grid = 
  #temp1 = BS_path(y, gam,1,T,0,NULL,NULL,1, N)  
   temp1 = new_WBS(y, gam,1,T,0,NULL,NULL,1,alpha,beta,N)
  #temp1 =  new_WBS(y[, (1+n/2):n], gam,1,T,0,NULL,NULL,1,alpha,beta,N)
  #n = max(N/2)
  #new_WBS(y1, gam=1,1,T,0,NULL,NULL,1, N)  
  Dval = temp1$Dval
  p1 =  parent(temp1)
  aux = sort(Dval,decreasing = TRUE)
  tau_grid = rev(aux[1:50]-10^{-4})
  tau_grid =  tau_grid[which(is.na(tau_grid)==FALSE)] ### *
  tau_grid = c(tau_grid,10) 
  
  S =  c() 
  for( j in 1:length(tau_grid))
  {
    aux = new_BS_threshold(temp1,tau_grid[j],p1)
    
    if(length(aux)==0)  ##*
      break;###*
    
    S[[j]] = sort(aux)
  }
  
  T= dim(y)[1]
  S = unique(S)
  score  = rep(0,length(S)+1)
  #n = dim(y)[2]
  
  for(i in 1:n)
  {
    for(j in 1:length(S))
    {
      
      score[j] = score[j]+ L_obj_v4(as.matrix(y[,i]),as.matrix(y[,i]),S[[j]])
    }
    score[j+1] = score[j+1]+ L_obj_v4(as.matrix(y[,i]),as.matrix(y[,i]),NULL)
    
  }
    
    
  best_ind = which.min(score)
  
  if(best_ind==length(score))
  {
    return(NULL)
  }
  return( S[[best_ind]])
}

L_obj_v4 = function(y1,y2,S)
{
  K =  length(S)
  # # 
  T =  dim(y1)[1]
  val = rep(0,length(grid)) 
  cost = 0
  sorted_y2 = sort(y2)
  # 
  for(k in 0:(K))
  {
    if(k== 0)
    {
      s =  1
      if(K>0)
      {
        e =  S[1]
      }
      if(K==0)
      {
        e = T
      }
    }
    ########
    ##3
    if(k   == K && k > 0)
    {
      s =  S[K]+1
      e =  T
    }
    ##############3
    if( 0 < k &&  k < K)
    {
      s= S[k]+1
      e= S[k+1]
    }
    ######################33
    
    aux =  as.vector(y2[s:e,])
    aux =  aux[which(is.na(aux)==FALSE)]
    
    Fhat = ecdf(  aux )
    F_hat_z = Fhat(sorted_y2)
    
    ny =  length(sorted_y2)
    a = 1:ny
    a = a*(ny-a)
    
    ind =  which(F_hat_z==0)
    if(length(ind)>0)
    {
      F_hat_z = F_hat_z[-ind]
      a = a[-ind] 
    }
    ind =  which(F_hat_z==1)
    if(length(ind)>0)
    {
      F_hat_z = F_hat_z[-ind]
      a = a[-ind] 
    }
    cost = cost +   ny*(e-s+1)*sum( (F_hat_z*log(F_hat_z)+(1-F_hat_z)*log(1-F_hat_z))/a     ) 
    
    
  }
  #0.5*(1/log(10))
  cost = -cost +  0.2*length(S)*(log(ny))^{2.1}
  return(cost   )
}

##############################################3

new_NWBS =  function(y, gam,s,e,flag,S,Dval,pos,alpha,beta,N)
{
  
  
  alpha_new =  pmax(alpha,s)
  beta_new = pmin(beta,e)
  # 
  ind = which( beta_new- alpha_new > 1+5)# 3*gam)
  #  alpha_new =  alpha_new[ind]h( beta_new- alpha_new >1)
  alpha_new =  alpha_new[ind] #+ gam
  beta_new=   beta_new[ind] #- gam
  M =  length(alpha_new)
  
  xi = 0#1/8
  alpha_new2 =alpha_new
  beta_new2  =  beta_new
  alpha_new =   ceiling((1-xi)*alpha_new2+  xi*beta_new2)
  beta_new =    ceiling((1-xi)*beta_new2 +  xi*alpha_new2)
  ind = which( beta_new- alpha_new >1)
  alpha_new =  alpha_new[ind]
  beta_new=   beta_new[ind]
  M =  length(alpha_new)
  
  # print(S)
  temp2 = cbind(alpha_new,beta_new)
  temp2 =unique(temp2)
  M  =  dim(temp2)[1]
  if(M>0)
  {
    alpha_new =  temp2[,1]
    beta_new =  temp2[,2]
  }
  
  
  #  beta_new=   beta_new[ind]
  # # 
  #  ind = whic
  if(M ==  0 ||  pos[length(pos)]>12)
  {
    S =  c(S,NULL)
    Dval =  c(Dval, NULL)
    pos  = NULL#c(pos,NULL)
    return(list(S=S,pos = pos,Dval = Dval))
  }
  
  b  =  rep(0,M)
  a =  rep(0,M)
  
  #print(beta_new)
  #print(alpha_new)
  #  print()
  for( m in 1:M  )
  {
    temp  =  rep(0,beta_new[m]-alpha_new[m]+1)
    for(t  in    (alpha_new[m]+1):(beta_new[m]-1 )  )
    {
      temp[t-(alpha_new[m]) ] =  Delta_se_t(y,alpha_new[m],beta_new[m],t,N)
    }
    best_ind  =  which.max(temp)
    a[m] =  alpha_new[m] +  best_ind
    b[m] =   temp[best_ind]
  }
  best_ind =  which.max(b)
  
  print(b[best_ind])
  #  print(b)
  #  if(b[best_ind]  <tau  )
  # {
  #  return(S)
  #}
  # print(S)
  #best_t =   s+ best_t 
  pos1 =  pos
  pos2 = pos
  pos1[length(pos)] = pos[length(pos)]+1
  pos2[length(pos)] = pos[length(pos)]+1
  
  #print(c(a[best_ind],b[best_ind]))
  temp2 = new_NWBS(y,gam,a[best_ind]+1,e,flag, S,Dval,pos2,alpha, beta,N)
  temp1 = new_NWBS(y,gam,s,a[best_ind]-1,flag, S,Dval,pos1,alpha, beta,N)
  S1 = temp1$S 
  Dval1 = temp1$Dval     
  pos1 = temp1$pos  
  S2 = temp2$S 
  Dval2 = temp2$Dval 
  pos2 = temp2$pos 
  
  S =  c(S,a[best_ind])
  Dval = c(Dval,b[best_ind])
  
  
  S  =   c(S,S1,S2)
  Dval =  c(Dval,Dval1,Dval2)
  pos = c(pos,pos1,pos2)
  
  return(list(S=S,Dval = Dval,pos=pos))
  # S  =   c(S,S1,S2)
  # 
  # 
  # return(S)
  
}