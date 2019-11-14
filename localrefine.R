
#usvt require symetric matrix A as input also tuning parameter tau
#usvt is good when p is large and rho is large
usvt=function(A, tau,tau1){
  result.mat=matrix(0,nrow=nrow(A),ncol=ncol(A))
  re=eigen(A,symmetric=T)
  kk=length(which(  re$values >tau)) 
  if(kk>0){
  for( i in 1:kk){
    result.mat=result.mat+  re$values[i]  *
       re$vectors[,i] %*% t( re$vectors[,i] )
  }}
  
  kk2=length(which( re$values< (-1)*tau)) 
  if (kk2>0){
  for( i in 1:kk2){
    result.mat=result.mat+  re$values[nrow(A)-i+1]  *
       re$vectors[,nrow(A)-i+1] %*% t( re$vectors[,nrow(A)-i+1] )
  }}
  result.mat[result.mat>tau1]=tau1
  result.mat[result.mat<-1*tau1]=-1*tau1
  return(result.mat)
}


#usvt apply to vector
usvt.vec.norm=function(diff.vec,p,tau){
  diff.mat=matrix(0,nrow=p  ,ncol=p )
  diff.mat[gen.lower.coordinate(p)]=diff.vec
  diff.mat=diff.mat+t(diff.mat)
  #tau=sqrt(p*rho)/2
  diff.mat.temp=usvt(diff.mat,tau,tau1=Inf)
  return(norm(diff.mat.temp,type="F"))
  
}




#local refinement
#data is the cutted vectors (cut the lower triangular) 
#data.long is the uncut version
#tau2 = rho*p*log(n)/2
local.refine=function(data,tau,nbs,p,rho){
  ll=length(nbs)
   result=nbs
    
  for ( ii in 2:ll) { #print(ii)
    if(length(result)>= (ii+1) ){
    record.temp= c(result[ii-1]: result[ii+1]) 
    
    diff.vec= sapply(record.temp, function(x) 
      compuate.cusum.vector(data,result[ii-1],result[ii+1],x))
    result[ii]=result[ii-1] -1+ which.max( 
     sapply(record.temp, function(x) usvt.vec.norm(diff.vec[ ,x-result[ii-1]+1 ],p,tau=tau) ) )
    #print(result[ii])
  result=unique(result)
    
  
 
  }
     }
   
  return(result[-c(1, length(result))]  ) 
  
  
}
  
  
  