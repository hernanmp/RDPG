# generate SBM given connetivity matrix and candidate vector
# matrix is the connectivity matrix, 
# can.vec is a p dimensional candidata matrix
generate.SBM.mean=function(matrix,can.vec,B,p){
  result= matrix(0,nrow=p,ncol=p)
  aa=p/B
  for ( i in 1: B){
    for ( j in 1: B){
      can.temp.1=rep(0,p)
      can.temp.1[can.vec[(1+(i-1)*aa):(i*aa)]]=1
      can.temp.2=rep(0,p)
      can.temp.2[can.vec[(1+(j-1)*aa):(j*aa)]]=1
      result=result+matrix[i,j]* can.temp.1%*% t(can.temp.2)
      
    }
    
  }
  
  return(result)
}

# generate unstationary network
generate.sec.net=function( mat.mean,n,p){
  mat.obs=matrix(0,nrow=p^2,ncol=n)
  mat.obs= t(sapply(mat.mean, function(x) rbinom(n, 1, x)))
   
  return(mat.obs)
}

#compuate cusum vector given data and (s,e], t
compuate.cusum.vector=function(the.data,s,e,t){
  result.vec=rep(0,nrow(the.data))
  if( t-s<3 | e-t<2){
    result.vec=rep(0,nrow(the.data))
  }else{
    result.vec= ( (e-t)/ ((e-s)*(t-s))  )^(0.5) *rowSums(the.data[,(s+1):t]) -
      ( (t-s)/ ((e-s)*(e-t) )  )^(0.5) *rowSums(the.data[,(t+1):e])
    
  }
  return(result.vec)
}

#compute the cusum norm using two data set
data.split.cusum.norm=function(data1,data2,s,e,t){
  return(sum(compuate.cusum.vector(data1,s,e,t)*compuate.cusum.vector(data2,s,e,t)))
}


# compute hausdorff distance
hausdorff.distance=function(v1,v2,n){
  p1=length(v1)
  p2=length(v2)
  if ( p1*p2==0){ dis=n} else{
  distance.mat=matrix(0,nrow=p1,ncol=p2)
  for( i in 1: p1){
    for( j in 1: p2){
      distance.mat[i,j]=abs(v1[i]-v2[j])
      
    }
  }
   dis=max(max( apply(distance.mat, 1, min)),max( apply(distance.mat, 2, min)))
  }
   return(dis)
  
}

#generate coordiantes of lower triangular matrix
gen.lower.coordinate=function(p){
  mat=matrix(1:p^2, nrow=p)
  
  return(mat[lower.tri(mat,diag = F)])  
}




