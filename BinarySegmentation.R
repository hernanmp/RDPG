#Binary Segmentation code 
#this function assume that cusum norm can be computed at every t\in (s,e]
#tau is the threshold
Binary.Segmentation=function(data1,data2,tau){
  start.set=c(1,ncol(data1))
  output.set=start.set
  while(length(start.set)>1){
    s=start.set[1]
    e=start.set[2]
    temp.cucum.value=seq(s,e,1)
    temp.cucum.value=sapply(temp.cucum.value,function(x) data.split.cusum.norm(data1,data2,s,e,x))
    break.point=which.max(temp.cucum.value)
    #print( temp.cucum.value [break.point]  )
    #print(break.point+s)
   # print(start.set)
    if(temp.cucum.value[break.point]>tau){
      start.set=sort(c(start.set,break.point+s-1))
      output.set =sort(c(output.set,break.point+s-1))
    }else{ start.set=start.set[-1]}
  }
  output.set[2:length(output.set)]=2*output.set[2:length(output.set)]
    return(output.set[c(-1,-length(output.set))])
}
  
#This version of BS is for local refinement only
local.Binary.Segmentation=function(data.long,s,e,diff.mat,p){
  lower.coor =gen.lower.coordinate(p)
    temp.cucum.value=seq( s, e,1)
    temp.cucum.value=sapply(temp.cucum.value,function(x) 
      sum(compuate.cusum.vector(data.long[lower.coor,],  s,  e,x)* diff.mat[lower.coor] ))
    
    break.point=which.max(temp.cucum.value)
    
  return(s+break.point-1)
}



