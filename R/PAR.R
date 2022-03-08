#===Calculate PAR====
PAR=function(x,y){
  out=mean(1*(rowSums(abs(x-y))==0))
  return(out)
}