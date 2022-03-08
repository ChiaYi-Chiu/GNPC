#===Calculate AAR====
AAR=function(x,y){
  aar=1-mean(abs(x-y))
  return(aar)
}
