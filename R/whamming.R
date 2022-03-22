#=====Estimation of Proficiency Class via weighted Hamming distance=====
whamming=function(Ideal,Y){
  M=nrow(Ideal)
  J=ncol(Ideal)
  N=nrow(Y)
  
  p.bar=apply(Y,2,mean)
  weight=1/(p.bar*(1-p.bar))
  WH.class=NULL
  ntie=0
  H=NULL
  for (i in 1:N)
  {
    ham=apply(matrix(rep(weight, M), M, J, byrow=TRUE)*abs(matrix(rep(Y[i,], M), M, J, byrow=TRUE)-Ideal), 1, sum)
    H=rbind(H,ham)
    min.ham=which(ham==min(ham))
    if (length(min.ham)!=1)
    {
      ntie=ntie+1
      min.ham=sample(min.ham,1,prob=rep(1/length(min.ham),length(min.ham))) ## temporarily fix ties
    }else ntie=ntie
    WH.class=c(WH.class, min.ham)
  }
  return(WH.class)
}
