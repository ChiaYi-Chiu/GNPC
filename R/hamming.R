#=====Estimation of Proficiency Class via Hamming distance=====
hamming = function(Ideal,Y){
  M=nrow(Ideal)
  J=ncol(Ideal)
  N=nrow(Y)
  
  H.class=NULL
  h.ntie=0
  for (i in 1:N)
  {
    ham=apply(abs(matrix(rep(Y[i,], M), M, J, byrow=TRUE)-Ideal), 1, sum)
    min.ham=which(ham==min(ham))
    if (length(min.ham)!=1)
    {
      h.ntie=h.ntie+1
      min.ham=sample(min.ham,1,prob=rep(1/length(min.ham),length(min.ham))) ## temporarily fix ties
    }else h.ntie=h.ntie
    H.class=c(H.class, min.ham)
  }
  return(H.class)
}
