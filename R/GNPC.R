# The first author of this R code is Yan Sun.
# Function NPC is required for GNPC
GNPC=function(Y,Q,  initial.gate = "AND", initial.dis = "Hamming"){
  N=dim(Y)[1]
  K=dim(Q)[2]
  J=dim(Q)[1]
  M=2^K
  
  pattern <- diag(K)
  for (l in 2:K){
    pattern <- rbind(pattern,t(apply(combn(K,l),2,function(x){apply(pattern[x,],2,sum)})))
  }
  pattern <- rbind(0,pattern)
  
  #============================
  # Conjunctive Ideal Response
  #============================
  Ideal=pattern%*%t(Q) #M*K %*% K*J = M * J
  Ideal.conj=1*(Ideal==(matrix(1,M)%*%t(rowSums(Q))))
  
  #============================
  # Disjunctive Ideal Response
  #============================
  Ideal.dis=1*(Ideal>=1)
  
  #=================================
  # Assigning Initial Weight
  # Initial Weighted Ideal Response
  #=================================
  weight=Ideal/matrix(rep(colSums(t(Q)),M),M,J,T)
  Ideal.mix=Ideal.conj+(Ideal.dis-Ideal.conj)*weight
  
  #===========================================
  # Classify Initial Latent Class
  #===========================================
  all.gate = c("AND", "OR", "Mix")
  all.weight = all.w.ideal = all.loss = all.class =NULL
  for (gg in 1:length(all.gate)){
    est = NPC(Y, Q, gate = all.gate[gg], method = initial.dis)
    fixed.weight = est$est.class
    initial.class=fixed.weight
    d=1
    time=0
    while(d>0.001)
    { #print(time)
      time=time+1
      #===================================================
      # Compute the general weights using the closed form
      #===================================================
      w=NULL
      temp=matrix(NA,M,1)
      Ideal.comb=Ideal.dis+Ideal.conj
      for (j in 1:J) {
        pQ=pattern*matrix(Q[j,]==1,M,K,T)
        upQ=unique(pQ)
        ng=nrow(upQ)
        for (g in 1:ng){
          match <- apply(pQ, 1, identical, upQ[g,])
          m=which(match)
          c= which(initial.class %in% m)
          c=c[!is.na(c)]
          if (length(c)==0){temp[m,]=0.01}
          else if (length(c)!=0){
            temp[m,]=sum(Y[c,j])/sum((1-Ideal.conj[initial.class[c],j])^2)}
        }
        temp[which(Ideal.comb[,j]==0)]=.01
        temp[which(Ideal.comb[,j]==2)]=.99
        w=cbind(w,temp)
      }
      Ideal.w=(1-w)*Ideal.conj+w*Ideal.dis
      w.class = NULL
      loss = 0
      for (i in 1:N)
      {
        ham=rowSums(((matrix(Y[i,], M, J, byrow=TRUE)-Ideal.w))^2)
        loss = loss + min(ham)
        min.ham=which(ham==min(ham))
        if (length(min.ham)!=1){
          min.ham=sample(min.ham,1,prob=rep(1/length(min.ham),length(min.ham))) ## temporarily fix ties
        }
        
        w.class=c(w.class, min.ham)
      }
      d=length(which(w.class-initial.class!=0))/N
      initial.class=w.class
    }
    all.weight[[gg]] = w
    all.w.ideal[[gg]] = Ideal.w
    all.loss = c(all.loss, loss)
    all.class[[gg]] = w.class
  }
  
  if (initial.gate == "Auto"){
    best = which(all.loss == min(all.loss))
    if (length(best)>1) best = sample(best, 1)
  }else best = which(all.gate == initial.gate)
  best
    #best.gate = paste0("Best gate is not available because initial.gate = '",initial.gate, "' is used.")
    
  best.gate = all.gate[best]
  best.weight = all.weight[[best]]
  best.w.ideal = all.w.ideal[[best]]
  best.class = all.class[[best]]
  best.att.est = pattern[unlist(best.class),]
  
  return(list(att.est = unlist(best.att.est), class = unlist(best.class), Ideal.w = unlist(best.w.ideal), 
              weight = unlist(best.weight), best.gate = unlist(best.gate), loss = all.loss))
}
