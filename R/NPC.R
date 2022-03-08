
###############################################################################
# AlphaNP:                                                                    #
#                                                                             #
# Estimate attribute profiles using the nonparametric approaches for the "AND #
# gate" or the "OR gate" cognitive diagnostic models.                         #
#                                                                             #
# Input:                                                                      #
# (1) Y: a matrix of binary responses (1=correct, 0=incorrect). Rows          #
#               represent persons and columns represent items.                #
# (2) Q: the Q-matrix of the test. Rows represent items and columns represent #
#        attributes.                                                          #
# (3) gate: "AND", the examinee needs to possess all related attributes       #
#                  in order to answer an item correctly                       #
#           "OR", the examinee needs to possess only one of the related       #
#                 attributes in order to answer an item correctly             #
# (4) method: "Hamming", the plain hamming distance method;                   #
#             "Weighted", the hamming distance weighted by inversed variance  #
#             "Penalized", the hamming distance weighted by inversed variance #
#                          and specified penalizing weights for guess and slip#
# Additional input for the "penalized" method:                                #
# (5) wg = weight assigned to guess                                           #
# (6) ws = weight assigned to slip                                            #
#                                                                             #
# Output:                                                                     #
# (1) alpha.est: estimated ability patterns                                   #
# (2) est.ideal: estimated ideal response to all the items for all the        #
#                 examinees                                                   #
# (3) n.ite: number of ties in the hamming distance. When we encounter ties,  # 
#            the best way is to radomly choose one. However, we want to       #
#           monitor the number of ties.                                       #
# (4) pattern: all attribute profiles in the search space.                    #
# (5) loss.matrix: The values for the loss function. Rows represent candidate #
#                  attribute profiles in the same order with the pattern      #
#                 matrix; Columns represent different examinees.              #
# For now, we report the estimated classes, estimated ideal response patterns #
# (mostly for the Q-refinement algorithm) and number of ties. We may want to  #
# report more later, such as the final RSS.                                   # 
#                                                                             #
###############################################################################
#library(prevalence) ## For CheckInput

NPC <- function(Y, Q, gate=c("AND", "OR", "Mix"), method=c("Hamming", "Weighted", "Penalized"), wg=1, ws=1) {
  
  #####
  # 1 #
  ##### Check dimension consistency and convert data to the right formats 
  
  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  #check <- NULL
  #check <- CheckInput(Y, Q)  
  #if (!is.null(check)) return(warning(check))

  gate <- match.arg(gate)
  method <- match.arg(method)
  
  #####
  # 2 #
  ##### Estimation
  
  nperson <- dim(Y)[1]
  nitem <- dim(Q)[1]
  natt <- dim(Q)[2]
  M <- 2 ^ natt
  
  # all possible alpha vectors
  
  pattern <- diag(natt)
  for (l in 2:natt){
    pattern <- rbind(pattern,t(apply(combn(natt,l),2,function(x){apply(pattern[x,],2,sum)})))
  }
  pattern <- rbind(0,pattern)
  
  # ideal responses for each possible alpha vector
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
  weight=Ideal/matrix(rep(colSums(t(Q)),M),M,nitem,T)
  Ideal.mix=Ideal.conj+(Ideal.dis-Ideal.conj)*weight
  
  if(gate == "AND"){
    Ideal=Ideal.conj} else if(gate == "OR") {
      Ideal=Ideal.dis} else if(gate == "Mix") {
        Ideal=Ideal.mix}
  
  if (method == "Hamming")
  {
    weight <- rep(1, nitem)
    ws <- wg <- 1
  } else if (method == "Weighted") 
  {
    p.bar <- apply(Y, 2, mean)
    weight <- 1 / (p.bar * (1 - p.bar))
    weight[weight > 1 / (0.95 * 0.05)] <- 1 / (0.95 * 0.05)
    ws <- wg <- 1
  } else if (method == "Penalized") 
  {
    p.bar <- apply(Y, 2, mean)
    weight <- 1 / (p.bar * (1 - p.bar))  
    weight[weight > 1 / (0.95 * 0.05)] <- 1 / (0.95 * 0.05)  
    if (ws == wg) warning("Penalzing weights for guess and slip are the same --> equivalent with the \"Weighted\" method.")
  } else 
  {
    return(warning("Method specification not valid.")) 
  }
  
  loss.matrix <- matrix(NA, nrow=M, ncol=nperson)
  est.class <- NULL
  est.pattern <- NULL
  n.tie <- rep(0, nperson)
  
  for (i in 1:nperson)
  {
    Y.matrix <- matrix(rep(Y[i,], M), M, nitem, byrow=TRUE)
    loss <- apply(matrix(rep(weight, M), M, nitem, byrow=TRUE) * (wg * abs(Y.matrix - Ideal) * Y.matrix + ws * abs(Y.matrix - Ideal) * (1 - Y.matrix)), 1, sum)
    loss.matrix[,i] <- loss
    
    min.loss <- which(loss == min(loss))

    if (length(min.loss) != 1) 
    {
      n.tie[i] <- length(min.loss)
      min.loss <- sample(min.loss, 1, prob=rep(1 / length(min.loss), length(min.loss))) 
    }
    
    est.class <- c(est.class, min.loss)
  }
  
  est.pattern <- pattern[est.class,]
  est.ideal <- Ideal[est.class,]
  output <- list(alpha.est=est.pattern, est.ideal=est.ideal, est.class=est.class, n.tie=n.tie, pattern=pattern, loss.matrix=loss.matrix, method=method, Q=Q, Y=Y)
  class(output) <- "AlphaNP"
  return(output)
}

