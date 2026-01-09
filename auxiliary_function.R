################################# Functions #############################################################
# Auxiliary function for Diversity (Solow-Polasky) measure
Div_SP <- function(X, Dis){ #X is the set of datapoints
  M <- as.matrix(exp(-(dist(as.data.frame(X), diag = T, upper = T, method = Dis))))
  #M_1 <- solve(M)
  M_1 <- ginv(M)
  #Div <- 0
  #for (i in 1:nrow(M_1)) {Div <- Div + sum(M_1[,i])/M_1[i,i]}
  Div <- sum(M_1)
  return(Div)
}

# Auxiliary function for diversity difference using Div_SP()
Div_SP_diff1 <- function(X, n){
  a <- Div_SP(X, Dis) - Div_SP(X[-n,], Dis)
  return(a)
}

# Auxiliary function for diversity difference with less computation
Div_SP_diff2 <- function(X, n, Dis){
  M <- as.matrix(exp(-(dist(as.data.frame(X), diag = T, upper = T, method = Dis))))
  M_1 <- solve(M)
  diff <- sum(M_1[,n])/M_1[n,n]
  return(diff)
}

# Auxiliary function for adaptive change of bound
BoundChange <- function(P, Pn, b, r, func){
  PPn <- as.data.frame(rbind(P, Pn))
  val <- ncol(PPn)+1
  PPn[,val] <- apply(PPn, 1, func) 
  P1 <- PPn[order(PPn[,val]),]
  P2 <- P1[1:r,1:(val-1)]
  b <- P1[r,val]
  return(list(Pop=P2, bound=b))
}

# Auxiliary function for crossover
crossover1 <- function(pr, L, U){
  Osp1 <- colMeans(pr)
  Osp2 <- colMeans(rbind(pr, L))
  Osp3 <- colMeans(rbind(pr, U))
  return(rbind(Osp1, Osp2, Osp3))
}

# Auxiliary function for mutation 
mutate <- function(p, L, U){
  mu <- if(sample(c(0, 1), 1)==0){L}else{U}
  dif <- mu - p
  teta <- runif(1, 0, 0.1)
  dir<- sample(c(0, 1), length(p), replace = T)*dif
  p1 <- p + dir*teta
  return(p1)
}

# Auxiliary function for ranking instances based diversity (Solow-Polasky) difference with less computation
Div_SP_diff_rank <- function(X, Dis){
  P_col <- ncol(X)+1
  M <- as.matrix(exp(-(dist(as.data.frame(X), diag = T, upper = T, method = Dis))))
  M_1 <- ginv(M)
  for (k in 1:nrow(X)) {
    X[k,P_col] <- sum(M_1[,k])/M_1[k,k]
  }
  P1 <- X[order(X[,P_col], decreasing = T),]
  return(P1)
}

# Auxiliary function for diversity selection
DivSelect_greedy <- function(P, r, Dis){
  P_col <- ncol(P)
  P_row <- nrow(P)
  P1 <- P
  i = 0
  while (i < (P_row-r)) {
    P1 <- Div_SP_diff_rank(P1, Dis)
    P1 <- P1[-nrow(P1),1:P_col]
    i <- i + 1
  }
  return(P1)
}

# Auxiliary function for diversity optimization
DivOpt <- function(P, r, b, c, L, U, func, Dis){
  i <- 0
  while(i < c){
    P1 <- NULL
    j <- 0
    while (j <= 2*r-nrow(P)) {
      bbb <- NULL
      aaa <- sample(nrow(P), 2)
      bbb <- crossover1(P[aaa,], L, U)
      bbb <- rbind(bbb, mutate(P[aaa[1],], L, U))
      bbb <- rbind(bbb, mutate(P[aaa[2],], L, U))
      bbb <- cbind(bbb, apply(bbb, 1, func))
      P1 <- rbind(P1, bbb[bbb[,ncol(bbb)] < b,-ncol(bbb)])
      j <- nrow(P1)
    }
    P2 <- rbind(P, P1)
    P4 <- DivSelect_greedy(P2, r, Dis)
    i <- i + 1
  }
  return(P)
}

## NOAH Function
# n = 200 #population size -> n = 2*r
# v = 0.85 #the higher the threshold, more diverse 
# g = 20 #number of geneartion iteration
# r = 50 #number of diverse population,the final output 
# c = 1
# L = c(0,0) # apply(df_b[,1:2], 2, min)
# U = c(1,1) # apply(df_b[,1:2], 2, max)
# func = ObjFunc_Eu
NOAH <- function(n, v, g, r, c, L, U, func, Dis) {
  b <- Inf
  P <- NULL
  i <- 0
  while (i < c) {
    # ObjOpt
    GA <- ga(type = "real-valued",
             fitness =  function(x) -func(x), 
             lower = L, upper = U, pcrossover = 0.5, pmutation = 0.5,
             popSize = n, maxiter = g, maxFitness = b)
    # BoundChange
    aa <- BoundChange(P, GA@population, b, r, func)
    P <- aa$Pop
    if((b - aa$bound) < v*b){i <- i+1}else{i <- 0}
    b <- aa$bound
    print(b)
    # DivOpt
    if(r > 1){P <- DivOpt(P, r, b, c, L, U, func, Dis)}
  }
  return(P)
}

# c = 1
# g = 5
NOAH2 <- function(n, v, g, r, c, L, U, func, Dis) {
  b <- Inf
  P <- NULL
  i <- 0
  while (i < c) {
    # ObjOpt
    GA <- ga(type = "real-valued",
             fitness =  function(x) -func(x), 
             lower = L, upper = U, pcrossover = 0.5, pmutation = 0.5,
             popSize = n, maxiter = g, maxFitness = b)
    # BoundChange
    aa <- BoundChange(P, GA@population, b, r, func)
    P <- aa$Pop
    # if((b - aa$bound) < v*b){i <- i+1}else{i <- 0}
    i <- i+1
    b <- aa$bound
    print(b)
    # DivOpt
    if(r > 1){P <- DivOpt(P, r, b, c, L, U, func, Dis)}
  }
  return(P)
}


###########################################################################
##### Optimization #####
# GA
ObjFunc_Eu <- function(x){sqrt(sum((x-p)^2))} # Objective function Euclidean distance
ObjFunc_M <- function(x){sum(abs(x-p))} # Objective function Manhattan distance
ObjFunc_C <- function(x){sum(abs(x-p)/(abs(x)+abs(p)))} # Objective function Canberra distance

ObjFunc1 <- function(x){mean(apply(df_b[,1:2], 1, FUN = function(y)(sqrt(sum((x-y)^2)))))} # Objective function sum of Euclidean distance
ObjFunc1 <- function(x){mean(apply(df_b[,1:2], 1, FUN = function(y)(abs(x-y))))} # Objective function sum of MAnhattan distance

ObjFunc_Eu_S <- function(x){mean(apply(df_b, 1, FUN = function(y)(sqrt(sum((x-y)^2)))))} # Objective function sum of Euclidean distance for real datasets
ObjFunc_M_S <- function(x){mean(apply(df_b, 1, FUN = function(y)(sum(abs(x-y)))))} # Objective function sum of Manhattan distance for real datasets
ObjFunc_C_S <- function(x){mean(apply(df_b, 1, FUN = function(y)(sum(abs(x-y)/(abs(x)+abs(y))))))} # Objective function sum of Canberra distance for real datasets

ObjFunc_Eu_S_dado <- function(x){mean(apply(df_b_cl_dado_sub, 1, FUN = function(y)(sqrt(sum((x-y)^2)))))} # Objective function sum of Euclidean distance for real datasets
ObjFunc_M_S_dado <- function(x){mean(apply(df_b_cl_dado_sub, 1, FUN = function(y)(sum(abs(x-y)))))} # Objective function sum of Manhattan distance for real datasets
ObjFunc_C_S_dado <- function(x){mean(apply(df_b_cl_dado_sub, 1, FUN = function(y)(sum(abs(x-y)/(abs(x)+abs(y))))))} # Objective function sum of Canberra distance for real datasets



