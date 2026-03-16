################################# Functions #############################################################
# Auxiliary function for Diversity (Solow-Polasky) measure
Div_SP <- function(X, Dis){ #X is the set of datapoints

  X_mat <- as.matrix(X)
  M <- as.matrix(exp(-(dist(X_mat, diag = TRUE, upper = TRUE, method = Dis))))
  
  # Try fast Cholesky inversion first, fallback to ginv if matrix is singular
  M_1 <- tryCatch({
    chol2inv(chol(M))
  }, error = function(e) {
    MASS::ginv(M)
  })
  
  return(sum(M_1))
}

# Auxiliary function for diversity difference using Div_SP()
Div_SP_diff1 <- function(X, n, Dis){
  a <- Div_SP(X, Dis) - Div_SP(X[-n,], Dis)
  return(a)
}

# Auxiliary function for diversity difference with less computation
Div_SP_diff2 <- function(X, n, Dis){
  X_mat <- as.matrix(X)
  M <- as.matrix(exp(-(dist(X_mat, diag = TRUE, upper = TRUE, method = Dis))))
  
  M_1 <- tryCatch({
    chol2inv(chol(M))
  }, error = function(e) {
    MASS::ginv(M)
  })
  
  diff <- sum(M_1[,n]) / M_1[n,n]
  return(diff)
}

# Auxiliary function for adaptive change of bound
BoundChange <- function(P, Pn, b, r, func){
  # Use matrix binding
  PPn_mat <- as.matrix(rbind(P, Pn))
  
  # Evaluate fitness
  vals <- apply(PPn_mat, 1, func) 
  
  # Bind and sort
  PPn_full <- cbind(PPn_mat, vals)
  val_col <- ncol(PPn_full)
  
  P1 <- PPn_full[order(PPn_full[, val_col]), ]
  P2 <- as.data.frame(P1[1:r, 1:(val_col-1), drop=FALSE])
  new_b <- P1[r, val_col]
  
  return(list(Pop = P2, bound = new_b))
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

# Auxiliary function for ranking instances based diversity difference
Div_SP_diff_rank <- function(X, Dis){
  X_mat <- as.matrix(X)
  P_col <- ncol(X_mat) + 1
  M <- as.matrix(exp(-(dist(X_mat, diag = TRUE, upper = TRUE, method = Dis))))
  
  M_1 <- tryCatch({
    chol2inv(chol(M))
  }, error = function(e) {
    MASS::ginv(M)
  })
  
  # Vectorized rank calculation instead of for loop
  rank_vals <- colSums(M_1) / diag(M_1)
  X_res <- cbind(X_mat, rank_vals)
  
  P1 <- X_res[order(X_res[, P_col], decreasing = TRUE), ]
  return(as.data.frame(P1))
}

# Auxiliary function for diversity selection
DivSelect_greedy <- function(P, r, Dis){
  P1 <- P
  P_row <- nrow(P)
  
  # Iteratively drop the least diverse rows
  for(i in seq_len(P_row - r)) {
    P1 <- Div_SP_diff_rank(P1, Dis)
    P1 <- P1[-nrow(P1), 1:ncol(P), drop=FALSE]
  }
  return(P1)
}

# Auxiliary function for diversity optimization
DivOpt <- function(P, r, b, c, L, U, func, Dis){
  i <- 0
  while(i < c){
    P1_list <- list()
    j <- 0
    target_rows <- 2 * r - nrow(P)
    
    # Store candidates in a list to prevent rbind memory fragmentation
    while (j <= target_rows) {
      aaa <- sample(nrow(P), 2)
      
      Osp <- crossover1(P[aaa, ], L, U)
      mut1 <- mutate(as.numeric(P[aaa[1], ]), L, U)
      mut2 <- mutate(as.numeric(P[aaa[2], ]), L, U)
      
      bbb <- rbind(Osp, mut1, mut2)
      bbb_vals <- apply(bbb, 1, func)
      
      valid_idx <- bbb_vals < b
      if (any(valid_idx)) {
        P1_list[[length(P1_list) + 1]] <- bbb[valid_idx, , drop=FALSE]
        j <- j + sum(valid_idx)
      }
    }
    
    P1 <- do.call(rbind, P1_list)
    P2 <- rbind(P, P1)
    P <- DivSelect_greedy(P2, r, Dis)
    i <- i + 1
  }
  return(P)
}

## NOAH Function
NOAH <- function(n, v, g, r, c, L, U, func, Dis) {
  b <- Inf
  P <- NULL
  i <- 0
  while (i < c) {
    GA <- ga(type = "real-valued",
             fitness =  function(x) -func(x), 
             lower = L, upper = U, pcrossover = 0.5, pmutation = 0.5,
             popSize = n, maxiter = g, maxFitness = b)
    
    aa <- BoundChange(P, GA@population, b, r, func)
    P <- aa$Pop
    
    if((b - aa$bound) < v*b){ i <- i+1 } else { i <- 0 }
    
    b <- aa$bound
    print(b)
    
    if(r > 1){ P <- DivOpt(P, r, b, c, L, U, func, Dis) }
  }
  return(P)
}

NOAH2 <- function(n, v, g, r, c, L, U, func, Dis) {
  b <- Inf
  P <- NULL
  i <- 0
  while (i < c) {
    GA <- ga(type = "real-valued",
             fitness =  function(x) -func(x), 
             lower = L, upper = U, pcrossover = 0.5, pmutation = 0.5,
             popSize = n, maxiter = g, maxFitness = b)
    
    aa <- BoundChange(P, GA@population, b, r, func)
    P <- aa$Pop
    
    i <- i+1
    b <- aa$bound
    print(b)
    
    if(r > 1){ P <- DivOpt(P, r, b, c, L, U, func, Dis) }
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



