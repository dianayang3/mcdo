# with majority cleaning

source(paste0(Path_Imb_Div,'auxiliary_function.r'))

source(paste0(Path_Imb_Div,'NBDOS2.r'))

clusterD_nbdos_multi_2 <- function(D, OSR, k1, k2, rTh, nTh) {
  
  
  ObjFunc_C <- function(x) { sum(abs(x - p_vec) / (abs(x) + abs(p_vec))) }
  ObjFunc_Eu_S_dado <- function(x) { mean(sqrt(rowSums(sweep(df_b_cl_dado_sub_mat, 2, x, "-")^2))) }
  
  ####### step 1 - NBDOS Cluster #######
  
  unique_classes <- unique(D$class)
  class_counts <- table(D$class)
  majority_count <- max(class_counts)
  
  minority_df_list <- list()
  synthetic_minority_df_list <- list()
  
  # 2. VECTORIZED ITERATIVE CLEANING
  iterative_cleaning <- function(df, k, rTh, max_iterations = 2) {
    cl = ncol(df)
    df$id <- seq_len(nrow(df))
    original_df <- df
    df$remove <- FALSE  
    iteration <- 0  
    
    repeat {
      iteration <- iteration + 1
      
      # Extract features as a matrix for faster FNN processing
      features <- as.matrix(df[, 1:(cl-1)])
      classes <- df$class
      
      # Compute k-nearest neighbors
      knn_result <- FNN::get.knn(features, k = k)
      
      # Vectorized evaluation of neighbor proportions
      neighbor_classes <- matrix(classes[knn_result$nn.index], nrow = nrow(features), ncol = k)
      proportion_different <- rowMeans(neighbor_classes != classes)
      removal_label <- proportion_different > rTh
      
      if (!any(removal_label & !df$remove) || iteration >= max_iterations) break
      
      df$remove <- removal_label
      df <- df[!df$remove, ]
    }
    
    original_df$remove <- !original_df$id %in% df$id
    return(original_df)
  }
  
  # Perform data cleaning
  majority_df <- iterative_cleaning(D, nTh, rTh, max_iterations = 3)
  majority_df <- subset(majority_df, class == 0)
  majority_df$class <- ifelse(majority_df$remove == T, paste(0, "Remove", sep = "_"), 0)
  
  minority_df <- subset(D, class != 0)
  majority_df_cleaned <- majority_df[majority_df$remove == FALSE, ]
  
  # Ensure clean column matching (replaces the hardcoded [1:cl] from original)
  D_cleaned <- rbind(majority_df_cleaned[, names(minority_df)], minority_df)
  
  for (class_name in unique_classes) { 
    if (class_name == 0) next
    
    sub_df <- D_cleaned
    sub_df$class <- ifelse(sub_df$class == class_name, class_name, 0) 
    
    minority_cl <- NULL 
    tryCatch({
      minority_cl = NBDOS2(sub_df, class_name, k1, k2, rTh, nTh) 
    }, error = function(e) {
      minority_cl <<- transform(sub_df[sub_df$class == class_name, ], Cl = 0)
    })
    
    minority_cl$Cl = ifelse(minority_cl$Cl == 0, paste(class_name, "Trapped", sep = "_"), paste(class_name, minority_cl$Cl, sep = "_"))
    minority_cl$id <- NULL
    minority_df_list[[paste0("minority_cl_", class_name)]] <- minority_cl 
  }
  
  combined_minority_df = do.call(rbind, minority_df_list)
  combined_minority_df$id <- NULL
  combined_minority_df$mth <- ifelse(grepl("Trapped", combined_minority_df$Cl), "diwo", "dado")
  
  ####### step 2 determine generation points ##############
  
  for (class_name in unique_classes) { 
    if (class_name == 0) next
    
    minority_df = combined_minority_df[combined_minority_df$class == class_name,]
    cluster_count <- as.data.frame(minority_df %>% dplyr::group_by(Cl) %>% dplyr::count())
    minority_cl_n <- merge(minority_df, cluster_count, by = "Cl", all.x = TRUE, all.y = FALSE)
    
    # 3. OPTIMIZED DUPLICATE CHECK
    # Bypasses the pipe-heavy find_duplicates() in favor of base R's duplicated()
    for (n_var in unique(minority_cl_n$Cl)) {
      idx <- minority_cl_n$mth == 'dado' & minority_cl_n$Cl == n_var
      dat_nvar <- minority_cl_n[idx, ]
      
      if (nrow(dat_nvar) > 0) {
        is_dup <- duplicated(dat_nvar) | duplicated(dat_nvar, fromLast = TRUE)
        if (all(is_dup)) {
          minority_cl_n$mth[idx] <- 'diwo'
          minority_cl_n$Cl[idx] <- 0
          minority_cl_n$n[minority_cl_n$mth == 'diwo'] <- sum(minority_cl_n$mth == 'diwo')
        }
      }
    }
    
    rmname <- c('Cl','n','density','mth', 'class')
    OS = OSR * majority_count - nrow(minority_cl_n)
    
    # DADO
    Snew_dado <- NULL
    df_b_cl_dado <- minority_cl_n[minority_cl_n$mth == 'dado',]
    
    if (length(unique(df_b_cl_dado$Cl)) >= 1) {
      df_b_cl_dado$density <- nrow(df_b_cl_dado)/nrow(minority_cl_n)
      
      for (i in unique(df_b_cl_dado$Cl)) {
        df_b_cl_dado_sub <- subset(df_b_cl_dado, Cl == i)
        density = unique(df_b_cl_dado_sub$density) * unique(df_b_cl_dado_sub$n) / nrow(df_b_cl_dado)
        
        df_b_cl_dado_sub_dat <- df_b_cl_dado_sub[, -which(names(df_b_cl_dado_sub) %in% rmname)]
        
        # Coerce to matrix globally for the current iteration so the GA objective function runs extremely fast
        df_b_cl_dado_sub_mat <<- as.matrix(df_b_cl_dado_sub_dat)
        
        P <- NOAH2(n = as.numeric(ceiling(OS*density))+1, v=0.5, g = 5, r = as.numeric(ceiling(OS*density)), c = 1, 
                   L = apply(df_b_cl_dado_sub_dat, 2, min)*0.9999, 
                   U = apply(df_b_cl_dado_sub_dat, 2, max)*1.0001,
                   func=ObjFunc_Eu_S_dado, Dis = "euclidean")
        
        Snew_dado <- rbind.data.frame(Snew_dado, P)
      }
    }
    
    # DIWO
    P <- NULL
    df_b_cl_diwo <- minority_cl_n[minority_cl_n$mth == 'diwo',]
    minority_cl_n_dat <- minority_cl_n[, -which(names(minority_cl_n) %in% rmname)]
    
    if (length(unique(df_b_cl_diwo$Cl)) >= 1) {
      df_b_cl_diwo$density <- nrow(df_b_cl_diwo)/nrow(minority_cl_n)
      r <- OS * unique(df_b_cl_diwo$density)
      df_b_cl_diwo_sub <- df_b_cl_diwo[, -which(names(df_b_cl_diwo) %in% rmname)]
      r1 <- if (r/nrow(df_b_cl_diwo_sub) < 1) { 1 } else { round(r/nrow(df_b_cl_diwo_sub), 0) }
      
      # Pull constants outside the loop
      min_bounds <- apply(minority_cl_n_dat, 2, min) * 0.999
      max_bounds <- apply(minority_cl_n_dat, 2, max) * 1.001
      
      for (ii in 1:nrow(df_b_cl_diwo_sub)) {
        # Cast to numeric vector to strip dataframe overhead from the objective function
        p_vec <<- as.numeric(df_b_cl_diwo_sub[ii,])
        
        P1 <- NOAH(r1+50, v = 0.5, g = 20, r = r1, c = 2, 
                   L = min_bounds, U = max_bounds, 
                   func = ObjFunc_C, Dis = "canberra")
        P <- rbind(P, P1)
      }
      Snew_diwo <- P
    }
    
    if (exists("Snew_dado") && exists("Snew_diwo") && length(Snew_dado) > 0 && length(Snew_diwo) > 0) {
      df_G <- merge(Snew_dado, Snew_diwo, all = TRUE)
    } else if (exists("Snew_dado") && length(Snew_dado) > 0) {
      df_G <- Snew_dado
    } else if (exists("Snew_diwo") && length(Snew_diwo) > 0) {
      df_G <- Snew_diwo
    }
    
    df_G$Class <- class_name
    synthetic_minority_df_list[[paste0("minority_cl_", class_name)]] <- df_G
  }
  
  synthetic_minority_df = do.call(rbind, synthetic_minority_df_list)
  return(synthetic_minority_df)
}


