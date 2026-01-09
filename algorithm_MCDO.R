# with majority cleaning

source(paste0(Path_Imb_Div,'auxiliary_function.r'))

source(paste0(Path_Imb_Div,'NBDOS2.r'))

clusterD_nbdos_multi_2 <- function(D, OSR,  k1, k2 , rTh , nTh  ) {
  
ObjFunc_C <- function(x){sum(abs(x-p)/(abs(x)+abs(p)))}
ObjFunc_Eu_S_dado <- function(x){mean(apply(df_b_cl_dado_sub_dat, 1, FUN = function(y)(sqrt(sum((x-y)^2)))))} # Objective function sum of Euclidean distance for real datasets

####### step 1 - NBDOS Cluster #######

unique_classes <- unique(D$class)

class_counts <- table(D$class)
majority_count <- max(class_counts)

minority_df_list <- list()

synthetic_minority_df_list <- list()


iterative_cleaning <- function(df, k, rTh, max_iterations = 2) {
  cl = ncol(df)
  
  df$id <- seq_len(nrow(df))
  
  original_df <- df


  df$remove <- FALSE  # Initialize removal column
  
  iteration <- 0  # Counter for iterations

  repeat {
    iteration <- iteration + 1
    
    # Extract features and classes

    features <- df[, 1:(cl-1)]
    classes <- df$class
    
    # Compute k-nearest neighbors
    knn_result <- get.knn(features, k = k)
    
    # Label points for removal
    removal_label <- sapply(1:nrow(df), function(i) {
      neighbor_indices <- knn_result$nn.index[i, ]
      neighbor_classes <- classes[neighbor_indices]
      proportion_different <- mean(neighbor_classes != classes[i])
      proportion_different > rTh
    })
    
    # Check if any new points are marked for removal
    if (!any(removal_label & !df$remove) || iteration >= max_iterations) break
    
    # Update removal labels
    df$remove <- removal_label
    
    # Remove points marked for removal
    df <- df[!df$remove, ]
  }
  
  # Add the final removal labels to the original dataset
  original_df$remove <- !original_df$id %in% df$id
  return(original_df)
}


# perform data cleaning

majority_df <- iterative_cleaning(D, nTh, rTh, max_iterations = 3)

majority_df =  subset(majority_df, class == 0)
majority_df$class = ifelse(majority_df$remove == T, paste(0, "Remove", sep = "_"), 0)
minority_df = subset(D, class != 0)
majority_df_cleaned = majority_df[majority_df$remove == FALSE, ]
D_cleaned = rbind(majority_df_cleaned[1:cl],minority_df )
#

for (class_name in unique_classes) { # for each class, apply NBDOS
  if (class_name == 0) {
    next
  }
  
# 
#   minority_df = subset(D, class != 0) 
#   
#   plot_data = rbind(majority_df %>% dplyr::select(x, y, class),minority_df )
#   
#   ggplot(plot_data, aes(x = x, y = y, color = class)) +
#     geom_point(size = 2) +  # Set point size 
#     geom_point(
#       data = subset(plot_data, class == "0_Remove"), # Filter for '0_Remove'
#       shape = 4, # 'X' mark shape
#       size = 3,  # Size of 'X' mark
#       color = "red" # Red color for 'X' mark
#     ) 
#   
  sub_df <- D_cleaned
#   
  #set class to 0
  sub_df$class <- ifelse(sub_df$class == class_name, class_name, 0) #update label of the rest of data points as majority
  
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

# Cl = 0 - isolated
# Cl = numbers other than 0 - clusters

combined_minority_df$mth <- ifelse(grepl("Trapped", combined_minority_df$Cl), "diwo", "dado")


#######  step 2 determine generation points ##############


# number of points within each cluster 
for (class_name in unique_classes) { # for each class, apply NBDOS
  if (class_name == 0) {
    next
  }
  minority_df = combined_minority_df[combined_minority_df$class == class_name,]
  
  cluster_count <- as.data.frame(minority_df %>% group_by(Cl) %>% count())
  minority_cl_n <- merge(minority_df, cluster_count, by.x = "Cl",
                   by.y = "Cl", all.x = TRUE, all.y = FALSE)

  
  # account for dups in dado
  for (n_var in unique(minority_cl_n$Cl)){
    dat_nvar <- minority_cl_n[minority_cl_n$mth == 'dado' & minority_cl_n$Cl == n_var,]
    dat_nvar_dup <- as.data.frame(dat_nvar %>% find_duplicates())
    if (nrow(dat_nvar) == nrow(dat_nvar_dup)&(nrow(dat_nvar) != 0 )){
      minority_cl_n[minority_cl_n$mth == 'dado' & minority_cl_n$Cl == n_var,]$mth <- 'diwo'
      minority_cl_n$Cl[minority_cl_n$mth == 'diwo' & minority_cl_n$Cl == n_var] <- 0
      minority_cl_n$n[minority_cl_n$mth == 'diwo'] <- as.integer(minority_cl_n[minority_cl_n$mth == 'diwo',] %>% count())
    }
  }
  
  rmname <- c('Cl','n','density','mth', 'class')
  
  #DADO
  Snew_dado <- NULL
  df_b_cl_dado <- minority_cl_n[minority_cl_n$mth == 'dado',]
  
  OS = OSR * majority_count - nrow(minority_cl_n)
  
  if(length(unique(df_b_cl_dado$Cl))>=1) {
    df_b_cl_dado$density <- nrow(df_b_cl_dado)/nrow(minority_cl_n)
    for (i in unique(df_b_cl_dado$Cl)){
      df_b_cl_dado_sub <- subset(df_b_cl_dado,df_b_cl_dado$Cl == i)
      density = unique(df_b_cl_dado_sub$density)*unique(df_b_cl_dado_sub$n)/nrow(df_b_cl_dado)
      
      df_b_cl_dado_sub_dat <- df_b_cl_dado_sub[, -which(names(df_b_cl_dado_sub) %in% rmname)]
      
      P <- NOAH2(n = as.numeric(ceiling(OS*density))+1, v=0.5, g = 5, r = as.numeric(ceiling(OS*density)), c = 1, 
                 L = apply(df_b_cl_dado_sub_dat, 2, min)*0.9999, 
                 U = apply(df_b_cl_dado_sub_dat, 2, max)*1.0001,
                 func=ObjFunc_Eu_S_dado,Dis = "euclidean")
      
      Snew_dado <- rbind.data.frame(Snew_dado,P)
      
    }
  }
  
  
  
  #DIWO
  
  # DIWO
  print("DIWO")
  P <- NULL
  df_b_cl_diwo <- minority_cl_n[minority_cl_n$mth == 'diwo',]
  OS = OSR * majority_count - nrow(minority_cl_n)
  
  minority_cl_n_dat <- minority_cl_n[, -which(names(minority_cl_n) %in% rmname)]
  
  if(length(unique(df_b_cl_diwo$Cl))>=1) {
    df_b_cl_diwo$density <- nrow(df_b_cl_diwo)/nrow(minority_cl_n)
    Dis = "canberra"
    r <- OS*unique(df_b_cl_diwo$density)
    df_b_cl_diwo_sub <- df_b_cl_diwo[, -which(names(df_b_cl_diwo) %in% rmname)]
    r1 <- if(r/nrow(df_b_cl_diwo_sub) < 1){1}else{round(r/nrow(df_b_cl_diwo_sub),0)}
    for (ii in 1:nrow(df_b_cl_diwo_sub)) {
      p <- df_b_cl_diwo_sub[ii,]
      p_dat <- p
      print(p)
      print(1)
      P1 <- NOAH(r1+50, v = 0.5, 20, r1, 2, apply(minority_cl_n_dat, 2, min)*0.999, apply(minority_cl_n_dat, 2, max)*1.001, ObjFunc_C, Dis = "canberra")
      print(2)
      P <- rbind(P, P1)
    }
    
    Snew_diwo <- P
  }
  
  
  if (exists("Snew_dado")&exists("Snew_diwo")&(length(Snew_dado)!=0)&(length(Snew_dado)!=0)){
    df_G <- merge(Snew_dado,Snew_diwo, all = TRUE)
  } else if (exists("Snew_dado")&(length(Snew_dado)!=0) ) {
    df_G <- Snew_dado
  } else if (exists("Snew_diwo")&(length(Snew_diwo)!=0)) {
    df_G <- Snew_diwo
  }
  df_G$Class <- class_name
  
  synthetic_minority_df_list[[paste0("minority_cl_", class_name)]] <- df_G

  }


synthetic_minority_df = do.call(rbind, synthetic_minority_df_list)

return(synthetic_minority_df)
  
  
}


