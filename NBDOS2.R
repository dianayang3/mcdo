

memory.limit(size = 56000)
options(scipen=10000)



## Auxiliary functions for SMOM
NBDOS <- function(Sc, k, ScNk, rTh, nTh){
  Nc <- ncol(Sc)-2
  sfSc <- NULL
  Hck <- NULL
  Sc$Cl <- 0
  ij = 0
  for (i in 1:nrow(Sc)) {
    Tem <- ScNk[i,ScNk[i,] %in% Sc$id]
    if(round(length(Tem)/k,3) >= rTh){
      ij = ij + 1
      sfSc <- rbind.data.frame(sfSc, Sc[i,])
      Hck[ij] <- list(Tem)
    }
  }
  for (ii in 1:nrow(sfSc)) {
    Tem <- get.knnx(Sc[,1:Nc], sfSc[ii,1:Nc], k=k2+1)
    Hck[[ii]] <- union(apply(Hck[[ii]],2,FUN=function(x)(x)), intersect(setdiff(Sc$id[Tem$nn.index],sfSc$id[ii]), sfSc$id))
  }
  curld <- 0
  for (iii in 1:nrow(sfSc)) {
    if(sfSc$Cl[iii]==0){
      curld <- curld + 1
      xi <- sfSc[iii,]
      xi$xj <- iii
      xi$Cl <- curld
      sfC <- xi
      Sc$Cl[Sc$id == iii] <- curld
      while(nrow(sfC)!=0){
        xj <- sfC[1,]
        for (l in 1:length(Hck[[xj$xj]])) {
          xl <- Sc[Sc$id == Hck[[xj$xj]][l],]
          if(xl$Cl==0){
            xl$Cl <- curld
            Sc$Cl[Sc$id == Hck[[xj$xj]][l]] <- curld
            if(xl$id %in% sfSc$id){
              xl$xj <- which(sfSc$id==xl$id)
              sfC <- rbind.data.frame(sfC, xl)
            }
          }
        }
        sfC <- sfC[-1,]
      }
    }
  }
  for (cr in 1:curld) {
    if(nrow(Sc[Sc$Cl==cr,]) < nTh){Sc$Cl[Sc$Cl==cr] <- 0}
  }
  return(Sc)
}

# --- 

NBDOS2 <- function(D, c, k1, k2, rTh, nTh){
  
SI <- NULL
D$id <- as.numeric(rownames(D))
Nc <- ncol(D)-2

Sc <- D[D$class==c,]
Sct <- D[D$class!=c,]

k3 <- max(k1, k2)
ScNk3 <- NULL
ScNk2 <- NULL
ScNk1 <- NULL
Scrk1 <- NULL
SctNk3 <- NULL
SFi <- as.data.frame(matrix(NA, nrow(Sc), nrow(D)))
SFd <- as.data.frame(matrix(NA, nrow(Sc), nrow(D)))

for (i in 1:nrow(Sc)) {
  temp <- get.knnx(Sc[,1:Nc], Sc[i,1:Nc], k=k3+1) #distance matrix of minority classes for the selected row
  ScNk3 <- rbind.data.frame(ScNk3, setdiff(temp$nn.index,i))
  dis_tem <- NULL
  for (ii in 1:ncol(ScNk3)) {
    dis_tem[ii] <- dist(rbind.data.frame(Sc[i,1:Nc], Sc[ScNk3[i,ii],1:Nc]))
  }
  Scrk1[i] <- dis_tem[order(dis_tem)][k1]
  dis_tem1 <- dis_tem[dis_tem <= Scrk1[i]]
  ScNk1 <- rbind.data.frame(ScNk1, ScNk3[i,dis_tem <= Scrk1[i]])
  for (jj in 1:length(dis_tem1)) {
    SFi[i,jj] <- Sc$id[ScNk1[i,jj]]
    SFd[i,jj] <- dis_tem1[jj]
  }
  
  temp <- get.knnx(Sct[,1:Nc], Sc[i,1:Nc], k=k3)
  SctNk3 <- rbind.data.frame(SctNk3, temp$nn.index)
  dis_tem <- NULL
  for (ii in 1:nrow(Sct)) {
    dis_tem[ii] <- dist(rbind.data.frame(Sc[i,1:Nc], Sct[ii,1:Nc]))
    if(dis_tem[ii] <= Scrk1[i]){
      SFi[i,jj+ii] <- Sct$id[ii]
      SFd[i,jj+ii] <- dis_tem[ii]
    }
  }
  union_tem <- rbind.data.frame(Sc[apply(ScNk3[i,],2,FUN=function(x)(x)),], Sct[apply(SctNk3[i,],2,FUN=function(x)(x)),])
  temp <- get.knnx(union_tem[,1:Nc], Sc[i,1:Nc], k=k2)
  ScNk2 <- rbind.data.frame(ScNk2, union_tem$id[temp$nn.index])
} #objective is to obtain a list of k-nearest neighbors of minority instances

ScCl <- NBDOS(Sc, k2, ScNk2, rTh, nTh)
ScO <- ScCl[ScCl$Cl != 0,] # outstanding instances by clusters
ScT <- ScCl[ScCl$Cl == 0,] # trapped instances

return(rbind(ScO,ScT))

}



