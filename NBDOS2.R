
## Auxiliary functions for NBDOS

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




