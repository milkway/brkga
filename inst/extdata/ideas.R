
library(tidyverse)

distances <- c(0, 1, 1, 1, 1, 4, 2, 1,
               1, 0, 3, 3, 3, 2, 1, 1, 
               1, 3, 0, 4, 2, 2, 1, 2,
               1, 3, 4, 0, 3, 3, 2, 1,
               1, 3, 2, 3, 0, 1, 2, 3,
               4, 2, 2, 3, 1, 0, 3, 1,
               2, 1, 1, 2, 2, 3, 0, 1,
               1, 1, 2, 1, 3, 1, 1, 0) %>% matrix(nrow = 8, byrow = TRUE)
colnames(distances) <- 1:8
rownames(distances) <- 1:8


#tour <- c(1, 2, 3)
tour <- c(2, 8, 6)
#tour <- sample(1:8, 3)
nontour <- setdiff(1:8, tour)

v <- matrix(c(rep(1,3), rep(0,5)), ncol = 1) 


M <- distances[c(tour, nontour), c(tour, nontour)]

#cat("Fitness: ", sum(distances[tour, tour])/2)
cat("Fitness: ", sum(M[1:3, 1:3])/2)

G <- M %*% v
R <- M[1:3, 4:8] %>% 
  apply(2, function(x){-x - G[1:3]}) %>% 
  apply(1, function(x){x + G[4:8]}) %>% t
#colnames(R) <- c(nontour)
#rownames(R) <- c(tour)


gtz <- as_tibble(which(R >= 0, arr.ind = TRUE)) %>% 
  mutate(delta = apply(., 1, function(x) {R[x[1], x[2]]}),
         prob = (delta+1)/sum(1+delta)) # +1 to allow zero delta a change


OutIn <- gtz %>% sample_n(1, weight = prob)
OutIn


aux <- tour[OutIn$row]
tour[OutIn$row] <- nontour[OutIn$col]
nontour[OutIn$col] <- aux

M2 <- distances[c(tour, nontour), c(tour, nontour)]

cat("Fitness: ", sum(M2[1:3, 1:3])/2)


G2 <- M2 %*% v
R2 <- M2[1:3, 4:8] %>% 
  apply(2, function(x){-x - G2[1:3]}) %>% 
  apply(1, function(x){x + G2[4:8]}) %>% t


gtz <- as_tibble(which(R2 >= 0, arr.ind = TRUE)) %>% 
  mutate(delta = apply(., 1, function(x) {R2[x[1], x[2]]}),
         prob = (delta+1)/sum(1+delta)) # +1 to allow zero delta a change


OutIn <- gtz %>% sample_n(1, weight = prob)
OutIn

aux <- tour[OutIn$row]
tour[OutIn$row] <- nontour[OutIn$col]
nontour[OutIn$col] <- aux


M3 <- distances[c(tour, nontour), c(tour, nontour)]
cat("Fitness: ", sum(M3[1:3, 1:3])/2)


#####


G3 <- M3 %*% v
R3 <- M3[1:3, 4:8] %>% 
  apply(2, function(x){-x - G3[1:3]}) %>% 
  apply(1, function(x){x + G3[4:8]}) %>% t


gtz <- as_tibble(which(R3 >= 0, arr.ind = TRUE)) %>% 
  mutate(delta = apply(., 1, function(x) {R3[x[1], x[2]]}),
         prob = (delta+1)/sum(1+delta)) # +1 to allow zero delta a change


OutIn <- gtz %>% sample_n(1, weight = prob)
OutIn

aux <- tour[OutIn$row]
tour[OutIn$row] <- nontour[OutIn$col]
nontour[OutIn$col] <- aux


M4 <- distances[c(tour, nontour), c(tour, nontour)]
cat("Fitness: ", sum(M4[1:3, 1:3])/2)

###

G4 <- M4 %*% v
R4 <- M4[1:3, 4:8] %>% 
  apply(2, function(x){-x - G4[1:3]}) %>% 
  apply(1, function(x){x + G4[4:8]}) %>% t


gtz <- as_tibble(which(R4 >= 0, arr.ind = TRUE)) %>% 
  mutate(delta = apply(., 1, function(x) {R4[x[1], x[2]]}),
         prob = (delta+1)/sum(1+delta)) # +1 to allow zero delta a change


OutIn <- gtz %>% sample_n(1, weight = prob)
OutIn

aux <- tour[OutIn$row]
tour[OutIn$row] <- nontour[OutIn$col]
nontour[OutIn$col] <- aux


M5 <- distances[c(tour, nontour), c(tour, nontour)]
cat("Fitness: ", sum(M5[1:3, 1:3])/2)

###


G5 <- M5 %*% v
R5 <- M5[1:3, 4:8] %>% 
  apply(2, function(x){-x - G5[1:3]}) %>% 
  apply(1, function(x){x + G5[4:8]}) %>% t


gtz <- as_tibble(which(R5 >= 0, arr.ind = TRUE)) %>% 
  mutate(delta = apply(., 1, function(x) {R5[x[1], x[2]]}),
         prob = (delta+1)/sum(1+delta)) # +1 to allow zero delta a change


OutIn <- gtz %>% sample_n(1, weight = prob)
OutIn

aux <- tour[OutIn$row]
tour[OutIn$row] <- nontour[OutIn$col]
nontour[OutIn$col] <- aux


M6 <- distances[c(tour, nontour), c(tour, nontour)]
cat("Fitness: ", sum(M6[1:3, 1:3])/2)


#######################################3
N <- 500
M <- 50
distances <- read_rds("inst/extdata/MDG.1.a.n500m50.rds")
colnames(distances) <- 1:N
rownames(distances) <- 1:N

#set.seed(1)
tour <- sample(1:N, M) #OutIn <- gtz %>% sample_n(1, weight = prob)
nontour <- setdiff(1:N, tour)


v <- matrix(c(rep(1,M), rep(0,N-M)), ncol = 1) 


distM <- distances[c(tour, nontour), c(tour, nontour)]

cat("Fitness: ", sum(distM[1:M, 1:M])/2)

for (j in 1:100){
  G <- distM %*% v
  R <- distM[1:M, (M+1):N] %>% 
    apply(2, function(x){-x - G[1:M]}) %>% 
    apply(1, function(x){x + G[(M+1):N]}) %>% t

  gtz <- as_tibble(which(R >= -100000, arr.ind = TRUE)) %>% 
    mutate(delta = apply(., 1, function(x) {R[x[1], x[2]]}),
           prob = exp(.5*delta)) # +1 to allow zero delta a change

  OutIn <- gtz %>% sample_n(1, weight = prob)
  #OutIn <- gtz %>% top_n(1, delta)
  #OutIn
  
  
  aux <- tour[OutIn$row]
  tour[OutIn$row] <- nontour[OutIn$col]
  nontour[OutIn$col] <- aux
  
  distM <- distances[c(tour, nontour), c(tour, nontour)]

    
  cat("\nFitness: ", round(sum(distM[1:M, 1:M])/2, 1), " Delta: ", OutIn$delta)
  
}

##################################


for (i in 1:15) {
  tour <- tibble(N = 1:500, Vetor = rst$vectors[,i]) %>% top_n(50, Vetor) %>% select(N) %>% unlist(use.names = FALSE)
  cat("\n i ", i, " Fitness: ", sum(distances[tour,tour])/2)
}


distances <- read_rds("inst/extdata/MDG.1.b.n500m50.rds")
DtM <- (distances - min(distances))/(max(distances) - min(distances))
rst <- eigen(DtM, symmetric = TRUE)
tour <- tibble(N = 1:500, Vetor = rst$vectors[,1]) %>% top_n(50) %>% select(N) %>% unlist(use.names = FALSE)
nontour <- setdiff(1:N, tour)


v <- matrix(c(rep(1,M), rep(0,N-M)), ncol = 1) 


distM <- distances[c(tour, nontour), c(tour, nontour)]

cat("Fitness: ", sum(distM[1:M, 1:M])/2)

for (j in 1:100){
  G <- distM %*% v
  R <- distM[1:M, (M+1):N] %>% 
    apply(2, function(x){-x - G[1:M]}) %>% 
    apply(1, function(x){x + G[(M+1):N]}) %>% t
  
  gtz <- as_tibble(which(R >= -500, arr.ind = TRUE)) %>% 
    mutate(delta = apply(., 1, function(x) {R[x[1], x[2]]}),
           prob = exp(delta/1000)) # +1 to allow zero delta a change
  
  OutIn <- gtz %>% sample_n(1, weight = prob)
  #OutIn <- gtz %>% top_n(10, delta)
  #OutIn
  
  
  aux <- tour[OutIn$row]
  #sair <- sample(1:50,1)
  #aux <- tour[sair]
  #tour[sair] <- nontour[OutIn$col]
  tour[OutIn$row] <- nontour[OutIn$col]
  nontour[OutIn$col] <- aux
  
  distM <- distances[c(tour, nontour), c(tour, nontour)]
  
  
  cat("\nFitness: ", round(sum(distM[1:M, 1:M])/2, 1), " Delta: ", OutIn$delta)
  
}

######
library(tidyverse)
N <- 2000
M <- 200
distances <- read_rds("inst/extdata/MDG.21.a.n2000m200.rds")

distances.d <- as.dist(distances)
cluster.mdp <- hclust(distances.d, method = "ward.D2")
tour <- 
  tibble(
    N = 1:N,
    Cluster = cutree(cluster.mdp, M) 
  ) %>% 
  group_by(Cluster) %>% 
  sample_n(1) %>% ungroup() %>% select(N) %>% unlist(use.names = FALSE) 

nontour <- setdiff(1:N, tour)


v <- matrix(c(rep(1,M), rep(0,N-M)), ncol = 1) 


distM <- distances[c(tour, nontour), c(tour, nontour)]

cat("Fitness: ", sum(distM[1:M, 1:M])/2)

start_ <- Sys.time()
end_ <- Sys.time()
MAX_TIME <- 60
Fitness <- 0
while((end_ - start_ ) <= MAX_TIME){
  G <- distM %*% v
  R <- distM[1:M, (M+1):N] %>% 
    apply(2, function(x){-x - G[1:M]}) %>% 
    apply(1, function(x){x + G[(M+1):N]}) %>% t
  
  gtz <- as_tibble(which(R >= 0, arr.ind = TRUE)) %>% 
    mutate(delta = apply(., 1, function(x) {R[x[1], x[2]]}),
           prob = exp(delta/1000)) #+1 to allow zero delta a change
  
  if (nrow(gtz) != 0) {
    OutIn <- gtz %>% sample_n(1, weight = prob)
    aux <- tour[OutIn$row]
    tour[OutIn$row] <- nontour[OutIn$col]
    nontour[OutIn$col] <- aux
  } 
  else{
    out_ <- sample(1:M,20)
    in_ <- sample(1:(N-M),20)
    aux <- tour[out_]
    tour[out_] <- nontour[in_]
    nontour[in_] <- aux
  }
  
  
  distM <- distances[c(tour, nontour), c(tour, nontour)]
  
  val <- round(sum(distM[1:M, 1:M])/2, 1)
  if (val >= Fitness) 
    Fitness <- val 
  cat("\nFitness: ", val, "Best: ", Fitness, " Delta: ", OutIn$delta)
  end_ <- Sys.time()
}
cat("\nBest Fitness: ", Fitness)
print(end_ - start_)
