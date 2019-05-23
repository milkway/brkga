
library(tidyverse)

distances <- c(0, 1, 1, 1, 1, 4, 2, 1,
               1, 0, 3, 3, 3, 2, 1, 1, 
               1, 3, 0, 4, 2, 2, 1, 2,
               1, 3, 4, 0, 3, 3, 2, 1,
               1, 3, 2, 3, 0, 1, 2, 3,
               4, 2, 2, 3, 1, 0, 3, 1,
               2, 1, 1, 2, 2, 3, 0, 1,
               1, 1, 2, 1, 3, 1, 1, 0) %>% matrix(nrow = 8, byrow = TRUE)

tour <- c(1, 2, 3)
nontour <- setdiff(1:8, tour)
v <- matrix(c(rep(1,3), rep(0,5)), ncol = 1) 


M <- distances[c(tour, nontour), c(tour, nontour)]
colnames(M) <- c(tour, nontour)
rownames(M) <- c(tour, nontour)
cat("Fitness: ", sum(distances[tour, tour])/2)
G <- M %*% v
R <- M[tour, nontour] %>% 
  apply(2, function(x){-x - G[1:3]}) %>% 
  apply(1, function(x){x + G[4:8]}) %>% t
colnames(R) <- c(nontour)
rownames(R) <- c(tour)


gtz <- as_tibble(which(R >= 0, arr.ind = TRUE)) %>% 
  mutate(delta = apply(., 1, function(x) {R[x[1], x[2]]}),
         prob = (delta+1)/sum(1+delta)) # +1 to allow zero delta a change


OutIn <- gtz %>% sample_n(1, weight = prob)
OutIn

aux <- tour[OutIn$row]
tour[OutIn$row] <- nontour[OutIn$col]
nontour[OutIn$col] <- aux

M2 <- M[c(tour, nontour), c(tour, nontour)]

cat("Fitness: ", sum(distances[tour, tour])/2)
#sum(M2[1:3,1:3])/2 M2[1:3, 1:3]

cat("Fitness: ", sum(distances[tour, tour])/2)


#colnames(M2) <- c(tour, nontour)
#rownames(M2) <- c(tour, nontour)
cat("Fitness: ", sum(distances[tour, tour])/2)
G2 <- M2 %*% v
R2 <- M2[1:3, 4:8] %>% 
  apply(2, function(x){-x - G[1:3]}) %>% 
  apply(1, function(x){x + G[4:8]}) %>% t


gtz <- as_tibble(which(R2 >= 0, arr.ind = TRUE)) %>% 
  mutate(delta = apply(., 1, function(x) {R2[x[1], x[2]]}),
         prob = (delta+1)/sum(1+delta)) # +1 to allow zero delta a change


OutIn <- gtz %>% sample_n(1, weight = prob)
OutIn

aux <- tour[OutIn$row]
tour[OutIn$row] <- nontour[OutIn$col]
nontour[OutIn$col] <- aux
