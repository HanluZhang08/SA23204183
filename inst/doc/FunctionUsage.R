## -----------------------------------------------------------------------------
library(SA23204183)
data(email)
head(email)

## ---- echo=FALSE--------------------------------------------------------------
pages <- data.frame(
  from = c(1, 1, 2, 3),
  to = c(2, 3, 3, 1),
  weight = c(1, 1, 1, 1)
)

n <- max(apply(pages,2,max))
A <- matrix(0,n,n)
for(i in 1 : nrow(pages)){
    A[pages[i, 2], pages[i, 1]] <- 1
}
adjM = A

## -----------------------------------------------------------------------------
library(microbenchmark)
promR = function(adjM, d){
  cs <- colSums(adjM)
  cs[cs==0] <- 1
  n <- nrow(adjM)
  delta <- (1-d)/n
  P <- matrix(NA, nrow(adjM), ncol(adjM))
  for (i in 1:n){
    P[i,] <- delta + d * adjM[i,]/cs
  }
  return(P)
}

tm1 <- microbenchmark(
  promC = PromC(adjM, d = 0.85),
  promR = promR(adjM, d = 0.85)
)
knitr::kable(summary(tm1)[,c(1,3,5,6)])

## ---- eval=FALSE--------------------------------------------------------------
#  pagerank = function(pages, d = 1){
#    n <- max(apply(pages[,c(1,2)],2,max))
#    A <- matrix(0,n,n)
#    for(i in 1 : nrow(pages)){
#      if (pages[i, 3] != 1) {
#        print('Error: pagerank() is a function for weight = 1.
#              Please use function itemrank() for Weight != 1')
#      } else {
#        A[pages[i, 2], pages[i, 1]] <- 1
#      }
#    }
#    adjM = A
#  
#    proM = PromC(adjM, d)
#  
#    eigenMatrix <- function(proM){
#      n <- nrow(proM)
#      x <- Re(eigen(proM)$vectors[,1])
#      return(x/sum(x))
#    }
#    rank = eigenMatrix(proM)
#  
#    return(rank)
#  }

## -----------------------------------------------------------------------------
library(igraph)
# prepare the data
edges <- data.frame(
  from = c(1, 1, 2, 3),
  to = c(2, 3, 3, 1),
  weight = c(1, 1, 1, 1)
)
graph <- graph_from_data_frame(edges, directed = TRUE, vertices = NULL)
E(graph)$weight <- edges$weight
# calculate
r1.0 <- page_rank(graph, weights = NA)$vector
r1.1 = pagerank(edges, d = 0.85)
print(r1.0)
print(r1.1)
# compare the effiency
tm2 <- microbenchmark(
  r1.0 = page_rank(graph, weights = NA),
  r1.1 = pagerank(edges, d = 0.85)
)
knitr::kable(summary(tm2)[,c(1,3,5,6)])

## -----------------------------------------------------------------------------
itemrank = function(pages, d = 1){
  n <- max(apply(pages[,c(1,2)],2,max))
  A <- matrix(0,n,n)
  for(i in 1 : nrow(pages)){
    A[pages[i, 2], pages[i, 1]] <- pages[i, 3]
  }
  if(all(pages[, 3]==1)){
    print('Note: All weight = 1, but itemrank() is a function for Weight != 1.
          A better choice is function pagerank().')
  }
  adjM = A
  
  proM = PromC(adjM, d)
  
  eigenMatrix <- function(proM){
    n <- nrow(proM)
    x <- Re(eigen(proM)$vectors[,1])
    return(x/sum(x))
  }
  rank = eigenMatrix(proM)
  
  return(rank)
}

## -----------------------------------------------------------------------------
ir = itemrank(email[,c(3,4,5)], d = 0.85)
head(ir)

## -----------------------------------------------------------------------------
personpagerank <- function(pages, d = 1, root, is_using_out = FALSE) {
  n <- max(apply(pages[,c(1,2)],2,max))
  nodes = c(1:n)
  A <- matrix(0,n,n)
  for(i in 1 : nrow(pages)){
    A[pages[i, 2], pages[i, 1]] <- pages[i, 3]
  }
  adjM = A
  if(is_using_out){
    degree <- function(adjM, nodes) {
      degrees = numeric(max(pages))
      for (i in unique(pages[,1])) {
        degrees[i] <- sum(adjM[ ,i])
      }
      for (i in unique(pages[,2])) {
        degrees[i] <- 10e6
      }
      return(degrees)
    }
    degrees <- 1/degree(adjM, nodes)
    adjM <- t(t(adjM) * degrees) 
  }
  matrix <- adjM
  
  r0 <- matrix(ifelse(nodes == root, 1, 0), nrow = length(nodes), ncol = 1)
  n <- nrow(matrix)
  A <- diag(n) - d * matrix
  b <- (1 - d) * r0
  r <- solve(A, b)
  rank <- vector()
  for (j in 1:n) {
    rank[nodes[j]] <- r[j, 1]
  }
  return(rank)
}

## -----------------------------------------------------------------------------
ppr = personpagerank(email[,c(3,4,5)], d = 0.85, root = 1, is_using_out = FALSE)
head(ppr)

