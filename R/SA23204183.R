#' @title A email dataset
#' @name email
#' @description A email dataset of Hillary Clinton used to calculate the ItemRank.
#' @examples
#' \dontrun{
#' data(email)
#' attach(email)
#' ir = itemrank(email, d = 0.85)
#' }
NULL

#' @import Rcpp
#' @import ggplot2
#' @import rpart
#' @import rpart.plot
#' @import labelled
#' @import knitr
#' @import bootstrap
#' @import boot
#' @import DAAG
#' @import coda
#' @import microbenchmark
#' @import igraph
#' @importFrom bootstrap bootstrap
#' @useDynLib SA23204183
NULL

#' @title PageRank algorithm using R
#' @description PageRank algorithm using R
#' @param pages Matrix of a graph, including the nodes with out-of-degree, nodes with in-of-degree and weight of edges.
#' @param d Damping factor, a positive value in (0,1), default 1.
#' @return The PageRank of the nodes in graph.
#' @examples
#' \dontrun{
#'   pagedata <- t(matrix(c(1,2,1,
#'   1,3,1,
#'   2,3,1
#'   ), nrow = 3))
#'   r1.1 = pagerank(pagedata, d = 1)
#'   r1.2 = pagerank(pagedata, d = 0.85)
#' }
#' @export
pagerank = function(pages, d = 1){
  n <- max(apply(pages[,c(1,2)],2,max))
  A <- matrix(0,n,n)
  for(i in 1 : nrow(pages)){
    if (pages[i, 3] != 1) {
      print('Error: pagerank() is a function for weight = 1.
            Please use function itemrank() for Weight != 1')
    } else {
      A[pages[i, 2], pages[i, 1]] <- 1
    }
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


#' @title ItemRank algorithm using R
#' @description ItemRank algorithm using R
#' @param pages Matrix of a graph, including the nodes with out-of-degree, nodes with in-of-degree and weight of edges.
#' @param d Damping factor, a positive value in (0,1), default 1.
#' @return The ItemRank of the nodes in graph.
#' @examples
#' \dontrun{
#'   weightdata <- t(matrix(c(1,2,0.3,
#'   1,3,0.7,
#'   2,3,1), nrow = 3))
#'   r2.1 = itemrank(weightdata, d = 1)
#'   r2.2 = itemrank(weightdata, d = 0.85)
#' }
#' @export
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








#' @title Personalized PageRank algorithm using R
#' @description Personalized PageRank algorithm using R
#' @param pages Matrix of a graph, including the nodes with out-of-degree, nodes with in-of-degree and weight of edges.
#' @param d Damping factor, a positive value in (0,1), default 1.
#' @param root The Personalized node with out-of-degree.
#' @param is_using_out Whether to use average out-of-degree as weight, default FALSE as not use.
#' @return The Personalized PageRank of the root node to nodes with in-of-degree in graph.
#' @examples
#' \dontrun{
#'   persondata <- t(matrix(c(1,5,1,
#'   1,6,0.1,
#'   2,5,0.2,
#'   3,6,0.4,
#'   4,5,0.7), nrow = 3))
#'   r3.1 = personpagerank(persondata, d = 0.85, root = 1, is_using_out=TRUE)
#'   r3.2 = personpagerank(persondata, d = 0.85, root = 1, is_using_out=FALSE)
#' }
#' @export

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