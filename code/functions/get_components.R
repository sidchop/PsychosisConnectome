get_components <- function(adj=NULL, return_comp=F, return_max=F) {
  library(sna)
  if (nrow(adj)!=ncol(adj)) {
    error('this adjacency matrix is not square')
  }
  
  diag(adj) <- 0
  
  #Dulmage-Mendelsohn decomposition ## not implemented in R
  # so using sna to docompose
  c <- component.largest(adj, return.as.edgelist = F, result = c("graph"))
  if(return_comp==F & return_max==T) {
    return(sum(c[upper.tri(c)]))
    stop()
  }
  compsize <- component.size.byvertex(adj)
  row.names(c) <- colnames(c) <-  which(compsize == max(compsize))
  edge_list <- reshape2::melt(c)
  edge_list <- edge_list[!edge_list$value==0,]
  out_comp <- matrix(0,nrow(adj), ncol(adj))
  
  out_comp[as.matrix(edge_list[,1:2])] <- 1
  if(return_comp==T & return_max==F) {
    return(out_comp)
  }
}
