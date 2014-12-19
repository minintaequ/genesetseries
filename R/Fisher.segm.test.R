#' A one sentence description of what your function does
#'
#' A more detailed description of what the function is and how
#' it works. It may be a paragraph that should not be separated
#' by any spaces.
#'
#' @param inputParameter1 A description of the input parameter \code{inputParameter1}
#' @param inputParameter2 A description of the input parameter \code{inputParameter2}
#'
#' @return output A description of the object the function outputs
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' R code here showing how your function works

Fisher.segm.test <- function(value.matrix,                               # matrix with values which define diff and not diff
                             groups,                                     # sets to be tested for enrichment
                             quantiles.up,                               # one or more quantiles to define up-regulation
                             quantiles.do)                               # one or more quantiles to define down-regulation
  {
    P.up <- NULL                                                         # create matrix of all p-values (G x quantiles)
    Time.list.up <-list()
    if(length(quantiles.up) > 0){
      for(t in seq(length=ncol(value.matrix))){ Time.list.up[[LETTERS[t]]] <- NULL}
      P <- sapply(quantiles.up, function(q){
                qs <- apply(value.matrix, 2, quantile, q)                # column-quantiles for q
                dec.matrix <- t(t(value.matrix) >= qs)                   # logical if more extreme than column quantile q
                n.int <- apply(dec.matrix, 2, sum)                       # number of interesting features in columns
                out <- sapply(groups, function(g){                       # calculate fisher.test p-value and tarone minimal p-values
                  int.in.group <- colSums(dec.matrix[g,])                # significantly up regulated group genes
                  p   <- phyper(int.in.group-1, n.int, nrow(dec.matrix)-n.int, length(g), lower.tail=FALSE) # p-value for up-regulation
                  return(p)})                                            # p-values return
                for(t in seq(length=ncol(dec.matrix))){                  # fill list for each time point with p.values corresponding to q
                  Time.list.up[[LETTERS[t]]] <<- cbind(Time.list.up[[LETTERS[t]]], out[t,])
                }})
                                        # matrix of p-values
      if(length(colnames(value.matrix)) > 0) names(Time.list.up) <- colnames(value.matrix)
      for(t in seq(length=ncol(value.matrix))) P.up <- cbind(P.up, Time.list.up[[t]])
      if(length(colnames(value.matrix)) > 0) colnames(P.up) <- rep(colnames(value.matrix), each=length(quantiles.up))
    }
    P.do <- NULL                                                         # create matrix of all p-values (G x quantiles)
    Time.list.do <-list()
    if(length(quantiles.do) > 0){
      for(t in seq(length=ncol(value.matrix))){ Time.list.do[[LETTERS[t]]] <- NULL}
      P <- sapply(quantiles.do, function(q){
                qs <- apply(value.matrix, 2, quantile, q)                # column-quantiles for q
                dec.matrix <- t(t(value.matrix) <= qs)                   # logical if more extreme than column quantile q
                n.int <- apply(dec.matrix, 2, sum)                       # number of interesting features in columns
                out <- sapply(groups, function(g){                       # calculate fisher.test p-value and tarone minimal p-values
                  int.in.group <- colSums(dec.matrix[g,])                # significantly up regulated group genes
                  p   <- phyper(int.in.group-1, n.int, nrow(dec.matrix)-n.int, length(g), lower.tail=FALSE) # p-value for up-regulation
                  return(p)})                                            # p-values return
                for(t in seq(length=ncol(dec.matrix))){                  # fill list for each time point with p.values corresponding to q
                  Time.list.do[[LETTERS[t]]] <<- cbind(Time.list.do[[LETTERS[t]]], out[t,])
                }})
                                        # matrix of p-values
      if(length(colnames(value.matrix)) > 0) names(Time.list.do) <- colnames(value.matrix)
      for(t in seq(length=ncol(value.matrix))) P.do <- cbind(P.do, Time.list.do[[t]])
      if(length(colnames(value.matrix)) > 0) colnames(P.do) <- rep(colnames(value.matrix), each=length(quantiles.do))
    }
    return(list(up=P.up, down=P.do))                                     # return the p-value matrix with groups in rows
  }
