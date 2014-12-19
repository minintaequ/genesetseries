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

## pval.adjust.FUN ... function with arguments:
##   pval.matrix   ... matrix of p.values in columns
##   method        ... adjustment method to keep significance niveau (Bonferroni, Holm, Hommel, Tarone) or False Discovery Rate (Hochberg)
##   signed        ... p-values signed (sign represents direction in profile up/down)
##   Neff          ... Neff-object from function effNumber() contains groups of same tests and number of effective tests
##   Tarone.matrix ... matrix of lowest possible significance niveau in current test (same dimension like pval.matrix)
##   sep.col       ... adjustment for columns separately, e.g. niveau holds only for columns separately
##   sep.sign      ... adjustment for each direction separately (sign splits into separately adjusted groups)
#
pval.adjust.FUN <- function(pval.matrix, method=c("bonferroni", "holm", "hommel", "tarone", "BH", "fdr", "TS-ABH")[1],
                            signed=FALSE, Neff=NULL){
  if(!require(multtest)){stop("R-package multtest seems to be not available.")}
  if(!is.numeric(pval.matrix)){stop("pval.matrix is supposed to be numeric")}
  if(min(pval.matrix) < -1 | max(pval.matrix > 1)){stop("pval.matrix seems not to include (signed) p-values")}
  if(!is.matrix(pval.matrix)){warning("pval.matrix was expected to be a matrix. Is transfered to a matrix with column dimension of 1")
                              rn <- names(pval.matrix)
                              pval.matrix <- matrix(pval.matrix)
                              rownames(pval.matrix) <- rn
                            }
  pval.matrix.adj <- matrix(NA,ncol=ncol(pval.matrix), nrow=nrow(pval.matrix), dimnames=dimnames(pval.matrix))
  if(!is.null(Neff)){                                               # object which determines the effective number of tests
    names(Neff[[2]]) <- paste("Neff_cl_", seq(length(Neff[[2]])), sep="")
    pval.names <- rownames(pval.matrix)
    cl.repr <- t(sapply(Neff[[2]], function(x){                     # cluster representants of multiple groups with identical genes
      y <- pval.matrix[x[1],]
      w <- which(rownames(pval.matrix) %in% x)
      pval.matrix <<- as.matrix(pval.matrix[-w,])
      pval.names <<- pval.names[-w]
      return(y)
    }))
    if(dim(cl.repr)[1] == 1) cl.repr <- t(cl.repr)
    rownames(cl.repr) <- names(Neff[[2]])
    pval.matrix <- rbind(pval.matrix, cl.repr)                      # pval.matrix includes representants of each cluster from Neff[[2]]
    rownames(pval.matrix) <- c(pval.names, rownames(cl.repr))
  }
  if(signed) s.matrix <- sign(pval.matrix)                          # keep sign
  if(method=="TS-ABH"){                                             # own implementation of FDR q-values two-stage step up BH-procedure
    P <- fdr.TST.ABH(abs(pval.matrix))
  } else {
    P <- p.adjust(abs(pval.matrix), method[1])}                     # p.value adjustment
  if(signed) P <- P*s.matrix                                        # restore sign
  pval.adj.matrix <- matrix(P, ncol=ncol(pval.matrix), dimnames=dimnames(pval.matrix))
  if(!is.null(Neff)){                                         # reconstruct groups with identical genes from Neff list
    equal.names <- rownames(pval.matrix.adj)[which(rownames(pval.matrix.adj) %in% rownames(pval.adj.matrix))]
    pval.matrix.adj[equal.names,] <- pval.adj.matrix[equal.names,,drop=FALSE]
    for(n in names(Neff[[2]])){
      fillIn <- t(replicate(length(Neff[[2]][[n]]), pval.adj.matrix[n,]))
      if(dim(fillIn)[1] == 1) fillIn <- t(fillIn)
      pval.matrix.adj[Neff[[2]][[n]],] <- fillIn
    }} else {
      pval.matrix.adj <- pval.adj.matrix
    }
  colnames(pval.matrix.adj) <- colnames(pval.matrix)
  return(pval.matrix.adj)
}
