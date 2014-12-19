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

Fisher.group.test <- function(value.matrix,                              # matrix with values which define diff and not diff
                              groups,                                    # sets to be tested for enrichment
                              threshold)                                 # one or more thresholds to define diff and not diff
  {
    dec.matrix <- value.matrix <= threshold                              # only smaller than threshold is interesting
    n.int <- apply(dec.matrix, 2, sum)
    Out <- sapply(groups, function(g){                                   # calculate fisher.test p-value and tarone minimal p-values
      int.in.group <- colSums(dec.matrix[g,])                            # significantly up regulated group genes
      p   <- phyper(int.in.group-1, n.int, nrow(dec.matrix)-n.int, length(g), lower.tail=FALSE) # p-value for up-regulation
      return(p)
    })
    return(t(Out))                             # return the p-value matrix with groups in rows
  }
