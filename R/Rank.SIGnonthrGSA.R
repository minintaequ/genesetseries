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

Rank.SIGnonthrGSA <- function(PROF=PROFILESdo.SIGnonthrGSA){
  q <- PROF$Quantiles
  P <- PROF$PVAL
  T <- PROF$GeneTmat
  A <- PROF$GOprofiles
  R <- t(apply(cbind(rownames(A), A), 1, function(prof){
    gr.name <- prof[1]
    lsig <- PROF$ListSIGgenes[[gr.name]]
    test <- T[rownames(lsig),]
    score <- sum(abs(lsig*test))/sum(lsig != 0)
    ngenes <- nrow(lsig)
    sgenes <- sum(apply(lsig, 1, function(x) any(x != 0)))
    profile <- PROF$simpleProfile[gr.name]
    term <- ifelse(grepl("GO:", gr.name), Term(gr.name), "KEGG pathway tba")
    out <- c(score, gr.name, profile, ngenes, sgenes, term)
    return(out)
  }))
  R <- R[order(as.numeric(R[,1]), decreasing=TRUE),]
  return(R)
}
