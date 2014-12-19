#' Function to get rows of name matrix matching a vector of arbitrarily IDs
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

ChipGeneNames <- function(v, N=hgu133plus2Names){
 N <- cbind(rownames(N), N)
 w <- numeric(length(v))        # vector of row numbers
 for(i in seq_along(v)){        # for all given IDs find row number in column
  a <- which(N[,1] %in% v[i]); b <- which(N[,2] %in% v[i]); c <- grep(v[i], N[,3], TRUE); d <- grep(v[i], N[,4], TRUE); e <- which(N[,5] %in% v[i])
  f <- which(N[,6] %in% v[i])
  if(length(d) > 1) d <- d[which(nchar(N[d,4]) == nchar(v[i]))]
  w[i] <- c(a,b,d,e,f,c)[1]     # set found row number in w
 }
 return(N[w,])                  # return corresponding rows
}
