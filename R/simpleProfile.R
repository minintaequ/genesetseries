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

simpleProfile <- function(STR.M){                                        # function to get symbolic profiles
    Prof.simple <- apply(STR.M, 1, function(x){                          # line per line change
     x <- ifelse(x==(-1), "-", ifelse(x==1, "+", ifelse(x=="2", 0, "o")))# 1,0,-1,2 -> +o-0
    return(paste(x, sep="", collapse=""))                                # combine symbols to string
   })
return(Prof.simple)}
