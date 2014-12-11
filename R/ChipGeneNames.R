## function to get rows of name matrix matching a vector of arbitrarily IDs
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