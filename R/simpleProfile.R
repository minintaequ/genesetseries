simpleProfile <- function(STR.M){                                        # function to get symbolic profiles
    Prof.simple <- apply(STR.M, 1, function(x){                          # line per line change
     x <- ifelse(x==(-1), "-", ifelse(x==1, "+", ifelse(x=="2", 0, "o")))# 1,0,-1,2 -> +o-0
    return(paste(x, sep="", collapse=""))                                # combine symbols to string
   })
return(Prof.simple)}
