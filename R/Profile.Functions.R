
#
## function to get rows of name matrix matching a vector of arbitrarily IDs
hgu133plus2names <- function(v, N=hgu133plus2Names){
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

#
## Filter genes according to their EntrezGene IDs, using Affymetrix-suffix-information for clusters with the same EntrezGene ID and
## the max-median-approach to choose the single probe set with the clearest signal in the study
#

## annotate genes to GO groups of a predefined ontology
# GO.anno(T)                               # list of gene group lists per ontology

KEGG.anno <- function(table, features){
    KEGG.sets <- NULL
    KEGG <- lapply(unique(table[,2]), function(set){
        ID <- paste("KEGG:", set, sep="")
        KEGG.sets[[ID]] <<- intersect(table[which(table[,2] == set),1], features)
    })
    return(KEGG.sets)
}


## PU <- pval.adjust.FUN(TTod.np.all$up.p.values, method="TS-ABH", signed=FALSE)
## apply(PU, 2, function(x) sum(x <= 0.05))
## 2643 2668 2932 2146 2760 2922
## PD <- pval.adjust.FUN(TTod.np.all$do.p.values, method="TS-ABH", signed=FALSE)
## apply(PD, 2, function(x) sum(x <= 0.05))
## 3320 3830 4079 3308 3683 3525

#
####################################################################################################
###                                                                                              ###
#     5. Estimating Gene Group Activation Profiles                                                 #
#                                                                                                  #
#       - threshold-GSA type with two Fisher-tests at every time point                             #
#       - non-threshold-GSA type with segmented Fisher tests at every time point                   #
#       - rotation-GSEA type with GSEA test at every time point with rotation testing              #
#       - STEM-algorithm clustering to obtain gene expression trajectory clusters                  #
#       - maSigFun-algorithm as linear model of gene expression values differences to reference    #
###                                                                                              ###
####################################################################################################

#
####################################################################################################
#       - threshold-GSA type with two Fisher-tests at every time point                             #
####################################################################################################
#

#
####################################################################################################
#       - non-threshold-GSA type with two Fisher-tests at every time point (1S-GSA)                #
####################################################################################################
#

#


####################################################################################################
#       - two threshold segmentation profile 2S-GSA                                                #
####################################################################################################
#

