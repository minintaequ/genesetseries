## Filter genes according to their EntrezGene IDs, using Affymetrix-suffix-information for clusters with the same EntrezGene ID and
## the max-median-approach to choose the single probe set with the clearest signal in the study
#
FilterGenes <- function(Names=hgu133plus2Names, clusters="EntrezGene.IDs", EM=PVAL.do$pvalues){
  NAMES <- Names[,clusters]
  FREQ  <- table(NAMES)
  g.NA  <- which(is.na(NAMES))                                     # features without grouping ID
  g.UN  <- which(NAMES %in% names(FREQ[which(FREQ==1)]))           # features with unique grouping ID
  f.MT  <- names(FREQ[which(FREQ > 1)])                            # grouping IDs with multiple features (2 to 11)
  L.MT  <- lapply(f.MT, function(x) return(names(NAMES[which(NAMES==x)])))
  F.MT  <- lapply(L.MT, function(f){                               # filtered list of representative features for EntrezGene ID
    suffixes <- sapply(strsplit(f, "_"), function(x) ifelse(length(x) > 2, x[2], NA))
    if(sum(is.na(suffixes)) > 0){ f <- f[is.na(suffixes)]} else{   # prefer XXX_at before
      if("a" %in% suffixes){ f <- f[which(suffixes == "a")]} else{ #        XXX_a_at before
      if("s" %in% suffixes){ f <- f[which(suffixes == "s")]} else{ #        XXX_s_at and before
        f <- f}}                                                   #        XXX_x_at!
    }
    if(length(f) == 1){ s <- f; return(s)}                         # give single gene back
    EM.f <- EM[f,]                                                 # expression of features in time
    s <- f[which.min(apply(EM.f, 1, min))]                         # minimum p.value at any time point
    #          apply(matrix(e, ncol=6), 2, median)}),2, max))]      # clearest signal = highest median expression value without reference
    return(s)                                                      # set or single representant for entrez gene / interesting gene
  })                                                               # end of groupingID cluster lapply
  l.un <- names(NAMES[g.UN]); names(l.un) <- NAMES[g.UN]
  L <- as.list(l.un)
  names(F.MT) <- f.MT
  L <- c(L, F.MT)                                                  # list named according grouping IDs and containing feature IDs
  ##                                                               # comments
  cat(paste("There are", length(NAMES), "features available.\n"))
  cat(paste("Thera are no", clusters, "available for", length(g.NA), "features.\n"))
  cat(paste("Thera are", length(g.UN), clusters, "available with unique chip features.\n"))
  cat(paste("Thera are", length(f.MT), clusters, "with multiple features.\n"))
  cat(paste("Thera are", sum(sapply(L, length)), "features in use for", length(L), clusters, "in whole.\n"))
  return(L)
}
