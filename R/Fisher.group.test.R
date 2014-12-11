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
