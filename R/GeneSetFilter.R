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

GeneSetFilter <- function(geneSets=ANNSETSdo, minSetSize=4, Rank=FALSE, features=unlist(EntrezGeneGroups.do), exclude="GO:0008150"){
  start <- proc.time()
  group.genes <- list()                                              # list of all gene sets
  if(is.list(geneSets[[1]])){
    for(g in names(geneSets)){                                       # for all sublists of a list of set-list-objects
      print(g)
      for(n in names(geneSets[[g]])){
        group.genes[[n]] <- geneSets[[g]][[n]]}                      # concatinate the lists into group.genes-list
        }} else {group.genes <- geneSets}
  geneSets <- group.genes                                            # rename the new list
  if(Rank){                                                          # Rank as number of different tests
     genes <- unique(unlist(geneSets))
     M_gene_group <- matrix(0, ncol=length(geneSets), nrow=length(genes))
     i <- 1                                                          # counter
     sapply(genes, function(g){                                      # binary matrix of genes (rows) in group (cols)
       a <- sapply(geneSets, function(GO){                           # for every group
          return(as.numeric(g %in% GO))})                            #  return vector of gene annotation
       M_gene_group[i,] <<- a                                        # update M_gene_group
       i <<- i+1         })
     Neff <- rank(M_gene_group)                                      # rank of binary matrix as effective number of tests
   } else{                                                           # only exclude identical sets from the FDR-calculation
     length.groups <- sapply(geneSets, length)                       # total of group genes
     Table.l.gr <- table(length.groups)                              # table of frequency of occurence
     if(minSetSize > 0){                                                # reduce the number of gene groups to those
       Table.l.gr <- Table.l.gr[-which(names(Table.l.gr) %in% seq_len(minSetSize-1))] # which fullfill a minSetSize of gene number
       geneSets <- geneSets[which(sapply(geneSets, length) >= minSetSize)]}
     cat(paste(sum(Table.l.gr), "GO groups fullfill the minsize requirement of", minSetSize, "genes.\n"))
     assign("Neff", sum(Table.l.gr == 1), envir=.GlobalEnv)          # number of different groups/tests
     assign("Index", NULL, envir=.GlobalEnv)                         # index of groups identical to previous groups
     assign("mult.GR", list(), envir=.GlobalEnv)                     # list of multiple groups
     for(j in names(Table.l.gr[which(Table.l.gr > 1)])){
       w <- names(length.groups)[which(length.groups==as.numeric(j))]   # those groups with total genes = j
       index <- rep(TRUE, length(w)); names(index) <- w
       test <- NULL
       for(i in w){
         test <- c(test,i)
         if(index[i]){
           inGroup.i <- geneSets[[i]]                                  # in current group
           identical <- sapply(geneSets[w[index]], function(x){        # which other groups with same gene number have all these genes?
             return(all(x %in% inGroup.i))
           })
           if(sum(identical)>1){
             mult.GR[[length(mult.GR)+1]] <- names(geneSets[w[index]])[identical] # update list of multiple groups
           }
           Neff <<- Neff + 1                                            # update number of effective different groups
           index[names(geneSets[w[index]])[identical]] <- FALSE        # set all identical columns to Index
         }
       }
     cat(paste(length(index), "groups with total genes", j, "had processed.\n"))
     }
   }
  duration <- (proc.time()-start)[3]
  cat("Neff calculation takes", round(duration/60), "minutes and", round(duration%%60), "seconds.",
      "Effective group number is Neff =", Neff, "and multiple groups split into", length(mult.GR), "clusters.\n")
  if(length(exclude) > 0){
    ALL.genes.without.excluded <- unique(unlist(geneSets[-which(names(geneSets)==exclude)]))
    cat(paste(length(features) - length(ALL.genes.without.excluded), "genes excluded, because of excluded sets (GO BP).\n"))
    for(e in exclude){
      geneSets[[e]] <- intersect(geneSets[[e]], ALL.genes.without.excluded)
    }
    features <- ALL.genes.without.excluded
  }
  time.excl <- (proc.time()-start)[3] - duration
  cat("Feature extraction takes", round(time.excl), "seconds and number of actually used genes is", length(features))
  return(list(effectiveNumber=Neff, idGroupList=mult.GR, groups=geneSets, genes=features))
}
