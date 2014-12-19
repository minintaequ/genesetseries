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

chipAnnot <- function(chip="hgu133plus2", geneIDs=rownames(PVAL.MeHg$up.q.values)){
    T <- toTable(get(paste(chip, "ENTREZID", sep="")))
    EntrezGene.IDs <- T[,2]; names(EntrezGene.IDs)=T[,1]
    T <- toTable(get(paste(chip, "GENENAME", sep="")))
    GeneNames  <- T[,2]; names(GeneNames)=T[,1]
    T <- toTable(get(paste(chip, "SYMBOL", sep="")))
    GeneSymbols <- T[,2]; names(GeneSymbols)=T[,1]
    T <- toTable(get(paste(chip, "UNIGENE", sep="")))
    UniGene.IDs <- T[,2]; names(UniGene.IDs)=T[,1]
    T <- toTable(get(paste(chip, "ACCNUM", sep="")))
    GenBank.IDs <- T[,2]; names(GenBank.IDs)=T[,1]
    makeGeneAnno <- function(set=c("EntrezGene.IDs", "GeneNames", "GeneSymbols", "UniGene.IDs", "GenBank.IDs"), Names=geneIDs){
        GeneIdentifier <- matrix(NA, ncol=length(set), nrow=length(Names))
        colnames(GeneIdentifier) <- set; rownames(GeneIdentifier) <- Names
        for(id in set){
            identifiers <- get(id)
            GeneIdentifier[names(identifiers), id] <- identifiers
        }
        return(GeneIdentifier)
    }
    Names <- makeGeneAnno()
}
