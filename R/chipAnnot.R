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