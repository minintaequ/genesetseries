annotate.KEGG <- function(table, features){
    KEGG.sets <- NULL
    KEGG <- lapply(unique(table[,2]), function(set){
        ID <- paste("KEGG:", set, sep="")
        KEGG.sets[[ID]] <<- intersect(table[which(table[,2] == set),1], features)
    })
    return(KEGG.sets)
}
