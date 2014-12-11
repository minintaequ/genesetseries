install.deps <- function(
    packages=c("Biobase", "affy", "limma", "graph", "GO.db",
        "AnnotationDbi", "DBI", "SparseM", "IRanges", "org.Hs.eg.db", "hgu133plus2cdf", "hgu133plus2.db"),
    lib=NULL, mirrorID=3){
    source("http://bioconductor.org/biocLite.R")
    chooseBioCmirror(ind=mirrorID)
    if(is.null(lib)) biocLite(packages) else biocLite(pkgs=packages, lib=lib)
}
