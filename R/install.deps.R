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

install.deps <- function(
    packages=c("Biobase", "affy", "limma", "graph", "GO.db",
        "AnnotationDbi", "DBI", "SparseM", "IRanges", "org.Hs.eg.db", "hgu133plus2cdf", "hgu133plus2.db"),
    lib=NULL, mirrorID=3){
    source("http://bioconductor.org/biocLite.R")
    chooseBioCmirror(ind=mirrorID)
    if(is.null(lib)) biocLite(packages) else biocLite(pkgs=packages, lib=lib)
}
