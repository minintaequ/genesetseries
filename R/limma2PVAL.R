limma2PVAL <- function(data, limma, pvadj="BH", ref=NULL){                 ## which column in design-matrix is reference
    PV <- matrix(p.adjust(limma$p.value, method="BH"), ncol=ncol(limma$t)) ## adjust p-value for multiple testing
    dimnames(PV) <- dimnames(limma$p.value)                                ## restore names after p-value adjustment
    dir <- ifelse(limma$t < 0, -1, 1)                                      ## direction with respect to controll
    PVu <- ifelse(limma$t > 0, PV, 0.999999)                               ## (adjusted) p-values for up-regulation
    PVd <- ifelse(limma$t < 0, PV, 0.999999)                               ## (adjusted) p-values for do-regulation
    calc.logFC <- function(E=data, des=limma$design, w=ref){
        if(!is.null(ref)){
            M = E[,des[,ref]]
            notref <- seq(ncol(des))[-ref]} else{
                M = matrix(0, ncol=1, nrow=nrow(E))
                notref <- seq(ncol(des))
            }
        O = NULL
        for(p in notref){
            O <- cbind(O, rowMeans(E[,des[,p]])-rowMeans(M))
        }
        rownames(O) <- rownames(E)
        colnames(O) <- colnames(des)[notref]
        return(O)
    }
    FCdiff <- calc.logFC()                                                 ## (log2) Fold-Changes for the difference to reference
    limma$statistic <- limma$t
    limma$q.value <- PV
    limma$up.p.values <- ifelse(limma$t > 0, limma$p.value, 0.99999)
    limma$do.p.values <- ifelse(limma$t < 0, limma$p.value, 0.99999)
    limma$up.q.values <- PVu
    limma$do.q.values <- PVd
    limma$logFC <- FCdiff
    limma$dir <- dir
    limma$adjust.method <- pvadj
    return(limma)                                                          ## return a limma-object with new elements for the gene set activation profiles
}