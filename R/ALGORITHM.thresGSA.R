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

ALGORITHM.thresGSA <- function(                                          # function for filtering diff. expressed genes from p-values matrix
   sets     = GROUPSod$groups,                                           # sets considered for thresGSA profile estimation
   genes    = GENESod,                                                   # genes considered in thresGSA profile estimation
   sigBord.genes=0.05,                                                   # threshold for considering genes as significantly differential expr
   adj.genes="unadjusted",                                               # p-value adjustment for gene-wise t-tests (function)
   sigBord.sets=0.01,                                                    # threshold for considering sets as significantly enriched
   adj.sets="unadjusted",                                                # p-value adjustment for gene-wise Fisher-tests (function)
   Neff=NULL,                                                            # objects including the effective number of distinguishable sets
   groupTest=Fisher.group.test,
   Pmatrix=TTod.np.all)
{
 start <- proc.time()[3]
 PM <- matrix(NA, ncol=ncol(Pmatrix$up.p.values)+ncol(Pmatrix$do.p.values), nrow=nrow(Pmatrix$do.p.values))
 PM[,seq(1, ncol(Pmatrix$up.p.values)*2, 2)] <- Pmatrix$up.p.values
 PM[,seq(2, ncol(Pmatrix$do.p.values)*2, 2)] <- Pmatrix$do.p.values
 rownames(PM) <- rownames(Pmatrix$up.p.values)
 colnames(PM) <- paste(rep(colnames(Pmatrix$do.p.values), each=2), rep(c("up", "do"), ncol(Pmatrix$do.p.values)), sep="-")
 Pmatrix <- PM[genes,]                                                   # use only annotated genes
 if(adj.genes %in% c("bonferroni", "holm", "hommel", "tarone", "BH", "fdr", "TS-ABH")){
   Pmatrix <- pval.adjust.FUN(Pmatrix, method=adj.genes, Neff=NULL)
 } else{warning("There is no p-value adjustment for gene to reference comparisons done!")}
 SIGgenes  <- apply(Pmatrix, 2, function(x) rownames(Pmatrix)[which(x <= sigBord.genes)])
 GO.pvalue <- groupTest(Pmatrix, sets, sigBord.genes)
 if(adj.sets %in% c("bonferroni", "holm", "hommel", "tarone", "BH", "fdr", "TS-ABH")){
   GO.pvalue <- pval.adjust.FUN(GO.pvalue, method=adj.sets, Neff=Neff)   # adjustment for multiple testing according to adj.sets function
 } else {warning("There is no p-value adjustment for group enrichment done!")}
 colnames(GO.pvalue) <- colnames(PM)
 SIGpvalGR <- GO.pvalue[apply(GO.pvalue, 1, function(x) any(x < sigBord.sets)),] # keep only sets, which are significant at least for one comparison
 SIGprofiles <- t(apply(SIGpvalGR, 1, function(p){                       # read profiles from p-values
   P <- matrix(p, nrow=2)                                                # create matrix with time points in columns
   A <- apply(P, 2, function(x){ sig <- (x < sigBord.sets)               # profile code up=1, down=-1, both=2, nothing=0
                    return(ifelse(!any(sig), 0, ifelse(all(sig), 2, sum(c(1, -1)*sig))))})
   return(A)}))
 SIGprofilesSimple <- simpleProfile(SIGprofiles)                         # create symbolic profile ++oo-- from 1100-1-1
 cat(paste("Algorithm thresGSA found", length(SIGprofilesSimple), "non-zero activation profiles in a run time of",
           round((proc.time()[3]-start)/60), "minutes.\n"))
                                                                         # give out summary of gene group enrichment
 return(list(PVAL=SIGpvalGR, GOprofiles=SIGprofiles, feasGenes=genes, ListSIGgenes=SIGgenes, n.SIGgenes=sapply(SIGgenes, length),
             GenesPerGroup=sets, groupsizes=sapply(sets, length), UsedGr=names(sets), p.val.adjust.sets.method=adj.sets,
             Alpha.gr=sigBord.sets, Alpha.ge=sigBord.genes, simpleProfile=SIGprofilesSimple, allPVAL=GO.pvalue,
             GenePVAL=Pmatrix, p.val.adjust.genes.method=adj.genes, p.val.adjust.sets.method=adj.sets))
}
