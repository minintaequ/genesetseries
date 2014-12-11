ALGORITHM.nonthrGSA <- function(                                         # function for filtering diff. expressed genes from p-values matrix
   sets     = GROUPSod$groups,                                           # sets considered for thresGSA profile estimation
   genes    = GENESod,                                                   # genes considered in thresGSA profile estimation
   Q.up = seq(0.75, 0.95, 0.05),                                         # quantiles defining the interesting genes in the segment test
   Q.do = rev(seq(0.05, 0.25, 0.05)),                                    # quantiles defining the interesting genes in the segment test
   sigBord.sets=0.01,                                                    # threshold for considering sets as significantly enriched
   adj.sets="unadjusted",                                                # p-value adjustment for gene-wise Fisher-tests (function)
   Neff=NULL,                                                            # objects including the effective number of distinguishable sets
   groupTest=Fisher.segm.test,
   Tmatrix=TTod.np.all$statistic)
{
 start <- proc.time()[3]
 Tmatrix <- Tmatrix[genes,]                                              # use only annotated genes
 SET.pvalues <- groupTest(Tmatrix, sets, Q.up, Q.do)
 time.segm <- proc.time()[3]
 cat(paste("Segment test calculation lasts ", (time.segm - start), "seconds.\n"))
 if(adj.sets %in% c("bonferroni", "holm", "hommel", "tarone", "BH", "fdr", "TS-ABH")){
    all.SET.pvalues <- pval.adjust.FUN(cbind(SET.pvalues$up, SET.pvalues$down), method=adj.sets, Neff=Neff)
    # adjustment for multiple testing according to adj.sets function
    time.adj <- proc.time()[3]
    cat(paste("Multiple testing adjustment lasts ", (time.adj - time.segm), "seconds.\n"))
 } else {warning("There is no p-value adjustment for group enrichment done!")}
 SIGquant <- NULL
 Pmatrix  <- t(apply(all.SET.pvalues, 1, function(q){
   M.up <- matrix(q[1:(length(Q.up)*ncol(Tmatrix))], ncol=ncol(Tmatrix))
   M.do <- matrix(q[(length(Q.up)*ncol(Tmatrix)+1):length(q)], ncol=ncol(Tmatrix))
   m.up <- apply(M.up, 2, min)
   m.do <- apply(M.do, 2, min)
   w.up <- apply(M.up, 2, function(x) which(x == min(x))[1])
   w.do <- apply(M.do, 2, function(x) which(x == min(x))[1])
   sigq <- numeric(length(m.up)+length(m.do))
   sigq[seq(1, length(sigq), 2)] <- w.up
   sigq[seq(2, length(sigq), 2)] <- w.do
   SIGquant <<- rbind(SIGquant, sigq)
   out  <- numeric(length(m.up)+length(m.do))
   out[seq(1, length(out), 2)] <- m.up
   out[seq(2, length(out), 2)] <- m.do
   return(out)
 }))
 rownames(SIGquant) <- rownames(Pmatrix)
 colnames(Pmatrix) <- colnames(SIGquant) <- paste(rep(colnames(Tmatrix), each=2), c("up", "down"), sep="-")
 SIGpvalGR <- Pmatrix[apply(Pmatrix, 1, function(x) any(x < sigBord.sets)),] # keep only sets, which are significant at least for one comparison
 SIGprofiles <- t(apply(SIGpvalGR, 1, function(p){                       # read profiles from p-values
   P <- matrix(p, nrow=2)                                                # create matrix with time points in columns
   A <- apply(P, 2, function(x){ sig <- (x < sigBord.sets)               # profile code up=1, down=-1, both=2, nothing=0
                    return(ifelse(!any(sig), 0, ifelse(all(sig), 2, sum(c(1, -1)*sig))))})
   return(A)}))
 SIGprofilesSimple <- simpleProfile(SIGprofiles)                         # create symbolic profile ++oo-- from 1100-1-1
 Q.numb.up <- ceiling(nrow(Tmatrix) - quantile(1:nrow(Tmatrix), Q.up))
 Q.numb.do <- floor(quantile(1:nrow(Tmatrix), Q.do))
 names(Q.numb.up) <- Q.up
 names(Q.numb.do) <- Q.do
 SIGquantsSet <- t(apply(cbind(SIGpvalGR, SIGquant[rownames(SIGpvalGR),]), 1, function(p){ # read profiles from p-values
   P <- matrix(p, nrow=2)                                                # create matrix with time points in columns
   Q <- rbind(P[,1:(ncol(P)/2)], P[,(ncol(P)/2 + 1):ncol(P)])            # 4 rows with q-values and numbers exceeding quantiles
   A <- apply(Q, 2, function(x){ sig <- (x[1:2] < sigBord.sets)          # profile code up=1, down=-1, both=2, nothing=0
                                 o <- numeric(2)
                                 o[sig] <- (x[3:4])[sig]                 # index of quantile leading to (up/down) in profile
                                 out <- c(ifelse(o[1] > 0, Q.numb.up[o[1]], o[1]), # relevant quantile for significant up
                                          ifelse(o[2] > 0, Q.numb.do[o[2]], o[2])) # relevant quantile for significant down
                    return(out)
   })
   Quants <- numeric(length(A))
   Quants[seq(1, length(Quants), 2)] <- A[1,]                            # vector of quantiles leading to significant for current group up direction
   Quants[seq(2, length(Quants), 2)] <- A[2,]                            # vector of quantiles leading to significant for current group down direction
   return(Quants)}))                                                     # output to SIGquantsSet
 SIGgenesPerQuantile <- list()
 for(t in colnames(Tmatrix)){
   for(q in Q.numb.up){
     SIGgenesPerQuantile[[t]][[paste("up", q, sep="_")]] <- names(sort(Tmatrix[,t], decreasing=TRUE))[1:q]
   }
   for(q in Q.numb.do){
     SIGgenesPerQuantile[[t]][[paste("do", q, sep="_")]] <- names(sort(Tmatrix[,t]))[1:q]
   }
 }
 SIGgenes <- list()
 for(s in names(sets)){
   genes <- sets[[s]]
   Out <- matrix(0,ncol=ncol(Tmatrix), nrow=length(genes), dimnames=list(Genes=genes, Times=colnames(Tmatrix)))
   if(s %in% rownames(SIGquantsSet)){
   for(i in seq_along(SIGquantsSet[s,])){
     if(SIGquantsSet[s,i] != 0){
       if(i %% 2 == 0){ # down
         wIn <- intersect(genes, SIGgenesPerQuantile[[ceiling(i/2)]][[paste("do", SIGquantsSet[s,i], sep="_")]])
         Out[wIn, ceiling(i/2)] <- -1
       }
       if(i %% 2 == 1){ # up
         wIn <- intersect(genes, SIGgenesPerQuantile[[ceiling(i/2)]][[paste("up", SIGquantsSet[s,i], sep="_")]])
         Out[wIn, ceiling(i/2)] <- 1
       }
     }
   }
   }
   SIGgenes[[s]] <- Out
 }
 time.prof <- proc.time()[3]
 cat(paste("Total time for non-threshold GSA profile algorithm:", round((time.prof - start)/60), "minutes.\n"))
 return(list(PVAL=SIGpvalGR, GOprofiles=SIGprofiles, feasGenes=genes, ListSIGgenes=SIGgenes,
             GenesPerGroup=sets, groupsizes=sapply(sets, length), UsedGr=names(sets),
             Alpha.gr=sigBord.sets, Quantiles=list(up=Q.up, down=Q.do), simpleProfile=SIGprofilesSimple, allPVAL=Pmatrix,
             GeneTmat=Tmatrix))
}
