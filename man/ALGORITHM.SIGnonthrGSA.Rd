\name{ALGORITHM.SIGnonthrGSA}
\alias{ALGORITHM.SIGnonthrGSA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ALGORITHM.SIGnonthrGSA(sets = GROUPSdo$groups, genes = GENESdo, sigBord.genes = 0.05, Q.up.l = 10, Q.do.l = 10, sigBord.sets = 0.01, adj.sets = "unadjusted", Neff = NULL, groupTest = Fisher.SIGsegm.test, Tobject = PVAL.do)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sets}{
%%     ~~Describe \code{sets} here~~
}
  \item{genes}{
%%     ~~Describe \code{genes} here~~
}
  \item{sigBord.genes}{
%%     ~~Describe \code{sigBord.genes} here~~
}
  \item{Q.up.l}{
%%     ~~Describe \code{Q.up.l} here~~
}
  \item{Q.do.l}{
%%     ~~Describe \code{Q.do.l} here~~
}
  \item{sigBord.sets}{
%%     ~~Describe \code{sigBord.sets} here~~
}
  \item{adj.sets}{
%%     ~~Describe \code{adj.sets} here~~
}
  \item{Neff}{
%%     ~~Describe \code{Neff} here~~
}
  \item{groupTest}{
%%     ~~Describe \code{groupTest} here~~
}
  \item{Tobject}{
%%     ~~Describe \code{Tobject} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (sets = GROUPSdo$groups, genes = GENESdo, sigBord.genes = 0.05, 
    Q.up.l = 10, Q.do.l = 10, sigBord.sets = 0.01, adj.sets = "unadjusted", 
    Neff = NULL, groupTest = Fisher.SIGsegm.test, Tobject = PVAL.do) 
{
    start <- proc.time()[3]
    q.up <- 1 - apply(Tobject$up.q.values, 2, function(x) sum(x <= 
        sigBord.genes))/length(genes)
    q.do <- apply(Tobject$do.q.values, 2, function(x) sum(x <= 
        sigBord.genes))/length(genes)
    Q.up <- sapply(q.up, function(x) seq(x, 0.99999, length.out = Q.up.l))
    Q.do <- sapply(q.do, function(x) rev(seq(1e-05, x, length.out = Q.do.l)))
    Tmatrix <- Tobject$statistic[genes, ]
    SET.pvalues <- groupTest(Tmatrix, sets, Q.up, Q.do)
    time.segm <- proc.time()[3]
    cat(paste("Segment test calculation lasts ", (time.segm - 
        start), "seconds.\n"))
    if (adj.sets \%in\% c("bonferroni", "holm", "hommel", "tarone", 
        "BH", "fdr", "TS-ABH")) {
        all.SET.pvalues <- pval.adjust.FUN(cbind(SET.pvalues$up, 
            SET.pvalues$down), method = adj.sets, Neff = Neff)
        time.adj <- proc.time()[3]
        cat(paste("Multiple testing adjustment lasts ", (time.adj - 
            time.segm), "seconds.\n"))
    }
    else {
        warning("There is no p-value adjustment for group enrichment done!")
    }
    SIGquant <- NULL
    Pmatrix <- t(apply(all.SET.pvalues, 1, function(q) {
        M.up <- matrix(q[1:(nrow(Q.up) * ncol(Tmatrix))], ncol = ncol(Tmatrix))
        M.do <- matrix(q[(nrow(Q.up) * ncol(Tmatrix) + 1):length(q)], 
            ncol = ncol(Tmatrix))
        m.up <- apply(M.up, 2, min)
        m.do <- apply(M.do, 2, min)
        w.up <- apply(M.up, 2, function(x) which(x == min(x))[1])
        w.do <- apply(M.do, 2, function(x) which(x == min(x))[1])
        sigq <- numeric(length(m.up) + length(m.do))
        sigq[seq(1, length(sigq), 2)] <- w.up
        sigq[seq(2, length(sigq), 2)] <- w.do
        SIGquant <<- rbind(SIGquant, sigq)
        out <- numeric(length(m.up) + length(m.do))
        out[seq(1, length(out), 2)] <- m.up
        out[seq(2, length(out), 2)] <- m.do
        return(out)
    }))
    rownames(SIGquant) <- rownames(Pmatrix)
    colnames(Pmatrix) <- colnames(SIGquant) <- paste(rep(colnames(Tmatrix), 
        each = 2), c("up", "down"), sep = "-")
    SIGpvalGR <- Pmatrix[apply(Pmatrix, 1, function(x) any(x < 
        sigBord.sets)), ]
    SIGprofiles <- t(apply(SIGpvalGR, 1, function(p) {
        P <- matrix(p, nrow = 2)
        A <- apply(P, 2, function(x) {
            sig <- (x < sigBord.sets)
            return(ifelse(!any(sig), 0, ifelse(all(sig), 2, sum(c(1, 
                -1) * sig))))
        })
        return(A)
    }))
    SIGprofilesSimple <- simpleProfile(SIGprofiles)
    Q.numb.up <- apply(Q.up, 2, function(x) {
        ceiling(nrow(Tmatrix) - quantile(1:nrow(Tmatrix), x))
    })
    Q.numb.do <- apply(Q.do, 2, function(x) {
        floor(quantile(1:nrow(Tmatrix), x))
    })
    rownames(Q.numb.up) <- 1:Q.up.l
    rownames(Q.numb.do) <- 1:Q.do.l
    SIGquantsSet <- t(apply(cbind(SIGpvalGR, SIGquant[rownames(SIGpvalGR), 
        ]), 1, function(p) {
        P <- matrix(p, nrow = 2)
        Q <- rbind(P[, 1:(ncol(P)/2)], P[, (ncol(P)/2 + 1):ncol(P)])
        Q <- rbind(Q, 1:ncol(Q))
        A <- apply(Q, 2, function(x) {
            sig <- (x[1:2] < sigBord.sets)
            o <- numeric(2)
            o[sig] <- (x[3:4])[sig]
            out <- c(ifelse(o[1] > 0, Q.numb.up[o[1], x[5]], 
                o[1]), ifelse(o[2] > 0, Q.numb.do[o[2], x[5]], 
                o[2]))
            return(out)
        })
        Quants <- numeric(length(A))
        Quants[seq(1, length(Quants), 2)] <- A[1, ]
        Quants[seq(2, length(Quants), 2)] <- A[2, ]
        return(Quants)
    }))
    SIGgenesPerQuantile <- list()
    for (t in colnames(Tmatrix)) {
        for (q in Q.numb.up[, t]) {
            SIGgenesPerQuantile[[t]][[paste("up", q, sep = "_")]] <- names(sort(Tmatrix[, 
                t], decreasing = TRUE))[1:q]
        }
        for (q in Q.numb.do[, t]) {
            SIGgenesPerQuantile[[t]][[paste("do", q, sep = "_")]] <- names(sort(Tmatrix[, 
                t]))[1:q]
        }
    }
    SIGgenes <- list()
    for (s in names(sets)) {
        genes <- sets[[s]]
        Out <- matrix(0, ncol = ncol(Tmatrix), nrow = length(genes), 
            dimnames = list(Genes = genes, Times = colnames(Tmatrix)))
        if (s \%in\% rownames(SIGquantsSet)) {
            for (i in seq_along(SIGquantsSet[s, ])) {
                if (SIGquantsSet[s, i] != 0) {
                  if (i\%\%2 == 0) {
                    wIn <- intersect(genes, SIGgenesPerQuantile[[ceiling(i/2)]][[paste("do", 
                      SIGquantsSet[s, i], sep = "_")]])
                    Out[wIn, ceiling(i/2)] <- -1
                  }
                  if (i\%\%2 == 1) {
                    wIn <- intersect(genes, SIGgenesPerQuantile[[ceiling(i/2)]][[paste("up", 
                      SIGquantsSet[s, i], sep = "_")]])
                    Out[wIn, ceiling(i/2)] <- 1
                  }
                }
            }
        }
        SIGgenes[[s]] <- Out
    }
    time.prof <- proc.time()[3]
    cat(paste("Total time for non-threshold GSA profile algorithm:", 
        round((time.prof - start)/60), "minutes.\n"))
    return(list(PVAL = SIGpvalGR, GOprofiles = SIGprofiles, feasGenes = genes, 
        ListSIGgenes = SIGgenes, GenesPerGroup = sets, groupsizes = sapply(sets, 
            length), UsedGr = names(sets), Alpha.gr = sigBord.sets, 
        Quantiles = list(up = Q.up, down = Q.do), simpleProfile = SIGprofilesSimple, 
        allPVAL = Pmatrix, GeneTmat = Tmatrix))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
