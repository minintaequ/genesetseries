\name{ALGORITHM.thresGSA}
\alias{ALGORITHM.thresGSA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ALGORITHM.thresGSA(sets = GROUPSod$groups, genes = GENESod, sigBord.genes = 0.05, adj.genes = "unadjusted", sigBord.sets = 0.01, adj.sets = "unadjusted", Neff = NULL, groupTest = Fisher.group.test, Pmatrix = TTod.np.all)
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
  \item{adj.genes}{
%%     ~~Describe \code{adj.genes} here~~
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
  \item{Pmatrix}{
%%     ~~Describe \code{Pmatrix} here~~
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
function (sets = GROUPSod$groups, genes = GENESod, sigBord.genes = 0.05, 
    adj.genes = "unadjusted", sigBord.sets = 0.01, adj.sets = "unadjusted", 
    Neff = NULL, groupTest = Fisher.group.test, Pmatrix = TTod.np.all) 
{
    start <- proc.time()[3]
    PM <- matrix(NA, ncol = ncol(Pmatrix$up.p.values) + ncol(Pmatrix$do.p.values), 
        nrow = nrow(Pmatrix$do.p.values))
    PM[, seq(1, ncol(Pmatrix$up.p.values) * 2, 2)] <- Pmatrix$up.p.values
    PM[, seq(2, ncol(Pmatrix$do.p.values) * 2, 2)] <- Pmatrix$do.p.values
    rownames(PM) <- rownames(Pmatrix$up.p.values)
    colnames(PM) <- paste(rep(colnames(Pmatrix$do.p.values), 
        each = 2), rep(c("up", "do"), ncol(Pmatrix$do.p.values)), 
        sep = "-")
    Pmatrix <- PM[genes, ]
    if (adj.genes \%in\% c("bonferroni", "holm", "hommel", "tarone", 
        "BH", "fdr", "TS-ABH")) {
        Pmatrix <- pval.adjust.FUN(Pmatrix, method = adj.genes, 
            Neff = NULL)
    }
    else {
        warning("There is no p-value adjustment for gene to reference comparisons done!")
    }
    SIGgenes <- apply(Pmatrix, 2, function(x) rownames(Pmatrix)[which(x <= 
        sigBord.genes)])
    GO.pvalue <- groupTest(Pmatrix, sets, sigBord.genes)
    if (adj.sets \%in\% c("bonferroni", "holm", "hommel", "tarone", 
        "BH", "fdr", "TS-ABH")) {
        GO.pvalue <- pval.adjust.FUN(GO.pvalue, method = adj.sets, 
            Neff = Neff)
    }
    else {
        warning("There is no p-value adjustment for group enrichment done!")
    }
    colnames(GO.pvalue) <- colnames(PM)
    SIGpvalGR <- GO.pvalue[apply(GO.pvalue, 1, function(x) any(x < 
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
    cat(paste("Algorithm thresGSA found", length(SIGprofilesSimple), 
        "non-zero activation profiles in a run time of", round((proc.time()[3] - 
            start)/60), "minutes.\n"))
    return(list(PVAL = SIGpvalGR, GOprofiles = SIGprofiles, feasGenes = genes, 
        ListSIGgenes = SIGgenes, n.SIGgenes = sapply(SIGgenes, 
            length), GenesPerGroup = sets, groupsizes = sapply(sets, 
            length), UsedGr = names(sets), p.val.adjust.sets.method = adj.sets, 
        Alpha.gr = sigBord.sets, Alpha.ge = sigBord.genes, simpleProfile = SIGprofilesSimple, 
        allPVAL = GO.pvalue, GenePVAL = Pmatrix, p.val.adjust.genes.method = adj.genes, 
        p.val.adjust.sets.method = adj.sets))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
