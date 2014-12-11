\name{GeneSetFilter}
\alias{GeneSetFilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GeneSetFilter(geneSets = ANNSETSdo, minSetSize = 4, Rank = FALSE, features = unlist(EntrezGeneGroups.do), exclude = "GO:0008150")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneSets}{
%%     ~~Describe \code{geneSets} here~~
}
  \item{minSetSize}{
%%     ~~Describe \code{minSetSize} here~~
}
  \item{Rank}{
%%     ~~Describe \code{Rank} here~~
}
  \item{features}{
%%     ~~Describe \code{features} here~~
}
  \item{exclude}{
%%     ~~Describe \code{exclude} here~~
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
function (geneSets = ANNSETSdo, minSetSize = 4, Rank = FALSE, 
    features = unlist(EntrezGeneGroups.do), exclude = "GO:0008150") 
{
    start <- proc.time()
    group.genes <- list()
    if (is.list(geneSets[[1]])) {
        for (g in names(geneSets)) {
            print(g)
            for (n in names(geneSets[[g]])) {
                group.genes[[n]] <- geneSets[[g]][[n]]
            }
        }
    }
    else {
        group.genes <- geneSets
    }
    geneSets <- group.genes
    if (Rank) {
        genes <- unique(unlist(geneSets))
        M_gene_group <- matrix(0, ncol = length(geneSets), nrow = length(genes))
        i <- 1
        sapply(genes, function(g) {
            a <- sapply(geneSets, function(GO) {
                return(as.numeric(g \%in\% GO))
            })
            M_gene_group[i, ] <<- a
            i <<- i + 1
        })
        Neff <- rank(M_gene_group)
    }
    else {
        length.groups <- sapply(geneSets, length)
        Table.l.gr <- table(length.groups)
        if (minSetSize > 0) {
            Table.l.gr <- Table.l.gr[-which(names(Table.l.gr) \%in\% 
                seq_len(minSetSize - 1))]
            geneSets <- geneSets[which(sapply(geneSets, length) >= 
                minSetSize)]
        }
        cat(paste(sum(Table.l.gr), "GO groups fullfill the minsize requirement of", 
            minSetSize, "genes.\n"))
        assign("Neff", sum(Table.l.gr == 1), envir = .GlobalEnv)
        assign("Index", NULL, envir = .GlobalEnv)
        assign("mult.GR", list(), envir = .GlobalEnv)
        for (j in names(Table.l.gr[which(Table.l.gr > 1)])) {
            w <- names(length.groups)[which(length.groups == 
                as.numeric(j))]
            index <- rep(TRUE, length(w))
            names(index) <- w
            test <- NULL
            for (i in w) {
                test <- c(test, i)
                if (index[i]) {
                  inGroup.i <- geneSets[[i]]
                  identical <- sapply(geneSets[w[index]], function(x) {
                    return(all(x \%in\% inGroup.i))
                  })
                  if (sum(identical) > 1) {
                    mult.GR[[length(mult.GR) + 1]] <- names(geneSets[w[index]])[identical]
                  }
                  Neff <<- Neff + 1
                  index[names(geneSets[w[index]])[identical]] <- FALSE
                }
            }
            cat(paste(length(index), "groups with total genes", 
                j, "had processed.\n"))
        }
    }
    duration <- (proc.time() - start)[3]
    cat("Neff calculation takes", round(duration/60), "minutes and", 
        round(duration\%\%60), "seconds.", "Effective group number is Neff =", 
        Neff, "and multiple groups split into", length(mult.GR), 
        "clusters.\n")
    if (length(exclude) > 0) {
        ALL.genes.without.excluded <- unique(unlist(geneSets[-which(names(geneSets) == 
            exclude)]))
        cat(paste(length(features) - length(ALL.genes.without.excluded), 
            "genes excluded, because of excluded sets (GO BP).\n"))
        for (e in exclude) {
            geneSets[[e]] <- intersect(geneSets[[e]], ALL.genes.without.excluded)
        }
        features <- ALL.genes.without.excluded
    }
    time.excl <- (proc.time() - start)[3] - duration
    cat("Feature extraction takes", round(time.excl), "seconds and number of actually used genes is", 
        length(features))
    return(list(effectiveNumber = Neff, idGroupList = mult.GR, 
        groups = geneSets, genes = features))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
