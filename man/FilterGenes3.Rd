\name{FilterGenes3}
\alias{FilterGenes3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FilterGenes3(Names = hgu133plus2Names, clusters = "EntrezGene.IDs", EM = PVAL.do$pvalues)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Names}{
%%     ~~Describe \code{Names} here~~
}
  \item{clusters}{
%%     ~~Describe \code{clusters} here~~
}
  \item{EM}{
%%     ~~Describe \code{EM} here~~
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
function (Names = hgu133plus2Names, clusters = "EntrezGene.IDs", 
    EM = PVAL.do$pvalues) 
{
    NAMES <- Names[, clusters]
    FREQ <- table(NAMES)
    g.NA <- which(is.na(NAMES))
    g.UN <- which(NAMES \%in\% names(FREQ[which(FREQ == 1)]))
    f.MT <- names(FREQ[which(FREQ > 1)])
    L.MT <- lapply(f.MT, function(x) return(names(NAMES[which(NAMES == 
        x)])))
    F.MT <- lapply(L.MT, function(f) {
        suffixes <- sapply(strsplit(f, "_"), function(x) ifelse(length(x) > 
            2, x[2], NA))
        if (sum(is.na(suffixes)) > 0) {
            f <- f[is.na(suffixes)]
        }
        else {
            if ("a" \%in\% suffixes) {
                f <- f[which(suffixes == "a")]
            }
            else {
                if ("s" \%in\% suffixes) {
                  f <- f[which(suffixes == "s")]
                }
                else {
                  f <- f
                }
            }
        }
        if (length(f) == 1) {
            s <- f
            return(s)
        }
        EM.f <- EM[f, ]
        s <- f[which.min(apply(EM.f, 1, min))]
        return(s)
    })
    l.un <- names(NAMES[g.UN])
    names(l.un) <- NAMES[g.UN]
    L <- as.list(l.un)
    names(F.MT) <- f.MT
    L <- c(L, F.MT)
    cat(paste("There are", length(NAMES), "features available.\n"))
    cat(paste("Thera are no", clusters, "available for", length(g.NA), 
        "features.\n"))
    cat(paste("Thera are", length(g.UN), clusters, "available with unique chip features.\n"))
    cat(paste("Thera are", length(f.MT), clusters, "with multiple features.\n"))
    cat(paste("Thera are", sum(sapply(L, length)), "features in use for", 
        length(L), clusters, "in whole.\n"))
    return(L)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
