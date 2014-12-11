\name{pval.adjust.FUN}
\alias{pval.adjust.FUN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pval.adjust.FUN(pval.matrix, method = c("bonferroni", "holm", "hommel", "tarone", "BH", "fdr", "TS-ABH")[1], signed = FALSE, Neff = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pval.matrix}{
%%     ~~Describe \code{pval.matrix} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{signed}{
%%     ~~Describe \code{signed} here~~
}
  \item{Neff}{
%%     ~~Describe \code{Neff} here~~
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
function (pval.matrix, method = c("bonferroni", "holm", "hommel", 
    "tarone", "BH", "fdr", "TS-ABH")[1], signed = FALSE, Neff = NULL) 
{
    if (!require(multtest)) {
        stop("R-package multtest seems to be not available.")
    }
    if (!is.numeric(pval.matrix)) {
        stop("pval.matrix is supposed to be numeric")
    }
    if (min(pval.matrix) < -1 | max(pval.matrix > 1)) {
        stop("pval.matrix seems not to include (signed) p-values")
    }
    if (!is.matrix(pval.matrix)) {
        warning("pval.matrix was expected to be a matrix. Is transfered to a matrix with column dimension of 1")
        rn <- names(pval.matrix)
        pval.matrix <- matrix(pval.matrix)
        rownames(pval.matrix) <- rn
    }
    pval.matrix.adj <- matrix(NA, ncol = ncol(pval.matrix), nrow = nrow(pval.matrix), 
        dimnames = dimnames(pval.matrix))
    if (!is.null(Neff)) {
        names(Neff[[2]]) <- paste("Neff_cl_", seq(length(Neff[[2]])), 
            sep = "")
        pval.names <- rownames(pval.matrix)
        cl.repr <- t(sapply(Neff[[2]], function(x) {
            y <- pval.matrix[x[1], ]
            w <- which(rownames(pval.matrix) \%in\% x)
            pval.matrix <<- as.matrix(pval.matrix[-w, ])
            pval.names <<- pval.names[-w]
            return(y)
        }))
        if (dim(cl.repr)[1] == 1) 
            cl.repr <- t(cl.repr)
        rownames(cl.repr) <- names(Neff[[2]])
        pval.matrix <- rbind(pval.matrix, cl.repr)
        rownames(pval.matrix) <- c(pval.names, rownames(cl.repr))
    }
    if (signed) 
        s.matrix <- sign(pval.matrix)
    if (method == "TS-ABH") {
        P <- fdr.TST.ABH(abs(pval.matrix))
    }
    else {
        P <- p.adjust(abs(pval.matrix), method[1])
    }
    if (signed) 
        P <- P * s.matrix
    pval.adj.matrix <- matrix(P, ncol = ncol(pval.matrix), dimnames = dimnames(pval.matrix))
    if (!is.null(Neff)) {
        equal.names <- rownames(pval.matrix.adj)[which(rownames(pval.matrix.adj) \%in\% 
            rownames(pval.adj.matrix))]
        pval.matrix.adj[equal.names, ] <- pval.adj.matrix[equal.names, 
            , drop = FALSE]
        for (n in names(Neff[[2]])) {
            fillIn <- t(replicate(length(Neff[[2]][[n]]), pval.adj.matrix[n, 
                ]))
            if (dim(fillIn)[1] == 1) 
                fillIn <- t(fillIn)
            pval.matrix.adj[Neff[[2]][[n]], ] <- fillIn
        }
    }
    else {
        pval.matrix.adj <- pval.adj.matrix
    }
    colnames(pval.matrix.adj) <- colnames(pval.matrix)
    return(pval.matrix.adj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
