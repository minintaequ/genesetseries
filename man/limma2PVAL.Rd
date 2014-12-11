\name{limma2PVAL}
\alias{limma2PVAL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
limma2PVAL(data, limma, pvadj = "BH", ref = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{limma}{
%%     ~~Describe \code{limma} here~~
}
  \item{pvadj}{
%%     ~~Describe \code{pvadj} here~~
}
  \item{ref}{
%%     ~~Describe \code{ref} here~~
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
function (data, limma, pvadj = "BH", ref = NULL) 
{
    PV <- matrix(p.adjust(limma$p.value, method = "BH"), ncol = ncol(limma$t))
    dimnames(PV) <- dimnames(limma$p.value)
    dir <- ifelse(limma$t < 0, -1, 1)
    PVu <- ifelse(limma$t > 0, PV, 0.999999)
    PVd <- ifelse(limma$t < 0, PV, 0.999999)
    calc.logFC <- function(E = data, des = limma$design, w = ref) {
        if (!is.null(ref)) {
            M = E[, des[, ref]]
            notref <- seq(ncol(des))[-ref]
        }
        else {
            M = matrix(0, ncol = 1, nrow = nrow(E))
            notref <- seq(ncol(des))
        }
        O = NULL
        for (p in notref) {
            O <- cbind(O, rowMeans(E[, des[, p]]) - rowMeans(M))
        }
        rownames(O) <- rownames(E)
        colnames(O) <- colnames(des)[notref]
        return(O)
    }
    FCdiff <- calc.logFC()
    limma$statistic <- limma$t
    limma$q.value <- PV
    limma$up.p.values <- ifelse(limma$t > 0, limma$p.value, 0.99999)
    limma$do.p.values <- ifelse(limma$t < 0, limma$p.value, 0.99999)
    limma$up.q.values <- PVu
    limma$do.q.values <- PVd
    limma$logFC <- FCdiff
    limma$dir <- dir
    limma$adjust.method <- pvadj
    return(limma)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
