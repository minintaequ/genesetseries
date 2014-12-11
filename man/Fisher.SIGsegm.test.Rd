\name{Fisher.SIGsegm.test}
\alias{Fisher.SIGsegm.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Fisher.SIGsegm.test(value.matrix, groups, quantiles.up, quantiles.do)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{value.matrix}{
%%     ~~Describe \code{value.matrix} here~~
}
  \item{groups}{
%%     ~~Describe \code{groups} here~~
}
  \item{quantiles.up}{
%%     ~~Describe \code{quantiles.up} here~~
}
  \item{quantiles.do}{
%%     ~~Describe \code{quantiles.do} here~~
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
function (value.matrix, groups, quantiles.up, quantiles.do) 
{
    P.up <- NULL
    Time.list.up <- list()
    if (length(quantiles.up) > 0) {
        for (t in seq(length = ncol(value.matrix))) {
            Time.list.up[[LETTERS[t]]] <- NULL
        }
        P <- apply(quantiles.up, 1, function(q) {
            qs <- apply(rbind(q, value.matrix), 2, function(t) return(quantile(t[-1], 
                t[1])))
            dec.matrix <- t(t(value.matrix) >= qs)
            n.int <- apply(dec.matrix, 2, sum)
            out <- sapply(groups, function(g) {
                int.in.group <- colSums(dec.matrix[g, ])
                p <- phyper(int.in.group - 1, n.int, nrow(dec.matrix) - 
                  n.int, length(g), lower.tail = FALSE)
                return(p)
            })
            for (t in seq(length = ncol(dec.matrix))) {
                Time.list.up[[LETTERS[t]]] <<- cbind(Time.list.up[[LETTERS[t]]], 
                  out[t, ])
            }
        })
        if (length(colnames(value.matrix)) > 0) 
            names(Time.list.up) <- colnames(value.matrix)
        for (t in seq(length = ncol(value.matrix))) P.up <- cbind(P.up, 
            Time.list.up[[t]])
        if (length(colnames(value.matrix)) > 0) 
            colnames(P.up) <- rep(colnames(value.matrix), each = nrow(quantiles.up))
    }
    P.do <- NULL
    Time.list.do <- list()
    if (length(quantiles.do) > 0) {
        for (t in seq(length = ncol(value.matrix))) {
            Time.list.do[[LETTERS[t]]] <- NULL
        }
        P <- apply(quantiles.do, 1, function(q) {
            qs <- apply(rbind(q, value.matrix), 2, function(t) return(quantile(t[-1], 
                t[1])))
            dec.matrix <- t(t(value.matrix) <= qs)
            n.int <- apply(dec.matrix, 2, sum)
            out <- sapply(groups, function(g) {
                int.in.group <- colSums(dec.matrix[g, ])
                p <- phyper(int.in.group - 1, n.int, nrow(dec.matrix) - 
                  n.int, length(g), lower.tail = FALSE)
                return(p)
            })
            for (t in seq(length = ncol(dec.matrix))) {
                Time.list.do[[LETTERS[t]]] <<- cbind(Time.list.do[[LETTERS[t]]], 
                  out[t, ])
            }
        })
        if (length(colnames(value.matrix)) > 0) 
            names(Time.list.do) <- colnames(value.matrix)
        for (t in seq(length = ncol(value.matrix))) P.do <- cbind(P.do, 
            Time.list.do[[t]])
        if (length(colnames(value.matrix)) > 0) 
            colnames(P.do) <- rep(colnames(value.matrix), each = nrow(quantiles.do))
    }
    return(list(up = P.up, down = P.do))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
