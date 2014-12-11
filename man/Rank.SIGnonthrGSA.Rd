\name{Rank.SIGnonthrGSA}
\alias{Rank.SIGnonthrGSA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Rank.SIGnonthrGSA(PROF = PROFILESdo.SIGnonthrGSA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PROF}{
%%     ~~Describe \code{PROF} here~~
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
function (PROF = PROFILESdo.SIGnonthrGSA) 
{
    q <- PROF$Quantiles
    P <- PROF$PVAL
    T <- PROF$GeneTmat
    A <- PROF$GOprofiles
    R <- t(apply(cbind(rownames(A), A), 1, function(prof) {
        gr.name <- prof[1]
        lsig <- PROF$ListSIGgenes[[gr.name]]
        test <- T[rownames(lsig), ]
        score <- sum(abs(lsig * test))/sum(lsig != 0)
        ngenes <- nrow(lsig)
        sgenes <- sum(apply(lsig, 1, function(x) any(x != 0)))
        profile <- PROF$simpleProfile[gr.name]
        term <- ifelse(grepl("GO:", gr.name), Term(gr.name), 
            "KEGG pathway tba")
        out <- c(score, gr.name, profile, ngenes, sgenes, term)
        return(out)
    }))
    R <- R[order(as.numeric(R[, 1]), decreasing = TRUE), ]
    return(R)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
