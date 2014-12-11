\name{Fisher.group.test}
\alias{Fisher.group.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Fisher.group.test(value.matrix, groups, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{value.matrix}{
%%     ~~Describe \code{value.matrix} here~~
}
  \item{groups}{
%%     ~~Describe \code{groups} here~~
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
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
function (value.matrix, groups, threshold) 
{
    dec.matrix <- value.matrix <= threshold
    n.int <- apply(dec.matrix, 2, sum)
    Out <- sapply(groups, function(g) {
        int.in.group <- colSums(dec.matrix[g, ])
        p <- phyper(int.in.group - 1, n.int, nrow(dec.matrix) - 
            n.int, length(g), lower.tail = FALSE)
        return(p)
    })
    return(t(Out))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
