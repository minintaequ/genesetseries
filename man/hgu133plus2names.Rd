\name{hgu133plus2names}
\alias{hgu133plus2names}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hgu133plus2names(v, N = hgu133plus2Names)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{
%%     ~~Describe \code{v} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
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
function (v, N = hgu133plus2Names) 
{
    N <- cbind(rownames(N), N)
    w <- numeric(length(v))
    for (i in seq_along(v)) {
        a <- which(N[, 1] \%in\% v[i])
        b <- which(N[, 2] \%in\% v[i])
        c <- grep(v[i], N[, 3], TRUE)
        d <- grep(v[i], N[, 4], TRUE)
        e <- which(N[, 5] \%in\% v[i])
        f <- which(N[, 6] \%in\% v[i])
        if (length(d) > 1) 
            d <- d[which(nchar(N[d, 4]) == nchar(v[i]))]
        w[i] <- c(a, b, d, e, f, c)[1]
    }
    return(N[w, ])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
