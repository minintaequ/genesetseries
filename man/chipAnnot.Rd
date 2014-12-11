\name{chipAnnot}
\alias{chipAnnot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
chipAnnot(chip = "hgu133plus2", geneIDs = rownames(PVAL.MeHg$up.q.values))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chip}{
%%     ~~Describe \code{chip} here~~
}
  \item{geneIDs}{
%%     ~~Describe \code{geneIDs} here~~
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
function (chip = "hgu133plus2", geneIDs = rownames(PVAL.MeHg$up.q.values)) 
{
    T <- toTable(get(paste(chip, "ENTREZID", sep = "")))
    EntrezGene.IDs <- T[, 2]
    names(EntrezGene.IDs) = T[, 1]
    T <- toTable(get(paste(chip, "GENENAME", sep = "")))
    GeneNames <- T[, 2]
    names(GeneNames) = T[, 1]
    T <- toTable(get(paste(chip, "SYMBOL", sep = "")))
    GeneSymbols <- T[, 2]
    names(GeneSymbols) = T[, 1]
    T <- toTable(get(paste(chip, "UNIGENE", sep = "")))
    UniGene.IDs <- T[, 2]
    names(UniGene.IDs) = T[, 1]
    T <- toTable(get(paste(chip, "ACCNUM", sep = "")))
    GenBank.IDs <- T[, 2]
    names(GenBank.IDs) = T[, 1]
    makeGeneAnno <- function(set = c("EntrezGene.IDs", "GeneNames", 
        "GeneSymbols", "UniGene.IDs", "GenBank.IDs"), Names = geneIDs) {
        GeneIdentifier <- matrix(NA, ncol = length(set), nrow = length(Names))
        colnames(GeneIdentifier) <- set
        rownames(GeneIdentifier) <- Names
        for (id in set) {
            identifiers <- get(id)
            GeneIdentifier[names(identifiers), id] <- identifiers
        }
        return(GeneIdentifier)
    }
    Names <- makeGeneAnno()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
