\name{annotate.GO}
\alias{annotate.GO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
annotate.GO(feasibleGenes = NULL, ontology = "BP", affyLib)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{feasibleGenes}{
%%     ~~Describe \code{feasibleGenes} here~~
}
  \item{ontology}{
%%     ~~Describe \code{ontology} here~~
}
  \item{affyLib}{
%%     ~~Describe \code{affyLib} here~~
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
function (feasibleGenes = NULL, ontology = "BP", affyLib) 
{
    affyLib <- paste(sub(".db$", "", affyLib), ".db", sep = "")
    require(affyLib, character.only = TRUE) || stop(paste("package", 
        affyLib, "is required", sep = " "))
    affyLib <- sub(".db$", "", affyLib)
    orgFile <- get(paste(get(paste(affyLib, "ORGPKG", sep = "")), 
        "_dbfile", sep = ""))
    try(dbGetQuery(get(paste(affyLib, "dbconn", sep = "_"))(), 
        paste("ATTACH '", orgFile(), "' as org;", sep = "")), 
        silent = TRUE)
    .sql <- paste("SELECT DISTINCT probe_id, go_id FROM probes INNER JOIN ", 
        "(SELECT * FROM org.genes INNER JOIN org.go_", tolower(ontology), 
        " USING('_id')) USING('gene_id');", sep = "")
    retVal <- dbGetQuery(get(paste(affyLib, "dbconn", sep = "_"))(), 
        .sql)
    if (!is.null(feasibleGenes)) 
        retVal <- retVal[retVal[["probe_id"]] \%in\% feasibleGenes, 
            ]
    genesInGroups <- split(retVal[["probe_id"]], retVal[["go_id"]])
    genes <- unique(unlist(genesInGroups))
    cat(paste(length(genes), "(of", length(feasibleGenes), "ca.", 
        round(100 * length(genes)/length(feasibleGenes), 1), 
        "\%) genes in", length(genesInGroups), "specific GO groups of the", 
        ontology, "ontology.\n"))
    groupNames <- names(genesInGroups)
    parentObject <- as.list(get(paste("GO", ontology, "PARENTS", 
        sep = "")))
    childrenObject <- as.list(get(paste("GO", ontology, "CHILDREN", 
        sep = "")))
    allgroups <- unique(c(unlist(parentObject), unlist(childrenObject)))
    nodesInFocus <- groupNames
    nodesWithoutChildren <- names(which(sapply(childrenObject[nodesInFocus], 
        function(x) {
            z <- all(is.na(x))
            y <- !any(x \%in\% nodesInFocus)
            return(z | y)
        })))
    nodesToConnect <- nodesWithoutChildren
    edgeL <- list()
    edgeData <- list()
    while (length(nodesToConnect) > 1) {
        PA <- parentObject[nodesToConnect]
        nodesToConnect <- NULL
        nPA <- names(PA)
        nPA <- nPA[!(is.na(nPA))]
        nothing <- sapply(nPA, function(x) {
            PA.x <- PA[[x]]
            nodesToConnect <<- unique(c(nodesToConnect, PA.x))
            edgeL[[x]] <<- list(edges = NULL, types = NULL, weights = NULL)
            if (length(PA.x) == 1 & PA.x[1] == "all") {
                return(invisible(NaN))
            }
            for (group in PA.x) {
                genesInGroups[[group]] <<- unique(c(genesInGroups[[group]], 
                  genesInGroups[[x]]))
                edgeL[[x]][["edges"]] <<- c(edgeL[[x]][["edges"]], 
                  group)
                edgeL[[x]][["types"]] <<- c(edgeL[[x]][["types"]], 
                  names(PA.x)[which(PA.x == group)])
                edgeL[[x]][["weights"]] <<- c(edgeL[[x]][["weights"]], 
                  -1)
                edgeData[[paste(x, group, sep = "|")]] <<- list(weights = -1, 
                  type = names(PA.x)[which(PA.x == group)])
            }
            return(invisible(NaN))
        })
    }
    cat(paste("Whole induced graph in", ontology, "ontology consists of", 
        length(edgeL), "nodes and", length(edgeData), "edges.\n"))
    groupNames <- names(edgeL)
    for (x in groupNames) {
        edgeL[[x]][["edges"]] <- which(groupNames \%in\% edgeL[[x]][["edges"]])
    }
    graph <- new("graphNEL", nodes = groupNames, edgeL = edgeL, 
        edgemode = "directed")
    return(list(graph = graph, edgeData = edgeData, groups = genesInGroups, 
        genes = genes))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
