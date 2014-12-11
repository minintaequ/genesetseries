## annotate genes to GO groups of a predefined ontology
#
annotate.GO <- function (feasibleGenes = NULL, ontology="BP", affyLib)
      {                                                                     # annFUN.db from topGO
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
            retVal <- dbGetQuery(get(paste(affyLib, "dbconn", sep = "_"))(), .sql)
            if (!is.null(feasibleGenes))
                      retVal <- retVal[retVal[["probe_id"]] %in% feasibleGenes,
                                                   ]
            genesInGroups <- split(retVal[["probe_id"]], retVal[["go_id"]]) # end annFUN.db from topGO
            genes <- unique(unlist(genesInGroups))
            cat(paste(length(genes), "(of", length(feasibleGenes), "ca.", round(100*length(genes)/length(feasibleGenes),1), "%) genes in", length(genesInGroups), "specific GO groups of the", ontology, "ontology.\n"))
            groupNames <- names(genesInGroups)
            parentObject <- as.list(get(paste("GO", ontology, "PARENTS", sep="")))
            childrenObject <- as.list(get(paste("GO", ontology, "CHILDREN", sep="")))
            allgroups <- unique(c(unlist(parentObject), unlist(childrenObject)))
            nodesInFocus <- groupNames
            nodesWithoutChildren <- names(which(sapply(childrenObject[nodesInFocus], function(x){
                                               z <- all(is.na(x))                    # only NA children
                                               y <- !any(x %in% nodesInFocus)        # only children within nodes not occupied with genes in matrix
                                               return(z | y)
                                               })))
            nodesToConnect <- nodesWithoutChildren
            edgeL <- list()
            edgeData <- list()
            while(length(nodesToConnect) > 1){                              # until root node do
              PA  <- parentObject[nodesToConnect]
              nodesToConnect <- NULL
              nPA <- names(PA)
              nPA <- nPA[!(is.na(nPA))]
            nothing <- sapply(nPA, function(x){
                PA.x <- PA[[x]]
                nodesToConnect <<- unique(c(nodesToConnect, PA.x))
                edgeL[[x]] <<- list(edges=NULL, types=NULL, weights=NULL)
                if(length(PA.x)==1 & PA.x[1] == "all"){return(invisible(NaN))}
                for(group in PA.x){
                  genesInGroups[[group]] <<- unique(c(genesInGroups[[group]], genesInGroups[[x]]))
                  edgeL[[x]][["edges"]] <<- c(edgeL[[x]][["edges"]], group)
                  edgeL[[x]][["types"]] <<- c(edgeL[[x]][["types"]], names(PA.x)[which(PA.x==group)])
                  edgeL[[x]][["weights"]] <<- c(edgeL[[x]][["weights"]], -1)
                  edgeData[[paste(x, group, sep="|")]] <<- list(weights=-1, type=names(PA.x)[which(PA.x==group)])
                }
                return(invisible(NaN))
              })
            }
            cat(paste("Whole induced graph in", ontology, "ontology consists of", length(edgeL), "nodes and", length(edgeData), "edges.\n"))
            groupNames <- names(edgeL)
            for(x in groupNames){ edgeL[[x]][["edges"]] <- which(groupNames %in% edgeL[[x]][["edges"]])}
            graph <- new("graphNEL", nodes=groupNames, edgeL=edgeL, edgemode="directed")
            return(list(graph=graph, edgeData=edgeData, groups=genesInGroups, genes=genes))
      }
