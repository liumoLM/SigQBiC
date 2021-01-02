#'Identify pathways that enriched with affected TFs by a given mutational signature
#'
#'@description A function gives the pathways enriched with affected TFs by a given mutational signature
#'
#'@param TF.list A list of affected TFs
#'
#'@param qvalue A numeric of qvalue as filter
#'
#'@param dbs.selected A character of database name in 'enrichR'. In our paper, we used 'Reactome_2016'
#'
#'@return A matrix of pathways enriched with given TFs
#'
#'@import enrichR
#'
#'@export
#'
EnrichrPathwayAnalysis <- function(TF.list, qvalue, dbs.selected) {
  enriched <- enrichr(TF.list, dbs.selected)
  enriched.matrix <- enriched[[dbs.selected]]
  enriched.matrix$TF.counts <- length(unique(TF.list))
  enriched.matrix$TFs <- apply(enriched.matrix.gain, 1, function(x) {
    x["TFs"] <- unlist(strsplit(x["Overlap"], "/"))[1]
  })
  enriched.matrix <- enriched.matrix[enriched.matrix$Adjusted.P.value < qvalue,
  ]
  enriched.matrix$New.Term <- apply(enriched.matrix, 1, function(x) {
    x["New.Term"] <- unlist(strsplit(x["Term"], " Homo", perl = T))[1]
  })
  return(enriched.matrix)
}

