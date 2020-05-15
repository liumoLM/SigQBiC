

#'@title Generate mutation spectrum from twelvemers
#'
#'@description This function generates a mutation spectrum of 96 mutation classes (mutations on pyrimidine centered trinucleotide) from a list of twelvemers
#'
#'@param twelvemers A list of twelvemers, with a 11mer centered at the mutation base and 1 mutated base appended
#'
#'@return A mutation spectrum with 96 mutation types
#'
#'@export
#'
Generate96ChannelSpectrumFromTwelvemers <- function(input.twelvemers.list) {

  mutation.type.list <- ICAMS::catalog.row.order$SBS96
  # Create 2 new columns that show the 3072 and 1536 mutation type
  context <- substring(input.twelvemers.list,5,7)
  mutations <- paste0(context,substring(input.twelvemers.list,12,12))
  # PyrPenta maps to strand-agnostic category
  # e.g. ATGCT>T "ATGCTT" maps to AGCAT>A, "AGCATA"
  pyr.mutation <- PyrPenta(mutations)
  mutation.spectrum <- data.frame(table(pyr.mutation))
  all.mutation.class <- data.frame(mutation.type.list)
  all.mutation.class[, 2] <- 0
  all.mutation.class[, 2] <- mutation.spectrum[match(all.mutation.class[, 1], mutation.spectrum[,
                                                                                                1]), 2]
  if (sum(is.na(all.mutation.class[, 2])) > 0) {
    all.mutation.class[which(is.na(all.mutation.class[, 2])), 2] <- 0
  }
  colnames(all.mutation.class) <- c("mutclass", "proportion")
  all.mutation.class$mutclass <- as.character(all.mutation.class$mutclass)
  return(all.mutation.class)
  # Create part of the 1536 catalog matrix but missing mutation
  # types have NA in the count column.
  # Create the 96 catalog matrix
  # if (is.null(trans.ranges)) {
  #   return(list(catSBS96 = mat96, catSBS1536 = mat1536))
  # }
  # One SBS mutation can be represented by more than 1 row in vcf2 if the mutation
  # position falls into the range of multiple transcripts. When creating the
  # 192 catalog, we only need to count these mutations once.
}

#'@title Convert mutations to pyrimidine centred.
#'
#'@description This function converts trinucleotide context mutations to pyrimidine centred. For example, CGA > CAA will be converted to TCG > TTG
#'
#'@param mutstring A string with for characters: proceeding base - reference base - following base - mutated base (For example, CGA > CCA was input as CGAC)
#'
#'@return A pyrimidine centred mutation
#'
#'@export
#'
PyrPenta <- function(mutstring) {
  stopifnot(nchar(mutstring) == rep(4, length(mutstring)))
  output <-
    ifelse(substr(mutstring, 2, 2) %in% c("A", "G"),
           paste0(revc(substr(mutstring, 1,3)),
                  revc(substr(mutstring, 4,4))),
           mutstring)
  return(output)
}





