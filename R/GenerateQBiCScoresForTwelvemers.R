#'@title Generate QBiC scores for all uPBMs for a given list of twelvemers
#'
#'@description This function generates QBiC scores for all uPBM experiments for given twelvemers
#'
#'@param twelvemers A list of twelvemers, with a 11mer centered at the mutation base and 1 mutated base appended
#'
#'@param uPBM_QBiC_scores All QBiC scores for a universal PBM
#'
#'@return A list of scores for the given twelvemers list
#'
#'@export
#'
GenerateQBiCScoresForTwelvemers <- function(uPBM_QBiC_scores, twelvemers) {
  twelvemers <- data.frame(twelvemers)
  seq_NA <- data.frame(setdiff(twelvemers[, 1], all.possible.twelvemers$seq))  ##select 12mers without mutation
  seq_mut <- data.frame(setdiff(twelvemers[, 1], seq_NA[, 1]))
  colnames(seq_NA) <- "seq"
  colnames(seq_mut) <- "seq"
  seq_mut[, 1] <- as.character(seq_mut[, 1])
  seq_scores <- seq_mut


  uPBM_QBiC_scores <- data.frame(uPBM_QBiC_scores)
  uPBM_QBiC_scores <- uPBM_QBiC_scores$z_score[which(!is.na(uPBM_QBiC_scores$z_score))]


  seq_scores[, 2] <- uPBM_QBiC_scores[match(seq_mut$seq, all.possible.twelvemers$seq)]
  if (dim(seq_NA)[1] > 0) {
    seq_NA[, 2] <- "NA"
    seq_scores <- rbind(seq_scores, seq_NA)
  }

  colnames(seq_scores) <- c("seq", "QBiC.scores")

  return(seq_scores)
}
