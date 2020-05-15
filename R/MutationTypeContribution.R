#'Compute the contribution
#'
#'@param uPBM_QBiC_scores All QBiC scores for a universal PBM
#'
#'@param spectrum A list of possibilities for each mutation type(mutational spectrum) or a character of name of mutational signature(mutational signature)
#'
#'@param p_values p_values for all twelvemers
#'
#'@return A matrix with contribution to the GR and LR for each of 96 mutation types
#'
#'@export
#'
MutationTypeContribution <- function(uPBM_QBiC_scores, p_values, spectrum) {
  summary.matrix <- data.frame(mutation.type.list)
  summary.matrix[, 2:5] <- 0
  summary.matrix[, 1] <- as.character(summary.matrix[, 1])
  uPBM_QBiC_scores <- data.frame(uPBM_QBiC_scores)
  uPBM_QBiC_scores <- uPBM_QBiC_scores$z_score[which(!is.na(uPBM_QBiC_scores$z_score))]
  PBM.scores <- data.frame(uPBM_QBiC_scores)

  row.names(PBM.scores) <- all.possible.twelvemers$seq
  PBM.scores$mutclass <- all.possible.twelvemers$mutclass
  number <- as.integer(max(PBM.scores[, 1])) + 2

  mutation.spectrum <- check.spectrum(spectrum)
  p_values <- data.frame(p_values)
  p_values <- p_values[!is.na(p_values[,1]),1]
  adjust_p_value <- p.adjust(p_values)
  PBM.scores <- PBM.scores[adjust_p_value<0.05,]

  k <- 0
  for (mutation.type in mutation.type.list) {
    k <- k + 1
    histogram_plot <- hist(PBM.scores[which(PBM.scores$mutclass ==
                                              mutation.type), 1], breaks = seq(-number, number, 0.5), plot = F)
    counts_signature.temp <- histogram_plot$counts * mutation.spectrum[mutation.type,
                                                                       2]
    neg_index <- which(histogram_plot$mids < 0)
    pos_index <- which(histogram_plot$mids > 0)
    summary.matrix[k, 2] <- sum(histogram_plot$mids[neg_index] * counts_signature.temp[neg_index])
    summary.matrix[k, 3] <- sum(histogram_plot$mids[pos_index] * counts_signature.temp[pos_index])
  }
  summary.matrix[, 4] <- as.numeric(summary.matrix[, 3])/sum(as.numeric(summary.matrix[,
                                                                                       3]))
  summary.matrix[, 5] <- as.numeric(summary.matrix[, 2])/sum(as.numeric(summary.matrix[,
                                                                                       2]))
  summary.matrix <- summary.matrix[, c(1, 4, 5)]  ##this gives a contribution of all final scores, in another word: all sums
  colnames(summary.matrix) <- c("Mutation.Class", "Contribute.to.Gain", "Contribute.to.Loss")
  return(summary.matrix)
}
