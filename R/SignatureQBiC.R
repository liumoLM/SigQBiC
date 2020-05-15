#'Generate GR and LR for a uPBM and a signature (Signature-QBiC)
#'
#'@description SignatureQBiC model. This function generates Gain Ratio and Loss Ratio for a TF (represented by universal PBM) with a mutational signature (or a mutation spectrum)
#'
#'@param uPBM_QBiC_scores All QBiC scores for a universal PBM
#'
#'@param p_values p_values for all twelvemers
#'
#'@param spectrum A list of possibilities for each mutation type(mutational spectrum) or a character of name of mutational signature(mutational signature)
#'
#'@return A list of GR and LR
#'

SignatureQBiC <- function(uPBM_QBiC_scores, p_values, spectrum) {
  uPBM_QBiC_scores <- data.frame(uPBM_QBiC_scores)
  uPBM_QBiC_scores <- uPBM_QBiC_scores$z_score[!is.na(uPBM_QBiC_scores$z_score)]
  PBM.scores <- data.frame(uPBM_QBiC_scores)

  row.names(PBM.scores) <- all.possible.twelvemers$seq
  PBM.scores$mutclass <- all.possible.twelvemers$mutclass
  number <- as.integer(max(PBM.scores[, 1])) + 2

  mutation.spectrum <- check.spectrum(spectrum)
  p_values <- data.frame(p_values)
  p_values <- p_values[!is.na(p_values[,1]),1]
  adjust_p_value <- p.adjust(p_values)
  PBM.scores <- PBM.scores[adjust_p_value<0.05,]

  counts_signature <- 0

  for (mutation.type in mutation.type.list) {
    histogram_plot_mutation <- graphics::hist(PBM.scores[which(PBM.scores$mutclass ==
                                                                 mutation.type), 1], breaks = seq(-number, number, 0.5), plot = F)
    counts_temp <- histogram_plot_mutation$counts * mutation.spectrum[mutation.type,
                                                                      2]
    counts_signature <- counts_signature + counts_temp
  }
  #######################################################################
  histogram_plot_raw <- hist(PBM.scores[, 1], breaks = seq(-number, number, 0.5), plot = F)
  neg_index <- which(histogram_plot_raw$mids < 0) ##index of positive scores
  pos_index <- which(histogram_plot_raw$mids > 0) ##index of negative scores
  neg.sum.raw <- sum(histogram_plot_raw$mids[neg_index] * histogram_plot_raw$counts[neg_index])/sum(histogram_plot_raw$counts[neg_index])
  pos.sum.raw <- sum(histogram_plot_raw$mids[pos_index] * histogram_plot_raw$counts[pos_index])/sum(histogram_plot_raw$counts[pos_index])
  normalized.counts_signature <- sum(histogram_plot_raw$counts) * counts_signature/sum(counts_signature)  #normalization
  histogram_plot_weighted <- histogram_plot_raw
  neg.sum.weighted <- sum(histogram_plot_weighted$mids[neg_index] * normalized.counts_signature[neg_index])/sum(normalized.counts_signature[neg_index])
  pos.sum.weighted <- sum(histogram_plot_weighted$mids[pos_index] * normalized.counts_signature[pos_index])/sum(normalized.counts_signature[pos_index])
  list <- c(pos.sum.weighted, neg.sum.weighted)/c(pos.sum.raw, neg.sum.raw)
  return(round(list, digits = 3))
}
