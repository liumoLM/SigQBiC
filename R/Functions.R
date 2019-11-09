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
Generate96ChannelSpectrumFromTwelvemers <- function(twelvemers) {
  twelvemers <- data.frame(twelvemers)
  twelvemers$ref <- substring(twelvemers[, 1], 6, 6)
  for (i in 1:nrow(twelvemers)) {
    if (twelvemers$ref[i] == "G" || twelvemers$ref[i] == "A") {
      twelvemers$mutclass[i] <- paste(spgs::reverseComplement(substring(twelvemers[i,
                                                                                             1], 5, 7), case = "upper"), spgs::reverseComplement(substring(twelvemers[i,
                                                                                                                                                                           1], 12, 12), case = "upper"), sep = "")
    } else {
      twelvemers$mutclass[i] <- paste(substring(twelvemers[i, 1],
                                                     5, 7), substring(twelvemers[i, 1], 12, 12), sep = "")
    }
  }

  mutation.spectrum <- data.frame(table(twelvemers$mutclass))
  mutation.spectrum[, 2] <- as.numeric(mutation.spectrum[, 2])/sum(as.numeric(mutation.spectrum[,
                                                                                                2]))

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

}


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


#'Plot histogram for each mutation class for a given signature and uPBM experiment name
#'
#'@description
#'
#'@param uPBM_QBiC_scores All QBiC scores for a universal PBM
#'
#'@param spectrum A list of GR/LR generated by Signature-QBiC from observed mutation spectrum
#'
#'@return A plot of histograms
#'
#'@export
#'
GenerateSignatureWeightedDistributionPlots <- function(uPBM_QBiC_scores, spectrum,
                                                            output_name) {

  uPBM_QBiC_scores <- data.frame(uPBM_QBiC_scores)
  uPBM_QBiC_scores <- uPBM_QBiC_scores$z_score[which(!is.na(uPBM_QBiC_scores$z_score))]
  PBM.scores <- data.frame(uPBM_QBiC_scores)

  row.names(PBM.scores) <- all.possible.twelvemers$seq
  PBM.scores$mutclass <- all.possible.twelvemers$mutclass
  number <- as.integer(max(PBM.scores[, 1])) + 2
  densities.weighted <- 0

  mutation.spectrum <- check.spectrum(spectrum)

  pdf(output_name, width = 16, height = 32)
  par(mfrow = c(8, 4))
  for (mutation.type in mutation.type.list) {
    histogram_plot <- hist(PBM.scores[which(PBM.scores$mutclass ==
                                                         mutation.type), 1], breaks = seq(-number, number, 0.5), plot = F)
    density_temp <- histogram_plot$density * mutation.spectrum[mutation.type, 2]
    densities.weighted <- densities.weighted + density_temp
    plot(histogram_plot, main = mutation.type, col = "black", freq = F, ylim = c(0,
                                                                                 ylimits))
    histogram_plot$density <- density_temp
    plot(histogram_plot, main = mutation.type, col = "grey", freq = F, ylim = c(0,
                                                                                ylimits))
  }
  histogram_plot <- hist(PBM.scores[, 1], breaks = seq(-number, number,
                                                                   0.5), plot = F)

  ylimits <- max(c(histogram_plot$density, densities.weighted)) + 0.05
  plot(histogram_plot, col = "black", freq = F, ylim = c(0, ylimits))
  histogram_plot$density <- densities.weighted
  plot(histogram_plot, col = "grey", freq = F, ylim = c(0, ylimits))
  dev.off()
  return(NULL)
}
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
#'@export
#'
Pathway_Selection <- function(TF.list, qvalue, dbs.selected) {
  enriched <- enrichr(TF.list, dbs.selected)
  enriched.matrix <- enriched[[dbs.selected]]
  enriched.matrix$TF.counts <- length(unique(TF.list))
  enriched.matrix$TFs <- apply(enriched.matrix.gain, 1, function(x) {
    x["TFs"] <- unlist(strsplit(x["Overlap"], "/"))[1]
  })
  enriched.matrix <- enriched.matrix[enriched.matrix$Adjusted.P.value < qvalue,
                                     ]
  enriched.matrix$New.Term <- apply(enriched.matrix, 1, function(x) {
    x["New.Term"] <- unlist(strsplit(x["Term"], "[_]", perl = T))[1]
  })
  return(enriched.matrix)
}



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


#'An internal function to transit spectrum into 96 channel
#'
#'@param spectrum A list of possibilities for each mutation type(mutational spectrum) with mutation types as rownames
#'
#'@return A matrix with contribution to the gain-binding and loss-binding for each of 96 mutation types
#'
#'@keywords internal

check.spectrum <- function(spectrum) {

  spectrum.96.channel <- data.frame(mutation.type.list)

  spectrum.96.channel$possibility <- 0

  spectrum.96.channel$possibility <- spectrum[match(row.names(spectrum), spectrum.96.channel[,
                                                                                             1]), 1]

  if (sum(is.na(spectrum.96.channel$possibility)) > 0) {
    spectrum.96.channel$possibility[which(is.na(spectrum.96.channel$possibility))] <- 0
  }

  row.names(spectrum.96.channel) <- mutation.type.list
  return(spectrum.96.channel)

}
