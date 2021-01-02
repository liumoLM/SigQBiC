#'Generate GR and LR for a TF-signature pair (Signature-QBiC)
#'
#'@description SignatureQBiC model. This function generates Gain Ratio and Loss Ratio
#'            for a TF (represented by universal PBM) with
#'            a mutational signature (or a mutation spectrum)
#'
#'@param QBiC_score_file_path All QBiC scores for a universal PBM. The file can be downloaded from
#'                         http://qbic.genome.duke.edu/downloads
#'
#'@param pvalue_file_path p_values for all twelvemers. The file can be downloaded from
#'                         http://qbic.genome.duke.edu/downloads
#'
#'@param sig A one-column matrix contains probability of each mutation types
#'           with row.names set to SigQBiC::mut_types.
#'
#'@return A list containing GR and LR.
#'


SignatureQBiC <- function(QBiC_score_file_path,
                          pvalue_file_path,
                          sig,
                          plot.path = NULL) {
  QBiC_scores_table <-
    data.table::fread(QBiC_score_file_path,
                      sep=" ", header=T, stringsAsFactors = F, fill = T)
  # This gives a data frame with colums diff and z_score
  # QBiC_scores_table contains NA for non-mutations, e.g AAAAAAAAAAA -> AAAAAAAAAAA



  pvalue <-
    scan(pvalue_file_path)
  # pvalue also contains NA for non-mutations
  pvalue <- pvalue[!is.na(pvalue)]

  QBiC_scores_matrix <-
    tibble(QBiC_mut = all.possible.twelvemers$seq,
           mut_type = all.possible.twelvemers$final_signature,
           scores   = QBiC_scores_table$z_score[!is.na(QBiC_scores_table$z_score)],
           p        = pvalue,
           q        = p.adjust(pvalue, method = "BH"))



  max.score <- as.integer(max(QBiC_scores_matrix$scores)) + 2


  summaryofscores <-data.frame(matrix(ncol=5,nrow=0))
  my.breaks <- seq(-max.score,max.score,0.001)

  if(!is.null(plot.path)){
    all.weighted.freq <- 0
    if (!dir.exists(plot.path)) {
      if (!dir.create(plot.path, recursive = T))
        stop("Cannot create plotting directory ", plot.path)
    }

    png(filename = paste0(plot.path, "/", "hist%03d.png"))
    par(mar = c(3,1,1,1))
    par(mfrow = c(8,4))

    # open.plot <- function(file.name.prefix) {
    #  png(paste0(plot.path, "/", file.name.prefix, ".png"))
    #  par(mar = c(3,1,1,1))
    #  par(mfrow = c(8,4))
    # }
  }

  for (mutation.type in mut.types) {

    stopifnot(mutation.type %in% QBiC_scores_matrix$mut_type)
    # Scores for the given mutation.type
    tmp.scores <-
      QBiC_scores_matrix$scores[QBiC_scores_matrix$mut_type==mutation.type] ##the scores were put into bins

    dist.hist <- hist(tmp.scores, breaks = my.breaks, plot=F)
    w.dist.hist <- dist.hist
    w.dist.hist$counts <- dist.hist$counts * sig[mutation.type, ]

    partial.summary <-
      data.frame(scores         = dist.hist$mids,
                 frequency      = dist.hist$counts,
                 mut_type       = mutation.type,
                 signature_freq = sig[mutation.type, ],
                 weighted.freq  = dist.hist$counts * sig[mutation.type, ])
    ##multiply the counts of each bin by the frequency of mutations in a signature

    if(!is.null(plot.path)){
      # Plot ...
      all.weighted.freq <- all.weighted.freq + dist.hist$counts * sig[mutation.type, ]
      # open.plot(mutation.type)
      TruncatedHist(QBiC_scores_matrix$scores,
                    original.scores = tmp.scores,
                    weighted.prop = sig[mutation.type, ],
                    mutation.type=mutation.type)
      # dev.off()
    }

    summaryofscores <- rbind(summaryofscores, partial.summary)
  }
  if(!is.null(plot.path)){
    all.scores <-  QBiC_scores_matrix$scores
    cut_off <- quantile( all.scores,seq(0,1,0.001))[1000] ##pile the 1% tail up
    all.scores[all.scores>cut_off] <- cut_off
    all.scores[all.scores<(-cut_off)] <- (-cut_off)
    weighted.hist <- original.hist <- hist(all.scores, breaks = my.breaks, plot=F)
    weighted.hist$density <-  all.weighted.freq
    # open.plot("summary")
    plot(original.hist,freq = F,ylim=c(0,max(original.hist$density)+0.05),main = "Original Distribution")
    plot(weighted.hist,freq=F,ylim = c(0,max(original.hist$density)+0.05),main = "Weighted Distribution")
    dev.off()
  }

  pos.sig.QBiC_scores_matrix <-
    QBiC_scores_matrix[QBiC_scores_matrix$q < 0.1 & QBiC_scores_matrix$scores>0,] #select Dpos


  qvalue.cutoff.score <- min(pos.sig.QBiC_scores_matrix$scores) ##get the cutoff of QBiC scores based on BH FDR

  summaryofscores$weighted.freq <-
    summaryofscores$weighted.freq *
    sum(summaryofscores$frequency)/sum(summaryofscores$weighted.freq)
  ##Normalize the weighted freqeuencies. After multiplying with signature probability, the weighted frequency is 96 times less than the original frequency. sum(freq) = 96*sum(weighted.freq)

  summaryofscores.Dpos <-
    summaryofscores[summaryofscores$scores>qvalue.cutoff.score, ] ##Select Dpos

  Dpos <- rep(summaryofscores.Dpos$scores,
              summaryofscores.Dpos$frequency)

  Dprimepos <- rep(summaryofscores.Dpos$scores,
                   round(summaryofscores.Dpos$weighted.freq, digits = 0))


  summaryofscores.Dneg <-
    summaryofscores[summaryofscores$scores<(-qvalue.cutoff.score), ] ##Select Dneg

  Dneg <- rep(summaryofscores.Dneg$scores,
              summaryofscores.Dneg$frequency)

  Dprimeneg <- rep(summaryofscores.Dneg$scores,
                   round(summaryofscores.Dneg$weighted.freq, digits = 0))

  GR = sum(Dprimepos)/sum(Dpos)
  LR = sum(Dprimeneg)/sum(Dneg)
  return(list(GR=GR,
              LR=LR))
}
