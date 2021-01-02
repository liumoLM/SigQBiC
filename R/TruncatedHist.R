#'Plot histograms with the long tails collected into a single bar
#'
#'@description This is an internal function called by SignatureQBiC. It takes
#'
#'@param all.QBiC.scores A numeric vector containing all QBiC scores for a universal PBM.
#'
#'@param original.scores A numeric vector containing all QBiC scores for the target mutation type
#'
#'@param weighted.prop A numeric variable. The probability of the target mutation type in the mutational signature
#'
#'@param mutation.type A character. E.g. "ACA_ATA"
#'
#'@return A list containing GR and LR.
#'


TruncatedHist <- function(all.QBiC.scores,original.scores,weighted.prop,mutation.type){
  cut_off <- quantile(all.QBiC.scores,seq(0,1,0.01))[100] ##pile the 1% tail up
  original.scores[original.scores>cut_off] <- cut_off
  original.scores[original.scores<(-cut_off)] <- (-cut_off)
  weighted.hist <- original.hist <- hist(original.scores, breaks = seq(-(cut_off+0.5),
                                                                       (cut_off+0.5),0.01), plot=F)
  weighted.hist$density <-  weighted.hist$density*weighted.prop
  plot(original.hist,freq = F,ylim=c(0,max(original.hist$density)+0.05),main=paste0("Original",mutation.type),xaxt="n",yaxt="n")
  plot(weighted.hist,freq=F,ylim = c(0,max(original.hist$density)+0.05),main=paste0("Weighted ",mutation.type),xaxt="n",yaxt="n")
}
