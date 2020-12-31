#'@title The resampling test generate 1000 flat signatures, to test if the GR or LR>1 is statistically significant.
#'
#'@description This function tests whether $D'_{Pos}$ ($D'_{Neg}$) are statistically > than $D_{Pos}$ ($D_{Neg}$)
#'
#'@param i Seednumber
#'
#'@return A list of frequencies of each mutation type with sum to 1
#'
#'@export
#'
##

ResampleMutationFrequency <- function(i){
  set.seed(i)
  resampling.of.mut.type <- table(sample(c(1:96),size=nrow(all.possible.twelvemers),replace=T))
  names(resampling.of.mut.type) <- mut.types
  resampling.of.mut.type <- resampling.of.mut.type/sum(resampling.of.mut.type)
  return(resampling.of.mut.type)
}

