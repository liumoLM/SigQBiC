library("gtools")
library(spgs)
DNAToBin <- function(DNA){
  temp <- DNA
  temp <- gsub("A", "00", temp)
  temp <- gsub("C", "01", temp)
  temp <- gsub("G", "10", temp)
  temp <- gsub("T", "11", temp)
  return(temp)
}

# Generate reference and mutant sequence #
dict <- c("A", "C", "G", "T")
n_half <- 5
seq_context <- permutations(4, 2*n_half, dict, repeats.allowed = T)

seq_A <- cbind(seq_context[,1:n_half], "A",
               seq_context[,(n_half+1):(2*n_half)])
seq_A <- apply(seq_A, 1, paste, collapse="")
seq_C <- cbind(seq_context[,1:n_half], "C",
               seq_context[,(n_half+1):(2*n_half)])
seq_C <- apply(seq_C, 1, paste, collapse="")
seq_G <- cbind(seq_context[,1:n_half], "G",
               seq_context[,(n_half+1):(2*n_half)])
seq_G <- apply(seq_G, 1, paste, collapse="")
seq_T <- cbind(seq_context[,1:n_half], "T",
               seq_context[,(n_half+1):(2*n_half)])
seq_T <- apply(seq_T, 1, paste, collapse="")

ref_A <- as.character(rep(seq_A, 4))
mut_A <- as.character(c(seq_A, seq_C, seq_G, seq_T))
ref_C <- as.character(rep(seq_C, 4))
mut_C <- mut_A
ref_G <- as.character(rep(seq_G, 4))
mut_G <- mut_A
ref_T <- as.character(rep(seq_T, 4))
mut_T <- mut_A

ref_all <- c(ref_A, ref_C, ref_G, ref_T)
mut_all <- c(mut_A, mut_C, mut_G, mut_T)

DNA_seq <- paste0(ref_all,
                  substr(mut_all, n_half+1, n_half+1))
DNA_bin <- DNAToBin(DNA_seq) #binary number
DNA_dec <- strtoi(DNA_bin, base = 2L) #decimal number

output <- data.frame(ref_all, mut_all, DNA_seq, DNA_dec)
colnames(output) <- c("ref", "mut", "seq", "num")
output <- output[order(output$num), ]

output$revc.mut <- apply(output,1,function(x){
  x["revc.mut"] <- reverseComplement(x["mut"],case="upper")
})
output$sum <- apply(output,1,function(x){
  x['sum'] <- sum(c(x["revc.mut"]==x["ref"],x["ref"]==x["mut"]))
})




output$final_signature <- apply(output,1,function(x){

  if(substring(x["seq"],6,6)=="A"||substring(x["seq"],6,6)=="G"){

    x["mutclass"] <- paste(reverseComplement(substring(x["seq"], 5, 7),case="upper"),
                           reverseComplement(substring(x["seq"], 12, 12),case="upper"),
                           sep = "")

    x["mutclass"] <- paste(reverseComplement(substring(x["seq"], 5, 7)),
                           paste0(reverseComplement(substring(x["seq"], 7, 7)),
                                  reverseComplement(substring(x["seq"], 12, 12)),
                                  reverseComplement(substring(x["seq"], 5, 5))),
                           sep="_")
  }else{

    x["mutclass"] <- paste(substring(x["seq"], 5, 7),
                           paste0(substring(x["seq"], 5, 5),
                                  substring(x["seq"], 12, 12),
                                  substring(x["seq"], 7, 7)),
                           sep="_")
  }

  return(x["mutclass"])
})

all.possible.twelvemers <- output[output$sum==0,c(3,7)]

mut.types <-
  c("ACA_AAA","ACC_AAC","ACG_AAG","ACT_AAT","CCA_CAA","CCC_CAC","CCG_CAG",
    "CCT_CAT","GCA_GAA","GCC_GAC","GCG_GAG","GCT_GAT","TCA_TAA","TCC_TAC",
    "TCG_TAG","TCT_TAT","ACA_AGA","ACC_AGC","ACG_AGG","ACT_AGT","CCA_CGA",
    "CCC_CGC","CCG_CGG","CCT_CGT","GCA_GGA","GCC_GGC","GCG_GGG","GCT_GGT",
    "TCA_TGA","TCC_TGC","TCG_TGG","TCT_TGT","ACA_ATA","ACC_ATC","ACG_ATG",
    "ACT_ATT","CCA_CTA","CCC_CTC","CCG_CTG","CCT_CTT","GCA_GTA","GCC_GTC",
    "GCG_GTG","GCT_GTT","TCA_TTA","TCC_TTC","TCG_TTG","TCT_TTT","ATA_AAA",
    "ATC_AAC","ATG_AAG","ATT_AAT","CTA_CAA","CTC_CAC","CTG_CAG","CTT_CAT",
    "GTA_GAA","GTC_GAC","GTG_GAG","GTT_GAT","TTA_TAA","TTC_TAC","TTG_TAG",
    "TTT_TAT","ATA_ACA","ATC_ACC","ATG_ACG","ATT_ACT","CTA_CCA","CTC_CCC",
    "CTG_CCG","CTT_CCT","GTA_GCA","GTC_GCC","GTG_GCG","GTT_GCT","TTA_TCA",
    "TTC_TCC","TTG_TCG","TTT_TCT","ATA_AGA","ATC_AGC","ATG_AGG","ATT_AGT",
    "CTA_CGA","CTC_CGC","CTG_CGG","CTT_CGT","GTA_GGA","GTC_GGC","GTG_GGG",
    "GTT_GGT","TTA_TGA","TTC_TGC","TTG_TGG","TTT_TGT")


usethis::use_data(all.possible.twelvemers,mut.types,internal=TRUE,overwrite = T)


