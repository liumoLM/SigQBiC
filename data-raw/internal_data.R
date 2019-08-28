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




output$mutclass <- apply(output,1,function(x){

  if(substring(x["seq"],6,6)=="A"||substring(x["seq"],6,6)=="G"){

    x["mutclass"] <- paste(reverseComplement(substring(x["seq"], 5, 7),case="upper"),
                           reverseComplement(substring(x["seq"], 12, 12),case="upper"),
                           sep = "")
  }else{

    x["mutclass"] <- paste(substring(x["seq"], 5, 7),
                           substring(x["seq"], 12, 12),
                           sep="")
  }

  return(x["mutclass"])
})

all.possible.twelvemers <- output[output$sum==0,c(3,7)]




mutation.type.list <-
  c("ACAA", "ACCA", "ACGA", "ACTA", "CCAA", "CCCA", "CCGA", "CCTA",
    "GCAA", "GCCA", "GCGA", "GCTA", "TCAA", "TCCA", "TCGA", "TCTA",
    "ACAG", "ACCG", "ACGG", "ACTG", "CCAG", "CCCG", "CCGG", "CCTG",
    "GCAG", "GCCG", "GCGG", "GCTG", "TCAG", "TCCG", "TCGG", "TCTG",
    "ACAT", "ACCT", "ACGT", "ACTT", "CCAT", "CCCT", "CCGT", "CCTT",
    "GCAT", "GCCT", "GCGT", "GCTT", "TCAT", "TCCT", "TCGT", "TCTT",
    "ATAA", "ATCA", "ATGA", "ATTA", "CTAA", "CTCA", "CTGA", "CTTA",
    "GTAA", "GTCA", "GTGA", "GTTA", "TTAA", "TTCA", "TTGA", "TTTA",
    "ATAC", "ATCC", "ATGC", "ATTC", "CTAC", "CTCC", "CTGC", "CTTC",
    "GTAC", "GTCC", "GTGC", "GTTC", "TTAC", "TTCC", "TTGC", "TTTC",
    "ATAG", "ATCG", "ATGG", "ATTG", "CTAG", "CTCG", "CTGG", "CTTG",
    "GTAG", "GTCG", "GTGG", "GTTG", "TTAG", "TTCG", "TTGG", "TTTG"
  )




usethis::use_data(all.possible.twelvemers,mutation.type.list,internal=TRUE,overwrite = T)


