all.possible.twelvemers <- readRDS("/Users/Steve Lab-Mo/Downloads/all.possible.twelvemers.rds")
devtools::use_data(all.possible.twelvemers,all.possible.twelvemers,internal=TRUE)


.canonical.96.row.order <-
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
devtools::use_data(.canonical.96.row.order,mutation.type.list,internal=TRUE)


