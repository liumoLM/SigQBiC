context("SignatureQBiC")

test_that("SignatureQBiC Model", {

  sigProfiler_SBS_sig <- fread("./testdata/sigProfiler_SBS_signatures_2019_05_22.csv")

  sigProfiler_SBS_sig <- data.frame(sigProfiler_SBS_sig)

  sigProfiler_SBS_sig$mutclass <- apply(sigProfiler_SBS_sig,1,function(x){

    x["mutclass"] <- paste(x["SubType"],substring(x["Type"],3,3),sep="")

  })

  spectrum <- data.frame(sigProfiler_SBS_sig$SBS7a)

  row.names(spectrum) <- mutation.type.list

  uPBM.scores <- readRDS("./testdata/prediction6mer.Mus_musculus_M01396_1.94d_Berger08_Arx_1738.2.rds")


  expect_equal(SignatureQBiC(uPBM.scores,spectrum), c(1.668,0.410))
})
