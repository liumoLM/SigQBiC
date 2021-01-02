context("SignatureQBiC")

test_that("SignatureQBiC Model", {

  load("testdata/SignatureQBiC.expected.Rdata")
  QBiC_score_file_path =
    "../../data-raw/prediction6mer.Homo_sapiens!M01252_1.94d!Barrera2016!HOXD13_I322L_R1.txt.gz"
  pvalue_file_path =
    "../../data-raw/pval6mer.Homo_sapiens!M01252_1.94d!Barrera2016!HOXD13_I322L_R1.csv.gz"
  PCAWG_subs_signature <- PCAWG7::signature$genome$SBS96
  row.names(PCAWG_subs_signature) <- mut.types
  sig <- PCAWG_subs_signature[,"SBS7a",drop=F]

  retval <- SignatureQBiC(QBiC_score_file_path,
                               pvalue_file_path,
                               sig,
                               plot.path = NULL)
  #save(test_retval,file="./tests/testthat/testdata/SignatureQBiC.expected.Rdata")

  expect_equal(retval, test_retval)
})
