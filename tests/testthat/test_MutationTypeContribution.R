context("MutationTypeContribution")

test_that("Compute Contribution of Each Mutation Class to GR and LR", {

  library(data.table)

  sigProfiler_SBS_sig <- fread("./testdata/sigProfiler_SBS_signatures_2019_05_22.csv")

  sigProfiler_SBS_sig <- data.frame(sigProfiler_SBS_sig)
  sigProfiler_SBS_sig$mutclass <- apply(sigProfiler_SBS_sig,1,function(x){

    x["mutclass"] <- paste(x["SubType"],substring(x["Type"],3,3),sep="")

  })

  spectrum <- data.frame(sigProfiler_SBS_sig$SBS7a)

  row.names(spectrum) <- sigProfiler_SBS_sig$mutclass

  example_contribution <- fread("./testdata/example_mutation_contribution_to_SBS7a.txt",header=T)

  example_contribution <- data.frame(example_contribution)

  uPBM.scores <- fread("./testdata/prediction6mer.Mus_musculus_M01396_1.94d_Berger08_Arx_1738.2.txt",fill=T)

  expect_equal(MutationTypeContribution(uPBM.scores,spectrum), example_contribution)
})
