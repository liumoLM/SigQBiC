context("GenerateQBiCScoresFromTwelvemers")

test_that("GenerateQBiCScoresFromTwelvemers", {
  uPBM.scores <- fread("./testdata/prediction6mer.Mus_musculus_M01396_1.94d_Berger08_Arx_1738.2.txt",fill=T)

  expected_scores <- fread("./testdata/example.twelvemer.QBiC.scores.txt",header=T)

  expected_scores <- data.frame(expected_scores)

  example_twelvemers <- unlist(fread("./testdata/example.twelvemer.txt",header=F))

  expect_equal(GenerateQBiCScoresFromTwelvemers(uPBM.scores,example_twelvemers), expected_scores)
})

