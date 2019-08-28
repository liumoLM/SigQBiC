
context("Generate96ChannelSpectrumFromTwelvemers")


test_that("Generate 96 Channel Spectrum From A Given Twelvemer List", {

  library(data.table)

  library(spgs)

  example_twelvemers <- unlist(fread("./testdata/example.twelvemer.txt",header=F))

  expected_spectrum <- fread("./testdata/expected_spectrum.txt",header=T)

  expected_spectrum <- data.frame(expected_spectrum)

  expect_equal(Generate96ChannelSpectrumFromTwelvemers(example_twelvemers), expected_spectrum)
})
