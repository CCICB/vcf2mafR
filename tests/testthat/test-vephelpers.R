test_that("cDNA position to PseudoTranscriptLength Works", {
  cDNA = c("100/1731", "100_100/2783", "100/1084", "100/1476", "100/1173", "100/1480")
  expected = c(1731, 2783, 1084, 1476, 1173, 1480)

  expect_equal(cdna_to_transcript_length(cDNA), expected)

  expect_equal(cdna_to_transcript_length(c("213/100", "")), c(100, 0))
  expect_equal(cdna_to_transcript_length(c("213/100", "")), c(100, 0))
  expect_equal(cdna_to_transcript_length(c("100", NA)), c(0, 0))
  expect_equal(cdna_to_transcript_length(c("STRING")), c(0))

})
