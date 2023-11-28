test_that("df2maf works with SO consequences", {
  df_maflike_so = read.csv(system.file(package = "vcf2mafR", "testfiles/test_so.tsv"), sep = "\t")

  # Throw no error with SO consequences
  df_maf = expect_error(df2maf(df_maflike_so, ref_genome = "hg38"), regexp = NA)

  # Throw no error with automatic consequence dictionary guessing
  expect_error(df2maf(df_maflike_so, ref_genome = "hg38", consequence_dictionary = "AUTO"), regexp = NA)

  # Returns a dataframe
  expect_s3_class(df_maf, class = "data.frame")

  # Number of rows should match input
  expect_equal(nrow(df_maf), nrow(df_maflike_so))

  # Throw error when forcing consequence dictionary search to be PAVE
  expect_error(df2maf(df_maflike_so, ref_genome = "hg38", consequence_dictionary = "PAVE"), regexp = "NOT a valid PAVE term")

})


test_that("df2maf works with PAVE consequences", {
  df_maflike_pave = read.csv(system.file(package = "vcf2mafR", "testfiles/test_pave.tsv"), sep = "\t")


  # Throw no error with PAVE consequences
  df_maf = expect_error(df2maf(df_maflike_pave, ref_genome = "hg38", consequence_dictionary = "PAVE"), regexp = NA)

  # Throw no error with automatic consequence dictionary guessing
  expect_error(df2maf(df_maflike_pave, ref_genome = "hg38", consequence_dictionary = "AUTO"), regexp = NA)

  # Returns a dataframe
  expect_s3_class(df_maf, class = "data.frame")

  # Number of rows should match input
  expect_equal(nrow(df_maf), nrow(df_maflike_pave))

  # Throw error when forcing consequence dictionary search to be SO
  expect_error(df2maf(df_maflike_pave, ref_genome = "hg38", consequence_dictionary = "SO"), regexp = "NOT valid SO terms")

})
