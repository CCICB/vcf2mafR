test_that("df2maf works with SO consequences", {
  df_maflike_so = read.csv(system.file(package = "vcf2mafR", "testfiles/test_so.tsv"), sep = "\t")

  # Throw no error with SO consequences
  df_maf = expect_error(df2maf(df_maflike_so, ref_genome = "hg38") |> suppressMessages(), regexp = NA)

  # Throw no error with automatic consequence dictionary guessing
  expect_error(df2maf(df_maflike_so, ref_genome = "hg38", consequence_dictionary = "AUTO") |> suppressMessages(), regexp = NA)

  # Returns a dataframe
  expect_s3_class(df_maf, class = "data.frame")

  # Number of rows should match input
  expect_equal(nrow(df_maf), nrow(df_maflike_so))

  # Throw error when forcing consequence dictionary search to be PAVE
  expect_error(df2maf(df_maflike_so, ref_genome = "hg38", consequence_dictionary = "PAVE") |> suppressMessages(), regexp = "NOT a valid PAVE term")

})


test_that("df2maf works with PAVE consequences", {
  df_maflike_pave = read.csv(system.file(package = "vcf2mafR", "testfiles/test_pave.tsv"), sep = "\t")

  # Throw no error with PAVE consequences
  df_maf = expect_error(df2maf(df_maflike_pave, ref_genome = "hg38", consequence_dictionary = "PAVE") |> suppressMessages(), regexp = NA)

  # Throw no error with automatic consequence dictionary guessing
  expect_error(df2maf(df_maflike_pave, ref_genome = "hg38", consequence_dictionary = "AUTO") |> suppressMessages(), regexp = NA)

  # Returns a dataframe
  expect_s3_class(df_maf, class = "data.frame")

  # Number of rows should match input
  expect_equal(nrow(df_maf), nrow(df_maflike_pave))

  # Throw error when forcing consequence dictionary search to be SO
  expect_error(df2maf(df_maflike_pave, ref_genome = "hg38", consequence_dictionary = "SO") |> suppressMessages(), regexp = "were NOT valid SO terms", fixed = TRUE)

})


test_that("df2maf missing_to_silent flag works as expected", {

  ### ---- PAVE ---- ###
  df_maflike_pave = read.csv(system.file(package = "vcf2mafR", "testfiles/test_pave.tsv"), sep = "\t")
  df_maflike_pave[['consequence']][1] <- ""

  # Throw error when empty string is present
  expect_error(df2maf(df_maflike_pave, ref_genome = "hg38"), "Found 1 variant with no mutation type")

  # Pass
  df <- expect_error(
    df2maf(df_maflike_pave, ref_genome = "hg38", missing_to_silent = TRUE, consequence_dictionary = "PAVE") |> suppressMessages(),
    NA
  )
  expect_equal(df[['Variant_Classification']][1], expected = 'Silent')

  # Repeat but for NA consequence instead of ""
  df_maflike_pave[['consequence']][1] <- NA_character_

  # Throw error when empty string is present
  expect_error(df2maf(df_maflike_pave, ref_genome = "hg38"), "no missing values! Found 1")

  # Pass
  df <- expect_error(
    df2maf(df_maflike_pave, ref_genome = "hg38", missing_to_silent = TRUE, consequence_dictionary = "PAVE") |> suppressMessages(),
    NA
  )
  expect_equal(df[['Variant_Classification']][1], expected = 'Silent')

  ### ---- SO ---- ###
  df_maflike_so = read.csv(system.file(package = "vcf2mafR", "testfiles/test_so.tsv"), sep = "\t")
  df_maflike_so[['consequence']][1] <- ""

  # Throw error when empty string is present
  expect_error(df2maf(df_maflike_so, ref_genome = "hg38"), "Found 1 variant with no mutation type")

  # Pass
  df <- expect_error(
    df2maf(df_maflike_so, ref_genome = "hg38", missing_to_silent = TRUE, consequence_dictionary = "SO") |> suppressMessages(),
    NA
  )
  expect_equal(df[['Variant_Classification']][1], expected = 'Silent')

  # Repeat but for NA consequence instead of ""
  df_maflike_so[['consequence']][1] <- NA_character_

  # Throw error when empty string is present
  expect_error(df2maf(df_maflike_so, ref_genome = "hg38"), "no missing values! Found 1")

  # Pass
  df <- expect_error(
    df2maf(df_maflike_so, ref_genome = "hg38", missing_to_silent = TRUE, consequence_dictionary = "SO") |> suppressMessages(),
    NA
  )
  expect_equal(df[['Variant_Classification']][1], expected = 'Silent')

})

