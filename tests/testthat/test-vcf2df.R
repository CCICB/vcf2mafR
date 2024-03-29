test_that("vcf2df works", {

  # Throw error if not vepped
  path_vcf_unannotated <- system.file("testfiles/test_b38.vcf", package = "vcf2mafR")
  expect_error(vcf2maf(vcf = path_vcf_unannotated) |> suppressMessages(), "Are you sure the VCF is VEP-annotated")

  # Throw error if transcript CANONICAL status not annotated by VEP
  path_vcf_annotated_no_canonical <- system.file("testfiles/test_b38.vepgui.canonical_not_included.vcf", package = "vcf2mafR")
  expect_error(vcf2maf(path_vcf_annotated_no_canonical) |> suppressMessages(), regexp = "Failed to find column [CANONICAL] in vep annotations. Please rerun VEP with `Identify Canonical Transcripts` option turned on", fixed = TRUE)


  # Throw error if transcript CANONICAL status not annotated by VEP
  path_vcf_vepped <- system.file("testfiles/test_b38.vepgui.vcf", package = "vcf2mafR")
  df <- expect_error(vcf2maf(path_vcf_vepped, ref_genome = 'hg38') |> suppressMessages(), regexp = NA)

  ## Check outputs make sense
  # Tumor_Sample_Barcode Makes Sense
  expect_equal(unique(df$Tumor_Sample_Barcode), "TUMOR")

  # Required Columns Present
  expect_error(assertions::assert_names_include(df, c('Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode', 'Hugo_Symbol', "Reference_Allele", 'Tumor_Seq_Allele2', 'Ref_Length', 'Alt_Length',  'Chromosome', 'Start_Position', 'End_Position')), regexp = NA)

  # Snapshot test
  expect_snapshot(df)
})

test_that("vcfs2maf works", {

  # Setup
  vcf_filepaths = dir(system.file(package='vcf2mafR', 'testfiles/cohort_of_vcfs/'), full.names = TRUE)

  # Works without error on valid inputs
  expect_error(vcfs2maf(vcf_filepaths, ref_genome = "hg38") |> suppressMessages(), NA)

  # Throws Error If Generic VCF sample names from different VCFs are at risk of being lumped into the same 'sample' in the MAF
  expect_error(vcfs2maf(vcfs = vcf_filepaths, tumor_id = c('a', 'b', 'c', 'd', 'd'), parse_tumor_id_from_filename = FALSE, ref_genome = "hg38") |> suppressMessages(), "ound duplicated")

  # (Same as above but due to duplicated filenames)
  expect_error(vcfs2maf(vcfs = c(vcf_filepaths, vcf_filepaths[1]), parse_tumor_id_from_filename = TRUE, ref_genome = "hg38") |> suppressMessages(), "Attempt to parse sample")


  # Works when vcf_tumor_id is manually specified as a vector
  vcf_tumor_ids <- rep('TUMOR', times=5)

  expect_error(
    vcfs2maf(
      vcfs = vcf_filepaths,
      parse_tumor_id_from_filename = TRUE,
      vcf_tumor_id = vcf_tumor_ids,
      ref_genome = "hg38") |> suppressMessages(),
    regexp = NA
    )

  expect_error(
    vcfs2maf(
      vcfs = vcf_filepaths,
      parse_tumor_id_from_filename = TRUE,
      vcf_tumor_id = vcf_tumor_ids[1:4],
      ref_genome = "hg38") |> suppressMessages(),
    regexp = "Mismatch between length of"
  )




})
