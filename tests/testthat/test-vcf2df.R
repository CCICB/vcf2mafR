test_that("vcf2df works", {

  # Throw error if not vepped
  path_vcf_unannotated <- system.file("testfiles/test_b38.vcf", package = "vcf2mafR")
  expect_error(vcf2maf(vcf = path_vcf_unannotated, verbose = FALSE), "Are you sure the VCF is VEP-annotated")

  # Throw error if transcript CANONICAL status not annotated by VEP
  path_vcf_annotated_no_canonical <- system.file("testfiles/test_b38.vepgui.canonical_not_included.vcf", package = "vcf2mafR")
  expect_error(vcf2maf(path_vcf_annotated_no_canonical, verbose = FALSE), regexp = "Failed to find column [CANONICAL] in vep annotations. Please rerun VEP with `Identify Canonical Transcripts` option turned on", fixed = TRUE)


  # Throw error if transcript CANONICAL status not annotated by VEP
  path_vcf_vepped <- system.file("testfiles/test_b38.vepgui.vcf", package = "vcf2mafR")
  df <- expect_error(vcf2maf(path_vcf_vepped, ref_genome = 'hg38', verbose = FALSE), regexp = NA)

  ## Check outputs make sense
  # Tumor_Sample_Barcode Makes Sense
  expect_equal(unique(df$Tumor_Sample_Barcode), "TUMOR")

  # Required Columns Present
  expect_error(assertions::assert_names_include(df, c('Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode', 'Hugo_Symbol', "Reference_Allele", 'Tumor_Seq_Allele2', 'Ref_Length', 'Alt_Length',  'Chromosome', 'Start_Position', 'End_Position')), regexp = NA)

  # Snapshot test
  expect_snapshot(df)
})
