# perl vcf2maf.pl --input-vcf WD4086.vcf --output-maf WD4086.maf --tumor-id WD4086 --normal-id NB4086

#' Convert a vepped VCF into a table ready to be passed into df2maf
#'
#'
#' @param vcf path to a VCF file
#' @param tumor_id desired value of Tumor_Sample_Barcode in maf (string)
#' @param normal_id desired value of Matched_Norm_Sample_Barcode in maf (string)
#'
#' @return data.frame
#' @export
#'
#' @examples
#' path_vcf <- system.file(package = "vcf2mafR", "testfiles/test_grch38.vep.vcf")
#' vcf2df(path_vcf)
vcf2df <- function(vcf, tumor_id = "Sample1Tumor", normal_id = "Sample1Normal",  debug_mode = TRUE){

  # Read VCF
  cli::cli_h1(text = "Reading VCF")
  vcfR <- vcfR::read.vcfR(vcf) # May need to increase default limit
  cli::cli_alert_success("VCF successfully read")

  # Test VCF is vep-annotated
  cli::cli_h1("Checking VCF is VEP-annotated")
  vep_in_meta <- any(grepl(x = vcfR@meta, pattern = "^##VEP="))

  cli::cli_progress_step(msg = "Looking for ##VEP entry in VCF header")
  assertions::assert(
    vep_in_meta,
    msg = "Failed to find VEP annotation step in VCF header (##VEP). Are you sure the VCF is VEP-annotated?"
    )

  cli::cli_progress_step(msg = "Checking at least some of our entries have CSQ field in INFO column")
  vep_csq <- vcfR::extract.info(x = vcfR, element = "CSQ")
  csq_info_field_present <- !all(is.na(vep_csq))

  assertions::assert(
    csq_info_field_present,
    msg = "Could not find CSQ field in INFO column of VCF. Are you sure the VCF is VEP-annotated?" # May need to extend to allow ANN fields from older VEP versions
  )



  # Convert to data.frame
  cli::cli_h1("Converting VCF to data.frame")
  cli::cli_progress_step("Converting to dataframe")

  ls_tidy_vcf <- vcfR::vcfR2tidy(vcfR, single_frame = TRUE)
  if(debug_mode) return(ls_tidy_vcf)
  df_vcf <- ls_tidy_vcf[["dat"]]

  df_meta <- ls_tidy_vcf[["meta"]]

  cli::cli_alert_success("{.strong vcf2df uccessful}")

  return(df_vcf)
}
