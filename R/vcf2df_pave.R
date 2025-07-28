vcf2df_pave <- function(vcf, tumor_id = vcf_tumor_id, normal_id = vcf_normal_id, vcf_tumor_id = "TUMOR", vcf_normal_id = "NORMAL", debug_mode = FALSE, verbose = TRUE) {

  # Assertions
  assertions::assert_string(vcf)
  assertions::assert_file_exists(vcf)
  assertions::assert_string(tumor_id)
  assertions::assert_string(normal_id)
  assertions::assert_string(vcf_tumor_id)
  assertions::assert_string(vcf_normal_id)
  assertions::assert_flag(debug_mode)
  assertions::assert_flag(verbose)

  # Read VCF
  if(verbose) cli::cli_h1(text = "Reading VCF")
  vcfR <- vcfR::read.vcfR(vcf, verbose = FALSE) # May need to increase default limit
  if(verbose) cli::cli_alert_success("VCF successfully read")

  # Test VCF is vep-annotated
  if(verbose) cli::cli_h1("Checking VCF is PAVE-annotated")
  pave_in_meta <- any(grepl(x = vcfR@meta, pattern = "^##PaveVersion"))

  if(verbose) cli::cli_progress_step(msg = "Looking for ##Pave entry in VCF header")

  assertions::assert(
    pave_in_meta,
    msg = "Failed to find PAVE annotation step in VCF header (##Pave). Are you sure the VCF is Pave-annotated?"
  )

  if(verbose) cli::cli_progress_step(msg = "Checking at least some of our entries have IMPACT field in INFO column")
  vep_impact <- vcfR::extract.info(x = vcfR, element = "IMPACT")
  impact_info_field_present <- !all(is.na(vep_impact))

  assertions::assert(
    impact_info_field_present,
    msg = "Could not find IMPACT field in INFO column of VCF. Are you sure the VCF is PAVE-annotated?"
  )

  # Convert to data.frame
  if(verbose)  cli::cli_h1("Converting VCF to data.frame")
  if(verbose)  cli::cli_progress_step("Converting to dataframe")

  ls_tidy_vcf <- vcfR::vcfR2tidy(vcfR, single_frame = TRUE, verbose = FALSE)
  if (debug_mode) {
    return(ls_tidy_vcf)
  }
  df_vcf <- ls_tidy_vcf[["dat"]]
  df_meta <- ls_tidy_vcf[["meta"]]


  # Rename columns to look more like expected input to df2maf
  data.table::setnames(df_vcf, old = c("CHROM", "POS", "REF", "ALT", "IMPACT", "Indiv"), new = c("chr", "pos", "ref", "alt", "consequence", "sample"))

  # Add variant ID column
  df_vcf[["variant_id"]] <- paste0(df_vcf[["chr"]], ":", df_vcf[["pos"]], " ", df_vcf[["ref"]], ">", df_vcf[["alt"]])

  # Add a column differentiating Tumours from Normals
  df_vcf <- df_vcf |>
    dplyr::mutate(sample_type = dplyr::case_when(
      sample == vcf_tumor_id ~ "Somatic",
      sample == vcf_normal_id ~ "Normal",
      .default = "ERROR"
    ))

  unexpected_sample_identifiers <- unique(df_vcf[["sample"]][df_vcf[["sample_type"]] == "ERROR"])
  assertions::assert(
    length(unexpected_sample_identifiers) == 0, msg =
      "The sample identifiers in your VCF/s [{unexpected_sample_identifiers}] do {.strong NOT} match the sample identifiers vcf2maf has been configured to look for ([{.strong {vcf_tumor_id}}] for tumor samples and [{.strong {vcf_normal_id}}] for normal samples).
    Fix this error by changing {.arg vcf_tumor_id} and {.arg vcf_normal_id} arguments the Tumor / Normal sample identifiers present in your VCF"
  )


  # Remove normal variants (technically we should pull out the normal sample ref/alt but in somatic MAFs we remove those anyway for privacy reasons, so we'll ignore
  df_vcf_somatic <- df_vcf |>
    dplyr::filter(sample_type == "Somatic")
  #dplyr::rename("Tumor_Seq_Allele1" = ref, "Tumor_Seq_Allele2" = alt, "Tumor_Sample_Barcode" = sample)

  # Add Tumor and Normal IDs
  df_vcf_somatic[['Matched_Norm_Sample_Barcode']] <- normal_id
  df_vcf_somatic[['Tumor_Sample_Barcode']] <- tumor_id


  # The final challenge is to filter multiple vep consequences -> a single most important consequence.
  # vcf2maf.pl takes the approach of ranking first by transcript biotype, then by severity, and then by longest transcript (https://github.com/mskcc/vcf2maf/blob/main/vcf2maf.pl).
  # Note Transcript_Length isn't separately reported, but can be parsed out from cDNA_position
  # Notes that multiple vep consequences based on different transcripts are separated by commas ','

  # Since pave doesn't report biotype, vcf2mafR will simply report variant from the 'IMPACT' field (Describes Canonical transcript affects)

  impact_meta <- unname(unlist(df_meta[df_meta[["ID"]]=='IMPACT','Description']))
  impact_fieldnames_string <- impact_meta |>
    sub(x=_, pattern = "^Variant Impact \\[", replacement = "") |>
    sub(x=_, pattern = "\\]$", replacement = "")

  impact_fieldnames <- unlist(strsplit(impact_fieldnames_string, split = ", "))

  # Check PAVE annotations include all the columns we need: SYMBOL, BIOTYPE, cDNA_position, Consequence, CANONICAL
  assertions::assert_includes(
    impact_fieldnames, "Gene", msg = 'Failed to find column [{.strong Gene}] in pave annotations.'
  )

  assertions::assert_includes(
    impact_fieldnames, "CanonicalEffect", msg = 'Failed to find column [{.strong CanonicalEffect}] in pave annotations'
  )

  # browser()


  # Fields we need to pull per variant to choose the most significant consequence: BIOTYPE, cDNA_position, Consequence
  ls_consequences = strsplit(df_vcf_somatic[['consequence']], split = ",")
  names(ls_consequences) <- df_vcf_somatic[['variant_id']]

  # df_consequences <- do.call(rbind, ls_consequences)

  ls_consequences_normalised <- lapply(ls_consequences, \(consequences){

    if(length(consequences == 1) & all(is.na(consequences))){
      consequences <- rep(NA_character_, times = length(impact_fieldnames))
    }

    names(consequences) <- impact_fieldnames
    return(consequences)
  })

  df_consequences <- as.data.frame(do.call("rbind", ls_consequences_normalised))
  df_consequences[["variant_id"]] <- rownames(df_consequences)
  rownames(df_consequences) <- NULL


  # Rename df_consequences to include maf_gene and maf_consequence
  colnames(df_consequences) <- ifelse(colnames(df_consequences) == "Gene", "maf_gene", colnames(df_consequences))
  colnames(df_consequences) <- ifelse(colnames(df_consequences) == "CanonicalEffect", "maf_effect", colnames(df_consequences))

  # If pave effects are NA or ""
  browser()
  #intergenic_variant
  df_maf <- df_vcf_somatic |>
    dplyr::left_join(df_consequences, by = "variant_id")

  if(verbose) cli::cli_alert_success("{.strong vcf2df successful}")

  return(df_maf)
}
