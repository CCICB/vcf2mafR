# perl vcf2maf.pl --input-vcf WD4086.vcf --output-maf WD4086.maf --tumor-id WD4086 --normal-id NB4086

#' Convert a vepped VCF into a table ready to be passed into df2maf
#'
#' Expects a vepped VCF. VEP must be run with a couple of options:
#' 1) Identifiers: Gene symbol & Transcript Version & (HGVS: optional)
#' 2) Transcript annotation: Transcript biotype & 'Identify canonical transcripts'
#' 3) (Optional) If using the commandline version of vep also make sure to use --total_length (helps in transcript choice)
#'
#' There are a couple of different types of inputs we supp
#'
#' @param vcf path to a VCF file
#' @param tumor_id desired value of Tumor_Sample_Barcode in maf (string)
#' @param normal_id desired value of Matched_Norm_Sample_Barcode in maf (string)
#' @param vcf_tumor_id the sample ID describing the tumor in the (string)
#' @param vcf_normal_id the sample ID describing the normal in the (string)
#'
#' @return data.frame
#' @export
#'
#' @examples
#' path_vcf <- system.file(package = "vcf2mafR", "testfiles/test_grch38.vep.vcf")
#' vcf2df(path_vcf)
vcf2df <- function(vcf, tumor_id = "TUMOR", normal_id = "NORMAL", vcf_tumor_id = tumor_id, vcf_normal_id = normal_id, debug_mode = FALSE) {

  # Read VCF
  cli::cli_h1(text = "Reading VCF")
  vcfR <- vcfR::read.vcfR(vcf, verbose = FALSE) # May need to increase default limit
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

  ls_tidy_vcf <- vcfR::vcfR2tidy(vcfR, single_frame = TRUE, verbose = FALSE)
  if (debug_mode) {
    return(ls_tidy_vcf)
  }
  df_vcf <- ls_tidy_vcf[["dat"]]
  df_meta <- ls_tidy_vcf[["meta"]]


  # Rename columns to look more like expected input to df2maf
  data.table::setnames(df_vcf, old = c("CHROM", "POS", "REF", "ALT", "CSQ", "Indiv"), new = c("chr", "pos", "ref", "alt", "consequence", "sample"))

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
  assertions::assert(length(unexpected_sample_identifiers) == 0, msg = "Unrecognised sample identifiers: [{unexpected_sample_identifiers}], All Sample Identifiers should be either {tumor_id} or {normal_id}")



  #df_vcf_normal <- df_vcf |>
   # dplyr::filter(sample_type == "Normal") |>
    #dplyr::rename("Matched_Norm_Sample_Barcode" = sample, "Match_Norm_Seq_Allele1" = ref, "Match_Norm_Seq_Allele2" = alt) |> # Clear Match_Norm_Seq_Allele1 and Match_Norm_Seq_Allele2 in somatic maf (could contain germline information)
    #dplyr::select(Matched_Norm_Sample_Barcode, Match_Norm_Seq_Allele1, Match_Norm_Seq_Allele2)

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

  csq_meta <- unlist(df_meta[df_meta[["ID"]]=='CSQ','Description'])
  csq_fieldnames_string <- gsub(x=csq_meta, '.*. Format: ', "")
  csq_fieldnames <- unlist(strsplit(csq_fieldnames_string, split = "\\|"))

  # Check VEP annotations include all the columns we need: SYMBOL, BIOTYPE, cDNA_position, Consequence, CANONICAL
  assertions::assert_includes(
    csq_fieldnames, "SYMBOL", msg = 'Failed to find column [{.strong SYMBOL}] in vep annotations. Please rerun VEP with the {.arg Identifier: Gene Symbol} option turned on!'
  )

  assertions::assert_includes(
    csq_fieldnames, "CANONICAL", msg = 'Failed to find column [{.strong CANONICAL}] in vep annotations. Please rerun VEP with {.arg Identify Canonical Transcripts} option turned on!'
  )
  assertions::assert_includes(
    csq_fieldnames, "BIOTYPE", msg = 'Failed to find column [{.strong BIOTYPE}] in vep annotations. Please rerun VEP with {.arg Annotate Transcript biotype} option turned on!'
  )
  assertions::assert_includes(
    csq_fieldnames, "cDNA_position", msg = 'Failed to find column [{.strong cDNA_position}] in vep annotations. Please rerun Reconfigure VEP to inlude annoation of cDNA_position {.arg Annotate Transcript biotype} option turned on!'
  )


  # Fields we need to pull per variant to choose the most significant consequence: BIOTYPE, cDNA_position, Consequence
  ls_consequences = strsplit(df_vcf_somatic[['consequence']], split = ",")
  names(ls_consequences) <- df_vcf_somatic[['variant_id']]

  ls_consequences_df <- lapply(ls_consequences, \(consequences){
    ls <- strsplit2(consequences, split = "|", fixed = TRUE)
    df <- data.frame(do.call(rbind, ls))
    names(df) <- csq_fieldnames
    df <- tibble::tibble(df)
    return(df)
    })

  # Rank the effect of each variant on each affected transcript first by transcript biotype, then by consequence severity, and then by longest transcript
  ls_consequences_df <- lapply(ls_consequences_df, \(df){

    # df is a data.frame describing all possible effects of a variant, with each row representing a different transcript
    df |>
      dplyr::mutate(
        transcript_lengths = cdna_to_transcript_length(cDNA_position),
        biotype_priority = vep_rank_biotypes(BIOTYPE),
        csq_priority = vep_rank_consequences(Consequence),
        has_gene_symbol = !is.na(SYMBOL) & nchar(SYMBOL) > 0
      ) |>

      # Sort effects first by transcript biotype, then by severity, and then by longest transcript
      dplyr::arrange(biotype_priority, csq_priority, dplyr::desc(transcript_lengths)) |>
      dplyr::mutate(rank = seq_len(dplyr::n()))
  })

  ls_res <- lapply(ls_consequences_df, \(df){

    # === Step 1: Find affected gene === #
    # Find the highest priority effect with a gene symbol i.e. the worst affected gene
    maf_gene = df |>
      dplyr::filter(has_gene_symbol) |>
      dplyr::pull(SYMBOL) |>
      head(n=1)

    # === Step 1: Choose which gene transcript to use === #

    # If that gene has a user-preferred isoform, report the effect on that isoform
    # <not implemented yet: see custom-enst option of vcf2maf.pl for implementation>

    # If that gene has no user-preferred isoform, then use the VEP-preferred (canonical) isoform
    maf_effect <- df |>
      dplyr::filter(.data[["SYMBOL"]] == maf_gene, .data[["CANONICAL"]] == TRUE) |>
      dplyr::pull(Consequence) |>
      head(n=1)

    # If that gene has no VEP-preferred isoform either, then choose the worst affected user-preferred isoform with a gene symbol
    # <not implemented yet>

    # If none of the isoforms are user-preferred, then choose the worst affected VEP-preferred isoform with a gene symbol
    if(length(maf_effect) == 0)
      df[['Consequence']][df$has_gene_symbol][1]

    # If none of the isoforms are user-preferred, then choose the worst affected VEP-preferred isoform with a gene symbol

    # If we still have nothing selected, then just report the worst effect
    if(length(maf_effect) == 0)
      maf_effect <- df[["Consequence"]][1]


    return(list('maf_gene' = maf_gene, 'maf_effect' = maf_effect))
  })

  maflike_df <- do.call(rbind.data.frame, ls_res)
  maflike_df[['variant_id']] <- rownames(maflike_df)

  # Create Final MAF data.frame
  # from a mix of the list and the original variant-level vcf data.frame
  df_final <- df_vcf_somatic |>
    dplyr::left_join(maflike_df, by = "variant_id")

  cli::cli_alert_success("{.strong vcf2df successful}")
  return(df_final)
}


vcf2maf <- function(vcf, ref_genome, tumor_id = "TUMOR", normal_id = "NORMAL", vcf_tumor_id = tumor_id, vcf_normal_id = normal_id, debug_mode = FALSE){
  df <- vcf2df(
    vcf = vcf,
    tumor_id = tumor_id,
    normal_id = normal_id,
    vcf_tumor_id = vcf_tumor_id,
    vcf_normal_id = vcf_normal_id,
    debug_mode = debug_mode
  )
  #browser()

  df2maf(
    data = df,
    ref_genome = ref_genome,
    keep_all = TRUE,
    col_chrom = "chr",
    col_pos = "pos",
    col_sample_identifier = "Tumor_Sample_Barcode",
    col_ref = "ref",
    col_alt = "alt",
    col_consequence = "maf_effect",
    col_gene = "maf_gene",
    col_center = NULL,
    col_entrez_gene_id = NULL,
    col_dbSNP_rsid = NULL,
    col_dbSNP_validation_status = NULL,
    col_matched_normal_sample_identifier = NULL,
    col_sequencer = NULL,
    col_sequence_source = NULL,
    col_mutation_status = NULL
  )
}


