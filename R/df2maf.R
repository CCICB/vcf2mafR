#' Convert VCF-like data.frame to MAF
#'
#' Converts a VEP-like dataframe with to a minimal MAF dataframe.
#' **Input dataframe should include the columns:**
#' * **chr**
#' * **pos** (1-based)
#' * **sample**
#' * **ref**
#' * **alt**
#' * **consequence** (Sequence Ontology terms)
#' * **gene**
#'
#' Alternate column names can be used, but require **col_chrom**, **col_pos**, **col_<etc>** arguments to be supplied.
#'
#' @param data a data.frame with 1 row per mutation in a cohort (data.frame)
#' @param ref_genome name of the reference genome used to call variants (string)
#' @param keep_all keep all columns in the original data.frame? (flag). If FALSE, only includes the minimal required columns and any column explicitly mapped using `col_<name>` arguments.
#' @param col_chrom name of column describing chromosome of the mutation (string)
#' @param col_pos name of column describing the 1-based position of the mutation (string)
#' @param col_ref name of column describing the reference allele (string)
#' @param col_alt name of column describing the alternate allele (string)
#' @param col_sample_identifier name of column describing the sample containing the mutation (string)
#' @param col_consequence name of column describing the consequence of the mutation (in SO ontology terms e.g. those that VEP would use) containing the mutation (string)
#' @param col_gene name of column containing Hugo_Symbol of the gene affected by the mutation (string)
#' @param col_entrez_gene_id name of column containing entrez gene IDs (string)
#' @param col_center name of column containing the genome sequencing center reporting the variant (string)
#' @param col_dbSNP_rsid name of column describing the dbSNP rsid of the variant, or "novel" if there is no dbSNP record (string)
#' @param col_dbSNP_validation_status name of column describing the validation status of the variant. Elements must be one of by1000genomes;by2Hit2Allele; byCluster; byFrequency; byHapMap; byOtherPop; bySubmitter; alternate_allele (string)
#' @param col_matched_normal_sample_identifier name of column describing the matched normal sample identifier
#' @param col_sequencer name of column describing the instrument used to produce data (string)
#' @param col_sequence_source name of column describing the molecular assay type used to produce the analytes used for sequencing (string). Elements are usually one of 'WGS', 'WGA', 'WXS', 'RNA-seq', etc
#' @return a maf-like data.frame (data.table)
#' @export
#'
#' @examples
#' path_mutations <- system.file("pedcbioportal_mutation_annotated.tsv", package = "vcf2mafR")
#' df_mutations <- read.csv(path_mutations, sep = "\t")
#' df_maf <- df2maf(df_mutations, ref_genome = "hg38")
df2maf <- function(
    data,
    ref_genome,
    keep_all = TRUE,
    col_chrom = "chr",
    col_pos = "pos" ,
    col_sample_identifier = "sample_id",
    col_ref = "ref",
    col_alt =  "alt",
    col_consequence = "consequence",
    col_gene = "gene",
    col_center = NULL,
    col_entrez_gene_id = NULL,
    col_dbSNP_rsid = NULL,
    col_dbSNP_validation_status = NULL,
    col_matched_normal_sample_identifier = NULL,
    col_sequencer = NULL,
    col_sequence_source = NULL
    ){

  # Assertions
  assertions::assert_dataframe(data)
  assertions::assert_string(col_chrom)
  assertions::assert_string(col_pos)
  assertions::assert_string(col_sample_identifier)
  assertions::assert_string(col_ref)
  assertions::assert_string(col_alt)
  assertions::assert_string(col_consequence)
  assertions::assert_string(col_gene)
  assertions::assert_names_include(data, names = c(col_chrom, col_pos, col_sample_identifier, col_ref, col_alt, col_consequence, col_gene))
  assertions::assert_string(ref_genome)
  assertions::assert_flag(keep_all)

  old_names <- c(col_chrom, col_pos, col_sample_identifier, col_ref, col_alt, col_consequence,  col_gene)
  new_names <- c("Chromosome", "Position_1based", "Tumor_Sample_Barcode", "Reference_Allele", "Tumor_Seq_Allele2", "Consequence", "Hugo_Symbol")

  # Eliminate no visible global binding
  Reference_Allele <- NULL
  Tumor_Seq_Allele2 <- NULL
  Ref_Length <- NULL
  Alt_Length <- NULL
  Position_1based <- NULL
  Consequence <- NULL
  Variant_Type <- NULL
  Inframe <- NULL

  # Conditional Assertions --------------------------------------------------
  #TODO: add content assertion for fields like Sequence_Source which have restricted language
  # Center
  if(!is.null(col_center)){
    assertions::assert_string(col_center)
    assertions::assert_names_include(data, col_center)
    old_names <- c(old_names, col_center)
    new_names <- c(new_names, "Center")
  }

  # Sequencer
  if(!is.null(col_sequencer)){
    assertions::assert_string(col_sequencer)
    assertions::assert_names_include(data, col_sequencer)
    old_names <- c(old_names, col_sequencer)
    new_names <- c(new_names, "Sequencer")
  }

  # Matched_Norm_Sample_Barcode
  if(!is.null(col_matched_normal_sample_identifier)){
    assertions::assert_string(col_matched_normal_sample_identifier)
    assertions::assert_names_include(data, col_matched_normal_sample_identifier)
    old_names <- c(old_names, col_matched_normal_sample_identifier)
    new_names <- c(new_names, "Matched_Norm_Sample_Barcode")
  }

  # dbSNP_RS
  if(!is.null(col_dbSNP_rsid)){
    assertions::assert_string(col_dbSNP_rsid)
    assertions::assert_names_include(data, col_dbSNP_rsid)
    old_names <- c(old_names, col_dbSNP_rsid)
    new_names <- c(new_names, "dbSNP_RS")
  }

  # dbSNP_Val_Status
  if(!is.null(col_dbSNP_validation_status)){
    assertions::assert_string(col_dbSNP_validation_status)
    assertions::assert_names_include(data, col_dbSNP_validation_status)
    old_names <- c(old_names, col_dbSNP_validation_status)
    new_names <- c(new_names, "dbSNP_Val_Status")
  }

  # Entrez_Gene_Id
  if(!is.null(col_entrez_gene_id)){
    assertions::assert_string(col_entrez_gene_id)
    assertions::assert_names_include(data, col_entrez_gene_id)
    old_names <- c(old_names, col_entrez_gene_id)
    new_names <- c(new_names, "Entrez_Gene_Id")
  }

  # Sequence_Source
  if(!is.null(col_sequence_source)){
    assertions::assert_string(col_sequence_source)
    assertions::assert_names_include(data, col_sequence_source)
    old_names <- c(old_names, col_sequence_source)
    new_names <- c(new_names, "Sequence_Source")
  }

  # Create Data.Table
  dt_maf <- data.table::data.table(data)

  #Rename Columns
  data.table::setnames(
    dt_maf,
    old = old_names,
    new = new_names
  )

  # Rename optional columns

  # Calculate Ref and Alt lengths
  dt_maf[,"Ref_Length" := nchar(Reference_Allele)]
  dt_maf[,"Alt_Length" := nchar(Tumor_Seq_Allele2)]

  # Standardise Ref and Alt allele representations, and fix Lengths & Positions
  df_fixed <-  fix_alleles(ref = dt_maf[["Reference_Allele"]], alt = dt_maf[["Tumor_Seq_Allele2"]])
  dt_maf[,"Ref_Length" := Ref_Length - df_fixed[["numdropped"]]]
  dt_maf[,"Alt_Length" := Alt_Length - df_fixed[["numdropped"]]]
  dt_maf[,"Position_1based" := Position_1based + df_fixed[["numdropped"]]]
  dt_maf[,"Reference_Allele" := df_fixed[["ref"]]]
  dt_maf[,"Tumor_Seq_Allele2" := df_fixed[["alt"]]]

  # Calculate Start_Position, End_Position, Variant_Types and Inframe status
  dt_maf[,"Start_Position" := data.table::fcase(
    # SNPs
    Ref_Length == Alt_Length, Position_1based,
    # Insertions (potentially with '-' Reference Alleles)
    Ref_Length < Alt_Length, ifelse(Reference_Allele == "-", yes = Position_1based - 1, Position_1based),
    # Deletions
    Ref_Length > Alt_Length, Position_1based,
    TRUE, stop("non-explicitly handled mutation type")
    )]
  dt_maf[,"End_Position" := data.table::fcase(
    # SNPs
    Ref_Length == Alt_Length, Position_1based + Alt_Length - 1,
    # Insertions (potentially with '-' Reference Alleles)
    Ref_Length < Alt_Length, ifelse(Reference_Allele == "-", yes = Position_1based, no = Position_1based + Ref_Length - 1),
    # Deletions
    Ref_Length > Alt_Length, Position_1based + Ref_Length - 1,
    TRUE, stop("non-explicitly handled mutation type")
  )]
  dt_maf[,"Inframe" := data.table::fcase(
    Ref_Length == Alt_Length, TRUE,
    abs(Ref_Length - Alt_Length) %% 3 == 0, TRUE,
    default = FALSE
    )]
  dt_maf[,"Variant_Type" := data.table::fcase(
    # SNPs
    Ref_Length == Alt_Length & Alt_Length == 1, "SNP",
    Ref_Length == Alt_Length & Alt_Length == 2, "DNP",
    Ref_Length == Alt_Length & Alt_Length == 3, "TNP",
    Ref_Length == Alt_Length & Alt_Length > 3, "ONP",
    # Insertions (potentially with '-' Reference Alleles)
    Ref_Length < Alt_Length, "INS",
    # Deletions
    Ref_Length > Alt_Length, "DEL",
    TRUE, stop("non-explicitly handled mutation type")
  )]


  # Convert SO to MAF mutation types
  dt_maf[, "Variant_Classification" := mutationtypes::mutation_types_convert_so_to_maf(so_mutation_types = Consequence, variant_type = Variant_Type, inframe = Inframe)]


  # Add reference genome
  dt_maf[, "NCBI_Build" := ref_genome]

  if(!keep_all){
   #cli::cli_alert_info("Dropping all non-essential columns without explicit mapping")
    colnames <- c(new_names, "NCBI_Build", "Start_Position", "End_Position", "Variant_Classification", "Variant_Type", "Inframe")
    colnames <- colnames[colnames != "Position_1based"]
    dt_maf <- dt_maf[, colnames, with=FALSE]
  }

  # Return
  return(dt_maf)

}

#Changing Ref:
# While the first Char of Reference_Allele and Tumor_Seq_Allele2 are the same, and both are non-empty
# go char by char and if the first char of Ref and Alt are the same, drop them. If empty - replace with a '-'
# Each time you drop a char, add 1 to Position_1based, and decrease 1 from --$ref_length; --$var_length;

# Need to add info about alleles
fix_alleles_scalar <- function(ref, alt) {
  assertions::assert_string(ref)
  assertions::assert_string(alt)

  numdropped = 0
  while(nchar(ref) != 0 & nchar(alt) != 0 & substr(ref, 1, 1) == substr(alt, 1, 1) & ref != alt){
    ref <- substr(ref, 2, nchar(ref))
    alt <- substr(alt, 2, nchar(alt))
    if(nchar(ref) == 0) ref <- "-"
    if(nchar(alt) == 0) alt <- "-"
    numdropped <- numdropped + 1
  }

  list(ref = ref, alt = alt, numdropped = numdropped)
}

# Returns a 3col dataframe.
# 'ref': new reference allele
# 'alt': new alt allele
# 'numdropped': number of chars dropped from the start of reference/alt.
# numdropped should be subtracted from Ref_Length & Alt_Length and added to Pos
fix_alleles <- function(ref, alt){
  assertions::assert_equal(length(ref), length(alt))

  ls_fixed_alleles <- lapply(seq_along(ref), FUN = function(i) {
    fix_alleles_scalar(
      ref = ref[i],
      alt = alt[i]
    )
  })

  df <- as.data.frame(do.call(rbind, ls_fixed_alleles), stringsAsFactors = FALSE)

  # Remove unnecessary list-ing of dataframe columns
  for (col in colnames(df)){ df[[col]] <- unlist(df[[col]]) }

  return(df)
}


