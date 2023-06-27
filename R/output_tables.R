#' Output clean variant table starting from combined output from VCF parsing
#'
#' @param combined_output
#'
#' @return
#' @export
#'
#' @examples
variant_table = function(combined_output){
  output = dplyr::select(combined_output, rowid, gene, transcript, chr,POS,coding,protein,exon,location, AF)
  return(output)
}

#' Output clean metric table starting from combined output from VCF parsing
#'
#' @param combined_output
#'
#' @return
#' @export
#'
#' @examples
metric_table = function(combined_output){
  output = dplyr::select(combined_output, rowid, gene, coding,protein, AF,QUAL)
  return(combined_output)
}


#' Output variant table to be used for variant interpretation
#' TSGinfo, cancerHotspots etc.
#' @param combined_output
#'
#' @return
#' @export
#'
#' @examples
interpretation_table = function(combined_output){
  output = dplyr::select(combined_output, rowid, gene, coding,protein,chr, POS, REF, ALT)
  return(output)
}


#' Creates output tables and directory by individual function calls
#'
#' @param combined_output
#' @param filepath
#'
#' @return
#' @export
#'
#' @examples
write_output_tables = function(combined_output, filepath){
  ## Create output dir first
  dirpath = dirname(filepath)
  dirpath = paste0(dirpath, '/watchdog_v2_output')
  dir.create(dirpath)
  ## write out parsed tables
  readr::write_tsv(interpretation_table(combined_output), file = "interpretation_table.tsv")
  readr::write_tsv(metric_table(combined_output), file = "variant_metrics_table.tsv")
  readr::write_tsv(variant_table(combined_output), file = "simplified_variant_table.tsv")
  }

