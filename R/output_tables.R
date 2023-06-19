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
  output = dplyr::select(combined_output, rowid, gene, coding,protein, AF,QUAL, GQ:VCFREF)
  return(output)
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



write_output_tables = function(combined_output, filepath){
  dirpath = dirname(filename)
  dirpath = paste0(dirpath, '/watchdog_v2_output')
  dir.create(dirpath)

  }



}
