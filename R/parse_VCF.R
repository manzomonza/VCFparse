#' Parse VCF -- create SNV
#'
#' @param vcf
#'
#' @return
#' @export
#'
#' @examples
parse_vcf_return_snv = function(vcf){
  # SNV table
  snv = dplyr::filter(vcf, variant_type != 'synonymous' & alt != "<CNV>")
  ## Remove zero AF entries, especially important for Genexus
  if("AF" %in% colnames(snv)){
    snv = dplyr::filter(snv, AF != 0)
  }
  snv = dplyr::relocate(snv, rowid)
  return(snv)
}


#' Parse VCF -- create CNV
#'
#' @param vcf
#'
#' @return
#' @export
#'
#' @examples
parse_vcf_return_cnv = function(vcf){
  cnv_rows = dplyr::filter(vcf, alt == "<CNV>")
  if(nrow(cnv_rows) > 0){
    cnv_rows = cnv_parse(cnv_rows)
    cnv_rows = dplyr::select(cnv_rows, -origPos)
    cnv_rows = dplyr::arrange(cnv_rows, desc(RAW_CN))
  }
    return(cnv_rows)
}

#' Parse VCF - return META
#'
#' @param vcf
#'
#' @return
#' @export
#'
#' @examples
parse_vcfpath_return_metainformation = function(vcfpath){
  cm = vcf_comment_section(vcfpath = vcfpath )
  metainf = aggregate_META_information(cm)
  return(metainf)
}

#' Write out parsed filepath information
#'
#' @param parsed_fp
#'
#' @return
#' @export
#'
#' @examples
write_parsed_fp_txt = function(parsed_fp, vcf_file, analysis_name){
  parsed_fp[1:5] = lapply(  parsed_fp[1:5], function(x) ifelse(file.exists(x), x, NA))
  filepath_df = data.frame(file = names(parsed_fp),
                           filepath  = unlist(unname(parsed_fp)))
  filepath_df = attach_ID(filepath_df, vcf_file = vcf_file, analysis_name = analysis_name)
  readr::write_tsv(filepath_df, file = parsed_fp$parsed_fp)
}
