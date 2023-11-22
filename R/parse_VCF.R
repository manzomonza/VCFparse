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
    cnv_rows = dplyr::distinct(cnv_rows)
    cnv_rows = tidyr::unnest(cnv_rows, CI)
    cnv_rows$CI = paste0('perc_',cnv_rows$CI)
    cnv_rows = tidyr::separate(cnv_rows, CI, into = c("percentile", "value"), sep = ":")
    cnv_rows = tidyr::pivot_wider(cnv_rows, names_from = "percentile", values_from = 'value')
    cnv_rows = dplyr::arrange(cnv_rows, desc(CN))
  }
    return(cnv_rows)
}

#' Return dataframe from percentile and value pair in one string
#'
#' @param percentile_value
#'
#' @return
#' @export
#'
#' @examples
cdf_extract = function(percentile_value){
  perc_val = stringr::str_split_1(percentile_value, pattern =":")
  perc_val = data.frame(percentile = perc_val[1], observed_val = perc_val[2])
  return(perc_val)
}


#' Parse CDF values from CDF_MAPD column to long and then wide data format.
#'
#' @param cnvrow_element
#'
#' @return
#' @export
#'
#' @examples
cnv_extract_cdf = function(cnvrow_element){
  # cnvrow_element == table
  # table needs to contain FUNC column
  cnvrow_element$gene = stringr::str_extract(cnvrow_element$FUNC, pattern = "(?<='gene':').*?(?=')")
  cdfs = lapply(cnvrow_element$CDF_MAPD[[1]], cdf_extract)
  cdfs = dplyr::bind_rows(cdfs)
  cdfs$percentile = paste0("percentile_", cdfs$percentile)
  cdfs = tidyr::pivot_wider(cdfs, names_from = percentile,
                            values_from = observed_val)
  cnvrow_element = cbind(cnvrow_element, cdfs)
  return(cnvrow_element)
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
