#' Extra#'#'#'racExtract character string from INFO field in VCF
#'
#' @param vcf_info
#'
#' @return
#' @export
#'
#' @examples
extract_variant_info = function(vcf_info){
  info = stringr::str_extract(vcf_info, pattern = '(?<=FUNC\\=\\[\\{).*(?=\\}\\]$)')
  return(info)
}

#' Convert character string to single row with multiple columns per variant info string
#'
#' @param vcf_info_string
#'
#' @return
#' @export
#'
#' @examples
variant_info_to_row = function(vcf_info_string){
  df = data.frame(entries = unlist(stringr::str_split(vcf_info_string, pattern = ";")))
  df$entries = gsub("'",'', df$entries)
  df = tidyr::separate(df, col = 'entries', into = c("parameter","value"), sep = ":")
  df = tidyr::pivot_wider(df, names_from = parameter, values_from = value)
  return(df)
}

vcf_to_info_row = function(filepath){
  vc = readr::read_tsv(filepath, comment = "##")
  infos = lapply(vc$INFO, extract_variant_info)
  df_row = lapply(infos, variant_info_to_row)
  return(df_row)
}





