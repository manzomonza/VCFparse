#' Extract commented section from VCF
#'
#' @param vcfpath
#'
#' @return
#' @export
#'
#' @examples
vcf_comment_section = function(vcfpath){
  comment_section = readr::read_tsv(vcfpath, col_names = FALSE)
  comment_section = dplyr::filter(comment_section, grepl("##", X1))
  comment_section$X1 = gsub("##", '', comment_section$X1)
  comment_section$X1 = gsub('\\"', '', comment_section$X1)
  comment_section = tidyr::separate(comment_section, col = X1, into = c("parameter", "value"), extra = "merge", sep = "=")
  path_df = tibble::tibble(parameter = "path_vcf",
                       value = vcfpath)
  comment_section = dplyr::bind_rows(path_df, comment_section)
  return(comment_section)
}

pull_comment_value = function(comment_section, stringoi){
  cm = dplyr::filter(comment_section, parameter == stringoi)
  value = cm$value
  return(value)
}
