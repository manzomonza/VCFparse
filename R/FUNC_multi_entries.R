## Challenge: Multiple entries in vcf FUNC field

## WORKAROUND
#' Check how many times 'gene' information is provided.
#' If more than once, 'combine_and_split' function is applied
#'
#' @param vcf_element
#'
#' @return
#' @export
#'
#' @examples
gene_No_check = function(FUNC_element){
  gene_No = length(grep("gene", unlist(FUNC_element), value = TRUE))
  return(gene_No)
}

#' Collapses string and splits based on curly brackets in string.
#' Returns single list with splits as elements
#' @param vcf_element_FUNC
#'
#' @return
#' @export
#'
#' @examples
combine_and_split = function(vcf_element_FUNC){
  multi_entries = paste(vcf_element_FUNC, collapse = ',')
  multi_entries = gsub("'", '', multi_entries)
  multi_entries = unlist(stringr::str_split(multi_entries, pattern = "\\},"))
  multi_entries = lapply(multi_entries, function(x) unlist(stringr::str_split(x, pattern = ",")))
  multi_entries = lapply(multi_entries, function(x) gsub("\\{|\\}|'",'', x))
  multi_entries = lapply(multi_entries, function(x) FUNC_extracts_to_df(x))
  return(multi_entries)
}

#' Extract variant information and combine with read metrics for 'data-complete' table
#'
#' @param vcf
#'
#' @return
#' @export
#'
#' @examples
combine_orig_with_FUNC_extracts = function(vcf){
  res_df = list()
  for (i in seq_along(vcf$FUNC)){
    number_of_gene_occurences = gene_No_check(vcf$FUNC[[i]])
    if(number_of_gene_occurences > 1){
      list_elem = combine_and_split(vcf$FUNC[[i]])
      variant_info_df = dplyr::bind_rows(list_elem)
    }else{
      variant_info_df = FUNC_extracts_to_df(vcf$FUNC[[i]])
    }
    variant_info_df$rowid = i
    res_df[[i]] = variant_info_df
  }
  res_df = dplyr::bind_rows(res_df)
  vcf_df = dplyr::full_join(res_df, vcf, by = 'rowid')
  return(vcf_df)
}


