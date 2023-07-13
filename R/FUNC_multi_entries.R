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
entry_No_check = function(vcf_element){
  gene_No = length(grep("gene", unlist(vcf_element$FUNC), value = TRUE))
  if(gene_No > 1){
    new_FUNC = combine_and_split(unlist(vcf_element$FUNC))
    vcf_element$FUNC = new_FUNC
  }
  return(vcf_element)
}

#' Collapses string and splits based on '}' in string.
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
  multi_entries = list(lapply(multi_entries, function(x) unlist(stringr::str_split(x, pattern = ","))))
  return(multi_entries)
}

