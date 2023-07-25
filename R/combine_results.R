## COMBINE RESULTS WITH ORIGINAL TABLE

#' Column bind vcf FUNC extracts with original vcf row
#'
#' @param vcf_element
#'
#' @return
#' @export
#'
#' @examples
combine_orig_with_FUNC_extracts = function(vcf_element){
  char_vec = unlist(vcf_element$FUNC)
  res_df = combine_FUNC_extracts(char_vec)
  combo_df = cbind(vcf_element, res_df)
  combo_df = tibble::as_tibble(combo_df)
  combo_df = dplyr::relocate(combo_df, rowid:origPos)
  return(combo_df)
}

#' Apply functions to extract information from FUNC to vcf
#'
#' @param vcffilepath
#'
#' @return
#' @export
#'
#' @examples
apply_FUNC_functions = function(vavcf){
  vavcf = dplyr::group_by(vavcf, rowid)
  vavcf = dplyr::group_split(vavcf)
  vavcf = lapply(vavcf, entry_No_check)
  vavcf = lapply(vavcf, function(x) from_multiple_FUNC_entries_to_multiple_vcf_rows(x))
  vavcf = dplyr::bind_rows(vavcf)
  vavcf$rowid = 1:nrow(vavcf)
  vavcf = dplyr::group_by(vavcf, rowid)
  vavcf = dplyr::group_split(vavcf)
  res_vcf = dplyr::bind_rows(lapply(vavcf, combine_orig_with_FUNC_extracts))
  return(res_vcf)
}


#' Parse multiple entries from FUNC field and convert into multiple rows.
#' Following functions parse FUNC field then uniformly.w
#'
#' @param vcf_element
#'
#' @return
#' @export
#'
#' @examples
from_multiple_FUNC_entries_to_multiple_vcf_rows = function(vcf_element){
  typevar = typeof(vcf_element$FUNC[[1]])
  if(typevar == 'list'){
    newvcf = list()
    for (i in seq_along(vcf_element$FUNC[[1]])){
      newfunc = dplyr::select(vcf_element, -FUNC)
      newfunc$FUNC = vcf_element$FUNC[[1]][i]
      newvcf[[i]] = newfunc
    }
    newvcf = dplyr::bind_rows(newvcf)
    return(newvcf)
  }else{
    return(vcf_element)
  }
}
