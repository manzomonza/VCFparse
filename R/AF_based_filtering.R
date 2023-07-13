# AF based filtering

#' Title
#'
#' @param AF_char_string
#'
#' @return
#' @export
#'
#' @examples
AF_based_index = function(AF_char_string){
  split_string = as.numeric(stringr::str_split_1(AF_char_string, pattern = ','))
  vec_index = which(split_string != 0)
  if( identical(vec_index, integer(0))){
    return(NA)
  }
  return(vec_index)
}


#' Check AF entries in vcf INFO entries
#'
#' @param vcfinfo
#'
#' @return
#' @export
#'
#' @examples
check_non_zero_AF = function(vcfinfo){
  # return true if AF is not only 0 entries
  afstring = stringr::str_extract(vcfinfo, pattern = "^AF.*?;")
  if(is.na(afstring)){
    return(TRUE)
  }
  afstring = gsub("AF=|;",'', afstring)
  af_vec = stringr::str_split_1(afstring, pattern = ',')
  equal_nonzero = !all(af_vec == "0")
  return(equal_nonzero)
}


#' Extract comma separated values with specified indeces
#'
#' @param char_string
#' @param vector_indeces
#'
#' @return
#' @export
#'
#' @examples
return_filtered_vectors = function(char_string, vector_indeces){
  if(!is.na(char_string)){
    split_string = stringr::str_split_1(char_string, pattern = ',')
    if(length(split_string) > 1){
      split_string = split_string[vector_indeces]
    }
    return(split_string)
  }
}

#' Index all column entries
#'
#' @param datadf
#'
#' @return
#' @export
#'
#' @examples
index_all_cols = function(datadf, index_pos){
  cols= colnames(df_list$dataf)
  res = dplyr::mutate(df_list$dataf, across(all_of(cols),
                                            function(x) return_filtered_vectors(char_string = x,
                                                                                vector_indeces = index_pos)))
  return(res)
}

