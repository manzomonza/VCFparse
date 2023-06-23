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
  return(vec_index)
}



#' Title
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


#' Reformat Info
#'
#' @param infostring
#'
#' @return
#' @export
#'
#' @examples
reformat_INFO = function(infostring){
  df = data.frame(entries = unlist(stringr::str_split(infostring, pattern = ";")))
  df$entries = gsub("'",'', df$entries)
  df_list = list(dataf = data.frame(entries =df[-grep("FUNC=", df$entries),]),
                 FUNC = df[grep("FUNC=", df$entries),])
  df_list$dataf = tidyr::separate(df_list$dataf, col = 'entries', into = c("parameter","value"), sep = "=")
  df_list$dataf = tidyr::pivot_wider(df_list$dataf, names_from = parameter, values_from = value)
  df_list$FUNC = convert_entries_to_table(df_list$FUNC)
  indeces = AF_based_index(df_list$dataf$AF)
  data_row = dplyr::bind_rows(lapply(indeces, function(x) index_all_cols(df_list$dataf, index_pos = x)))
  data_reshaped = cbind(df_list$FUNC, data_row)
  return(data_reshaped)
}