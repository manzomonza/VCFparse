### problematic genes and locations

#' Extract all shortest possible strings within curly brackets
#'
#' @param func_string
#'
#' @return
#' @export
#'
#' @examples {gene:string}
gene_extract = function(info_string){
  func_string = stringr::str_extract_all(info_string, pattern = '(?<=\\{).*?(?=\\})')
  func_string = lapply(func_string, function(x) gsub("'", "",x))
  return(func_string)
}

#' Extract values from FUNC string
#'
#' @param info_string
#'
#' @return
#' @export
#'
#' @examples
extract_gene_parameters = function(gene_string){
  FUNC = stringr::str_split(gene_string, pattern = ",")[[1]]
  df = data.frame(parameter = FUNC, value = NA)
  df = tidyr::separate(df, col = parameter, into = c("parameter", 'value'), sep = ":")
  df = tidyr::pivot_wider(df, names_from = parameter, values_from = value)
  return(df)
}

#' Apply
#'
#' @param info_string
#'
#' @return
#' @export
#'
#' @examples
convert_entries_to_table = function(info_string){
  gene_strings =  gene_extract(info_string)
  roi = dplyr::bind_rows(lapply(gene_strings, function(x) lapply(x, extract_gene_parameters)))
  roi = as.data.frame(roi)
  return(roi)
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
  data_reshaped = dplyr::relocate(data_reshaped, gene, coding, protein, location, contains("Ref"),
                                  contains("Alt"), contains("orig"), contains("normalized"))
  return(data_reshaped)
}




