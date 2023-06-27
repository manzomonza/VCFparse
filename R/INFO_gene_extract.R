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

#' Title
#'
#' @param info_string
#'
#' @return
#' @export
#'
#' @examples
extract_gene_parameters = function(info_string){
  FUNC = stringr::str_split(info_string, pattern = ",")[[1]]
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



