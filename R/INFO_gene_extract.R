### problematic genes and locations

#' Extract st
#'
#' @param func_string
#'
#' @return
#' @export
#'
#' @examples
gene_extract = function(func_string){
  func_string = stringr::str_extract_all(func_string, pattern = '(?<=\\{).*?(?=\\})')
  func_string = lapply(func_string, function(x) gsub("'", "",x))
  return(func_string)
}

extract_FUNC = function(info_string){
  FUNC = stringr::str_split(info_string, pattern = ",")[[1]]
  df = data.frame(parameter = FUNC, value = NA)
  df = tidyr::separate(df, col = parameter, into = c("parameter", 'value'), sep = ":")
  df = tidyr::pivot_wider(df, names_from = parameter, values_from = value)
  return(df)
}

parse_FUNC = function(vcfInfo){
  func_extracts =  gene_extract(vcfInfo)
  roi = dplyr::bind_rows(lapply(func_extracts, function(x) lapply(x, extract_FUNC))) |> as.data.frame()
  return(roi)
}

i = 4
teststring = vcfs[[8]][[i]]$INFO
parse_FUNC(teststring)
lapply(vcfs[[1]], function(x) parse_FUNC(x$INFO))

vcfs
