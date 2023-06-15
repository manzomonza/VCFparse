### problematic genes and locations

func_extract = function(func_string){
  func_string = stringr::str_extract(func_string, pattern = '(?<=FUNC=\\[).*(?=\\])')
  func_string = gsub("\\[|\\]", '',func_string)
  roi = stringr::str_split(func_string, pattern = ",\\{")[[1]]
  return(roi)
}

parse_FUNC = function(vcfInfo){
  func_extracts =  func_extract(vcfInfo)
  roi = dplyr::bind_rows(lapply(func_extracts, function(x) lapply(x, extract_FUNC))) |> as.data.frame()
  roi = dplyr::relocate(roi, contains('orig'), contains('normalized'))
  return(roi)
}




