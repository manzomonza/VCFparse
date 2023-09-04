#' Introduce function to remove buggy output
#'
#' @param tableofinterest
#'
#' @return
#' @export
#'
#' @examples
remove_multientry_leftover_characters = function(tableofinterest){
  check_cols = tableofinterest
  for (col in check_cols){
    tableofinterest[[col]] = gsub("\\[|\\{|\\}|\\]", '', tableofinterest[[col]])
  }
  return(tableofinterest)
}
