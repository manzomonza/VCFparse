### CNV PARSING


#' Parse CNV rows
#'
#' @param cnv_rows
#'
#' @return
#' @export
#'
#' @examples
cnv_parse = function(cnv_rows){
  cnv_rows = dplyr::select(cnv_rows, seqnames, gene, NUMTILES, CI, CN, PVAL)
  cnv_rows = dplyr::distinct(cnv_rows)
  cnv_rows = tidyr::unnest(cnv_rows, CI)
  cnv_rows$CI = paste0('perc_',cnv_rows$CI)
  cnv_rows = tidyr::separate(cnv_rows, CI, into = c("percentile", "value"), sep = ":")
  cnv_rows = tidyr::pivot_wider(cnv_rows, names_from = "percentile", values_from = 'value')
  return(cnv_rows)
}
