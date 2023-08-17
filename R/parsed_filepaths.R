# filepaths created from dirname

#' Generate paths for output files based on vcfpath
#'
#' @param vcfpath
#'
#' @return
#' @export
#'
#' @examples
parsed_filepaths = function(output_path){
  path_parsed_output = paste0(output_path, "/parsed_output")
  if(!dir.exists(path_parsed_output)){
    dir.create(path_parsed_output)
  }
  parsed_info = paste0(path_parsed_output, '/parsed_vcf_info.yaml')
  parsed_snv = paste0(path_parsed_output, '/parsed_snv.tsv')
  parsed_complete = paste0(path_parsed_output, '/parsed_complete.tsv')
  parsed_cnv = paste0(path_parsed_output, '/parsed_cnv.tsv')
  parsed_fusions = paste0(path_parsed_output, '/parsed_fusions.tsv')
  parsed_fp = paste0(path_parsed_output, '/parsed_filepaths.tsv')
  return(list(parsed_info = parsed_info,
              parsed_snv = parsed_snv,
              parsed_cnv = parsed_cnv,
              parsed_fusions = parsed_fusions,
              parsed_complete = parsed_complete,
              parsed_fp = parsed_fp))
}

