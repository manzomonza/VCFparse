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
  path_file_info = paste0(path_parsed_output, '/parsed_vcf_info.yaml')
  path_file_snv = paste0(path_parsed_output, '/parsed_snv.tsv')
  path_file_complete = paste0(path_parsed_output, '/parsed_complete.tsv')
  path_file_cnv = paste0(path_parsed_output, '/parsed_cnv.tsv')
  path_file_fusions = paste0(path_parsed_output, '/parsed_fusions.tsv')
  filepaths = paste0(path_parsed_output, '/parsed_filepaths.tsv')
  return(list(path_file_info = path_file_info,
              path_file_snv = path_file_snv,
              path_file_cnv = path_file_cnv,
              path_file_fusions = path_file_fusions,
              path_file_complete = path_file_complete,
              path_filepath = filepaths))
}

