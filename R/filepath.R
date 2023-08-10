# filepaths created from dirname

#' Generate paths for output files based on vcfpath
#'
#' @param vcfpath
#'
#' @return
#' @export
#'
#' @examples
generate_paths = function(vcfpath){
  analysis_name = aggregate_META_information(vcf_comment_section(vcfpath))
  analysis_name = analysis_name$IonReporter$AnalysisName
  analysis_name = stringr::str_remove(analysis_name, pattern=" ")
  vcfdir = dirname(vcfpath)
  path_parsed_output = paste0(vcfdir, "/", analysis_name, "_parsed_output")
  if(!dir.exists(path_parsed_output)){
    dir.create(path_parsed_output)
  }
  path_file_info = paste0(path_parsed_output, '/file_info.yaml')
  path_file_snv = paste0(path_parsed_output, '/snv.txt')
  path_file_complete = paste0(path_parsed_output, '/complete.txt')
  path_file_cnv = paste0(path_parsed_output, '/cnv.txt')
  path_file_fusions = paste0(path_parsed_output, '/fusions.txt')
  filepaths = paste0(path_parsed_output, '/filepaths.txt')
  return(list(path_file_info = path_file_info,
              path_file_snv = path_file_snv,
              path_file_cnv = path_file_cnv,
              path_file_fusions = path_file_fusions,
              path_file_complete = path_file_complete,
              path_filepath = filepaths))
}

