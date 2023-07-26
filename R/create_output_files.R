## Create outputs

#' Title
#'
#' @param path_vcf
#'
#' @return
#' @export
#'
#' @examples
create_output_folder = function(path_vcf){
  output_dir = dirname(path_vcf)
  sampleName = as.character(unique(readVCF(path_vcf)$sampleNames))
  output_dirname = paste0(output_dir,"/", sampleName, '_djinni')
  if(!dir.exists(output_dirname)){
    dir.create(output_dirname)
  }
  return(output_dirname)
}


