
# read_vcf = function(vcf_filepath){
#   vcf = readr::read_tsv(vcf_filepath, comment = "##")
#   vcf = vcf[which(unlist(lapply(vcf$INFO, check_non_zero_AF))),]
#   vcf = dplyr::filter(vcf, FILTER != 'FAIL')
#   if(nrow(vcf)){
#     vcf$rowid = 1:nrow(vcf)
#     vcf = dplyr::group_by(vcf, rowid)
#     vcf = dplyr::group_split(vcf)
#     return(vcf)
#   }
# }

#' Uses VariantAnnotation package to read in VCF as GRanges object
#'
#' @param vcffilepath
#'
#' @return
#' @export
#'
#' @examples
readVCF = function(vcffilepath){
  vavcf = VariantAnnotation::readVcfAsVRanges(vcffilepath)
  vavcf = tibble::as_tibble(vavcf)
  # filter_column = readr::read_tsv(vcffilepath, comment = "##")$FILTER
  # vavcf$FILTER = filter_column
  #vavcf = dplyr::filter(vavcf, AF !=0 | alt == "<CNV>")
  vavcf$rowid = 1:nrow(vavcf)
  return(vavcf)
}


