#' Read VCF file skipping ## commented entries.
#' Also skips entries with 'AF=0;' as these are likely only derived by hotspot annotation and not by variant calling.
#'
#' @param vcf_filepath
#'
#' @return
#' @export
#'
#' @examples
read_vcf = function(vcf_filepath){
  vcf = readr::read_tsv(vcf_filepath, comment = "##")
  vcf = check_non_zero_AF(vcf)
  vcf = dplyr::filter(vcf,FILTER != 'FAIL')
  if(nrow(vcf)){
    vcf$rowid = 1:nrow(vcf)
    vcf = dplyr::group_by(vcf, rowid)
    vcf = dplyr::group_split(vcf)
    return(vcf)
  }
}


