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
  vcf = dplyr::filter(vcf,!grepl('AF=0;', INFO) & FILTER != 'FAIL')
  vcf$rowid = 1:nrow(vcf)
  vcf = dplyr::group_by(vcf, rowid)
  vcf = dplyr::group_split(vcf)
  return(vcf)
}
