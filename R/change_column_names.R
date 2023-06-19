#' Change column names from VCF
#'
#' @param vcf
#'
#' @return
#' @export
#'
#' @examples
change_col_name = function(vcf){
  chrom_col_index = which(colnames(vcf) == '#CHROM')
  colnames(vcf)[chrom_col_index] = 'chr'
  return(vcf)
  }
