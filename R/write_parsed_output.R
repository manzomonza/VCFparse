#' Write out VCF parse information
#'
#' @param vcf
#' @param vcf_path
#' @param parsed_fp
#'
#' @return
#' @export
#'
#' @examples
write_parsed_output = function(vcf, vcf_path, parsed_fp){
  metainf = parse_vcfpath_return_metainformation(vcfpath)
  analysis_name = metainf$IonReporter$AnalysisName
  ## Output files
  if(nrow(vcf) > 0){
    vcf = attach_ID(vcf, vcf_file = vcf_file, analysis_name = analysis_name)
    readr::write_tsv(vcf, file = parsed_fp$parsed_complete)
  }
  snv = parse_vcf_return_snv(vcf)
  if(nrow(snv) > 0){
    snv = attach_ID(snv, vcf_file = vcf_file, analysis_name = analysis_name)
    readr::write_tsv(snv, file = parsed_fp$parsed_snv)
  }
  cnv = parse_vcf_return_cnv(vcf)
  if(nrow(cnv) > 0){
    cnv = attach_ID(cnv, vcf_file = vcf_file, analysis_name = analysis_name)
    readr::write_tsv(cnv, file = parsed_fp$parsed_cnv)
  }

  ## Meta information
  write_out_META_information(metainf, filename = parsed_fp$parsed_info)

  ## Parsed filepath info
  write_parsed_fp_txt(parsed_fp)
}
