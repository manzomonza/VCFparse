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
  metainf = parse_vcfpath_return_metainformation(vcf_path)
  analysis_name = metainf$IonReporter$AnalysisName
  vcf_file = basename(vcf_path)
  ## Output files
  if(nrow(vcf) > 0){
    vcf_w = attach_ID(vcf, vcf_file = vcf_file, analysis_name = analysis_name)
    vcf_w = remove_multientry_leftover_characters(vcf_w)
    readr::write_tsv(vcf_w, file = parsed_fp$parsed_complete)
  }
  snv = parse_vcf_return_snv(vcf)
  snv = remove_multientry_leftover_characters(snv)
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
  write_parsed_fp_txt(parsed_fp, vcf_file = vcf_file, analysis_name = analysis_name)
}
