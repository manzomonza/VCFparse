## VCFparse pipeline


################################################# OPTPARSE
# library(optparse)
#
# option_list = list(
#   make_option(c("-f", "--file"), type="character", default=NULL,
#               help="vcf file (path)", metavar="character"))
#
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser)
#
# print(opt$file)

#################################################
##
path_parsed_output = paste0(analysis_dir, "/parsed_output") # analysis_dir gets inherited from genie.R
if(!dir.exists(path_parsed_output)){
  dir.create(path_parsed_output)
}

# parsed_fusions = paste0(path_parsed_output, '/parsed_fusions.tsv') # not implemented yet




metainf = parse_vcfpath_return_metainformation(vcfpath)
analysis_name = metainf$IonReporter$AnalysisName
vcf_file = basename(vcfpath)

## Output files
if(nrow(vcf) > 0){
  vcf_w = attach_ID(vcf, vcf_file = vcf_file, analysis_name = analysis_name)
  vcf_w = remove_multientry_leftover_characters(vcf_w)
  readr::write_tsv(vcf_w, file = paste0(path_parsed_output, '/parsed_complete.tsv'))
}
snv = parse_vcf_return_snv(vcf)
snv = remove_multientry_leftover_characters(snv)
if(nrow(snv) > 0){
  snv = attach_ID(snv, vcf_file = vcf_file, analysis_name = analysis_name)
  readr::write_tsv(snv, file = paste0(path_parsed_output, '/parsed_snv.tsv'))
}
cnv = parse_vcf_return_cnv(vcf)
if(nrow(cnv) > 0){
  cnv = attach_ID(cnv, vcf_file = vcf_file, analysis_name = analysis_name)
  readr::write_tsv(cnv, file = paste0(path_parsed_output, '/parsed_cnv.tsv'))
}

## Meta information
write_out_META_information(metainf, filename = paste0(path_parsed_output, '/parsed_vcf_info.yaml'))


