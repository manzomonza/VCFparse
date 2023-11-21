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
path_snv = paste0(path_parsed_output, '/parsed_snv.tsv')
path_complete = paste0(path_parsed_output, '/parsed_complete.tsv')
path_fusions = paste0(path_parsed_output, '/parsed_fusions.tsv') # not implemented yet
path_yaml = paste0(path_parsed_output, '/parsed_vcf_info.yaml')

file_nrw_write = function(table_obj, filestring){
  if(nrow(table_obj) > 0){
    readr::write_tsv(table_obj, file = filestring)
  }
}

cm = vcf_comment_section(vcfpath = vcfpath )
metainf = aggregate_META_information(cm)
analysis_name = metainf$IonReporter$AnalysisName
vcf_file = basename(vcfpath)

## Output files
vcf_w = attach_ID(vcf, vcf_file = vcf_file, analysis_name = analysis_name)
vcf_w = remove_multientry_leftover_characters(vcf_w)
file_nrw_write(vcf_w, filestring = path_complete)

snv = parse_vcf_return_snv(vcf)
snv = remove_multientry_leftover_characters(snv)
if(nrow(snv) > 0){
  snv = attach_ID(snv, vcf_file = vcf_file, analysis_name = analysis_name)
}
file_nrw_write(snv, filestring = path_snv)


cnv = parse_vcf_return_cnv(vcf)
if(nrow(cnv) > 0){
  cnv = attach_ID(cnv, vcf_file = vcf_file, analysis_name = analysis_name)
}
file_nrw_write(cnv, filestring = path_cnv)

## Meta information
out_yaml = as.yaml(metainf)
yaml::write_yaml(out_yaml, file = path_yaml)

