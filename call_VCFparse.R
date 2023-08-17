## VCFparse pipeline
library(VCFparse)
library(VariantAnnotation)
library(yaml)

print("Packages: loaded")


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


#   dplyr::mutate_all(as.character) |>
# tibble::as_tibble(vcf) |>
#   tidyr::pivot_longer(-gene)|>
#   dplyr::filter(grepl("\\[", value)) |>
#   dplyr::filter(name != 'FUNC')




# COMMENT SECTION

#write_parsed_output(vcf = vcf, vcf_path = vcfpath, parsed_fp = parsed_fp  )

print("TROUBLESHOOT")

metainf = parse_vcfpath_return_metainformation(vcfpath)
print(paste0("metainf of vcf:", metainf))

analysis_name = metainf$IonReporter$AnalysisName
print(paste0("base name of analysis_name:", analysis_name))

vcf_file = basename(vcfpath)
print(paste0("base name of vcfpath:", vcf_file))

print("complete file")
if(nrow(vcf) > 0){
  vcf_w = attach_ID(vcf, vcf_file = vcf_file, analysis_name = analysis_name)
  readr::write_tsv(vcf_w, file = parsed_fp$parsed_complete)
}
print("snv file")

snv = parse_vcf_return_snv(vcf)
if(nrow(snv) > 0){
  snv = attach_ID(snv, vcf_file = vcf_file, analysis_name = analysis_name)
  readr::write_tsv(snv, file = parsed_fp$parsed_snv)
}
print("cnv file")

cnv = parse_vcf_return_cnv(vcf)
if(nrow(cnv) > 0){
  cnv = attach_ID(cnv, vcf_file = vcf_file, analysis_name = analysis_name)
  readr::write_tsv(cnv, file = parsed_fp$parsed_cnv)
}

## Meta information
print("meta file")
write_out_META_information(metainf, filename = parsed_fp$parsed_info)

print("filepaths file")
## Parsed filepath info
write_parsed_fp_txt(parsed_fp, vcf_file = vcf_file, analysis_name = analysis_name)





