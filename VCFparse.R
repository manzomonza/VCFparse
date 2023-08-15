## VCFparse pipeline
library(VCFparse)
library(VariantAnnotation)
library(yaml)

print("Packages: loaded")


################################################# OPTPARSE
library(optparse)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="vcf file (path)", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

print(opt$file)

#################################################

vcfpath = opt$file

## Retrieve dirname and generate filepaths
filepaths = parsed_filepaths(outputpath)


## Read in file
vcf = VCFparse::readVCF(vcfpath)
# Parse FUNC
if(nrow(vcf) > 0){
  vcf = combine_orig_with_FUNC_extracts(vcf)
}
  ## remove double columns
vcf =  dplyr::select(vcf, -contains(".1"))
vcf =  dplyr::select(vcf, -FUNC)

# Generate tables
vcf$variant_type = gsub("[^[:alnum:] ]", "", vcf$variant_type)
vcf$protein = gsub("\\[|\\]", "", vcf$protein)
vcf$transcript = gsub("\\[|\\]", "", vcf$transcript)
vcf = dplyr::relocate(vcf, rowid)


# SNV table
snv = dplyr::filter(vcf, variant_type != 'synonymous' & alt != "<CNV>")

## Remove zero AF entries, especially important for Genexus
if("AF" %in% colnames(snv)){
  snv = dplyr::filter(snv, AF != 0)
}
snv = dplyr::relocate(snv, rowid)
#   dplyr::mutate_all(as.character) |>
# tibble::as_tibble(vcf) |>
#   tidyr::pivot_longer(-gene)|>
#   dplyr::filter(grepl("\\[", value)) |>
#   dplyr::filter(name != 'FUNC')


# CNV entries only
cnv_rows = dplyr::filter(vcf, alt == "<CNV>")
if(nrow(cnv_rows) > 0){
  cnv_rows = cnv_parse(cnv_rows)
  cnv_rows = dplyr::select(cnv_rows, -origPos)
  cnv_rows = dplyr::arrange(cnv_rows, desc(RAW_CN))
}

# COMMENT SECTION
cm = vcf_comment_section(vcfpath = vcfpath )
metainf = aggregate_META_information(cm)


vcf_file = metainf$vcf_path
analysis_name = metainf$IonReporter$AnalysisName
## Output files
if(nrow(vcf) > 0){
  vcf = attach_ID(vcf, vcf_file = vcf_file, analysis_name = analysis_name)
  readr::write_tsv(vcf, file = filepaths$path_file_complete)
}

if(nrow(snv) > 0){
  snv = attach_ID(snv, vcf_file = vcf_file, analysis_name = analysis_name)
  readr::write_tsv(snv, file = filepaths$path_file_snv)
}

if(nrow(cnv_rows) > 0){
  cnv_rows = attach_ID(cnv_rows, vcf_file = vcf_file, analysis_name = analysis_name)
  readr::write_tsv(cnv_rows, file = filepaths$path_file_cnv)
}

write_out_META_information(metainf, filename = filepaths$path_file_info)

filepaths[1:5] = lapply(filepaths[1:5], function(x) ifelse(file.exists(x), x, NA))


filepath_df = data.frame(file = names(filepaths),
           filepath  = unlist(unname(filepaths)))
filepath_df = attach_ID(filepath_df, vcf_file = vcf_file, analysis_name = analysis_name)

readr::write_tsv(filepath_df, filepaths$path_filepath)
