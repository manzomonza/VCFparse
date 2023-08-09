## VCFparse pipeline
library(VCFparse)
library(VariantAnnotation)
library(yaml)

print("Packages: loaded")


################################################# OPTPARSE #################################################
library(optparse)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="vcf file (path)", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

print(opt$file)

#################################################

vcfpath = opt$file

## Retrieve dirname
filepaths = generate_paths(vcfpath = vcfpath)
## Filepaths

## Read in file
vcf = VCFparse::readVCF(vcfpath)
# Parse FUNC
if(nrow(vcf) > 0){
  vcf = combine_orig_with_FUNC_extracts(vcf)
}
  ## remove double columns
vcf =  dplyr::select(vcf, -contains(".1"))
# Generate complete table
vcf$variant_type = gsub("[^[:alnum:] ]", "", vcf$variant_type)
vcf$protein = gsub("\\[|\\]", "", vcf$protein)
vcf$transcript = gsub("\\[|\\]", "", vcf$transcript)
table(vcf$protein)
complete_file = dplyr::filter(vcf, variant_type != 'synonymous' & alt != "<CNV>")

## Remove zero AF entries, especially important for Genexus
if("AF" %in% colnames(complete_file)){
  complete_file = dplyr::filter(complete_file, AF != 0)
}

# tibble::as_tibble(vcf) |>
#   dplyr::mutate_all(as.character) |>
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

## Output files
if(nrow(complete_file) > 0){
  readr::write_tsv(complete_file, file = filepaths$path_file_complete)
}
if(nrow(cnv_rows) > 0){
  readr::write_tsv(cnv_rows, file = filepaths$path_file_cnv)
}

write_out_META_information(metainf, filename = filepaths$path_file_info)



