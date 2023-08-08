## VCFparse pipeline
library(VCFparse)
library(VariantAnnotation)
library(yaml)

print("Packages: loaded")


################################################# OPTPARSE #################################################
library(optparse)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name", metavar="character"))

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
vcf = combine_orig_with_FUNC_extracts(vcf)
## remove double columns
vcf =  dplyr::select(vcf, -contains(".1"))
# Generate complete table
complete_file = dplyr::filter(vcf, variant_type != 'synonymous' & alt != "<CNV>")
# CNV entries only
cnv_rows = dplyr::filter(vcf, alt == "<CNV>")
cnv_rows = cnv_parse(cnv_rows)

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



