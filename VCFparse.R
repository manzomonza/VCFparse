## VCFparse pipeline
library(VCFparse)
library(VariantAnnotation)


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
vcf = combine_orig_with_FUNC_extracts(vcf)
print(colnames(vcf))
## p
complete_file = dplyr::filter(vcf, variant_type != 'synonymous' & alt != "<CNV>")
cnv_rows = dplyr::filter(vcf, alt == "<CNV>")
cnv_rows = cnv_parse(cnv_rows)
## Output files
readr::write_tsv(complete_file, file = filepaths$path_file_complete)
readr::write_tsv(cnv_rows, file = filepaths$path_file_cnv)



