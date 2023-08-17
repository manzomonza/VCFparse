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

write_parsed_output(vcf = vcf, vcf_path = vcfpath, parsed_fp = parsed_fp  )


