

vcffiles = list.files(path = "VCF", pattern = '*.vcf', recursive = TRUE, full.names = TRUE)
vcfs = lapply(vcffiles[1:5], read_vcf)
parsed_vcf = lapply(vcfs, function(x) dplyr::bind_rows(lapply(x, combine_function_calls)))
parsed_vcf = parsed_vcf[which(unlist(lapply(parsed_vcf, function(x) nrow(x))) > 0)]
parsed_vcf = lapply(parsed_vcf, change_col_name)
interprets = lapply(parsed_vcf[4:5], interpretation_table)

i = 14
interprets[[i]]$protein

sapply(interprets[[i]]$protein, amino_acid_conversion_three_to_one)

