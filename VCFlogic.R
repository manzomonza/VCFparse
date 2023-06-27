vcffiles = list.files(path = "VCF", pattern = '*.vcf', recursive = TRUE, full.names = TRUE)
vcffiles = grep("/Variants/", vcffiles, value = TRUE, invert = TRUE)

vcfs = lapply(vcffiles, read_vcf)
parsed_vcf = lapply(vcfs[90:95], function(x) dplyr::bind_rows(lapply(x, combine_function_calls)))

vcfs[93]


parsed_vcf = parsed_vcf[which(unlist(lapply(parsed_vcf, function(x) nrow(x))) > 0)]
parsed_vcf = lapply(parsed_vcf, change_col_name)
interprets = lapply(parsed_vcf[4:5], interpretation_table)

sapply(interprets[[i]]$protein, amino_acid_conversion_three_to_one)

write_output_tables(parsed_vcf[[4]], "test/testpath")

vcfs[[93]][[2]]$Y341_1_B2023.5338

interp = interpretation_table(parsed_vcf[[4]])
met_tab = metric_table(parsed_vcf[[4]])
var_tab = variant_table(parsed_vcf[[4]])
reformat = dplyr::bind_rows(lapply(vcfs[[93]], function(x) reformat_INFO(x$INFO)))
excel_file = list(var_tab, met_tab, interp, reformat)
names(excel_file) <- c("Variants", "Metadata", "Interpretation", "Metrics")
writexl::write_xlsx(excel_file, path = "Y341_example_output.xlsx")






reformatted = colnames(reformat_INFO(vcfs[[93]][[2]]$INFO))
parsed = colnames(parsed_vcf[[4]])

parsed[!parsed %in% reformatted]

lapply(vcfs[[i]], function(x) x$INFO)

all_em_infos = unlist(lapply(vcfs, function(x) lapply(x, function(z) collect_AF(z$INFO))))

all_em_infos = readRDS('all_em_infos.RDS')
all_em_values = lapply(all_em_infos, extract_AF_values)
all_qs = (lapply(all_em_values, function(x) extract_AF_quantiles(x$df$values)))
#all_qs = all_qs[!unlist(lapply(all_qs, function(x) is.null(x)))]

all_qs

library(ggplot2)
wrong_numbers = which(unlist(lapply(all_qs, function(x) unname(x[3] > 0.4))))

all_qs[wrong_numbers]

unname(all_qs[[1]][3] > 0.5)

dplyr::bind_rows(all_qs) |>
  ggplot(aes(`50%`, `100%`)) +
  geom_point()

dplyr::bind_rows(all_qs) |>
  dplyr::filter(`50%` > .4)


dplyr::bind_rows(all_qs) |>
  dplyr::filter(`50%` > 0.4)

lapply(vcfs[[7]], function(x) x$INFO)
vcfinfo = vcfs[[93]][[2]]$INFO
dplyr::distinct(reformat_INFO(vcfinfo))
indexed = tidyr::pivot_longer(res, -pivot_col) |> as.data.frame()

indexed$pivot_col = 'pivot'
indexed$OMAPALT
tidyr::pivot_longer(indexed, -pivot_col)

