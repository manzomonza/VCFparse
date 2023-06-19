## VCF parsing logic

#' Read VCF file in skipping ## commented entries
#'
#' @param vcf_filepath
#'
#' @return
#' @export
#'
#' @examples
vcf_read = function(vcf_filepath){
  vcf = readr::read_tsv(vcf_filepath, comment = "##")
  vcf = dplyr::filter(vcf,!grepl('AF=0;', INFO))
  vcf$rowid = 1:nrow(vcf)
  vcf = dplyr::group_by(vcf, rowid)
  vcf = dplyr::group_split(vcf)
  return(vcf)
}

vcfs = lapply(vcffiles[180:195],vcf_read)
lengths(vcfs)
lapply(vcfs[1:10], parse_FUNC)
lapply(vcfs, parse_FUNC)

lapply(vcfs[[1]], parse_FUNC)
variant_info_to_row(vcfs[[1]][[3]]$INFO)



combine_metric_gene_info = function(vcfinfo){
  df = data.frame(entries = unlist(stringr::str_split(vcfinfo, pattern = ";")))
  df$entries = gsub("'",'', df$entries)
  df_list = list(dataf = data.frame(entries =df[-grep("FUNC=", df$entries),]),
                 FUNC = df[grep("FUNC=", df$entries),])
  df_list$dataf = tidyr::separate(df_list$dataf, col = 'entries', into = c("parameter","value"), sep = "=")
  df_list$dataf = tidyr::pivot_wider(df_list$dataf, names_from = parameter, values_from = value)
  df_list$FUNC = parse_FUNC(df_list$FUNC)

  res = cbind(df_list$dataf, df_list$FUNC)
  return(res)
}

dplyr::bind_rows(lapply(vcfs[[8]], function(x) combine_metric_gene_info(x$INFO)))
lapply(vcfs[[1]], function(x) grepl('TYPE=snp', x$INFO))

parse_selector(vcfs[[1]][[1]]$INFO)

vcfs[[1]][[5]]$INFO


combine_rowid = function(vcf_row){
  gene_metric = combine_metric_gene_info(vcf_row$INFO)
  gene_metric$rowid = vcf_row$rowid
  return(gene_metric)
}

standard_columns = function(vcf){
  vcf = dplyr::select(vcf, -INFO)
  vcf$sample_id = colnames(vcf)[9]
  vcf = dplyr::select(vcf, -9,-FORMAT)
  return(vcf)
  }


extract_format_values = function(vcf_row){
  format_parameters = vcf_row$FORMAT
  format_parameters = stringr::str_split(format_parameters, pattern = ":")[[1]]
  format_values = vcf_row[,ncol(vcf_row)-1]
  format_values = stringr::str_split(format_values, pattern = ":")[[1]]
  df = data.frame(parameter = format_parameters, value = format_values)
  df = tidyr::pivot_wider(df, names_from = parameter, values_from = value)
  df$rowid = vcf_row$rowid
  return(df)
}

combine_function_calls = function(vcfrow){
  if(!grepl('AF=0;', vcfrow$INFO)){
   if(!grepl('SVTYPE=CNV', vcfrow$INFO)){
    format_vals = extract_format_values(vcfrow)
    std_cols = standard_columns(vcfrow)
    info_vals = combine_rowid(vcfrow)
    first_join = dplyr::full_join(std_cols, format_vals)
    second_join = dplyr::full_join(first_join, info_vals)
    return(second_join)
    }
  }
}


colnames(complete)

vcfrow = vcfs[[1]][[1]]
format_vals = extract_format_values(vcfrow)
std_cols = standard_columns(vcfrow)
info_vals = combine_rowid(vcfrow)
first_join = dplyr::full_join(std_cols, format_vals)
second_join = dplyr::full_join(first_join, info_vals)
complete = dplyr::bind_rows(lapply(vcfs[[3]], combine_function_calls)) |> dplyr::relocate(gene)

parsed_vcf = lapply(vcfs, function(x) dplyr::bind_rows(lapply(x, combine_function_calls)))
vcffiles[194:200]
colnames(complete)[2] = 'chr'
parsed_vcf = lapply(parsed_vcf, change_col_name)


vcfs[[7]][[1]]$INFO

lapply(vcfs[[7]], extract_format_values)
lapply(vcfs[[7]], standard_columns)
lapply(vcfs[[7]], combine_rowid)

library(NGSannotation)

colnames(parsed_vcf[[i]])
vcfs[[12]]
length(parsed_vcf)
parsed_vcf[12]
lapply(parsed_vcf[11:12], variant_table)
metric_table(complete) |> dplyr::arrange(desc(as.numeric(AF)))

parsed_vcf[4]
parsed_vcf = parsed_vcf[which(unlist(lapply(parsed_vcf, function(x) nrow(x))) > 0)]
interprets = lapply(parsed_vcf, interpretation_table)

i = 14
interprets[[i]]$protein

sapply(interprets[[i]]$protein, amino_acid_conversion_three_to_one)

