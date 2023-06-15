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
  vcf$rowid = 1:nrow(vcf)
  vcf = dplyr::group_by(vcf, rowid)
  vcf = dplyr::group_split(vcf)
  return(vcf)
}

(probs[2])
vcfs = lapply(vcfs[1:10],vcf_read)

lapply(vcfs[1:10], parse_FUNC)
lapply(vcfs, parse_FUNC)

lapply(vcfs[[1]], parse_FUNC)
variant_info_to_row(vcfs[[1]][[3]]$INFO)



parse_selector = function(infostring){
  if(!grepl('SVTYPE=CNV', infostring)){
    result = combine_metric_gene_info(infostring)
    return(result)
  }
}


combine_metric_gene_info = function(vcfinfo){
  df = data.frame(entries = unlist(stringr::str_split(vcfinfo, pattern = ";")))
  df$entries = gsub("'",'', df$entries)
  df_list = list(dataf = data.frame(entries =df[-grep("^FUNC=", df$entries),]),
                 FUNC = df[grep("^FUNC=", df$entries),])
  df_list$dataf = tidyr::separate(df_list$dataf, col = 'entries', into = c("parameter","value"), sep = "=")
  df_list$dataf = tidyr::pivot_wider(df_list$dataf, names_from = parameter, values_from = value)
  df_list$FUNC = parse_FUNC(df_list$FUNC)

  res = cbind(df_list$dataf, df_list$FUNC)
  return(res)
}

dplyr::bind_rows(lapply(vcfs[[2]], function(x) parse_selector(x$INFO)))
lapply(vcfs[[1]], function(x) grepl('TYPE=snp', x$INFO))

parse_selector(vcfs[[1]][[1]]$INFO)

vcfs[[1]][[5]]$INFO


combine_rowid = function(vcf_row){
  gene_metric = parse_selector(vcf_row$INFO)
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
  if(!grepl('SVTYPE=CNV', vcfrow$INFO)){
    format_vals = extract_format_values(vcfrow)
    std_cols = standard_columns(vcfrow)
    info_vals = combine_rowid(vcfrow)
    first_join = dplyr::full_join(std_cols, format_vals)
    second_join = dplyr::full_join(first_join, info_vals)
    return(second_join)
  }
}


extract_format_values(vcfs[[1]][[13]])
standard_columns(vcfs[[1]][[13]])
combine_rowid(vcfs[[1]][[13]])


colnames(complete)
variant_display = function(combined_output){
  output = dplyr::select(combined_output, gene, transcript, chr, coding,protein,exon,location, AF)
  return(output)
}



complete = dplyr::bind_rows(lapply(vcfs[[10]], combine_function_calls)) |> dplyr::relocate(gene)
colnames(complete)[2] = 'chr'
variant_display(complete)


