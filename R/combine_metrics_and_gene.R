#' Title
#'
#' @param vcfinfo
#'
#' @return
#' @export
#'
#' @examples
combine_metric_gene_info = function(vcfinfo){
  df = data.frame(entries = unlist(stringr::str_split(vcfinfo, pattern = ";")))
  df$entries = gsub("'",'', df$entries)
  df_list = list(dataf = data.frame(entries =df[-grep("FUNC=", df$entries),]),
                 FUNC = df[grep("FUNC=", df$entries),])
  df_list$dataf = tidyr::separate(df_list$dataf, col = 'entries', into = c("parameter","value"), sep = "=")
  df_list$dataf = tidyr::pivot_wider(df_list$dataf, names_from = parameter, values_from = value)
  df_list$FUNC = convert_entries_to_table(df_list$FUNC)

  res = cbind(df_list$dataf, df_list$FUNC)
  return(res)
}




#' Add rowid to gene metric table
#'
#' @param vcf_row
#'
#' @return
#' @export
#'
#' @examples
combine_rowid = function(vcf_row){
  gene_metric = combine_metric_gene_info(vcf_row$INFO)
  gene_metric$rowid = vcf_row$rowid
  return(gene_metric)
}

#' Extract VCF standard columns. Additionally extract sample name from column 9 and remove column
#'
#' @param vcf
#'
#' @return
#' @export
#'
#' @examples
standard_columns = function(vcf){
  vcf = dplyr::select(vcf, -INFO)
  vcf$sample_id = colnames(vcf)[9]
  vcf = dplyr::select(vcf, -9,-FORMAT)
  return(vcf)
}


#' Extract all values from FORMAT column
#'
#' @param vcf_row
#'
#' @return
#' @export
#'
#' @examples
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

#' Combines function call
#'
#' @param vcfrow
#'
#' @return
#' @export
#'
#' @examples
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
