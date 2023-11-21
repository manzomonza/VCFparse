#' Extract commented section from VCF
#'
#' @param vcfpath
#'
#' @return
#' @export
#'
#' @examples
vcf_comment_section = function(vcfpath){
  comment_section = readr::read_tsv(vcfpath, col_names = FALSE)
  comment_section = dplyr::filter(comment_section, grepl("##", X1))
  comment_section$X1 = gsub("##", '', comment_section$X1)
  comment_section$X1 = gsub('\\"', '', comment_section$X1)
  comment_section = tidyr::separate(comment_section, col = X1, into = c("parameter", "value"), extra = "merge", sep = "=")
  path_df = tibble::tibble(parameter = "path_vcf",
                       value = vcfpath)
  comment_section = dplyr::bind_rows(path_df, comment_section)
  return(comment_section)
}

#' Retrieve
#'
#' @param comment_section
#' @param stringoi
#'
#' @return
#' @export
#'
#' @examples
pull_comment_value = function(comment_section, stringoi){
  cm = dplyr::filter(comment_section, parameter == stringoi)
  if(nrow(cm) == 1){
    value = cm$value
  }else{
    value = NA
  }
  return(value)
}


#' Aggregate required information from VCF comment section
#'
#' @param comment_section
#'
#' @return
#' @export
#'
#' @examples
aggregate_META_information = function(comment_section){
  analysis_name = pull_comment_value(comment_section, stringoi = 'IonReporterAnalysisName')
  vcf_path = pull_comment_value(comment_section, stringoi = 'path_vcf')
  vcf_path = basename(vcf_path)
  meta_info = list(vcf_path = vcf_path,
                   IonReporter = list(AnalysisName = analysis_name,
                                      fileDate = pull_comment_value(comment_section,
                                                                    stringoi = 'fileDate'),
                                      fileUTCtime = pull_comment_value(comment_section,
                                                                       stringoi = 'fileUTCtime'),
                                      Software_version = pull_comment_value(comment_section,
                                                                            stringoi = 'IonReporterSoftwareVersion'),
                                      Export_version = pull_comment_value(comment_section,
                                                                          stringoi = 'IonReporterExportVersion'),
                                      reference = pull_comment_value(comment_section,
                                                                            stringoi = 'reference'),
                                      Workflow_version = pull_comment_value(comment_section,
                                                                            stringoi = 'IonReporterWorkflowVersion'),
                                      WorkflowName = pull_comment_value(comment_section,
                                                                        stringoi = 'IonReporterWorkflowName')),

                   SampleInfo = list(disease = pull_comment_value(comment_section,
                                                                  stringoi = 'sampleDiseaseType'),
                                     tumor_cellularity_manual = pull_comment_value(comment_section,
                                                                                   stringoi = 'manually_input_percent_tumor_cellularity'),
                                     tumor_cellularity_calculated = pull_comment_value(comment_section,
                                                                                       stringoi = 'calculated_tumor_cellularity')),

                   SequencingQC = list(percent_aligned_reads = pull_comment_value(comment_section,
                                                                                  stringoi = 'percent_aligned_reads'),
                                percent_non_zero_amplicons = pull_comment_value(comment_section,
                                                                                stringoi = 'percent_non_zero_amplicons'),
                                total_read_count = pull_comment_value(comment_section, stringoi = 'total_read_count'),
                                median_reads_per_amplicon = pull_comment_value(comment_section,
                                                                               stringoi = 'median_reads_per_amplicon'),
                                MAPD = pull_comment_value(comment_section,
                                                          stringoi = 'mapd'),
                                deamination_metric = pull_comment_value(comment_section,
                                                                        stringoi = 'deamination_metric'))
                   )
  return(meta_info)
}

#' Attach sample ID from original vcf path and analysis name to output dataframes.
#' These entries serve as unique identifiers across tables
#' @param dataframe_oi
#' @param sample_id
#' @param analysis_name
#'
#' @return
#' @export
#'
#' @examples
attach_ID = function(dataframe_oi, vcf_file, analysis_name){
  dataframe_oi$vcf_file = vcf_file
  dataframe_oi$analysis_name = analysis_name
  return(dataframe_oi)
}

