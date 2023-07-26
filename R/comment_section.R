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
  if(nrow(cm == 1)){
    value = cm$value
  }else{
    value = NA
  }
  return(value)
}

## IR_info

## sample info


## seq QC
pull_comment_value(cm, stringoi = 'percent_aligned_reads')
pull_comment_value(cm, stringoi = 'percent_non_zero_amplicons')
pull_comment_value(cm, stringoi = 'total_read_count')
pull_comment_value(cm, stringoi = 'median_reads_per_amplicon')
pull_comment_value(cm, stringoi = 'mapd')

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
  meta_info = list(IonReporter = list(AnalysisName = analysis_name,
                                      Software_version = pull_comment_value(comment_section, stringoi = 'IonReporterSoftwareVersion'),
                                      Export_version = pull_comment_value(comment_section, stringoi = 'IonReporterExportVersion'),
                                      Workflow_version = pull_comment_value(comment_section, stringoi = 'IonReporterWorkflowVersion'),
                                      WorkflowName = pull_comment_value(comment_section, stringoi = 'IonReporterWorkflowName')),

                   SampleInfo = list(disease = pull_comment_value(comment_section, stringoi = 'sampleDiseaseType'),
                                     tumor_cellularity_manual = pull_comment_value(comment_section, stringoi = 'manually_input_percent_tumor_cellularity'),
                                     tumor_cellularity_calculated = pull_comment_value(comment_section, stringoi = 'calculated_tumor_cellularity')))
  return(meta_info)
}


write_out_META_information = function(meta_information, basedir){
  filestring = paste0(basedir,'/parsed_output/meta_information.yaml')
  out = as.yaml(meta_information)
  yaml::write_yaml(out, file = filestring)
}

