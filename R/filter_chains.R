# Scope of functions is to generate specific filter chain based
# on inputs available from meta_information.yaml

#' Title
#'
#' @param meta_info
#'
#' @return
#' @export
#'
#' @examples
panel_detection = function(meta_info){
  panel_name = meta_info$IonReporter$WorkflowName
  # INTEGRATE function that assigns panel name based on Look up table
  if(panel_name == "RDX Ampliseq Germline Template HD"){
    panel = 'Genexus DNA'
  }else if(panel_name == "RDX Ampliseq Fusions Template"){
    panel = 'Genexus RNA'
  }else{
    panel = "generic"
  }
  return(panel)
}


# filterChain_gnxsDNA = function(vcf){
#   vcf
# }
# vavcf = VariantAnnotation::readVcfAsVRanges(testvcf)
# vavcf = tibble::as_tibble(vavcf)
# filter_column = readr::read_tsv(testvcf, comment = "##")$FILTER
# vavcf$FILTER = filter_column


