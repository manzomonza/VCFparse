# Scope of functions is to generate specific filter chain based
# on inputs available from meta_information.yaml

panel_detection = function(meta_info){
  panel_name = meta_info$IonReporter$WorkflowName
  # INTEGRATE function that assigns panel name based on Look up table
  if(panel_name == "RDX Ampliseq Germline Template HD"){
    panel = 'Genexus'
  }else{
    panel = "generic"
  }
  return(panel)
}


meta_info = aggregate_META_information(cm)
panel_detection(meta_info)
