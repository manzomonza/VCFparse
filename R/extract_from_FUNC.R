### character vector should be 'cleaned' already  via
### character_vector = gsub("'","",character_vector)


#' Remove unwanted characters from vcf FUNC character vector
#'
#' @param character_vector
#'
#' @return
#' @export
#'
#' @examples
clean_FUNC_character = function(character_vector){
  character_vector = gsub("\\[\\{|\\}\\]|'",'',character_vector)
  return(character_vector)
}

#' Extract Gene symbol from FUNC vcf field
#'
#' @param character_vector
#'
#' @return
#' @export
#'
#' @examples
extract_FUNC_gene = function(character_vector){
  string_oi = grep("gene:", character_vector, value = TRUE)
  if(identical(string_oi,character(0))){
    return(NA)
  }else{
  #string_oi = gsub("'","",string_oi)
  string_oi = stringr::str_remove(string_oi, pattern = "gene:")
  return(string_oi)
  }
}


#' Extract coding information from FUNC vcf field
#'
#' @param character_vector
#'
#' @return
#' @export
#'
#' @examples
extract_FUNC_coding = function(character_vector){
  string_oi = grep("coding:", character_vector, value = TRUE)
  if(identical(string_oi,character(0))){
    return(NA)
  }else{
    #string_oi = gsub("'","",string_oi)
    string_oi = stringr::str_remove(string_oi, pattern = "coding:")
    return(string_oi)
  }
}



#' Extract transcript information from FUNC vcf field
#'
#' @param character_vector
#'
#' @return
#' @export
#'
#' @examples
extract_FUNC_transcript = function(character_vector){
  string_oi = grep("transcript:", character_vector, value = TRUE)
  if(identical(string_oi,character(0))){
    return(NA)
  }else{
    #string_oi = gsub("'","",string_oi)
    string_oi = stringr::str_remove(string_oi, pattern = "transcript:")
    return(string_oi)
  }
}



#' Extract protein information from FUNC vcf field
#'
#' @param character_vector
#'
#' @return
#' @export
#'
#' @examples
extract_FUNC_protein = function(character_vector){
  string_oi = grep("protein:", character_vector, value = TRUE)
  if(identical(string_oi,character(0))){
    return(NA)
  }else{
    #string_oi = gsub("'","",string_oi)
    string_oi = stringr::str_remove(string_oi, pattern = "protein:")
    return(string_oi)
  }
}



#' Extract origPos information from FUNC vcf field
#'
#' @param character_vector
#'
#' @return
#' @export
#'
#' @examples
extract_FUNC_origPos = function(character_vector){
  string_oi = grep("origPos:", character_vector, value = TRUE)
  if(identical(string_oi,character(0))){
    return(NA)
  }else{
    #string_oi = gsub("'","",string_oi)
    string_oi = stringr::str_remove(string_oi, pattern = "origPos:")
    return(string_oi)
  }
}

#' Extract origRef information from FUNC vcf field
#'
#' @param character_vector
#'
#' @return
#' @export
#'
#' @examples
extract_FUNC_origRef = function(character_vector){
  string_oi = grep("origRef:", character_vector, value = TRUE)
  if(identical(string_oi, character(0))){
    return(NA)
  }else{
    #string_oi = gsub("'","",string_oi)
    string_oi = stringr::str_remove(string_oi, pattern = "origRef:")
    return(string_oi)
  }
}

#' Extract origAlt information from FUNC vcf field
#'
#' @param character_vector
#'
#' @return
#' @export
#'
#' @examples
extract_FUNC_origAlt = function(character_vector){
  string_oi = grep("origAlt:", character_vector, value = TRUE)
  if(identical(string_oi, character(0))){
    return(NA)
  }else{
    #string_oi = gsub("'","",string_oi)
    string_oi = stringr::str_remove(string_oi, pattern = "origAlt:")
    return(string_oi)
  }
}





#' Extract location information from FUNC vcf field
#'
#' @param character_vector
#'
#' @return
#' @export
#'
#' @examples
extract_FUNC_location = function(character_vector){
  string_oi = grep("location:", character_vector, value = TRUE)
  if(identical(string_oi,character(0))){
    return(NA)
  }else{
    #string_oi = gsub("'","",string_oi)
    string_oi = stringr::str_remove(string_oi, pattern = "location:")
    return(string_oi)
  }
}

#' Extract location information from FUNC vcf field
#'
#' @param character_vector
#'
#' @return
#' @export
#'
#' @examples
extract_FUNC_exon = function(character_vector){
  string_oi = grep("exon:", character_vector, value = TRUE)
  if(identical(string_oi,character(0))){
    return(NA)
  }else{
    #string_oi = gsub("'","",string_oi)
    string_oi = stringr::str_remove(string_oi, pattern = "exon:")
    return(string_oi)
  }
}


#' Extract location information from FUNC vcf field
#'
#' @param character_vector
#'
#' @return
#' @export
#'
#' @examples
extract_FUNC_funct = function(character_vector){
  string_oi = grep("function:", character_vector, value = TRUE)
  if(identical(string_oi,character(0))){
    return(NA)
  }else{
    #string_oi = gsub("'","",string_oi)
    string_oi = stringr::str_remove(string_oi, pattern = "function:")
    return(string_oi)
  }
}

#' Call all vcf FUNC extract functions
#'
#' @param character_vector
#'
#' @return
#' @export
#'
#' @examples
FUNC_extracts_to_df = function(character_vector){
  character_vector = clean_FUNC_character(character_vector)
  gene = extract_FUNC_gene(character_vector)
  transcript = extract_FUNC_transcript(character_vector)
  coding = extract_FUNC_coding(character_vector)
  protein = extract_FUNC_protein(character_vector)
  location = extract_FUNC_location(character_vector)
  exon = extract_FUNC_exon(character_vector)
  origPos = extract_FUNC_origPos(character_vector)
  origRef = extract_FUNC_origRef(character_vector)
  origAlt = extract_FUNC_origAlt(character_vector)
  funct = extract_FUNC_funct(character_vector)
  res_df = data.frame(gene = gene,
                transcript = transcript,
                variant_type = funct,
                coding = coding,
                protein = protein,
                location = location,
                exon = exon,
                origPos = origPos,
                origRef = origRef,
                origAlt = origAlt)
  return(res_df)
}



