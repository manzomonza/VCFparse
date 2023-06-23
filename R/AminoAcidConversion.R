aminocode  <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr",
                "Trp", "Tyr", "Val", "Sec", "Pyl", "Asx", "Xle", "Glx", "Xaa")
names(aminocode) <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T",
                      "W", "Y", "V", "U", "O", "B", "J", "Z","X")
aminocode['\\*']  <- "Ter"
aminocode['dup']  <- "dup"


#' ###### Change three to one letter code
#'
#' @param amino_acid_change_entry
#'
#' @return
#' @export
#'
#' @examples
amino_acid_conversion_three_to_one <- function(amino_acid_change_entry){
  amino_acid_change_parse = amino_acid_change_entry
  if(!is.na(amino_acid_change_entry) & !grepl("p\\.\\?", amino_acid_change_entry) ){
    for(i in seq_along(aminocode)){
      amino_acid_change_parse <- gsub(aminocode[i], names(aminocode)[i], amino_acid_change_parse, ignore.case = FALSE)
    }
    if(!grepl("p.", amino_acid_change_parse)){
      amino_acid_change_parse <- paste0("p.", amino_acid_change_parse)
    }
    return(amino_acid_change_parse)
  }else{
    return(amino_acid_change_entry)
  }
}


