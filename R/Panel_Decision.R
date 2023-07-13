### panel decision

vcf_comment_section = function(vcfpath){
  comment_section = readr::read_tsv(vcfpath, col_names = FALSE)
  comment_section = dplyr::filter(comment_section, grepl("##", X1))
  comment_section$X1 = gsub("##", '', comment_section$X1)
  comment_section$X1 = gsub('\\"', '', comment_section$X1)
  comment_section = tidyr::separate(comment_section, col = X1, into = c("parameter", "value"), sep = "=")
  return(comment_section)
}
glob2rx('\\"tvc 5')
vcf_comment_section(testfiles[1])

panel_string = function(vcf_comments)
