## FAFO

vcffiles = list.files(path = "VCF", recursive = TRUE, pattern = '.*.vcf', full.names = TRUE)
vcffiles = grep("SmallVariants", vcffiles, invert = TRUE, value = TRUE)
vcffiles = grep("Fusions", vcffiles, invert = TRUE, value = TRUE)
vcffiles = grep("RNA", vcffiles, invert = TRUE, value = TRUE)
vcffiles = grep("Non-Filtered", vcffiles, invert = TRUE, value = TRUE)
vcffiles = grep("Filtered", vcffiles, invert = FALSE, value = TRUE)


testfile = grep("Y676",vcfs, value = TRUE)

library(tidyverse)
split_alt_type = function(vcf){
  vcf_snv = dplyr::filter(vcf, ALT != "<CNV>")
  vcf_cnv = dplyr::filter(vcf, ALT == "<CNV>")
  return(list(snv = vcf_snv,
              cnv = vcf_cnv))
}


info_parse = function(info_string){
  info = stringr::str_split(info_string, pattern = ";")[[1]]
  return(info)
}

element_sep = function(element_string){
  element = stringr::str_split(element_string, pattern = "\\=")
  return(element)
}


extract_info_parse = function(info_string){
  infos = info_parse(info_string)
  infos = grep('FUNC=\\[\\{',infos, invert = TRUE, value = TRUE)
  df = data.frame(parameter = infos, value = NA)
  df = tidyr::separate(df, col = parameter, into = c("parameter", 'value'), sep = "=", fill = 'right')
  df = tidyr::pivot_wider(df, names_from = parameter, values_from = value)
  return(df)
}


extract_format_values = function(vcf_row){
  format_parameters = vcf_row$FORMAT
  format_parameters = stringr::str_split(format_parameters, pattern = ":")[[1]]
  format_values = vcf_row[,ncol(vcf_row)]
  format_values = stringr::str_split(format_values, pattern = ":")[[1]]
  df = data.frame(parameter = format_parameters, value = format_values)
  df = tidyr::pivot_wider(df, names_from = parameter, values_from = value)
  return(df)
}

####
i = 123

combine_extracts = function(vcf_row){
  vcf = vcf_row[,1:7]
  #fv = extract_format_values(vcf_row)
  ip = extract_info_parse(vcf_row$INFO)
  func = extract_FUNC(vcf_row$INFO)
  vcf = cbind(vcf, ip, func)
  return(vcf)
}


integrate_extracts = function(vcf_path){
  vc = readr::read_tsv(vcf_path, comment = "##")
  vc = dplyr::filter(vc, !grepl("AF=0;", INFO))
  vc = dplyr::filter(vc, FILTER != "FAIL")
  vc = split_alt_type(vc)
  snv = vc$snv
  if(nrow(snv) >0){
    snv$rowid = 1:nrow(snv)
    if(length(snv$rowid) > 50){
      vc_table = lapply(snv$rowid, function(x) combine_extracts(snv[x,]))
      }else{
        vc_table = lapply(snv$rowid, function(x) combine_extracts(snv[x,]))
        }
    vc_table = lapply(vc_table, function(x) x %>% dplyr::mutate_all(as.character))
  return(vc_table)
  }
}


start = Sys.time()
res2 = lapply(vcfs[c(41,46)], function(x) dplyr::bind_rows(integrate_extracts(x)))
end = Sys.time()
end - start

res2 = res2[which(unlist(lapply(res2, function(x) nrow(x) > 0)))]
all_col_names = lapply(res2, colnames)
(use_colnames = Reduce(intersect, all_col_names))


ress = lapply(res2, function(x) x |> dplyr::select(all_of(use_colnames)))
ress = dplyr::bind_rows(ress)

table(ress$location)

saveRDS(res, "Intersect_vcf_res.RDS")
res=readRDS("Intersect_vcf_res.RDS")

##############

res |>
  dplyr::group_by(`#CHROM`, POS,ALT,`function`)|>
  dplyr::count(sort = TRUE)

filtered_res = res |>
  dplyr::filter(`#CHROM` == "chr4" & POS == '1807894' |
                `#CHROM` == "chr3" & POS == '41266100')

filtered_res = res |>
  dplyr::filter(`#CHROM` == "chr7" & POS == '55228053' )

quantile(as.numeric(filtered_res$AF))

res |>
  dplyr::filter(FILTER == "NOCALL") |>
  ggplot(aes(as.numeric(AF))) +
  geom_histogram()


filtered_res = filtered_res |> dplyr::select(1:5,7,QUAL, AF:FXX,MLLD, QD:STBP, VARB) |>
  dplyr::mutate(across(.cols = QUAL:VARB, .fns = as.numeric))


filtered_mat = dplyr::select(filtered_res, FILTER,where(is.numeric)) |>
  dplyr::select(-SSEN, -SSEP, -SPD, -FR,-FDVR)

filtered_na = na.omit(filtered_mat)

toi = filtered_na[sample(1:nrow(filtered_na), 4000),]
toi = filtered_na

colnames(toi)
metadata = toi
metadata = dplyr::select(metadata, FILTER,QUAL, AF, QD,FSAF,FSAR,FSRF,FSRR, FRO, FAO) |> as.data.frame()
rownames(metadata) = paste0('v', 1:nrow(metadata))
metadata$artefact = ifelse(metadata$AF < 0.001, TRUE, FALSE)
metadata$artefact = ifelse(metadata$FILTER == "NOCALL", TRUE, FALSE)

######## CLUSTER
toi = filtered_na[sample(1:nrow(filtered_na), 4000),]
toi = dplyr::select(toi, -FILTER)

skimr::skim(toi)

table(metadata$artefact)
toi = scale(toi) |> t() |> as.matrix()

colnames(toi) = paste0('v', 1:ncol(toi))
########## PCA ##########
########## ########## ##########
skimr::skim(toi)

head(toi)

pcap = PCAtools::pca(toi, scale = FALSE, metadata = metadata)
PCAtools::biplot(pcap, showLoadings = TRUE, colby = 'artefact', lab = NULL, alpha = .3)
########## ########## ##########
colnames(res)
teststring = "AF=0.994426;AO=1775;DP=1795;FAO=1784;FDP=1794;FDVR=10;FR=.;FRO=10;FSAF=717;FSAR=1067;FSRF=7;FSRR=3;FWDB=-0.0137457;FXX=5.57103E-4;HRUN=4;HS_ONLY=0;LEN=1;MLLD=548.977;OALT=G;OID=.;OMAPALT=G;OPOS=2488153;OREF=A;PB=0.5;PBP=1;PPD=0;QD=71.9427;RBI=0.0140999;REFB=-0.00501424;REVB=0.0031407;RO=2;SAF=708;SAR=1067;SPD=0;SRF=2;SRR=0;SSEN=0;SSEP=0;SSSB=-0.0013484;STB=0.501727;STBP=0.062;TYPE=snp;VARB=4.25283E-6;FUNC=[{'transcript':'NR_037844.2','gene':'TNFRSF14-AS1','location':'intronic_nc','CLNACC1':'135349','CLNREVSTAT1':'no_assertion_provided','CLNSIG1':'not_provided','CLNID1':'rs4870'},{'origPos':'2488153','origRef':'A','normalizedRef':'A','gene':'TNFRSF14','normalizedPos':'2488153','normalizedAlt':'G','polyphen':'0.388','gt':'pos','codon':'AGA','coding':'c.50A>G','sift':'0.01','grantham':'26.0','transcript':'NM_003820.3','function':'missense','protein':'p.Lys17Arg','location':'exonic','origAlt':'G','exon':'1','CLNACC1':'135349','CLNREVSTAT1':'no_assertion_provided','CLNSIG1':'not_provided','CLNID1':'rs4870'}]"
stringr::str_extract(teststring, pattern = "(?<=FUNC=\\[\\{).*(?=\\}\\])")
########## PLOTS ##########
############################
library(ggplot2)
metadata$allele_count_ratio = (metadata$FAO+1) / (metadata$FRO+1)
metadata$strand_alt = (metadata$FSAF+1) / (metadata$FSAR+1)
metadata$strand_ref = (metadata$FSRF+1) / (metadata$FSRR+1)
metadata$strand_ratios = (metadata$strand_ref) / (metadata$strand_alt)
metadata |>
  ggplot(aes(artefact, QD)) +
  geom_violin() +
  geom_jitter(height = 0, width = .2,alpha = .2) +
  scale_y_log10()

metadata |>
  ggplot(aes(artefact, FE)) +
  geom_violin() +
  geom_jitter(height = 0, width = .2,alpha = .2) +
  scale_y_log10()

metadata |>
  ggplot(aes(FRO, FAO,col = artefact)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10()

metadata |>
  ggplot(aes(allele_count_ratio, strand_ratios,col = artefact)) +
  geom_point(alpha  = .3)  +
  scale_y_log10() +
  scale_x_log10()

############################

library(rpart)
library(rpart.plot)
library(rattle)

metadata$FE = metadata$QUAL * metadata$AF

nytree = rpart(artefact ~ ., data = dplyr::select(metadata, -c(AF, FILTER, allele_count_ratio,FAO)),
               method = 'class')
fancyRpartPlot(nytree)

############## LOG MODEL ##############
############## ############## ##############
log_model = glm(artefact~QD + strand_ratios, family = 'binomial', data = metadata)
summary(log_model)
probs = predict(log_model, metadata, type = 'response')
preds = ifelse(probs > 0.5, TRUE, FALSE)

metadata$pred = preds
metadata$probs = probs

metadata |>
  dplyr::filter(artefact != pred)

head(metadata)

############## UNDER CONSTRUCTION ##############
############## ############## ##############


extract_FUNC(vc[3,]$INFO)

vc[1,]
fifo = info_parse(vc$INFO[1])
element_sep(fifo[10])

vc |>
  dplyr::filter(grepl('coding', INFO))
  dplyr::pull(INFO)


vc |>
  dplyr::filter(grepl('^SVTYPE', INFO)) |>
  dplyr::filter(FILTER != "FAIL")
df = data.frame(entries = unlist(stringr::str_split(infos, pattern = ",")))
df$entries = gsub("'",'', df$entries)
df = tidyr::separate(df, col = 'entries', into = c("parameter","value"), sep = ":", fill = 'right')
df = tidyr::pivot_wider(df, names_from = parameter, values_from = value)
df = dplyr::relocate(df, contains("gene"), transcript,exon, coding, protein, contains("orig"), contains("normalized"), polyphen,sift )




