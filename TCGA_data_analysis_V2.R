#TCGAbiolinks(TCGA RNAseq data download)-----------------------------------

#BiocManager::install("TCGAbiolinks")
#browseVignettes("TCGAbiolinks")
#BiocManager::install("SummarizedExperiment")
#install.packages("DT")

setwd('C:/Users/dannyj/Documents/TCGA_data_analysis/LGG_RNAseq')

library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(DT)

query <- GDCquery(project = "TCGA-LGG",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts",
                  experimental.strategy = "RNA-Seq",
                  legacy = FALSE)
GDCdownload(query, method = "api", files.per.chunk = 10)
RNAseq_data <- GDCprepare(query, save = TRUE, save.filename = "LGG_RNAseq_counts.rda")


# datatable(as.data.frame(colData(data)),
#           options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
#           rownames = FALSE)
# datatable(assay(data)[1:100,],
#           options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
#           rownames = TRUE)
#rowRanges(data)


#deseq2 RowData prepare----------------------------------------------------

#load("LGG_RNAseq_counts.rda")
RNAseq_data_matrix= as.data.frame(assay(RNAseq_data))
write.csv(RNAseq_data_matrix, 
          file="LGG_RNAseq_data_matrix.csv",
          row.names = TRUE,
          col.names = TRUE,
          quote=FALSE)

colnames(RNAseq_data_matrix)=
  gsub("-\\S\\S\\S-\\S\\S\\S-\\S\\S\\S\\S-\\S\\S$", 
       "",
       colnames(RNAseq_data_matrix))

write.csv(RNAseq_data_matrix, 
          file="TCGA_LGG_RowData.csv", 
          row.names = TRUE,
          quote=FALSE)



#TCGAbiolinks(TCGA mutation data download)------
query_mut <- GDCquery(project = "TCGA-LGG",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts",
                      experimental.strategy = "RNA-Seq",
                      legacy = FALSE)
GDCdownload(query_mut, method = "api", files.per.chunk = 10)

LGG.muse.maf <- GDCquery_Maf("LGG", 
                             pipelines = "muse",
                             save.csv=FALSE)
write.csv(LGG.muse.maf,
          file = "LGG.muse.maf.csv",
          quote = FALSE,
          row.names = FALSE)
# LGG.varscan2.maf <- GDCquery_Maf("LGG", 
#                                  pipelines = "varscan2",
#                                  save.csv=TRUE)
# LGG.somaticsniper.maf <- GDCquery_Maf("LGG", 
#                                       pipelines = "somaticsniper",
#                                       save.csv=TRUE)
# LGG.mutect.maf <- GDCquery_Maf("LGG", 
#                                pipelines = "mutect2",
#                                save.csv=TRUE)

#setwd('C:/Users/dannyj/Documents/TCGA_data_analysis/LGG_RNAseq/GDCdata')
#LGG.muse.maf=read.csv("TCGA.LGG.muse.59a84472-27d4-497c-8f37-8bc447ff9374.DR-10.0.somatic.maf.csv",  header=FALSE, stringsAsFactors = FALSE)

head(LGG.muse.maf)
str(LGG.muse.maf)
summary(LGG.muse.maf)

LGG.muse.maf_2=LGG.muse.maf%>%
  select(Hugo_Symbol,
         Entrez_Gene_Id,
         Variant_Classification,
         dbSNP_RS,
         dbSNP_Val_Status,
         Tumor_Sample_Barcode,
         Matched_Norm_Sample_Barcode,
         Tumor_Sample_UUID,
         Matched_Norm_Sample_UUID,
         HGVSc,
         HGVSp,
         HGVSp_Short,
         Transcript_ID,
         Exon_Number,
         Gene,
         One_Consequence,
         Consequence,
         Protein_position,
         Amino_acids,
         Codons,
         SWISSPROT,
         SIFT,
         PolyPhen,
         DOMAINS,
         IMPACT,
         COSMIC)
LGG.muse.maf_TP53=LGG.muse.maf_2%>%
  filter(Hugo_Symbol=="TP53")

LGG.muse.maf_TP53_clean=LGG.muse.maf_TP53%>%
  mutate(Tumor_Case_Barcode=gsub("-\\S\\S\\S-\\S\\S\\S-\\S\\S\\S\\S-\\S\\S$", 
                                 "",
                                 Tumor_Sample_Barcode))%>%
  select(Tumor_Case_Barcode,
         Protein_position,
         HGVSp_Short,
         Exon_Number,
         Variant_Classification,
         Consequence,
         IMPACT,
         SIFT,
         PolyPhen,
  )

write.csv(LGG.muse.maf_TP53_clean, 
          file="LGG.muse.maf_TP53_clean.csv", 
          row.names = FALSE,
          col.names = TRUE, 
          quote=FALSE)

#create list of all case IDs with mutation data-------------------------------
LGG.muse.maf_case_list=LGG.muse.maf_2%>%
  mutate(Tumor_Case_Barcode=gsub("-\\S\\S\\S-\\S\\S\\S-\\S\\S\\S\\S-\\S\\S$", 
                                 "",
                                 Tumor_Sample_Barcode))%>%
  select(Tumor_Case_Barcode)
LGG.muse.maf_case_list=unique(LGG.muse.maf_case_list)

write.csv(LGG.muse.maf_case_list, 
          file="LGG.muse.maf_case_list.csv", 
          row.names = FALSE,
          col.names = TRUE, 
          quote=FALSE)

#create the Coldata manually

#create the RowData that matches the ColData--------------------------------------------------------

RowData=read.csv("TCGA_LGG_RowData.csv", 
                 header=TRUE,
                 check.names=FALSE,
                 stringsAsFactors = FALSE)
names(RowData)[1] <- "ENSEMBL"
RowData_match= RowData%>%
  select(c(1),"TCGA-TQ-A7RF",
         "TCGA-DU-6407",
         "TCGA-FG-8191",
         "TCGA-DU-A7T8",
         "TCGA-FG-8182",
         "TCGA-S9-A7R4",
         "TCGA-DU-7007",
         "TCGA-DU-6399",
         "TCGA-E1-5304",
         "TCGA-QH-A6CW",
         "TCGA-DH-A7UT",
         "TCGA-DU-A5TU",
         "TCGA-DU-7301",
         "TCGA-E1-A7YD",
         "TCGA-RY-A845",
         "TCGA-P5-A5F1",
         "TCGA-HT-8563",
         "TCGA-FG-A4MX",
         "TCGA-QH-A6X9",
         "TCGA-TM-A84J",
         "TCGA-FG-5963",
         "TCGA-DU-8165",
         "TCGA-DU-A6S7",
         "TCGA-HT-8111",
         "TCGA-E1-5322",
         "TCGA-P5-A780",
         "TCGA-VV-A86M",
         "TCGA-S9-A7J0",
         "TCGA-TQ-A7RH",
         "TCGA-HT-7479",
         "TCGA-DU-8167",
         "TCGA-DU-7019",
         "TCGA-HT-7858",
         "TCGA-WY-A85C",
         "TCGA-RY-A83Z",
         "TCGA-DH-A66D",
         "TCGA-TQ-A7RM",
         "TCGA-HT-7902",
         "TCGA-FG-A4MY",
         "TCGA-P5-A733",
         "TCGA-CS-6290",
         "TCGA-HT-7477",
         "TCGA-DU-5872",
         "TCGA-HT-7676",
         "TCGA-HT-7688",
         "TCGA-DU-7010",
         "TCGA-DH-A66G",
         "TCGA-P5-A5EU",
         "TCGA-P5-A72X",
         "TCGA-HT-7879",
         "TCGA-DB-A4XF",
         "TCGA-DH-5142",
         "TCGA-DB-A64S",
         "TCGA-DU-7015",
         "TCGA-P5-A736",
         "TCGA-TQ-A8XE",
         "TCGA-HT-A5RB",
         "TCGA-HT-A74J",
         "TCGA-S9-A7QW",
         "TCGA-TM-A84L",
         "TCGA-TM-A7CF",
         "TCGA-HT-7485",
         "TCGA-DU-7013",
         "TCGA-WY-A85E",
         "TCGA-DU-6542",
         "TCGA-VM-A8C8",
         "TCGA-S9-A6TU",
         "TCGA-CS-4943",
         "TCGA-HT-A61B",
         "TCGA-DB-A4XD",
         "TCGA-S9-A7QX",
         "TCGA-HW-A5KL",
         "TCGA-HT-7690",
         "TCGA-DB-5277",
         "TCGA-S9-A6U8",
         "TCGA-TQ-A7RR",
         "TCGA-P5-A72W",
         "TCGA-DB-A4XE",
         "TCGA-HT-7483",
         "TCGA-HW-8321",
         "TCGA-HT-7603",
         "TCGA-S9-A6U9",
         "TCGA-E1-A7YW",
         "TCGA-TQ-A7RK",
         "TCGA-TQ-A7RV",
         "TCGA-DB-5281",
         "TCGA-S9-A7R7",
         "TCGA-TM-A84M",
         "TCGA-CS-4938",
         "TCGA-DU-7309",
         "TCGA-WY-A859",
         "TCGA-DU-7306",
         "TCGA-DB-A4XB",
         "TCGA-TM-A84Q",
         "TCGA-CS-5393",
         "TCGA-DB-5275",
         "TCGA-HT-7472",
         "TCGA-HT-7474",
         "TCGA-TM-A84F",
         "TCGA-E1-5305",
         "TCGA-DU-7008",
         "TCGA-S9-A7R8",
         "TCGA-IK-7675",
         "TCGA-HT-A5R5",
         "TCGA-S9-A7R3",
         "TCGA-QH-A6XA",
         "TCGA-DU-A5TR",
         "TCGA-HT-7478",
         "TCGA-S9-A6U6",
         "TCGA-TM-A84T",
         "TCGA-DU-5871",
         "TCGA-E1-A7YK",
         "TCGA-HT-7693",
         "TCGA-S9-A89Z",
         "TCGA-P5-A5F4",
         "TCGA-WY-A85B",
         "TCGA-TQ-A7RW",
         "TCGA-DU-6396",
         "TCGA-FG-5965",
         "TCGA-DH-A7UU",
         "TCGA-HT-7686",
         "TCGA-P5-A5EW",
         "TCGA-HT-7689",
         "TCGA-DH-A7UV",
         "TCGA-HT-A4DS",
         "TCGA-HT-7473",
         "TCGA-S9-A6TS",
         "TCGA-FG-A711",
         "TCGA-DU-5853",
         "TCGA-DU-5855",
         "TCGA-VM-A8CF",
         "TCGA-QH-A870",
         "TCGA-HW-7490",
         "TCGA-DU-7304",
         "TCGA-CS-5396",
         "TCGA-E1-5307",
         "TCGA-DU-A7TC",
         "TCGA-HT-A614",
         "TCGA-P5-A5F2",
         "TCGA-DU-8163",
         "TCGA-HT-7606",
         "TCGA-HW-A5KM",
         "TCGA-DB-5273",
         "TCGA-HT-8105",
         "TCGA-HT-7470",
         "TCGA-HT-7476",
         "TCGA-HT-7884",
         "TCGA-HT-7601",
         "TCGA-DU-8166",
         "TCGA-FG-7636",
         "TCGA-HT-8018",
         "TCGA-TM-A7C4",
         "TCGA-DB-A75M",
         "TCGA-CS-4942",
         "TCGA-E1-A7YH",
         "TCGA-TQ-A7RN",
         "TCGA-DU-7299",
         "TCGA-FN-7833",
         "TCGA-E1-5303",
         "TCGA-DU-A5TP",
         "TCGA-HT-A618",
         "TCGA-DU-7011",
         "TCGA-E1-A7YU",
         "TCGA-DU-7300",
         "TCGA-HT-7611",
         "TCGA-TQ-A7RJ",
         "TCGA-E1-A7YI",
         "TCGA-DU-A5TW",
         "TCGA-DU-6401",
         "TCGA-TM-A7CA",
         "TCGA-S9-A6WO",
         "TCGA-DB-A4XC",
         "TCGA-CS-4941",
         "TCGA-CS-4944",
         "TCGA-CS-5390",
         "TCGA-CS-5394",
         "TCGA-CS-5395",
         "TCGA-CS-5397",
         "TCGA-CS-6186",
         "TCGA-CS-6188",
         "TCGA-CS-6668",
         "TCGA-CS-6670",
         "TCGA-DB-5274",
         "TCGA-DB-5278",
         "TCGA-DB-5279",
         "TCGA-DB-A4X9",
         "TCGA-DB-A4XA",
         "TCGA-DB-A4XG",
         "TCGA-DB-A4XH",
         "TCGA-DB-A64L",
         "TCGA-DB-A64O",
         "TCGA-DB-A64P",
         "TCGA-DB-A64Q",
         "TCGA-DB-A64R",
         "TCGA-DB-A64U",
         "TCGA-DB-A64V",
         "TCGA-DB-A64W",
         "TCGA-DB-A75K",
         "TCGA-DB-A75L",
         "TCGA-DB-A75O",
         "TCGA-DB-A75P",
         "TCGA-DH-5141",
         "TCGA-DH-5144",
         "TCGA-DH-A669",
         "TCGA-DH-A66B",
         "TCGA-DH-A66F",
         "TCGA-DH-A7UR",
         "TCGA-DH-A7US",
         "TCGA-DU-5847",
         "TCGA-DU-5849",
         "TCGA-DU-5852",
         "TCGA-DU-5854",
         "TCGA-DU-5870",
         "TCGA-DU-5874",
         "TCGA-DU-6393",
         "TCGA-DU-6394",
         "TCGA-DU-6397",
         "TCGA-DU-6400",
         "TCGA-DU-6402",
         "TCGA-DU-6403",
         "TCGA-DU-6404",
         "TCGA-DU-6405",
         "TCGA-DU-6406",
         "TCGA-DU-6410",
         "TCGA-DU-7006",
         "TCGA-DU-7009",
         "TCGA-DU-7012",
         "TCGA-DU-7018",
         "TCGA-DU-7290",
         "TCGA-DU-7292",
         "TCGA-DU-7294",
         "TCGA-DU-7302",
         "TCGA-DU-8158",
         "TCGA-DU-8161",
         "TCGA-DU-8162",
         "TCGA-DU-8164",
         "TCGA-DU-8168",
         "TCGA-DU-A5TS",
         "TCGA-DU-A5TT",
         "TCGA-DU-A5TY",
         "TCGA-DU-A6S2",
         "TCGA-DU-A6S3",
         "TCGA-DU-A6S6",
         "TCGA-DU-A6S8",
         "TCGA-DU-A76K",
         "TCGA-DU-A76L",
         "TCGA-DU-A76R",
         "TCGA-DU-A7T6",
         "TCGA-DU-A7TA",
         "TCGA-DU-A7TB",
         "TCGA-DU-A7TD",
         "TCGA-DU-A7TG",
         "TCGA-DU-A7TJ",
         "TCGA-E1-5302",
         "TCGA-E1-5311",
         "TCGA-E1-5318",
         "TCGA-E1-5319",
         "TCGA-E1-A7YJ",
         "TCGA-E1-A7YM",
         "TCGA-E1-A7YN",
         "TCGA-E1-A7YO",
         "TCGA-E1-A7YQ",
         "TCGA-E1-A7YS",
         "TCGA-E1-A7YV",
         "TCGA-E1-A7YY",
         "TCGA-E1-A7Z2",
         "TCGA-E1-A7Z3",
         "TCGA-E1-A7Z4",
         "TCGA-E1-A7Z6",
         "TCGA-EZ-7264",
         "TCGA-F6-A8O3",
         "TCGA-FG-5962",
         "TCGA-FG-5964",
         "TCGA-FG-6688",
         "TCGA-FG-6689",
         "TCGA-FG-6691",
         "TCGA-FG-6692",
         "TCGA-FG-7634",
         "TCGA-FG-7637",
         "TCGA-FG-7638",
         "TCGA-FG-7641",
         "TCGA-FG-7643",
         "TCGA-FG-8186",
         "TCGA-FG-8187",
         "TCGA-FG-8189",
         "TCGA-FG-A4MU",
         "TCGA-FG-A4MW",
         "TCGA-FG-A60K",
         "TCGA-FG-A60L",
         "TCGA-FG-A6IZ",
         "TCGA-FG-A6J1",
         "TCGA-FG-A6J3",
         "TCGA-FG-A70Y",
         "TCGA-FG-A70Z",
         "TCGA-FG-A710",
         "TCGA-FG-A713",
         "TCGA-FG-A87N",
         "TCGA-FG-A87Q",
         "TCGA-HT-7467",
         "TCGA-HT-7468",
         "TCGA-HT-7471",
         "TCGA-HT-7480",
         "TCGA-HT-7481",
         "TCGA-HT-7602",
         "TCGA-HT-7605",
         "TCGA-HT-7607",
         "TCGA-HT-7608",
         "TCGA-HT-7616",
         "TCGA-HT-7620",
         "TCGA-HT-7677",
         "TCGA-HT-7680",
         "TCGA-HT-7681",
         "TCGA-HT-7684",
         "TCGA-HT-7687",
         "TCGA-HT-7691",
         "TCGA-HT-7692",
         "TCGA-HT-7694",
         "TCGA-HT-7695",
         "TCGA-HT-7854",
         "TCGA-HT-7855",
         "TCGA-HT-7856",
         "TCGA-HT-7857",
         "TCGA-HT-7860",
         "TCGA-HT-7873",
         "TCGA-HT-7874",
         "TCGA-HT-7875",
         "TCGA-HT-7877",
         "TCGA-HT-7880",
         "TCGA-HT-7881",
         "TCGA-HT-7882",
         "TCGA-HT-8010",
         "TCGA-HT-8011",
         "TCGA-HT-8012",
         "TCGA-HT-8104",
         "TCGA-HT-8107",
         "TCGA-HT-8108",
         "TCGA-HT-8109",
         "TCGA-HT-8110",
         "TCGA-HT-8113",
         "TCGA-HT-A4DV",
         "TCGA-HT-A5R7",
         "TCGA-HT-A5R9",
         "TCGA-HT-A5RA",
         "TCGA-HT-A5RC",
         "TCGA-HT-A615",
         "TCGA-HT-A616",
         "TCGA-HT-A617",
         "TCGA-HT-A619",
         "TCGA-HT-A61A",
         "TCGA-HT-A61C",
         "TCGA-HT-A74H",
         "TCGA-HT-A74K",
         "TCGA-HT-A74L",
         "TCGA-HT-A74O",
         "TCGA-HW-7486",
         "TCGA-HW-7487",
         "TCGA-HW-7489",
         "TCGA-HW-7491",
         "TCGA-HW-7495",
         "TCGA-HW-8320",
         "TCGA-HW-8322",
         "TCGA-HW-A5KJ",
         "TCGA-HW-A5KK",
         "TCGA-IK-8125",
         "TCGA-KT-A74X",
         "TCGA-KT-A7W1",
         "TCGA-P5-A5ET",
         "TCGA-P5-A5EV",
         "TCGA-P5-A5EX",
         "TCGA-P5-A5EY",
         "TCGA-P5-A5EZ",
         "TCGA-P5-A5F0",
         "TCGA-P5-A5F6",
         "TCGA-P5-A72U",
         "TCGA-P5-A72Z",
         "TCGA-P5-A730",
         "TCGA-P5-A731",
         "TCGA-P5-A737",
         "TCGA-P5-A77W",
         "TCGA-P5-A77X",
         "TCGA-P5-A781",
         "TCGA-QH-A65R",
         "TCGA-QH-A65S",
         "TCGA-QH-A65V",
         "TCGA-QH-A65X",
         "TCGA-QH-A65Z",
         "TCGA-QH-A6CS",
         "TCGA-QH-A6CU",
         "TCGA-QH-A6CV",
         "TCGA-QH-A6CX",
         "TCGA-QH-A6CY",
         "TCGA-QH-A6CZ",
         "TCGA-QH-A6X3",
         "TCGA-QH-A6X4",
         "TCGA-QH-A6X5",
         "TCGA-QH-A6X8",
         "TCGA-QH-A6XC",
         "TCGA-QH-A86X",
         "TCGA-R8-A6MK",
         "TCGA-R8-A6ML",
         "TCGA-R8-A6MO",
         "TCGA-R8-A6YH",
         "TCGA-R8-A73M",
         "TCGA-RY-A83X",
         "TCGA-RY-A83Y",
         "TCGA-RY-A840",
         "TCGA-RY-A843",
         "TCGA-RY-A847",
         "TCGA-S9-A6TV",
         "TCGA-S9-A6TW",
         "TCGA-S9-A6TX",
         "TCGA-S9-A6TY",
         "TCGA-S9-A6TZ",
         "TCGA-S9-A6U0",
         "TCGA-S9-A6U1",
         "TCGA-S9-A6U2",
         "TCGA-S9-A6U5",
         "TCGA-S9-A6UA",
         "TCGA-S9-A6UB",
         "TCGA-S9-A6WE",
         "TCGA-S9-A6WH",
         "TCGA-S9-A6WI",
         "TCGA-S9-A6WL",
         "TCGA-S9-A6WM",
         "TCGA-S9-A6WN",
         "TCGA-S9-A6WP",
         "TCGA-S9-A7IQ",
         "TCGA-S9-A7IS",
         "TCGA-S9-A7IX",
         "TCGA-S9-A7IY",
         "TCGA-S9-A7IZ",
         "TCGA-S9-A7J1",
         "TCGA-S9-A7J2",
         "TCGA-S9-A7J3",
         "TCGA-S9-A7QY",
         "TCGA-S9-A7QZ",
         "TCGA-S9-A7R1",
         "TCGA-S9-A7R2",
         "TCGA-S9-A89V",
         "TCGA-TM-A7C3",
         "TCGA-TM-A7C5",
         "TCGA-TM-A84B",
         "TCGA-TM-A84C",
         "TCGA-TM-A84G",
         "TCGA-TM-A84H",
         "TCGA-TM-A84O",
         "TCGA-TM-A84R",
         "TCGA-TM-A84S",
         "TCGA-TQ-A7RG",
         "TCGA-TQ-A7RI",
         "TCGA-TQ-A7RO",
         "TCGA-TQ-A7RP",
         "TCGA-TQ-A7RQ",
         "TCGA-TQ-A7RS",
         "TCGA-VM-A8C9",
         "TCGA-VM-A8CA",
         "TCGA-VM-A8CB",
         "TCGA-VM-A8CD",
         "TCGA-VM-A8CE",
         "TCGA-VM-A8CH",
         "TCGA-VV-A829",
         "TCGA-VW-A7QS",
         "TCGA-W9-A837",
         "TCGA-WY-A85A")
write.csv(RowData_match, 
          file="TCGA_LGG_RowData_match.csv", 
          row.names = FALSE,
          quote=FALSE)


#launch DESeq2 analysis----------------------------------------------
library("DESeq2")
#cts <- RowData_match
cts= read.csv('TCGA_LGG_RowData_match.csv',
              header=TRUE,
              stringsAsFactors = FALSE)
cts2 <- cts[,-1]
rownames(cts2) <- cts[,1]
coldata= read.table('TCGA_LGG_ColData.txt', 
                    header=TRUE, 
                    sep="\t", 
                    stringsAsFactors = FALSE)


dds <- DESeqDataSetFromMatrix(countData = cts2,
                              colData = coldata,
                              design= ~ IMPACT)
dds$IMPACT <- relevel(dds$IMPACT, ref = "WT")
dds <- DESeq(dds)
dds
resultsNames(dds) # lists the coefficients
res1 <- results(dds, name="IMPACT_HIGH_vs_WT")

write.csv(as.data.frame(res1), 
          file="IMPACT_HIGH_vs_WT.csv")

res2 <- results(dds, name="IMPACT_LOW_vs_WT")

write.csv(as.data.frame(res2), 
          file="IMPACT_LOW_vs_WT.csv")

res3 <- results(dds, name="IMPACT_MODERATE_vs_WT")

write.csv(as.data.frame(res3),
          file="IMPACT_MODERATE_vs_WT.csv")

#launch GSEA analysis----------------------------------------------------

library(dplyr)
library(AnnotationHub)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(msigdbr)
setwd('C:/Users/dannyj/Documents/TCGA_data_analysis/LGG_RNAseq')
D <- read.csv("IMPACT_MODERATE_vs_WT.csv")
head(D)
str(D)
summary(D)

names(D)[1] <- "ENSEMBL_ID"
D2=D%>%
  dplyr::filter(ENSEMBL_ID !="__no_feature"&
                  ENSEMBL_ID !="__ambiguous"& 
                  ENSEMBL_ID !="__too_low_aQual"&
                  ENSEMBL_ID !="__not_aligned"&
                  ENSEMBL_ID !="__alignment_not_unique")%>%
  dplyr::filter(is.na(log2FoldChange)==FALSE)%>%
  dplyr::filter(is.na(padj)==FALSE)


D2_ENTREZID<- bitr(D2[,1], fromType = "ENSEMBL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Hs.eg.db,
                   drop = TRUE)

D3=right_join(D2,D2_ENTREZID, by=c("ENSEMBL_ID"="ENSEMBL"))
head(D3)
str(D3)
summary(D3)

## feature 1: numeric vector
geneList_ENTREZ <- D3[,3]
## feature 2: named vector
names(geneList_ENTREZ) <- D3[,8]
## feature 3: decreasing order
geneList_ENTREZ <- sort(geneList_ENTREZ, decreasing = TRUE)
head(geneList_ENTREZ)
str(geneList_ENTREZ)

keytypes(org.Hs.eg.db)

m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

C6=GSEA(
  geneList_ENTREZ,
  exponent = 1,
  minGSSize = 15,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.25,
  pAdjustMethod = "BH",
  TERM2GENE= m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea"
)
write.table(as.data.frame(C6@result), 
            file="LGG_IMPACT_MODERATE_vs_WT_ENTERZ_padjNAcut_C6_25.txt",
            sep="\t", 
            row.names = FALSE,
            col.names = TRUE,
            quote=FALSE)


# launch DEG analysis-----------------------------------------------------
library(AnnotationHub)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)

d <- read.csv("IMPACT_HIGH_vs_WT.csv")
head(d)
str(d)
summary(d)

names(d)[1] <- "ENSEMBL_ID"

#positive
d2=d%>%
  dplyr::filter(ENSEMBL_ID !="__no_feature"&
                  ENSEMBL_ID !="__ambiguous"& 
                  ENSEMBL_ID !="__too_low_aQual"&
                  ENSEMBL_ID !="__not_aligned"&
                  ENSEMBL_ID !="__alignment_not_unique")%>%
  dplyr::filter(is.na(log2FoldChange)==FALSE)%>%
  dplyr::filter(is.na(padj)==FALSE)%>%
  dplyr::filter(log2FoldChange>=0)%>%
  dplyr::filter(padj<0.05)
d2_ENTREZID<- bitr(d2[,1], fromType = "ENSEMBL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Hs.eg.db,
                   drop = TRUE)
d2_SYMBOL<- bitr(d2[,1], fromType = "ENSEMBL",
                 toType = c("SYMBOL"),
                 OrgDb = org.Hs.eg.db,
                 drop = TRUE)
d3=left_join(d2,d2_ENTREZID, by=c("ENSEMBL_ID"="ENSEMBL"))
d4=left_join(d3,d2_SYMBOL, by=c("ENSEMBL_ID"="ENSEMBL"))
d4<- d4[ , c(1,8, 9, 2:7)]
d5=d4%>%
  arrange(desc(log2FoldChange))
write.csv(d5,file="LGG_IMPACT_HIGH_vs_WT_DEG_pos_unfilter.csv",
          row.names = FALSE,
          quote=FALSE)

#negative
d12=d%>%
  dplyr::filter(ENSEMBL_ID !="__no_feature"&
                  ENSEMBL_ID !="__ambiguous"& 
                  ENSEMBL_ID !="__too_low_aQual"&
                  ENSEMBL_ID !="__not_aligned"&
                  ENSEMBL_ID !="__alignment_not_unique")%>%
  dplyr::filter(is.na(log2FoldChange)==FALSE)%>%
  dplyr::filter(is.na(padj)==FALSE)%>%
  dplyr::filter(log2FoldChange<=0)%>%
  dplyr::filter(padj<0.05)
d12_ENTREZID<- bitr(d12[,1], fromType = "ENSEMBL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db,
                    drop = TRUE)
d12_SYMBOL<- bitr(d12[,1], fromType = "ENSEMBL",
                  toType = c("SYMBOL"),
                  OrgDb = org.Hs.eg.db,
                  drop = TRUE)
d13=left_join(d12,d12_ENTREZID, by=c("ENSEMBL_ID"="ENSEMBL"))
d14=left_join(d13,d12_SYMBOL, by=c("ENSEMBL_ID"="ENSEMBL"))
d14<- d14[ , c(1,8, 9, 2:7)]
d15=d14%>%
  arrange(log2FoldChange)
write.csv(d15,file="LGG_IMPACT_HIGH_vs_WT_DEG_neg_unfilter.csv",
          row.names = FALSE,
          quote=FALSE)
