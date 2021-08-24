#TCGAbiolinks(TCGA RNAseq data download)-----------------------------------

#BiocManager::install("TCGAbiolinks")
#browseVignettes("TCGAbiolinks")
#BiocManager::install("SummarizedExperiment")
#install.packages("DT")

library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(DT)

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts",
                  experimental.strategy = "RNA-Seq",
                  legacy = FALSE)
GDCdownload(query,
            method = "api",
            files.per.chunk = 10)

RNAseq_data <- GDCprepare(query,
                          save = TRUE,
                          save.filename = "BRCA_RNAseq_counts.rda")


# datatable(as.data.frame(colData(data)),
#           options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
#           rownames = FALSE)
# datatable(assay(data)[1:100,],
#           options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
#           rownames = TRUE)
#rowRanges(data)


#deseq2 RowData prepare----------------------------------------------------

load("BRCA_RNAseq_counts.rda")
RNAseq_data <- data
RNAseq_data_matrix= as.data.frame(assay(RNAseq_data))
save(RNAseq_data_matrix,
     file = "BRCA_RNAseq_data_matrix.rdata")
write.csv(RNAseq_data_matrix, 
          file="BRCA_RNAseq_data_matrix.csv",
          row.names = TRUE,
          col.names = TRUE,
          quote=FALSE)
load("BRCA_RNAseq_data_matrix.rdata")

RNAseq_data_matrix_sample <- RNAseq_data_matrix
colnames(RNAseq_data_matrix_sample)=
  gsub("\\S-\\S\\S\\S-\\S\\S\\S\\S-\\S\\S$", 
       "",
       colnames(RNAseq_data_matrix_sample))


write.csv(RNAseq_data_matrix_sample, 
          file="RNAseq_data_matrix_sample.csv",
          quote=FALSE)



#TCGAbiolinks(TCGA mutation data download)------

BRCA.muse.maf <- GDCquery_Maf("BRCA", 
                             pipelines = "muse",
                             save.csv=TRUE)
save(BRCA.muse.maf,
     file = "BRCA.muse.maf.rdata")


load("BRCA.muse.maf.rdata")
#BRCA.muse.maf=read.csv("TCGA.BRCA.muse.59a84472-27d4-497c-8f37-8bc447ff9374.DR-10.0.somatic.maf.csv",  header=FALSE, stringsAsFactors = FALSE)


BRCA.muse.maf_TP53_clean=BRCA.muse.maf %>% 
  filter(Hugo_Symbol=="TP53") %>% 
  mutate(Tumor_Case_Barcode=gsub("\\S-\\S\\S\\S-\\S\\S\\S\\S-\\S\\S$", 
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

write.csv(BRCA.muse.maf_TP53_clean, 
          file="BRCA.muse.maf_TP53_clean.csv", 
          quote=FALSE)

#create list of all case IDs with mutation data-------------------------------
BRCA.muse.maf_case_list=BRCA.muse.maf%>%
  mutate(Tumor_Case_Barcode=gsub("\\S-\\S\\S\\S-\\S\\S\\S\\S-\\S\\S$", 
                                 "",
                                 Tumor_Sample_Barcode))%>%
  select(Tumor_Case_Barcode)
write.csv(BRCA.muse.maf_case_list,
          file = "BRCA.muse.maf_case_list.csv")

BRCA.muse.maf_case_list_unique <- BRCA.muse.maf_case_list
BRCA.muse.maf_case_list_unique=unique(BRCA.muse.maf_case_list_unique)
#write.csv(BRCA.muse.maf_case_list_unique, file = "BRCA.muse.maf_case_list_unique.csv")

write.csv(BRCA.muse.maf_case_list, 
          file="BRCA.muse.maf_case_list.csv", 
          row.names = FALSE,
          col.names = TRUE, 
          quote=FALSE)

#create the Coldata manually



#create the RowData that matches the ColData--------------------------------------------------------

#delete duplicated RNAseq data
RNAseq_data_matrix_sample_delete_duplicated <- 
  RNAseq_data_matrix_sample[, !duplicated(colnames(RNAseq_data_matrix_sample))]

mutation_list <- BRCA.muse.maf_case_list_unique[,1, drop = TRUE]

#select samples that have mutation data
RowData <- 
  RNAseq_data_matrix_sample_delete_duplicated[, colnames(RNAseq_data_matrix_sample_delete_duplicated)%in%
                                         (mutation_list)]
write.csv(RowData, file = "RowData.csv")

#launch DESeq2 analysis----------------------------------------------
library("DESeq2")
#cts <- RowData_match
cts= read.csv('TCGA_BRCA_RowData_match.csv',
              header=TRUE,
              stringsAsFactors = FALSE)
cts2 <- cts[,-1]
rownames(cts2) <- cts[,1]
coldata= read.table('TCGA_BRCA_ColData.txt', 
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
setwd('C:/Users/dannyj/Documents/TCGA_data_analysis/BRCA_RNAseq')
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
            file="BRCA_IMPACT_MODERATE_vs_WT_ENTERZ_padjNAcut_C6_25.txt",
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
write.csv(d5,file="BRCA_IMPACT_HIGH_vs_WT_DEG_pos_unfilter.csv",
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
write.csv(d15,file="BRCA_IMPACT_HIGH_vs_WT_DEG_neg_unfilter.csv",
          row.names = FALSE,
          quote=FALSE)
