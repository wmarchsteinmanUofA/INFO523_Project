#Imports
library("Rsamtools")
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library("GenomicAlignments")
library("BiocParallel")
library("Rsubread")
library("pheatmap")
library("annotate")
library(dplyr) 
library(clusterProfiler)
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(apeglm)
library(CBNplot)
library(enrichplot)
library(msigdbr)
library(fgsea)
library(GSEAplot)


#Data Loading and preprocessing -- data not included in github due to size.
#filenames <- list.files(path="./SampleData/All_data", pattern=".bam", all.files=FALSE,
#                        full.names=TRUE)
#filenames_bam <- filenames[!grepl(".bai", filenames)]
#filenames_bai <- filenames[grepl(".bai", filenames)]
#filenames_bam_50 <- filenames_bam[grepl("50", filenames_bam)]
#
#fc.fifty <- featureCounts(filenames_bam_50, annot.inbuilt = "hg38", isPairedEnd = TRUE)
#mam.samples_50 <- read.table("./SampleData/All_data/SampleFile50.txt",header=T)
#mam.samples_50$Treatment <-factor(mam.samples_50$Treatment)
#mam.samples_50$Group <-factor(mam.samples_50$Group)
#
#colnames(fc.fifty$counts) <- mam.samples_50$Sample 
#colnames(fc.ut$counts) <- mam.samples_UT$Sample 
#
#
#retain <- rowSums(fc.fifty$counts) >= 10
#counts_filtered.fifty <- fc.fifty$counts[retain,]
#saveRDS(counts_filtered.fifty, file = "C:/Users/wmarchsteinman/Desktop/INFO523_Project/fiftyCounts.RDS")
#saveRDS(mam.samples_50, file = "C:/Users/wmarchsteinman/Desktop/INFO523_Project/mamsamples_50.RDS")

#load preprocessed data
counts_filtered.fifty <- readRDS(file = "./fiftyCounts.RDS")
mam.samples_50 <- readRDS(file = "./mamsamples_50.RDS")

diff_fc.fifty <- DESeqDataSetFromMatrix(countData = counts_filtered.fifty, colData = DataFrame(mam.samples_50), design = ~Treatment)

keep <- rowSums(counts(diff_fc.fifty)) >= 10
diff_fc.fifty <- diff_fc.fifty[keep,]

ap.dso_fifty <- diff_fc.fifty
ap.samples_fifty <- mam.samples_50
ap.vst_fifty <- assay(vst(ap.dso_fifty, blind = FALSE))

#Differential Expression Results for SFRX vs J14 50
ap.dso_fifty$Group <- relevel(ap.dso_fifty$Group, ref = "J14")
ap.deseq2_fifty <-DESeq(ap.dso_fifty)
resultsNames(ap.deseq2_fifty)
J14.results <- results(ap.deseq2_fifty, name = "Treatment_SFRX_vs_J14", alpha = 0.05)
J14.results <- J14.results %>% as.data.frame %>% arrange(pvalue) %>% filter(!is.na(padj))
#differentially-expressed genes with largest fold changes
J14.test <- J14.results %>% as.data.frame %>% filter(padj < 0.05)

#enrichmentAnalysis
KeggPathway <- clusterProfiler::enrichKEGG(rownames(J14.test))
GoPathway <- clusterProfiler::enrichGO(rownames(J14.test), ont = "BP", OrgDb = org.Hs.eg.db)
go_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
go_gene_sets
g_set <- go_gene_sets %>% dplyr::select(gs_name, entrez_gene)
g_set

gsea_input <- J14.test$stat
names(gsea_input) <- rownames(J14.test)
gsea_input <- gsea_input[order(gsea_input,decreasing=T)]
gseaPathways <- clusterProfiler::GSEA(gsea_input, TERM2GENE = g_set)

lfc <- J14.test[order(J14.test$log2FoldChange, decreasing = TRUE),]
lfc
geneList <- lfc$log2FoldChange
names(geneList) <- rownames(lfc)
pway <- ReactomePA::enrichPathway(gene = rownames(J14.test))
pway <- setReadable(pway,OrgDb = org.Hs.eg.db)
pway <- enrichplot::pairwise_termsim(pway)
pwayGSE <- ReactomePA::gsePathway(geneList)

sigPathway <- subset(pway@result, p.adjust<0.05)
enrichplot::gseaplot2(gseaPathways, geneSetID = 1, title = gseaPathways@result$ID[1], pvalue_table = T)

#Data Formatting
gseaPathways@organism <- "human"
gseaPathways@keytype <- "ENTREZID"
ensm <- clusterProfiler::bitr(rownames(ap.vst_fifty), fromType="ENTREZID", toType="ENSEMBL", OrgDb=org.Hs.eg.db)
ap.vst_fifty_ensm <- ap.vst_fifty[ensm$ENTREZID, ]
rownames(ap.vst_fifty_ensm) <- ensm$ENSEMBL
bngeneplot(results = gseaPathways@result, exp = ap.vst_fifty_ensm, convertSymbol = TRUE, pathNum = 5)
