#-------------------------------------------------------------------------------
# 08_GEP9.R
#-------------------------------------------------------------------------------
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(cluster)
library(ComplexHeatmap)
library(circlize)
library(viridisLite)
library(NMF)
library(ggrastr)
library(ggrepel)
"%ni%" <- Negate("%in%")
se <- readRDS("rds/04_se.rds")
theme_cb <- function(color="black"){theme_classic() + theme(axis.text = element_text(color=color), axis.ticks = element_line(color=color))}
SampleInfo <- read.csv("output/Tables/00_SampleListSummary_Curated.csv", row.names = 1)

cancertype_colors <- ArchR::ArchRPalettes$stallion2[c(1:17)]
names(cancertype_colors) <- table(SampleInfo$CancerType2) %>% sort(., decreasing = T) %>%names()
ecclass_colors <- c("ecDNA_p" = "#E7298A", "ecDNA_n" = "gray")
ampclass_colors <- c("ecDNA" = "#E7298A", "BFB"= "#7570B3", "Complex-non-cyclic" = "#D95F02", "Linear" = "#1B9E77")

#---------- Sample classification ----------#
AA_result_T <- read.csv("output/Tables/02_AA_Summary_Tumor_filtered.csv")
tmp <- SampleInfo$StudyID
names(tmp) <- SampleInfo$SampleID_T
AA_result_T$StudyID <- tmp[AA_result_T$Sample.name]
AA_result_T$CancerType <- SampleInfo[AA_result_T$StudyID,]$CancerType2

#ecDNA samples
ecDNA_p <- AA_result_T$StudyID[AA_result_T$Classification == "ecDNA"] %>% unique() #277 sample w/ ecDNA
BFBCnC_p <- AA_result_T$StudyID[AA_result_T$Classification %in% c("BFB","Complex-non-cyclic")] %>% unique()
intersect(ecDNA_p, BFBCnC_p)
BFBCnC_p <- BFBCnC_p[BFBCnC_p %ni% ecDNA_p] #139 samples

se <- se[,se$StudyID %ni% BFBCnC_p]
se$category1 <- "ecDNA_n"
se$category1[which(se$StudyID %in% ecDNA_p)] <- "ecDNA_p"

###check###
ecDNA_n <- se$StudyID[se$category1 == "ecDNA_n"]
AA_result_T[AA_result_T$StudyID %in% ecDNA_n, ]$Classification
table(se$category1)
###end###

#------------ Diff Test ------------#
#Diff test ecDNA+ vs ecDNA-
GEP8_genes <- read.table("output/Tables/07_geneList_GEP8.txt")[,1]
p_values <- sapply(GEP8_genes, function(x){
  df_tmp <- data.frame(expr = assays(se)$log2tpm[rowData(se)$gene_name == x,], ecDNA = se$category1)
  out <- wilcox.test(df_tmp$expr[df_tmp$ecDNA == "ecDNA_p"], df_tmp$expr[df_tmp$ecDNA == "ecDNA_n"])
  return(out$p.value)
}) 
fdr_values <- p.adjust(p_values, method = "fdr")
diff_values <- sapply(GEP8_genes, function(x){
  df_tmp <- data.frame(expr = assays(se)$log2tpm[rowData(se)$gene_name == x,], ecDNA = se$category1)
  out <- mean(df_tmp$expr[df_tmp$ecDNA == "ecDNA_p"]) - mean(df_tmp$expr[df_tmp$ecDNA == "ecDNA_n"])
  return(out)
}) 

df_difftest_ec <- data.frame(expr_diff = diff_values, FDR = -log10(fdr_values), name = names(diff_values))
df_difftest_ec$sig <- "N"
df_difftest_ec$sig[which(df_difftest_ec$FDR > 2 & df_difftest_ec$expr_diff > 0)] <- "Y"

p1 <- ggplot(df_difftest_ec, aes(x = expr_diff, y = FDR, label = name, color = sig)) + 
  geom_point_rast(alpha = 0.5, size = 2.5) + theme_cb() +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 2, lty = "dotted") +
  scale_color_manual(values = c("Y" = "red", "N" = "darkgray")) +
  geom_text_repel() + labs(x = "Mean Diff", y = "-log10(FDR)")
p2 <- ggplot(df_difftest_ec, aes(x = expr_diff, y = FDR, color = sig)) + 
  geom_point_rast(alpha = 0.5, size = 2.5) + theme_cb() +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 2, lty = "dotted") +
  scale_color_manual(values = c("Y" = "red", "N" = "darkgray")) + labs(x = "Mean Diff", y = "-log10(FDR)")
pdf("output/Plots/08_ScatterPlot_DiffTest_Labeled.pdf", width = 18, height = 13)
p1
dev.off()
pdf("output/Plots/08_ScatterPlot_DiffTest_Labeled2.pdf", width = 8, height = 4)
p1
dev.off()
pdf("output/Plots/08_ScatterPlot_DiffTest_NonLabeled.pdf", width = 6, height = 3)
p2
dev.off()

write.table(rownames(df_difftest_ec)[df_difftest_ec$sig == "Y"], "output/Tables/08_ecDNA_geneSig.txt", row.names = F, col.names = F, quote = F)

#------------ overlapping ecDNA amplification ------------#
ecDNA_geneSig <- read.table("output/Tables/08_ecDNA_geneSig.txt")[,1]

###amplified genes?
#gene annotation: gencode.v37.annotation.gtf
library(GenomicFeatures)
gff_anno <- rtracklayer::import.gff("ref/gencode.v37.annotation.gtf")
gff_geneanno <- gff_anno[which(gff_anno$type == "gene")]

#ecDNA granges
gr_T_filtered <- readRDS("rds/02_gr_T_filtered.rds")
gr_T_filtered <- unlist(gr_T_filtered)
gr_ec <- gr_T_filtered[gr_T_filtered$Classification == "ecDNA"]
ecAmplifiedGenes <- gff_geneanno$gene_name[queryHits(findOverlaps(gff_geneanno,gr_ec))] %>% unique

ovlpGenes_ecAmp <- intersect(ecDNA_geneSig, ecAmplifiedGenes) # 45, 20.1%

p3 <- lapply(ovlpGenes_ecAmp, function(x){
  print(x)
  mtx_tmp <- data.frame(expr = assays(se)$log2tpm[which(rowData(se)$gene_name == x),], ecDNA = se$category1)
  mtx_tmp$ecDNAgene <- "No"
  fo1 <- findOverlaps(gr_ec, gff_geneanno[gff_geneanno$gene_name == x])
  ecAmp_Case <- AA_result_T$StudyID[AA_result_T$Feature.ID %in% gr_ec$FeatureID[queryHits(fo1)]] %>% unique
  mtx_tmp[ecAmp_Case,]$ecDNAgene <- "Yes"
  mtx_tmp <- mtx_tmp[order(mtx_tmp$ecDNAgene),]
  out <- ggplot(mtx_tmp, aes(x=ecDNA, y=expr, fill=ecDNA)) + 
    geom_jitter_rast(height = 0, width = 0.2, alpha = 0.5, aes(color = ecDNAgene)) + 
    geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() + ggtitle(x) +
    scale_color_manual(values = c("Yes" = "red", "No" = "darkgray")) + 
    labs(x = "", y = "log2(TPM+1)") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    geom_signif(comparisons = list(c("ecDNA_n", "ecDNA_p")), textsize = 2) + theme(legend.key.size = unit(0.3, 'cm'))
  return(out)
})
pdf("output/Plots/08_Boxplot_expr_ecDNAgeneSig_overlapAmpgene.pdf", width = 2.5, height = 3)
p3
dev.off()

ecAmpOverlapGene_pvalue <- sapply(ovlpGenes_ecAmp, function(x){
  mtx_tmp <- data.frame(expr = assays(se)$log2tpm[which(rowData(se)$gene_name == x),], ecDNA = se$category1)
  mtx_tmp$ecDNAgene <- "No"
  fo1 <- findOverlaps(gr_ec, gff_geneanno[gff_geneanno$gene_name == x])
  ecAmp_Case <- AA_result_T$StudyID[AA_result_T$Feature.ID %in% gr_ec$FeatureID[queryHits(fo1)]] %>% unique
  mtx_tmp[ecAmp_Case,]$ecDNAgene <- "Yes"
  res1 <- wilcox.test(mtx_tmp$expr~mtx_tmp$ecDNA)
  res2 <- wilcox.test(mtx_tmp[mtx_tmp$ecDNAgene == "No",]$expr~mtx_tmp[mtx_tmp$ecDNAgene == "No",]$ecDNA)
  return(c(Include = -log10(res1$p.value), Exclude = -log10(res2$p.value)))
}) %>% as.data.frame
ecAmpOverlapGene_pvalue <- t(ecAmpOverlapGene_pvalue) %>% as.data.frame()
ecAmpOverlapGene_pvalue <- ecAmpOverlapGene_pvalue[order(ecAmpOverlapGene_pvalue$Include),]
ecAmpOverlapGene_pvalue$gene <- rownames(ecAmpOverlapGene_pvalue)
ecAmpOverlapGene_pvalue <- reshape2::melt(ecAmpOverlapGene_pvalue)
p4 <- ggplot(ecAmpOverlapGene_pvalue, aes(x = reorder(gene, value), y = value, color = variable)) + 
  geom_point(alpha = 0.5, size = 2.5) + scale_color_manual(values = c(Include = "red", Exclude = "blue")) +
  scale_y_continuous(limits = c(0,27), breaks=c(0,5,10,15,20,25)) + 
  theme_cb() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "Gene", y = "-log10(P-value)")
pdf("output/Plots/08_Dotplot_Pvalue_ecDNAgeneSig_overlapAmpgene.pdf", width = 6.5, height = 3)
p4
dev.off()

#------------ overlapping eLife paper ------------#
CorEx_genesDF <- read.csv("ref/elife_corEx.csv")[,c(1,3)]
CorEx_genesDF <- data.frame(gene = stringr::str_split(CorEx_genesDF$X.Gene.GeneId[c(1:643)], "\\|", simplify = T)[,1],
                            direct = CorEx_genesDF$direction.for_geneset_enrichment.[c(1:643)])
CorEx_genesAL <- CorEx_genesDF$gene %>% sort 
CorEx_genesUP <- CorEx_genesDF$gene[CorEx_genesDF$direct == "UP"] %>% sort 
CorEx_genesDN <- CorEx_genesDF$gene[CorEx_genesDF$direct == "DOWN"] %>% sort

ovlpGenes_CorEx_genesAL <- intersect(ecDNA_geneSig, CorEx_genesAL) # 47
ovlpGenes_CorEx_genesUP <- intersect(ecDNA_geneSig, CorEx_genesUP) # 43
ovlpGenes_CorEx_genesDN <- intersect(ecDNA_geneSig, CorEx_genesDN) # 2


df5 <- data.frame(geneSet = factor(c("CorEx_all","CorEx_UP","CorEx_DN"), levels = c("CorEx_all","CorEx_UP","CorEx_DN")),
                  Overlap = c(rep("Overlap",3),rep("NonOverlap",3)),
                  NumGenes = c(c(47, 43, 2), 224 - c(47, 43, 2)))
df6 <- data.frame(geneSet = factor(c("CorEx_all","CorEx_UP","CorEx_DN"), levels = c("CorEx_all","CorEx_UP","CorEx_DN")),
                  Overlap = c(rep("Overlap",3),rep("NonOverlap",3)),
                  NumGenes = c(c(47, 43, 2), c(643-46, 262-43, 271-1)))

p5 <- ggplot(df5, aes(x = geneSet, y = NumGenes, fill = Overlap)) + 
  geom_bar(stat = "identity", position = "fill") + scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  theme_cb() + scale_fill_manual(values = c("NonOverlap" = "darkgray", "Overlap" = "orange")) + 
  labs(x = "gene set", y = "% genes in ecDNA Transcriptional Signature") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p6 <- ggplot(df6, aes(x = geneSet, y = NumGenes, fill = Overlap)) + 
  geom_bar(stat = "identity", position = "fill") + scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  theme_cb() + scale_fill_manual(values = c("NonOverlap" = "darkgray", "Overlap" = "deepskyblue")) +
  labs(x = "gene set", y = "% genes in CorEx gene sets") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf("output/Plots/08_Barplot_Overlap_ecDNASig_CorEx.pdf", width =3, height = 3)
p5
p6
dev.off()

#Summary sheet
df6 <- data.frame(gene = ecDNA_geneSig, CorEx_AL = "N", CorEx_UP = "N", CorEx_DN = "N")
df6$CorEx_AL[df6$gene %in% ovlpGenes_CorEx_genesAL] <- "Y"
df6$CorEx_UP[df6$gene %in% ovlpGenes_CorEx_genesUP] <- "Y"
df6$CorEx_DN[df6$gene %in% ovlpGenes_CorEx_genesDN] <- "Y"

CorEx_genesDF$overlap_ecDNATS <- "N"
CorEx_genesDF$overlap_ecDNATS[CorEx_genesDF$gene %in% ecDNA_geneSig] <- "Y"

write.csv(df6, "output/Tables/08_OverlapTable_ecDNATS_CorEx.csv")
write.csv(CorEx_genesDF, "output/Tables/08_OverlapTable_CorEx_ecDNATS.csv")

#------------ gene enrichment by Enrichr ------------#
res_enrichr <- fread("output/Tables/MSigDB_Hallmark_2020_table.txt")[,c(1,4)] %>% as.data.frame()
res_enrichr$log <- -log10(res_enrichr$`Adjusted P-value`)

p7 <- ggplot(res_enrichr, aes(x = log, y = reorder(Term, log))) + geom_bar(stat = "identity") + theme_cb() + 
  labs(x = "-log10(P adjusted)", y = "") + scale_x_continuous(expand = c(0,0))
p8 <- ggplot(res_enrichr[which(res_enrichr$log > 2),], aes(x = log, y = reorder(Term, log))) + geom_bar(stat = "identity") + theme_cb() + 
  labs(x = "-log10(P adjusted)", y = "") + scale_x_continuous(expand = c(0,0))

pdf("output/Plots/08_Barplot_GE_Hallmark.pdf", width = 7, height = 6)
p7
p8
dev.off()

#------------ ssGSEA ------------#
res_ssGSEA <- fread("ref/result_ssGSEA_1066samples_20250806.txt")
res_ssGSEA <- as.data.frame(res_ssGSEA)
rownames(res_ssGSEA) <- res_ssGSEA$gene_set
res_ssGSEA <- res_ssGSEA[,-1]
res_ssGSEA <- res_ssGSEA[,SampleInfo$SampleID_RNA]

tmp <- SampleInfo$StudyID
names(tmp) <- SampleInfo$SampleID_RNA

colnames(res_ssGSEA) <- tmp[colnames(res_ssGSEA)]
res_ssGSEA <- res_ssGSEA[, se$StudyID]

difftest_ssGSEA <- mclapply(rownames(res_ssGSEA), function(x){
  df_tmp <- data.frame(value = as.numeric(res_ssGSEA[x,]), ecDNA = se$category1)
  meandiff <- mean(df_tmp$value[which(df_tmp$ecDNA=="ecDNA_p")]) - mean(df_tmp$value[which(df_tmp$ecDNA=="ecDNA_n")])
  p_val <- wilcox.test(df_tmp$value~df_tmp$ecDNA)$p.value
  return(c(meandiff,p_val))
}, mc.cores = 24)

difftest_ssGSEA <- do.call(rbind,difftest_ssGSEA)
difftest_ssGSEA <- as.data.frame(difftest_ssGSEA)
rownames(difftest_ssGSEA) <- rownames(res_ssGSEA)
colnames(difftest_ssGSEA) <- c("MeanDiff","pvalue")
difftest_ssGSEA$FDR <- p.adjust(difftest_ssGSEA$pvalue, method = "fdr")
difftest_ssGSEA$logFDR <- -log10(difftest_ssGSEA$FDR)

write.csv(difftest_ssGSEA, "output/Tables/08_DiffTest_ssGSEA_ecDNA.csv")

difftest_ssGSEA_HM <- difftest_ssGSEA[grep("HALLMARK_",rownames(difftest_ssGSEA)),]
difftest_ssGSEA_HM$logP <- -log10(difftest_ssGSEA_HM$pvalue)
ggplot(difftest_ssGSEA_HM, aes(x=MeanDiff, y=logP)) + geom_point(alpha=0.25) + theme_cb() +
  labs(x="MeanDiff",y="Wilcoxon P-value [-log10]")
res_ssGSEA_HM <- res_ssGSEA[grep("HALLMARK_",rownames(res_ssGSEA)),]

p9 <- lapply(rownames(res_ssGSEA_HM), function(x){
  df_tmp <- data.frame(value = as.numeric(res_ssGSEA_HM[x,]), ecDNA = factor(se$category1, levels = c("ecDNA_n", "ecDNA_p")))
  out <- ggplot(df_tmp, aes(x=ecDNA, y=value, fill=ecDNA)) + 
    geom_jitter_rast(height = 0, width = 0.2, alpha = 0.5, aes(color = ecDNA)) + 
    geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() +
    scale_color_manual(values = ecclass_colors) +
    labs(x = "", y = "ssGSEA score") + theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 6)) + 
    geom_signif(comparisons = list(c("ecDNA_n", "ecDNA_p")), textsize = 2, test = "wilcox.test") + ggtitle(x) + theme(legend.key.size = unit(0.3, 'cm'))
  return(out)
})
pdf("output/Plots/08_Boxplot_ssGSEAscore_ecDNA.pdf", width = 2.5, height = 3)
p9
dev.off()

#------------ Module score ------------#
norm_counts <- assays(se)$log2tpm
gene_ls <- list(CorExUP=intersect(CorEx_genesUP, rowData(se)$gene_name), 
                CorExDN=intersect(CorEx_genesDN, rowData(se)$gene_name), 
                ecDNA_geneSig=ecDNA_geneSig)

addGeneModuleScore <- function(
  expr_matrix,
  gene_list,
  nbin = 24,
  ctrl = 100,
  seed = 123
){
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  
  # gene average and ranking
  genes.avg <- Matrix::rowMeans(x = expr_matrix)
  genes.avg <- sort(genes.avg)
  genes.cut <- ggplot2::cut_number(x = genes.avg + rnorm(n = length(genes.avg))/1e+30, n = nbin)
  names(genes.cut) <- names(genes.avg)
  
  # picking up ctrl genes
  ctrl.genes.list <- vector(mode = "list", length = length(gene_list))
  for (i in 1:length(gene_list)) {
    gene_set <- gene_list[[i]]
    for (j in 1:length(gene_set)) {
      ctrl.genes <- names(sample(genes.cut[which(genes.cut==genes.cut[gene_set[j]])],size = ctrl,replace = FALSE))
      ctrl.genes.list[[i]] <- c(ctrl.genes.list[[i]], ctrl.genes)
    }
  }
  ctrl.genes.list <- lapply(ctrl.genes.list, unique)
  
  # score for ctrl genes
  ctrl.scores <- matrix(data = 1,
                        nrow = length(ctrl.genes.list),
                        ncol = ncol(expr_matrix))
  for (i in 1:length(ctrl.genes.list)){
    ctrl.genes <- ctrl.genes.list[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(expr_matrix[ctrl.genes,])
  }
  
  # score for features of interest
  gene.set.scores <- matrix(data = 1, 
                            nrow = length(gene_list), 
                            ncol = ncol(expr_matrix))
  for (i in 1:length(gene_list)) {
    gene_set <- gene_list[[i]]  
    gene.set.scores[i, ] <- Matrix::colMeans(expr_matrix[gene_set, ,drop = FALSE])
  }
  
  # gene score subtracted by ctrl score
  gene.set.scores.final <- gene.set.scores - ctrl.scores
  
  rownames(gene.set.scores.final) <- names(gene_list)
  colnames(gene.set.scores.final) <- colnames(expr_matrix)
  gene.set.scores.final <- data.frame(t(gene.set.scores.final))
  
  return(gene.set.scores.final)
}

mat <- norm_counts
rownames(mat) <- rowData(se)[rownames(mat),]$gene_name
moduleScoreDF <- addGeneModuleScore(expr_matrix = mat, gene_list = gene_ls)
moduleScoreDF$CancerType <- se$CancerType2
moduleScoreDF$ecDNA <- se$category1

p10 <- lapply(names(gene_ls), function(x){
  mtx_tmp <- data.frame(ModuleScore = moduleScoreDF[,x])
  mtx_tmp$ecDNA <- factor(moduleScoreDF$ecDNA, levels = c("ecDNA_n", "ecDNA_p"))
  out <- ggplot(mtx_tmp, aes(x=ecDNA, y=ModuleScore, fill=ecDNA)) + 
    geom_jitter_rast(height = 0, width = 0.2, alpha = 0.5, aes(color = ecDNA)) + 
    geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() +
    scale_color_manual(values = ecclass_colors) +
    labs(x = "", y = paste0("Module score")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    geom_signif(comparisons = list(c("ecDNA_n", "ecDNA_p")), textsize = 2) + ggtitle(x) + theme(legend.key.size = unit(0.3, 'cm'))
  return(out)
})
pdf("output/Plots/08_Boxplot_ModuleScore_ecDNASig_CorEx.pdf", width = 2.5, height = 3)
p10
dev.off()

#------------ Module score vs CN ------------#
idx1 <- rownames(moduleScoreDF)[which(moduleScoreDF$ecDNA == "ecDNA_p")]
MaxCN_value <- sapply(idx1, function(x) max(AA_result_T$Feature.maximum.copy.number[which(AA_result_T$StudyID == x & AA_result_T$Classification == "ecDNA")]))

df11 <- data.frame(moduleScore = moduleScoreDF[idx1,]$ecDNA_geneSig, CN = MaxCN_value)
p11 <- ggplot(df11, aes(x = CN, y = moduleScore)) + geom_point_rast(alpha = 0.5, size = 2, color ="#E7298A") + theme_cb() + 
  scale_x_continuous(trans = "log2", breaks = c(0,2,4,8,16,32,64,128,256)) + 
  labs(x = "WGS-derived CN", y = "ecDNA geneSignature module Score") +
  geom_smooth(method = "lm", se = F, col = "black", lty = "dashed")
pdf("output/Plots/08_Scatter_ModuleScore_ecDNASig_MaxCN.pdf", width = 3, height = 3)
p11
dev.off()

cor.test(df11$moduleScore, df11$CN)
# Pearson's product-moment correlation
# 
# data:  df11$moduleScore and df11$CN
# t = -0.21548, df = 275, p-value = 0.8296
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1306480  0.1050237
# sample estimates:
#         cor 
# -0.01299257 

#------------ Module score ~ oncogene ------------#
idx_oncogene <- grep("\\[\\]", AA_result_T$Oncogenes, invert = T)
sample_oncogene <- unique(AA_result_T$StudyID[idx_oncogene])

moduleScoreDF$oncogene <- "nononcogene"
moduleScoreDF$oncogene[rownames(moduleScoreDF) %in% sample_oncogene] <- "oncogene"
moduleScoreDF$oncogene <- paste0(moduleScoreDF$ecDNA, "#", moduleScoreDF$oncogene)
moduleScoreDF$oncogene[moduleScoreDF$ecDNA == "ecDNA_n"] <- "ecDNA_n"
moduleScoreDF$oncogene <- factor(moduleScoreDF$oncogene, levels = c("ecDNA_n","ecDNA_p#nononcogene","ecDNA_p#oncogene"))

p12 <- ggplot(moduleScoreDF, aes(x = oncogene, y = ecDNA_geneSig)) +
  geom_jitter_rast(height = 0, width = 0.2, alpha = 0.5, aes(color = oncogene)) + 
  geom_boxplot(outlier.shape = NA, alpha=0) +
  theme_cb() + labs(x = "", y = "Module score") + 
  scale_color_manual(values = c("ecDNA_p#oncogene" = "red", "ecDNA_p#nononcogene" = "#E7298A" ,"ecDNA_n"= "gray")) +
  geom_signif(comparisons = list(c("ecDNA_p#oncogene", "ecDNA_p#nononcogene"), c("ecDNA_p#oncogene", "ecDNA_n"), c("ecDNA_p#nononcogene","ecDNA_n")), test = "wilcox.test",
              y_position = c(3,2.7,2.4))
pdf("output/Plots/08_Boxplot_moduleScore_ecDNASig_oncogene.pdf", width = 5, height = 4)
p12
dev.off()

#------------ Module score correlation to CorEx ------------#
df13 <- data.frame(ecDNA_geneSig = moduleScoreDF[,"ecDNA_geneSig"], 
                    CorExUP = moduleScoreDF[,"CorExUP"],
                    CorExDN = moduleScoreDF[,"CorExDN"])
df13$ecDNA <- factor(moduleScoreDF$ecDNA, levels = c("ecDNA_n", "ecDNA_p"))
df13 <- df13[order(df13$ecDNA),]

p13.1 <- ggplot(df13, aes(x=ecDNA_geneSig, y=CorExUP, color=ecDNA)) + geom_point_rast(alpha = 0.5, size = 2) + 
  scale_color_manual(values = ecclass_colors) + labs(x = "ecDNA geneSig", y = "CorExUP") + theme_cb() +
  geom_smooth(method = "lm", se = F)
p13.2 <- ggplot(df13, aes(x=ecDNA_geneSig, y=CorExDN, color=ecDNA)) + geom_point_rast(alpha = 0.5, size = 2) + 
  scale_color_manual(values = ecclass_colors) + labs(x = "ecDNA geneSig", y = "CorExDN") + theme_cb() +
  geom_smooth(method = "lm", se = F)
p13.3 <- ggplot(df13, aes(x=CorExUP, y=CorExDN, color=ecDNA)) + geom_point_rast(alpha = 0.5, size = 2) + 
  scale_color_manual(values = ecclass_colors) + labs(x = "CorExUP", y = "CorExDN") + theme_cb() +
  geom_smooth(method = "lm", se = F)

pdf("output/Plots/08_ScatterPlot_Correlation_MS_ecDNASig_CorEx.pdf", width = 4, height = 3)
p13.1
p13.2
p13.3
dev.off()

cor.test(df13$ecDNA_geneSig[which(df13$ecDNA == "ecDNA_p")], df13$CorExUP[which(df13$ecDNA == "ecDNA_p")])
# Pearson's product-moment correlation
# 
# data:  df13$ecDNA_geneSig[which(df13$ecDNA == "ecDNA_p")] and df13$CorExUP[which(df13$ecDNA == "ecDNA_p")]
# t = 28.229, df = 275, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.8285734 0.8896781
# sample estimates:
#       cor 
# 0.8622304 

#------------ TCGA RNA-seq ------------#
# TCGA_SampleInfo <- read.csv("/mnt/host_mnt/Users/nextganken/Analysis/S24-2_ecDNA-RE/ref/elife_sampleInfo.csv")
# setwd("/mnt/host_mnt/Volumes/Shared2/Kume/Analysis/TCGAData/RNA/")
# 
# TCGA_SampleType <- names(TCGA_SampleInfo$Tumor %>% table) %>% sort
# 
# # TCGA RNA-seq Counts data
# library(TCGAbiolinks)
# query_ls <- list()
# for(x in TCGA_SampleType){
#   print(x)
#   SampleBarcode <- TCGA_SampleInfo$full_barcode[which(TCGA_SampleInfo$Tumor == x)]
#   query <- GDCquery(project=paste0("TCGA-", x),
#                     data.category="Transcriptome Profiling",
#                     data.type="Gene Expression Quantification",
#                     workflow.type="STAR - Counts", 
#                     barcode = SampleBarcode)
#   query_ls[[x]] <- query
# }
# 
# TCGA_RNAExp <- list()
# for(x in names(query_ls)){
#   print(x)
#   TCGA_RNAExp[[x]] <- GDCprepare(query_ls[[x]])
# }
# 
# se_TCGA <- lapply(names(TCGA_RNAExp), function(x){
#   out <- TCGA_RNAExp[[x]]
#   colData(out) <- colData(out)[,c("barcode","patient","sample","sample_type")]
#   colData(out)$CancerType <- x
#   assays(out)$unstranded <- assays(TCGA_RNAExp[[x]])$unstranded
#   assays(out)$log2tpm <- log2(assays(TCGA_RNAExp[[x]])$tpm_unstrand+1)
#   assays(out) <- assays(out)[which(names(assays(out)) %in% c("unstranded", "log2tpm"))]
#   return(out)
# })
# se_TCGA <- do.call(cbind, se_TCGA) #868 samples
# 
# TCGA_SampleInfo[which(TCGA_SampleInfo$full_barcode %ni% colData(se_TCGA)$barcode),]
# # ecDNA_status Tumor  Sample_barcode          SampleType                 full_barcode Center platform
# # 103     ecDNA(+)  STAD TCGA-CG-4449-01 Primary Solid Tumor TCGA-CG-4449-01A-01R-1157-13  BCGSC       GA
# # 576     ecDNA(-)  STAD TCGA-CG-4474-01 Primary Solid Tumor TCGA-CG-4474-01A-02R-1157-13  BCGSC       GA
# 
# setwd("/mnt/host_mnt/Volumes/NEXTSSD2/Analysis/S24-2_ecDNA-RE/v4/")
# 
# rownames(TCGA_SampleInfo) <- TCGA_SampleInfo$full_barcode
# se_TCGA$ecDNA <- TCGA_SampleInfo[colnames(se_TCGA),]$ecDNA_status
# 
# saveRDS(se_TCGA, "rds/08_se_TCGARNA.rds")

se_TCGA <- readRDS("/mnt/host_mnt/Volumes/NEXTSSD2/Analysis/S24-2_ecDNA-RE/v4/rds/08_se_TCGARNA.rds")
rownames(TCGA_SampleInfo) <- TCGA_SampleInfo$full_barcode
se_TCGA$ecDNA <- TCGA_SampleInfo[colnames(se_TCGA),]$ecDNA_status
table(se_TCGA$ecDNA)
# ecDNA(-) ecDNA(+) 
# 635      233

#module score
gene_ls_TCGA <- list(CorExUP=intersect(CorEx_genesUP, rowData(se_TCGA)$gene_name), 
                     CorExDN=intersect(CorEx_genesDN, rowData(se_TCGA)$gene_name), 
                     ecDNA_geneSig=intersect(ecDNA_geneSig, rowData(se_TCGA)$gene_name))

mat <- assays(se_TCGA)$log2tpm
rownames(mat) <- rowData(se_TCGA)[rownames(mat),]$gene_name
moduleScoreDF_TCGA <- addGeneModuleScore(expr_matrix = mat, gene_list = gene_ls_TCGA)
moduleScoreDF_TCGA$CancerType <- se_TCGA$CancerType
moduleScoreDF_TCGA$ecDNA <- se_TCGA$ecDNA
moduleScoreDF_TCGA$ecDNA[moduleScoreDF_TCGA$ecDNA == "ecDNA(-)"] <- "ecDNA_n"
moduleScoreDF_TCGA$ecDNA[moduleScoreDF_TCGA$ecDNA == "ecDNA(+)"] <- "ecDNA_p"

p14 <- lapply(names(gene_ls_TCGA), function(x){
  mtx_tmp <- data.frame(ModuleScore = moduleScoreDF_TCGA[,x])
  mtx_tmp$ecDNA <- factor(moduleScoreDF_TCGA$ecDNA, levels = c("ecDNA_n", "ecDNA_p"))
  out <- ggplot(mtx_tmp, aes(x=ecDNA, y=ModuleScore, fill=ecDNA)) + 
    geom_jitter_rast(height = 0, width = 0.2, alpha = 0.5, aes(color = ecDNA)) + 
    geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() +
    scale_color_manual(values = ecclass_colors) +
    labs(x = "", y = "Module score") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    geom_signif(comparisons = list(c("ecDNA_n", "ecDNA_p")), textsize = 2) + ggtitle(x) + theme(legend.key.size = unit(0.3, 'cm'))
  return(out)
})
pdf("output/Plots/08_Boxplot_ModuleScore_ecDNASig_CorEx_TCGA.pdf", width = 2.5, height = 3)
p14
dev.off()

df15 <- data.frame(ecDNA_geneSig = moduleScoreDF_TCGA[,"ecDNA_geneSig"], 
                   CorExUP = moduleScoreDF_TCGA[,"CorExUP"],
                   CorExDN = moduleScoreDF_TCGA[,"CorExDN"])
df15$ecDNA <- factor(moduleScoreDF_TCGA$ecDNA, levels = c("ecDNA_n", "ecDNA_p"))
df15 <- df15[order(df15$ecDNA),]

p15.1 <- ggplot(df15, aes(x=ecDNA_geneSig, y=CorExUP, color=ecDNA)) + geom_point_rast(alpha = 0.5, size = 2) + 
  scale_color_manual(values = ecclass_colors) + labs(x = "ecDNA geneSig", y = "CorExUP") + theme_cb() +
  geom_smooth(method = "lm", se = F)
p15.2 <- ggplot(df15, aes(x=ecDNA_geneSig, y=CorExDN, color=ecDNA)) + geom_point_rast(alpha = 0.5, size = 2) + 
  scale_color_manual(values = ecclass_colors) + labs(x = "ecDNA geneSig", y = "CorExDN") + theme_cb() +
  geom_smooth(method = "lm", se = F)
p15.3 <- ggplot(df15, aes(x=CorExUP, y=CorExDN, color=ecDNA)) + geom_point_rast(alpha = 0.5, size = 2) + 
  scale_color_manual(values = ecclass_colors) + labs(x = "CorExUP", y = "CorExDN") + theme_cb() +
  geom_smooth(method = "lm", se = F)

pdf("output/Plots/08_ScatterPlot_Correlation_MS_ecDNASig_CorEx_TCGA.pdf", width = 4, height = 3)
p15.1
p15.2
p15.3
dev.off()

cor.test(df15$ecDNA_geneSig[which(df15$ecDNA == "ecDNA_p")], df15$CorExUP[which(df15$ecDNA == "ecDNA_p")])
# Pearson's product-moment correlation
# 
# data:  df15$ecDNA_geneSig[which(df15$ecDNA == "ecDNA_p")] and df15$CorExUP[which(df15$ecDNA == "ecDNA_p")]
# t = 13.612, df = 231, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.5891408 0.7328339
# sample estimates:
#       cor
# 0.6671477
