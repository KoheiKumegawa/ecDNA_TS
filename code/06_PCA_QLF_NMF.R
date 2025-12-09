#-------------------------------------------------------------------------------
# 06_PCA_QLF_NMF.R
#-------------------------------------------------------------------------------
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggrastr)
library(ggrepel)
"%ni%" <- Negate("%in%")
theme_cb <- function(color="black"){theme_classic() + theme(axis.text = element_text(color=color), axis.ticks = element_line(color=color))}

se <- readRDS("rds/04_se.rds")
SampleInfo <- read.csv("output/Tables/00_SampleListSummary_Curated.csv", row.names = 1)
cancertype_colors <- ArchR::ArchRPalettes$stallion2[c(1:17)]
names(cancertype_colors) <- table(SampleInfo$CancerType2) %>% sort(., decreasing = T) %>%names()
ampclass_colors <- c("ecDNA" = "#E7298A", "BFB"= "#7570B3", "Complex-non-cyclic" = "#D95F02", "Linear" = "#1B9E77")

#---------- Sample classification ----------#
AA_result_T <- read.csv("output/Tables/02_AA_Summary_Tumor_filtered.csv")
tmp <- SampleInfo$StudyID
names(tmp) <- SampleInfo$SampleID_T
AA_result_T$StudyID <- tmp[AA_result_T$Sample.name]
AA_result_T$CancerType <- SampleInfo[AA_result_T$StudyID,]$CancerType2

#ecDNA samples
ecDNA_p <- AA_result_T$StudyID[AA_result_T$Classification == "ecDNA"] %>% unique() #277 sample w/ ecDNA
BFBCnC_p <- AA_result_T$StudyID[AA_result_T$Classification %in% c("BFB","Complex-non-cyclic")] %>% unique() #288 sample w/ ecDNA
intersect(ecDNA_p, BFBCnC_p) %>% length #149 overlaps
BFBCnC_p <- BFBCnC_p[BFBCnC_p %ni% ecDNA_p] #139 samples

se <- se[,se$StudyID %ni% BFBCnC_p]
se$category <- "ecDNA_n"
se$category[which(se$StudyID %in% ecDNA_p)] <- "ecDNA_p"

###check###
ecDNA_n <- se$StudyID[se$category == "ecDNA_n"]
AA_result_T[AA_result_T$StudyID %in% ecDNA_n, ]$Classification
table(se$category)
# ecDNA_n ecDNA_p 
# 463     277
###end###

#---------- Principal Component Analysis ----------#
se_pc <- se[rowData(se)$gene_type == "protein_coding",]
saveRDS(se_pc, "rds/06_se_pc.rds")

row_vars1 <- rowVars(assays(se_pc)$log2tpm)
upper25perc <- which(row_vars1 > quantile(row_vars1, 0.75)) # 4842 genes
pca1 <- prcomp(t(assays(se_pc)$log2tpm[upper25perc,]))

### variance explained ###
pca1_summ <- summary(pca1)

pca1_bstick <- vegan::bstick(pca1) #baseline variance (broken stick)
pca1_bstick <- pca1_bstick / sum(pca1_summ$importance[1,]^2) #baseline variance explained
pca1_bstick <- pca1_bstick*100 #percent
pca1_bstick <- pca1_bstick[c(1:30)]

df1 <- data.frame(PC = c(1:30), 
                  observed = pca1_summ$importance[2, c(1:30)]*100, 
                  baseline = pca1_bstick)
df1 <- reshape2::melt(df1, id = "PC")
p1 <- ggplot(df1, aes(x = PC, y = value, color = variable, group = variable)) + geom_line() + geom_point() + 
  scale_y_log10() + theme_cb() + scale_color_manual(values = c("baseline" = "darkgray", "observed" = "darkblue")) +
  labs(x = "PC", y = "% variance explained")
pdf("output/Plots/061_Lineplot_PCA_VE.pdf", width = 6, height = 2.5)
p1
dev.off()

### Average within sample type ###
pca1_avrg <- lapply(names(cancertype_colors), function(x){
  idx <- colnames(se_pc)[se_pc$CancerType2 == x]
  out <- lapply(c(1:22), function(n) mean(pca1$x[idx,n])) %>% unlist
  return(out)
})
pca1_avrg <- do.call(cbind,pca1_avrg)
colnames(pca1_avrg) <- names(cancertype_colors)
rownames(pca1_avrg) <- c(1:22)

df2 <- reshape2::melt(pca1_avrg)
p2 <- ggplot(df2, aes(x = Var1, y = value, color = Var2)) + geom_point(alpha = 0.8, size = 2.5) + theme_cb() + 
  scale_color_manual(values = cancertype_colors) + geom_hline(yintercept = 0, lty = "dotted") +
  labs(x = "PC", y = "Average PC values for each SampleType") +
  theme(legend.key.size = unit(0.2, 'cm'))
pdf("output/Plots/06_Dotplot_PCA_Avrg_SampleType.pdf", width = 7, height = 2.5)
p2
dev.off()
write.csv(pca1_avrg, "output/Tables/06_PCA_Avrg_SampleType.csv")

### standard deviation within sample type ###
df_PC <- pca1$x[, c(1:22)] %>% data.frame
df_PC$CancerType <- factor(se_pc$CancerType2, levels = names(cancertype_colors))

pca1_sd <- lapply(names(cancertype_colors), function(x){
  idx <- colnames(se_pc)[se_pc$CancerType2 == x]
  out <- lapply(c(1:22), function(n) sd(pca1$x[idx,n])) %>% unlist
  return(out)
})
pca1_sd <- do.call(cbind,pca1_sd)
colnames(pca1_sd) <- names(cancertype_colors)
rownames(pca1_sd) <- c(1:22)

df3 <- reshape2::melt(pca1_sd)
p3 <- ggplot(df3, aes(x = Var1, y = value, color = Var2)) + geom_point(alpha = 0.6, size = 2.5) + theme_cb() + 
  scale_color_manual(values = cancertype_colors) + labs(x = "PC", y = "SD within SampleType") +
  theme(legend.key.size = unit(0.2, 'cm'))
pdf("output/Plots/06_Dotplot_PCA_SD_SampleType.pdf", width = 7, height = 2.5)
p3
dev.off()
write.csv(pca1_sd, "output/Tables/06_PCA_SD_SampleType.csv")

### PC visualization ###
p4 <- ggplot(df_PC, aes(x = PC1, y = PC2, color = CancerType)) + geom_point(alpha = 0.6, size = 2)  + 
  geom_vline(xintercept = 0, lty = "dotted", size = 0.5) + geom_hline(yintercept = 0, lty = "dotted", size = 0.5) +
  theme_cb() + scale_color_manual(values = cancertype_colors) +
  theme(legend.key.size = unit(0.2, 'cm'))
p5 <- ggplot(df_PC, aes(x = PC11, y = PC18, color = CancerType)) + geom_point(alpha = 0.6, size = 2)  + 
  geom_vline(xintercept = 0, lty = "dotted", size = 0.5) + geom_hline(yintercept = 0, lty = "dotted", size = 0.5) +
  theme_cb() + scale_color_manual(values = cancertype_colors) +
  theme(legend.key.size = unit(0.2, 'cm'))
pdf("output/Plots/06_PCAplot_select.pdf", width = 4, height = 2.5)
p4
p5
dev.off()

### PC loading ###
pca1_loading <- t(t(pca1$rotation)*pca1$sdev)[,c(1:22)] %>% as.data.frame()
pca1_loading$gene_name <- rowData(se_pc[rownames(pca1_loading),])$gene_name
write.csv(pca1_loading, "output/Tables/06_PCA_loading.csv")

#---------- QLF test for DEG between ecDNA status ----------#
source("code/edgeR_PairwiseFunction.R")
DiffTest <- edgeR_pairwise(se_pc, compareCol = "category", topGroup = "ecDNA_p", bottomGroup = "ecDNA_n")
DiffTestDF <- as.data.frame(assay(DiffTest))
DiffTestDF$gene_name <- rowData(se_pc)$gene_name
DiffTestDF$gene_name2 <- DiffTestDF$gene_name
DiffTestDF$gene_name2[which(DiffTestDF$FDR >= 0.01)] <- ""
write.csv(DiffTestDF, "output/Tables/06_DiffTest_QLF_posneg_ecDNA.csv")

p6 <- ggplot(DiffTestDF, aes(x = log2FC, y = -log10(FDR))) + 
  geom_point_rast(shape = 1) + 
  geom_vline(xintercept = c(-1,1), lty = "dotted") + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 2, lty = "dotted") +
  theme_cb() + labs(x = "Log2FC", y = "-log10(FDR)")

ecDNA_UP <- DiffTestDF$gene_name[which(DiffTestDF$log2FC > 1 & DiffTestDF$FDR < 0.01)]

expr_ecDNAUP <- colMeans(assays(se_pc)$log2tpm[which(rowData(se_pc)$gene_name %in% ecDNA_UP), ]) %>% as.data.frame() %>% `colnames<-`("expr")
expr_ecDNAUP$CancerType <- factor(se_pc$CancerType2, levels = names(cancertype_colors))
expr_ecDNAUP$category <- factor(se_pc$category, levels = c("ecDNA_n", "ecDNA_p"))

ecclass_colors <- c("ecDNA_p" = "#E7298A", "ecDNA_n" = "gray")
p7.1 <- ggplot(expr_ecDNAUP, aes(x=category,y=expr)) +
  geom_jitter(height = 0, width = 0.2, shape = 19, alpha = 1, aes(color = CancerType)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() +
  scale_color_manual(values = cancertype_colors) +
  labs(x = "ecDNA status", y = "Average Expression") + theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  geom_signif(comparisons = list(c("ecDNA_p", "ecDNA_n")), test = "wilcox.test")
p7.2 <- ggplot(expr_ecDNAUP, aes(x=CancerType,y=expr,fill=category)) +
  geom_jitter(shape = 19, alpha = 1, aes(color = category),
              position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() +
  scale_color_manual(values = ecclass_colors) +
  labs(x = "Cancer type", y = "Average Expression") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

df8 <- table(se_pc$CancerType2, se_pc$category) %>% reshape2::melt(.)
df8$Var1 <- factor(df8$Var1, levels = names(cancertype_colors))
p8 <- ggplot(df8, aes(x = Var2, y = value, fill = Var1)) + geom_bar(stat = "identity", position = "fill", color = "black", size = 0.25) +
  theme_cb() + scale_y_continuous(expand = c(0,0)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = cancertype_colors) +
  labs(x = "", y = "%Sample")

pdf("output/Plots/06_VolcanoPlot_QLF_ecDNAstatus.pdf", width = 5, height = 5)
p6
dev.off()
pdf("output/Plots/06_Boxplot_QLF_ecDNAstatus.pdf", width = 4, height = 5)
p7.1
dev.off()
pdf("output/Plots/06_Boxplot_QLF_ecDNAstatus_CancerType.pdf", width = 8, height = 5)
p7.2
dev.off()

pdf("output/Plots/06_Barplot_ecDNAstatus_CancerType.pdf", width = 3, height = 5)
p8
dev.off()

#---------- NMF input mat ----------#
mat1 <- assays(se_pc)$log2tpm[upper25perc,]
mat1[c(1:5), c(1:5)]
write.csv(mat1, "output/Tables/RNMF_InputMatrix.csv")

write.csv(data.frame(StudyID = se_pc$StudyID, ecDNA = se_pc$category, CancerType = se_pc$CancerType2), 
          "output/Tables/06_Table_StudyID_ecDNAstatus.csv")
