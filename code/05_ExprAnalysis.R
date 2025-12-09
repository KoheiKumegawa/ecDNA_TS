#-------------------------------------------------------------------------------
# 05_ExprAnalysis.R
#-------------------------------------------------------------------------------
library(SummarizedExperiment)
library(parallel)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggrastr)
library(ggrepel)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(viridisLite)
library(stringr)
library(scales)
"%ni%" <- Negate("%in%")
theme_cb <- function(color="black"){theme_classic() + theme(axis.text = element_text(color=color), axis.ticks = element_line(color=color))}

se <- readRDS("rds/04_se.rds")
SampleInfo <- read.csv("output/Tables/00_SampleListSummary_Curated.csv", row.names = 1)
gr_amp <- readRDS("rds/02_gr_T_filtered.rds")

cancertype_colors <- ArchR::ArchRPalettes$stallion2[c(1:17)]
names(cancertype_colors) <- table(SampleInfo$CancerType2) %>% sort(., decreasing = T) %>%names()
ecclass_colors <- c("ecDNA" = "#E7298A", "ChromAmp" = "yellow3", "No-fSCNA" = "gray")
ampclass_colors <- c("ecDNA" = "#E7298A", "BFB"= "#7570B3", "Complex-non-cyclic" = "#D95F02", "Linear" = "#1B9E77")

#---------- Constructing a correlation matrix of GEx and amplicons ----------#
fo1 <- findOverlaps(gr_amp, rowRanges(se))

AA_result_T <- read.csv("output/Tables/02_AA_Summary_Tumor_filtered.csv")
tmp <- SampleInfo$StudyID
names(tmp) <- SampleInfo$SampleID_T
AA_result_T$StudyID <- tmp[AA_result_T$Sample.name]
AA_result_T$CancerType <- SampleInfo[AA_result_T$StudyID,]$CancerType2

AA_Expr_DF <- mclapply(seq_along(fo1), function(x){
  ampID <- queryHits(fo1)[x]
  genID <- subjectHits(fo1)[x]
  out <- data.frame(Sample = AA_result_T$Sample.name[ampID],
                    ID = AA_result_T$Feature.ID[ampID],
                    Class = AA_result_T$Classification[ampID], 
                    CN = AA_result_T$Feature.maximum.copy.number[ampID],
                    EnsemblID = rownames(se)[genID],
                    gene_name = rowData(se)$gene_name[genID],
                    gene_type = rowData(se)$gene_type[genID],
                    StudyID = AA_result_T$StudyID[ampID],
                    CancerType = AA_result_T$CancerType[ampID])
  return(out)
}, mc.cores = 24)
AA_Expr_DF <- do.call(rbind, AA_Expr_DF)

#check#
head(AA_Expr_DF)
#check-end#

## TPM values ##
tpm_value <- mclapply(c(1:nrow(AA_Expr_DF)), function(i){
  genID <- which(rownames(se) == AA_Expr_DF$EnsemblID[i])
  samID <- AA_Expr_DF$StudyID[i]
  out <- assays(se)$log2tpm[genID,samID]
  out <- 2^out-1
  return(out)
}, mc.cores = 24)
tpm_value <- unlist(tpm_value)
AA_Expr_DF$TPM <- tpm_value

## No-fSCNA TPM values ##
tpm_value_nofsca <- mclapply(c(1:nrow(AA_Expr_DF)), function(i){
  genID <- which(rownames(se) == AA_Expr_DF$EnsemblID[i])
  samID_amp <- AA_Expr_DF$StudyID[which(AA_Expr_DF$EnsemblID == AA_Expr_DF$EnsemblID[i])]
  samID <- se$StudyID[which(se$StudyID %ni% samID_amp & se$CancerType2 == AA_Expr_DF$CancerType[i])]
  out <- assays(se)$log2tpm[genID,samID]
  out <- mean(2^out-1)
  return(out)
}, mc.cores = 24)
tpm_value_nofsca <- unlist(tpm_value_nofsca)
AA_Expr_DF$TPM_nofsca <- tpm_value_nofsca

## Fold change ##
AA_Expr_DF$FC <- (AA_Expr_DF$TPM+1) / (AA_Expr_DF$TPM_nofsca+1)
AA_Expr_DF$FC2 <- AA_Expr_DF$TPM / AA_Expr_DF$TPM_nofsca
AA_Expr_DF$Class <- factor(AA_Expr_DF$Class, levels = c("ecDNA", "BFB", "Complex-non-cyclic", "Linear"))

plot(log2(AA_Expr_DF$FC), log2(AA_Expr_DF$FC2))

## oncogene annotation ##
AA_oncogenes <- mclapply(c(1:nrow(AA_result_T)), function(i){
  input_string <- AA_result_T[i,]$Oncogenes
  clean_string <- gsub("\\[", "", input_string)
  clean_string <- gsub("\\]", "", clean_string)
  clean_string <- gsub("'", "", clean_string)
  clean_string <- gsub(" ", "", clean_string)
  result <- unlist(strsplit(clean_string, ","))
  return(result)
}, mc.cores = 20)
AA_oncogenes_UQ <- sort(unique(unlist(AA_oncogenes))) #522 unique genes

AA_Expr_DF$oncogene <- "Nononcogene"
AA_Expr_DF$oncogene[which(AA_Expr_DF$gene_name %in% AA_oncogenes)] <- "oncogene"
write.csv(AA_Expr_DF, "output/Tables/05_AA_Expr_DF_preCurated.csv")

#---------- Top amplified oncogene GEx ----------#
AA_Expr_DF <- read.csv("output/Tables/05_AA_Expr_DF_preCurated.csv", row.names = 1)
AA_Expr_DF$Class <- factor(AA_Expr_DF$Class, levels = c("ecDNA", "BFB", "Complex-non-cyclic", "Linear"))

oncogene_Amp_df <- read.csv("output/Tables/03_Number_oncogene_amplicon.csv",row.names = 1)
geneOI <- sapply(names(cancertype_colors), function(x) colSums(oncogene_Amp_df[which(oncogene_Amp_df$CancerType == x & oncogene_Amp_df$AmpClass == "ecDNA"),c(1:522)]) %>% sort(., decreasing = T) %>% names() %>% .[c(1:5)])
geneOI <- reshape2::melt(geneOI)
geneOI <- geneOI[-c(7,10),]
#geneOI <- geneOI[-c(7,10,37:40,45,68:70,76:90),]

p1 <- lapply(c(1:nrow(geneOI)), function(i){
  print(i)
  ampID <- AA_Expr_DF$ID[AA_Expr_DF$gene_name == geneOI$value[i] & AA_Expr_DF$CancerType == geneOI$Var2[i]] %>% unique
  
  df1 <- data.frame(row.names = colnames(se), expr = assays(se)$log2tpm[rowRanges(se)$gene_name == geneOI$value[i]], CancerType = se$CancerType2)
  df1 <- df1[which(df1$CancerType == geneOI$Var2[i]),]
  
  tmp <- data.frame(row.names = AA_result_T$StudyID[which(AA_result_T$Feature.ID %in% ampID)],
                    Classification = AA_result_T$Classification[which(AA_result_T$Feature.ID %in% ampID)],
                    CN = AA_result_T$Feature.maximum.copy.number[which(AA_result_T$Feature.ID %in% ampID)])
  df1$AmpClass <- tmp[rownames(df1),"Classification"]
  df1$AmpClass[is.na(df1$AmpClass)] <- "No-fSCA"
  df1$AmpClass <- factor(df1$AmpClass, levels = c("ecDNA", "BFB", "Complex-non-cyclic", "Linear", "No-fSCA"))
  
  compList <- as.character(sort(unique(df1$AmpClass)))
  compList <- compList[which(compList %ni% c("ecDNA", "No-fSCA"))]
  compList <- lapply(compList, function(n) c("ecDNA", n))
  compList <- c(compList, list(c("ecDNA","No-fSCA")))
  
  out <- ggplot(df1, aes(x = AmpClass, y = expr)) + 
    geom_jitter(height = 0, width = 0.1, shape = 3, alpha = 1, aes(color = AmpClass)) + 
    stat_summary(fun.y = "median", geom = "crossbar", width = 0.5) + 
    theme_cb() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_color_manual(values = c(ampclass_colors, "No-fSCA" = "gray40")) +
    labs(x = "Amplicon Classification", y = "log2(TPM+1)") + ggtitle(paste0(geneOI$value[i], ": ", geneOI$Var2[i])) +
    geom_signif(comparisons = compList, test = "wilcox.test")
  return(out)
})

pdf("output/Plots/05_Jitterplot_Log2TPM_freqAmpOncogenes.pdf", width = 4.5, height = 6)
p1
dev.off()

# CN vs Expr
p2 <- lapply(c(1:nrow(geneOI)), function(i){
  print(i)
  ampID <- AA_Expr_DF$ID[AA_Expr_DF$gene_name == geneOI$value[i] & AA_Expr_DF$CancerType == geneOI$Var2[i]] %>% unique
  
  df2 <- data.frame(row.names = colnames(se), expr = assays(se)$log2tpm[rowRanges(se)$gene_name == geneOI$value[i]], CancerType = se$CancerType2)
  df2 <- df2[which(df2$CancerType == geneOI$Var2[i]),]
  
  tmp <- data.frame(row.names = AA_result_T$StudyID[which(AA_result_T$Feature.ID %in% ampID)],
                    Classification = AA_result_T$Classification[which(AA_result_T$Feature.ID %in% ampID)],
                    CN = AA_result_T$Feature.maximum.copy.number[which(AA_result_T$Feature.ID %in% ampID)])
  df2$AmpClass <- tmp[rownames(df2),"Classification"]
  df2 <- df2[which(is.na(df2$AmpClass) == F),]
  df2$AmpClass <- factor(df2$AmpClass, levels = c("ecDNA", "BFB", "Complex-non-cyclic", "Linear"))
  df2$CN <- tmp[rownames(df2),"CN"]
  
  out <- ggplot(df2, aes(x = CN, y = expr, color = AmpClass)) + geom_point(shape = 3) + theme_cb() +
    scale_color_manual(values = ampclass_colors) +
    scale_x_continuous(trans = "log2", breaks = c(0,2,4,8,16,32,64,128,256)) +
    geom_smooth(method = "lm", se = F) + ggtitle(paste0(geneOI$value[i], ": ", geneOI$Var2[i])) +
    labs(x = "WGS-derived CN (Maximum segment)", y = "Expression log2(TPM+1)")
  return(out)
})
pdf("output/Plots/05_Scatterplot_Log2TPM_CN_freqAmpOncogenes.pdf", width = 6, height = 4)
p2
dev.off()

#---------- lncRNA near CCND1 in headneck ----------#
# gene expression per amplicon class, overlapping genes (the same intervals on CCND1)
ampID <- AA_Expr_DF$ID[AA_Expr_DF$gene_name == "AP000439.3" & AA_Expr_DF$CancerType == "Headneck"] %>% unique

df3 <- data.frame(row.names = colnames(se), expr = assays(se)$log2tpm[rowRanges(se)$gene_name == "AP000439.3"], CancerType = se$CancerType2)
df3 <- df3[which(df3$CancerType == "Headneck"),]

tmp <- data.frame(row.names = AA_result_T$StudyID[which(AA_result_T$Feature.ID %in% ampID)],
                  Classification = AA_result_T$Classification[which(AA_result_T$Feature.ID %in% ampID)],
                  CN = AA_result_T$Feature.maximum.copy.number[which(AA_result_T$Feature.ID %in% ampID)])
df3$AmpClass <- tmp[rownames(df3),"Classification"]
df3$AmpClass[is.na(df3$AmpClass)] <- "No-fSCA"
df3$AmpClass <- factor(df3$AmpClass, levels = c("ecDNA", "BFB", "Complex-non-cyclic", "Linear", "No-fSCA"))

compList <- as.character(sort(unique(df3$AmpClass)))
compList <- compList[which(compList %ni% c("ecDNA", "No-fSCA"))]
compList <- lapply(compList, function(n) c("ecDNA", n))
compList <- c(compList, list(c("ecDNA","No-fSCA")))

p3 <- ggplot(df3, aes(x = AmpClass, y = expr)) + 
  geom_jitter(height = 0, width = 0.1, shape = 3, alpha = 1, aes(color = AmpClass)) + 
  stat_summary(fun.y = "median", geom = "crossbar", width = 0.5) + 
  theme_cb() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values = c(ampclass_colors, "No-fSCA" = "gray40")) +
  labs(x = "Amplicon Classification", y = "log2(TPM+1)") + ggtitle("AP000439.3: Headneck") +
  geom_signif(comparisons = compList, test = "wilcox.test")
pdf("output/Plots/05_Jitterplot_Log2TPM_AP000439.3.pdf", width = 4.5, height = 6)
p3
dev.off()

#---------- Relative expression of each amplicon per each gene type ----------#
# removing genes with sample TPM ≤ 5 & FC over mean+5SD
max(AA_Expr_DF$FC) # 1787.689
cutoff <- mean(AA_Expr_DF$FC) + sd(AA_Expr_DF$FC) * 5
AA_Expr_DFV <- AA_Expr_DF[which(AA_Expr_DF$TPM > 5 & AA_Expr_DF$FC < cutoff),] #8563 genes

AA_Expr_DFV$oncogene2 <- AA_Expr_DFV$oncogene
AA_Expr_DFV$oncogene2[AA_Expr_DFV$gene_type != "protein_coding"] <- "noncoding"
AA_Expr_DFV$oncogene2 <- factor(AA_Expr_DFV$oncogene2, levels = c("oncogene", "Nononcogene", "noncoding"))

p4 <- ggplot(AA_Expr_DFV, aes(x = Class, y = log2(FC+1))) + 
  geom_jitter_rast(height = 0, width = 0.2, shape = 3, alpha = 0.2, aes(color = Class)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() +
  scale_color_manual(values = ampclass_colors) +
  labs(x = "Amplicon Type", y = "log2FC") + theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  geom_signif(comparisons = list(c("ecDNA", "BFB"), c("ecDNA", "Complex-non-cyclic"), c("ecDNA", "Linear")), test = "wilcox.test", y_position = c(6.5,6.8,7.1)) + 
  facet_grid(~oncogene2) 
pdf("output/Plots/05_Boxplot_RelativeExpr_eachClassGenetype.pdf", width = 8, height = 6)
p4
dev.off()

#---------- Correlation of Expr & CN ----------#
#All cancer types
p5.1 <- ggplot(AA_Expr_DFV, aes(x = CN, y = log2(FC+1), color = Class)) + 
  geom_point_rast(shape = 3, alpha = 0.2) + theme_cb() +
  scale_color_manual(values = ampclass_colors) +
  scale_x_continuous(trans = "log2", breaks = c(0,2,4,8,16,32,64,128,256)) +
  geom_smooth(method = "glm", se = F) + ggtitle("All genes, all cancer types") +
  labs(x = "WGS-derived CN (Maximum segment)", y = "log2(FC+1)")
p5.2 <- ggplot(AA_Expr_DFV[which(AA_Expr_DFV$oncogene2 == "oncogene"),], aes(x = CN, y = log2(FC+1), color = Class)) + 
  geom_point_rast(shape = 3, alpha = 0.2) + theme_cb() +
  scale_color_manual(values = ampclass_colors) +
  scale_x_continuous(trans = "log2", breaks = c(0,2,4,8,16,32,64,128,256)) +
  geom_smooth(method = "glm", se = F) + ggtitle("Oncogenes, all cancer types") +
  labs(x = "WGS-derived CN (Maximum segment)", y = "log2(FC+1)")
p5.3 <- ggplot(AA_Expr_DFV[which(AA_Expr_DFV$oncogene2 == "Nononcogene"),], aes(x = CN, y = log2(FC+1), color = Class)) + 
  geom_point_rast(shape = 3, alpha = 0.2) + theme_cb() +
  scale_color_manual(values = ampclass_colors) +
  scale_x_continuous(trans = "log2", breaks = c(0,2,4,8,16,32,64,128,256)) +
  geom_smooth(method = "glm", se = F) + ggtitle("Non-oncogenes, all cancer types") +
  labs(x = "WGS-derived CN (Maximum segment)", y = "log2(FC+1)")
p5.4 <- ggplot(AA_Expr_DFV[which(AA_Expr_DFV$oncogene2 == "noncoding"),], aes(x = CN, y = log2(FC+1), color = Class)) + 
  geom_point_rast(shape = 3, alpha = 0.2) + theme_cb() +
  scale_color_manual(values = ampclass_colors) +
  scale_x_continuous(trans = "log2", breaks = c(0,2,4,8,16,32,64,128,256)) +
  geom_smooth(method = "glm", se = F) + ggtitle("Non-coding genes, all cancer types") +
  labs(x = "WGS-derived CN (Maximum segment)", y = "log2(FC+1)")

pdf("output/Plots/05_Scatter_CN_Expr_AmpClassAllCancertype.pdf", width = 6, height = 5)
p5.1
p5.2
p5.3
p5.4
dev.off()

#each cancer types
p6 <- lapply(names(cancertype_colors), function(x){
  df6 <- AA_Expr_DFV[which(AA_Expr_DFV$CancerType == x),]
  res1 <- ggplot(df6, aes(x = CN, y = log2(FC+1), color = Class)) + 
    geom_point_rast(shape = 3, alpha = 0.2) + theme_cb() +
    scale_color_manual(values = ampclass_colors) +
    scale_x_continuous(trans = "log2", breaks = c(0,2,4,8,16,32,64,128,256)) +
    geom_smooth(method = "glm", se = F) + ggtitle(paste0("All genes, ", x)) +
    labs(x = "WGS-derived CN (Maximum segment)", y = "log2(FC+1)")
  
  res2 <- lapply(c("oncogene", "Nononcogene", "noncoding"), function(y){
    ggplot(df6[which(df6$oncogene2 == y),], aes(x = CN, y = log2(FC+1), color = Class)) + 
      geom_point_rast(shape = 3, alpha = 0.2) + theme_cb() +
      scale_color_manual(values = ampclass_colors) +
      scale_x_continuous(trans = "log2", breaks = c(0,2,4,8,16,32,64,128,256)) +
      geom_smooth(method = "glm", se = F) + ggtitle(paste0(y, ": ", x)) +
      labs(x = "WGS-derived CN (Maximum segment)", y = "log2(FC+1)")
  })
  out <- c(list(res1), res2)
  return(out)
})
pdf("output/Plots/05_Scatter_CN_Expr_AmpClassEachCancertype.pdf", width = 6, height = 5)
p6
dev.off()
