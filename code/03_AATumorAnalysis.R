#-------------------------------------------------------------------------------
# 03_AATumorAnalysis.R
#-------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(viridisLite)
library(ggridges)
library(ggsignif)
library(scales)
library(ggrastr)
"%ni%" <- Negate("%in%")
theme_cb <- function(color="black"){theme_classic() + theme(axis.text = element_text(color=color), axis.ticks = element_line(color=color))}

SampleInfo <- read.csv("output/Tables/00_SampleListSummary_Curated.csv", row.names = 1)

#---------- # of amplicon & sample ----------#
cancertype_colors <- ArchR::ArchRPalettes$stallion2[c(1:17)]
names(cancertype_colors) <- table(SampleInfo$CancerType2) %>% sort(., decreasing = T) %>% names()
ecclass_colors <- c("ecDNA" = "#E7298A", "ChromAmp" = "yellow3", "No-fSCNA" = "gray")
ampclass_colors <- c("ecDNA" = "#E7298A", "BFB"= "#7570B3", "Complex-non-cyclic" = "#D95F02", "Linear" = "#1B9E77")

AA_result_T <- read.csv("output/Tables/02_AA_Summary_Tumor_filtered.csv")

tmp <- SampleInfo$StudyID
names(tmp) <- SampleInfo$SampleID_T
AA_result_T$StudyID <- tmp[AA_result_T$Sample.name]
AA_result_T$CancerType <- SampleInfo[AA_result_T$StudyID,]$CancerType2

write.csv(AA_result_T, "output/Tables/02_AA_Summary_Tumor_filtered_v2.csv")

# number of samples
nrow(SampleInfo) #879 samples
ecDNA_p <- AA_result_T$StudyID[AA_result_T$Classification == "ecDNA"] %>% unique() #277 sample w/ ecDNA
ecDNA_n <- AA_result_T$StudyID[AA_result_T$StudyID %ni% ecDNA_p] %>% unique() #207 sample w/ ecDNA

SampleInfo$category <- "No-fSCNA"
SampleInfo$category[SampleInfo$StudyID %in% ecDNA_n] <- "ChromAmp"
SampleInfo$category[SampleInfo$StudyID %in% ecDNA_p] <- "ecDNA"

tmp <- lapply(SampleInfo$StudyID, function(x){
  res <- AA_result_T$Classification[AA_result_T$StudyID == x]
  out <- c("ecDNA" = length(which(res == "ecDNA")),
           "BFB" = length(which(res == "BFB")),
           "Complex-non-cyclic" = length(which(res == "Complex-non-cyclic")),
           "Linear" = length(which(res == "Linear")))
  return(out)
})
AmpClassSample <- do.call(rbind, tmp)
rownames(AmpClassSample) <- SampleInfo$StudyID
max(AmpClassSample) #12

mtx1 <- t(AmpClassSample)
mtx1 <- lapply(names(cancertype_colors), function(x){
  out <- mtx1[,which(SampleInfo$CancerType2 %in% x)]
  out <- out[,order(out[4,],decreasing = T)]
  out <- out[,order(out[3,],decreasing = T)]
  out <- out[,order(out[2,],decreasing = T)]
  out <- out[,order(out[1,],decreasing = T)]
  return(out)
})
mtx1 <- do.call(cbind,mtx1)

col_fun1 <- colorRamp2(c(0,1,12), c("lightgray","orange","red"))
fh <- function(x) hclust(dist(x), method = "ward.D2")

column_ha1 = HeatmapAnnotation(NumAmp = anno_barplot(colSums(mtx1), gp=gpar(col ="gray30",fill="gray30")),
                               NumEc = anno_barplot(mtx1[1,], gp=gpar(col ="#E7298A",fill="#E7298A")))
column_ha2 = HeatmapAnnotation(CancerType = SampleInfo[colnames(mtx1),]$CancerType2, AmpCategory = SampleInfo[colnames(mtx1),]$category,
                               col = list(CancerType=cancertype_colors, AmpCategory=ecclass_colors))
ht1 <- Heatmap(mtx1, name = "Amp.occurence", col = col_fun1, cluster_columns = F, cluster_rows = F,
               show_column_names = F, top_annotation = c(column_ha1,column_ha2),
               column_split = factor(SampleInfo[colnames(mtx1),]$CancerType2, levels = names(cancertype_colors)),
               heatmap_legend_param = list(col_fun = col_fun1, at = c(0, 1, 5, 10, 15)))
p1 <- draw(ht1)
pdf("output/Plots/03_Heatmap_ComplexAmpSample.pdf", width = 18, height = 3)
p1
dev.off()

df2 <- data.frame(AmpNo = colSums(mtx1), category = factor(SampleInfo[colnames(mtx1),]$category, names(ecclass_colors)),
                  CancerType = factor(SampleInfo[colnames(mtx1),]$CancerType2, names(cancertype_colors)))
df2 <- df2[which(df2$category != "No-fSCNA"),]
p2 <- ggplot(df2, aes(x=category, y=AmpNo, fill=category)) + 
  geom_jitter_rast(height = 0, width = 0.2, alpha = 0.5, aes(color = category)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() +
  scale_color_manual(values = ecclass_colors) + guides(colour=F, fill = F) +
  labs(x = "", y = "# Amplicon") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_signif(comparisons = list(c("ecDNA", "ChromAmp")), textsize = 4) + theme(legend.key.size = unit(0.3, 'cm'))
p3 <- ggplot(df2, aes(x=category, y=AmpNo, fill=category)) + 
  geom_jitter_rast(height = 0, width = 0.2, alpha = 0.5, aes(color = CancerType)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() +
  scale_color_manual(values = cancertype_colors) +
  labs(x = "", y = "# Amplicon") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_signif(comparisons = list(c("ecDNA", "ChromAmp")), textsize = 4) + theme(legend.key.size = unit(0.3, 'cm'))

pdf("output/Plots/03_Boxplot_NoAmpInSample_ecStatus.pdf", width = 2, height = 4)
p2
dev.off()
pdf("output/Plots/03_Boxplot_NoAmpInSample_ecStatus_TypeColored.pdf", width = 4, height = 4)
p3
dev.off()

write.csv(SampleInfo, "output/Tables/03_SampleListSummary_Curated_addClass.csv")


df4 <- table(SampleInfo$CancerType2, SampleInfo$category) %>% as.data.frame.matrix()
df4 <- df4[names(cancertype_colors),c("ecDNA", "ChromAmp", "No-fSCNA")]
write.csv(df4, "output/Tables/03_cM_CancerType_ecDNAamp.csv")

df4.2 <- table(SampleInfo$CancerType2, SampleInfo$category) %>% as.data.frame()
colnames(df4.2) <- c("CancerType", "AmpClass", "N")
df4.2$CancerType <- factor(df4.2$CancerType, levels = c(names(cancertype_colors)))
df4.2$AmpClass <- factor(df4.2$AmpClass, levels = c("No-fSCNA", "ChromAmp", "ecDNA"))

p4 <- ggplot(df4.2, aes(x = CancerType, y = N, fill = AmpClass)) + geom_bar(stat = "identity", color = "black", size = 0.25) +
  theme_cb() + scale_y_continuous(expand = c(0,0)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = ecclass_colors) +
  labs(x = "Cancer Type", y = "# Samples")

tmp <- lapply(c("ecDNA", "ChromAmp", "No-fSCNA"), function(x) sum(df4.2[which(df4.2$AmpClass == x), "N"])) %>% unlist()
df5 <- data.frame(CancerType = "All", AmpClass = factor(c("ecDNA", "ChromAmp", "No-fSCNA"), levels = c("ecDNA", "ChromAmp", "No-fSCNA")), N = tmp)

p5 <- ggplot(df5, aes(x="", y=N, fill=AmpClass)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) +
  theme_cb() + scale_fill_manual(values = ecclass_colors)

pdf("output/Plots/03_Barplot_Sample.pdf", width = 8, height = 6)
p4
dev.off()
pdf("output/Plots/03_Piechart_AllSample.pdf", width = 6, height = 5)
p5
dev.off()

# number of amplicon
df6 <- table(AA_result_T$CancerType, AA_result_T$Classification) %>% as.data.frame.matrix()
df6 <- df6[names(cancertype_colors),c("Linear", "Complex-non-cyclic", "BFB", "ecDNA")]
write.csv(df6, "output/Tables/03_cM_CancerType_AmpClass.csv")

df6.2 <- table(AA_result_T$Classification, AA_result_T$CancerType) %>% as.data.frame
colnames(df6.2) <- c("AmpClass", "CancerType", "N")

df6.2$AmpClass <- factor(df6.2$AmpClass, levels = c("Linear", "Complex-non-cyclic", "BFB", "ecDNA"))
df6.2$CancerType <- factor(df6.2$CancerType, levels = c(names(cancertype_colors)))
df_tmp <- data.frame(AmpClass = c("Linear", "Complex-non-cyclic", "BFB", "ecDNA"), CancerType = rep("Kidney",4), N = c(0,0,0,0))
df6.2 <- rbind(df6.2, df_tmp)

p6 <- ggplot(df6.2, aes(x = CancerType, y = N, fill = AmpClass)) + geom_bar(stat = "identity", color = "black", size = 0.25) +
  theme_cb() + scale_y_continuous(expand = c(0,0)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = ampclass_colors) + labs(x = "Cancer Type", y = "# detected amplicons")

tmp <- lapply(c("Linear", "Complex-non-cyclic", "BFB", "ecDNA"), function(x) sum(df6.2[which(df6.2$AmpClass == x), "N"])) %>% unlist()
df7 <- data.frame(CancerType = "All", AmpClass = factor(c("Linear", "Complex-non-cyclic", "BFB", "ecDNA"), levels = c("Linear", "Complex-non-cyclic", "BFB", "ecDNA")), N = tmp)

p7 <- ggplot(df7, aes(x="", y=N, fill=AmpClass)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) +
  theme_cb() + scale_fill_manual(values = ampclass_colors)
pdf("output/Plots/03_Barplot_Amplicon.pdf", width = 8, height = 6)
p6
dev.off()
pdf("output/Plots/03_Piechart_AllAmplicon.pdf", width = 6, height = 5)
p7
dev.off()

#---------- Histogram of Amplified regions ----------#
source("code/func_makeWindows.R")
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
blacklist <- import.bed("ref/hg38-blacklist.v2.bed")

#1MB window(no overlapping)
windows <- makeWindows(genome = BSgenome.Hsapiens.UCSC.hg38, blacklist = blacklist, windowSize = 1e6, slidingSize = 1e6)

#gene overlaps
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gene_gr <- genes(txdb)
names(gene_gr) <- NULL
gene_gr$SYMBOL <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_gr$gene_id, keytype = "ENTREZID", columns = "SYMBOL")[,2]

gene_ls <- parallel::mclapply(seq_along(windows), function(x) paste(gene_gr$SYMBOL[queryHits(findOverlaps(gene_gr, windows[x]))], collapse = ";"), mc.cores = 26)
gene_ls <- unlist(gene_ls)

windows$gene <- gene_ls
saveRDS(windows, "rds/03_windows_1Mb_noOverlapping.rds")

#window overlaps
gr_T <- readRDS("rds/02_gr_T_filtered.rds")
gr_T <- unlist(gr_T)
gr_T$FeatureID2 <- stringr::str_split(gr_T$FeatureID, "_amplicon", simplify = T)[,1]

#calculate frequency of each class across cancertypes
tmp <- AA_result_T$CancerType
names(tmp) <- AA_result_T$Feature.ID
gr_T$CancerType <- tmp[gr_T$FeatureID]
df8 <- lapply(names(cancertype_colors), function(y){
  out <- lapply(c("Linear", "Complex-non-cyclic", "BFB", "ecDNA"), function(x){
    gr <- gr_T[gr_T$Classification == x & gr_T$CancerType == y]
    overlapDF <- DataFrame(findOverlaps(windows, gr, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any"))
    overlapDF$name <- mcols(gr)[overlapDF[, 2], "FeatureID2"]
    overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
    sparseM <- Matrix::sparseMatrix(
      i = overlapTDF[, 1], 
      j = overlapTDF[, 4],
      x = rep(1, nrow(overlapTDF)), 
      dims = c(length(windows), length(unique(overlapDF$name))))
    
    sparseM <- as.matrix(sparseM)
    colnames(sparseM) <- unique(overlapDF$name)
    sparseM[sparseM > 0] <- 1
    sparseM2 <- rowSums(sparseM)
    
    df_out <- data.frame(freq = sparseM2)
    return(df_out)
  })
  out <- do.call(cbind, out)
  colnames(out) <- c("Linear", "Complex-non-cyclic", "BFB", "ecDNA")
  return(out)
})
names(df8) <- names(cancertype_colors)

idx1 <- sapply(paste0("chr", c(1:22, "X", "Y")), function(x) which(windows$wSeq == x)[1])

p8 <- lapply(names(cancertype_colors), function(y){
  max_ind <- max(df8[[y]])
  p <- lapply(c("Linear", "Complex-non-cyclic", "BFB", "ecDNA"), function(x){
    df_tmp <- df8[[y]]
    df <- data.frame(index = c(1:nrow(df_tmp)), freq = df_tmp[,x])
    color_tmp <- c(ampclass_colors)[x]
    
    out <- ggplot(df, aes(x = index, y = freq)) + geom_bar(stat = "identity", color = as.character(color_tmp)) + 
      geom_vline(xintercept = c(idx1, length(windows)), lty = "dotted", color = "black") +
      scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0), limits = c(0,max_ind)) + theme_cb() +
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),axis.title.x = element_blank()) +
      labs(x = "Genomic position", y = "# Samples") + ggtitle(y)
    return(out)
  })
  return(p)
})
pdf("output/Plots/03_Barplot_NumSample_AmpClass.pdf", width = 9, height = 2)
p8
dev.off()

df8.2 <- do.call(cbind, df8)
mcols(windows)[,colnames(df8.2)] <- df8.2
names(windows) <- NULL
write.csv(data.frame(windows), "output/Tables/03_windows_1Mb_noOverlapping_NumSample.csv")

#---------- # Oncogene cargo ----------#
oncogene_Amp <- lapply(c(1:nrow(AA_result_T)), function(i){
  input_string <- AA_result_T[i,]$Oncogenes
  clean_string <- gsub("\\[", "", input_string)
  clean_string <- gsub("\\]", "", clean_string)
  clean_string <- gsub("'", "", clean_string)
  clean_string <- gsub(" ", "", clean_string)
  result <- unlist(strsplit(clean_string, ","))
  return(result)
})
oncogene_Amp_uq <- sort(unique(unlist(oncogene_Amp)))

oncogene_Amp_df <- lapply(oncogene_Amp_uq, function(x){
  out <- lapply(oncogene_Amp, function(ls){
    length(which(x %in% ls))
  }) %>% unlist
  return(out)
})
oncogene_Amp_df <- do.call(cbind, oncogene_Amp_df) %>% as.data.frame()
colnames(oncogene_Amp_df) <- oncogene_Amp_uq
rownames(oncogene_Amp_df) <- AA_result_T$Feature.ID
oncogene_Amp_df$AmpClass <- AA_result_T$Classification
oncogene_Amp_df$CancerType <- AA_result_T$CancerType
write.csv(oncogene_Amp_df, "output/Tables/03_Number_oncogene_amplicon.csv")

df9 <- reshape2::melt(oncogene_Amp_df)
p9 <- lapply(names(cancertype_colors), function(x){
  df_tmp <- df9[which(df9$CancerType == x),]
  idx <- colSums(oncogene_Amp_df[which(oncogene_Amp_df$CancerType == x & oncogene_Amp_df$AmpClass == "ecDNA"),c(1:522)]) %>% sort(., decreasing = T) %>% names() %>% .[c(1:30)]
  df_tmp <- df_tmp[df_tmp$variable %in% idx,]
  df_tmp$variable <- factor(df_tmp$variable, levels = idx)
  df_tmp$AmpClass <- factor(df_tmp$AmpClass, levels = c("Linear", "Complex-non-cyclic", "BFB", "ecDNA"))
  out <-  ggplot(df_tmp, aes(x = variable, y = value, fill = AmpClass)) +
    geom_bar(stat = "identity") +
    theme_cb() + scale_y_continuous(expand = c(0,0)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = ampclass_colors) +
    labs(x = "Oncogene", y = "# amplicon") + ggtitle(x)
  return(out)
})
pdf("output/Plots/03_Barplot_Oncogene_AmpClass_CancerType.pdf", width = 8, height = 3)
p9
dev.off()

#---------- Amplicon CN, length, # genes ----------#
AA_result_T_geneNo <- lapply(c(1:nrow(AA_result_T)), function(i){
  input_string <- AA_result_T[i,]$All.genes
  clean_string <- gsub("\\[", "", input_string)
  clean_string <- gsub("\\]", "", clean_string)
  clean_string <- gsub("'", "", clean_string)
  clean_string <- gsub(" ", "", clean_string)
  result <- unlist(strsplit(clean_string, ","))
  result <- length(result)
  return(result)
}) %>% unlist

df10 <- data.frame(amplicon = AA_result_T$Feature.ID,
                  class = factor(AA_result_T$Classification, levels = rev(c("Linear", "Complex-non-cyclic", "BFB", "ecDNA"))),
                  CN = AA_result_T$Feature.maximum.copy.number,
                  length = AA_result_T$Captured.interval.length,
                  NumGenes = AA_result_T_geneNo,
                  CancerType = factor(AA_result_T$CancerType, levels = names(cancertype_colors)))

p10.1 <- ggplot(df10, aes(x = class, y = CN)) + geom_jitter(height = 0, width = 0.25, shape = 3, alpha = 0.5, aes(color = class)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() + 
  scale_color_manual(values = ampclass_colors) + labs(x = "Amplicon Classification", y = "WGS-derived CN (Maximum segment)") + 
  scale_y_continuous(trans = "log2", breaks = c(0,2,4,8,16,32,64,128,256)) + 
  geom_signif(comparisons = list(c("ecDNA", "BFB"), c("ecDNA", "Complex-non-cyclic"), c("ecDNA", "Linear")), test = "wilcox.test", y_position = c(8,8.3,8.6))
p10.2 <- ggplot(df10, aes(x = class, y = length)) + geom_jitter(height = 0, width = 0.25, shape = 3, alpha = 0.5, aes(color = class)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() + scale_color_manual(values = ampclass_colors)  +
  scale_y_log10(breaks=10^(4:8), labels=trans_format("log10",math_format(10^.x))) +
  labs(x = "Amplicon Classification", y = "Amplicon length (bp)") +
  geom_signif(comparisons = list(c("ecDNA", "BFB"), c("ecDNA", "Complex-non-cyclic"), c("ecDNA", "Linear")), test = "wilcox.test", y_position = c(7.8,8,8.2))
p10.3 <- ggplot(df10, aes(x = class, y = NumGenes)) + geom_jitter(height = 0, width = 0.25, shape = 3, alpha = 0.5, aes(color = class)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() + scale_color_manual(values = ampclass_colors)  +
  labs(x = "Amplicon Classification", y = "# genes on amplicons") +
  geom_signif(comparisons = list(c("ecDNA", "BFB"), c("ecDNA", "Complex-non-cyclic"), c("ecDNA", "Linear")), test = "wilcox.test", y_position = c(380,400,420))

pdf("output/Plots/03_Boxplot_Properties_AmpliconClass.pdf", width = 6, height = 5)
p10.1
p10.2
p10.3
dev.off()

p10.4 <- ggplot(df10, aes(x = class, y = CN)) + geom_jitter(height = 0, width = 0.25, shape = 3, alpha = 0.5, aes(color = class)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values = ampclass_colors) + labs(x = "Amplicon Classification", y = "WGS-derived CN (Maximum segment)") + 
  scale_y_continuous(trans = "log2", breaks = c(0,2,4,8,16,32,64,128,256)) + 
  facet_grid(~CancerType)
p10.5 <- ggplot(df10, aes(x = class, y = length)) + geom_jitter(height = 0, width = 0.25, shape = 3, alpha = 0.5, aes(color = class)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() + 
  scale_y_log10(breaks=10^(4:8), labels=trans_format("log10",math_format(10^.x))) + scale_color_manual(values = ampclass_colors) +
  labs(x = "Amplicon Classification", y = "Amplicon length (bp)") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(~CancerType)
p10.6 <- ggplot(df10, aes(x = class, y = NumGenes)) + geom_jitter(height = 0, width = 0.25, shape = 3, alpha = 0.5, aes(color = class)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() + scale_color_manual(values = ampclass_colors)  +
  labs(x = "Amplicon Classification", y = "# genes on amplicons") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(~CancerType)

pdf("output/Plots/03_Boxplot_Properties_AmpliconClass_CancerType.pdf", width = 15, height = 5)
p10.4
p10.5
p10.6
dev.off()

df11 <- df10[which(df10$class == "ecDNA"),]
p11.1 <- ggplot(df11, aes(x = CancerType, y = CN)) + geom_jitter(height = 0, width = 0.25, shape = 3, alpha = 0.5, aes(color = class)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values = ampclass_colors) + labs(x = "Cancer Type", y = "WGS-derived CN (Maximum segment)") + 
  scale_y_continuous(trans = "log2", breaks = c(0,2,4,8,16,32,64,128,256)) +
  geom_signif(comparisons = list(c("Headneck", "Sarcoma")), test = "wilcox.test")
p11.2 <- ggplot(df11, aes(x = CancerType, y = length)) + geom_jitter(height = 0, width = 0.25, shape = 3, alpha = 0.5, aes(color = class)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() + scale_color_manual(values = ampclass_colors) +
  scale_y_log10(breaks=10^(4:8), labels=trans_format("log10",math_format(10^.x))) +
  labs(x = "Amplicon Classification", y = "Amplicon length (bp)") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_signif(comparisons = list(c("Headneck", "Sarcoma")), test = "wilcox.test")
p11.3 <- ggplot(df11, aes(x = CancerType, y = NumGenes)) + geom_jitter(height = 0, width = 0.25, shape = 3, alpha = 0.5, aes(color = class)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() + scale_color_manual(values = ampclass_colors)  +
  labs(x = "Amplicon Classification", y = "# genes on amplicons") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_signif(comparisons = list(c("Headneck", "Sarcoma")), test = "wilcox.test")
pdf("output/Plots/03_Boxplot_Properties_ecDNA_CancerType.pdf", width = 6, height = 4)
p11.1
p11.2
p11.3
dev.off()

df12 <- lapply(names(cancertype_colors), function(x){
  out <- lapply(c("Linear", "Complex-non-cyclic", "BFB", "ecDNA"), function(y){
    res <- data.frame(CN = mean(df10$CN[which(df10$CancerType==x & df10$class==y)]),
                      length = mean(df10$length[which(df10$CancerType==x & df10$class==y)]),
                      NumGenes = mean(df10$NumGenes[which(df10$CancerType==x & df10$class==y)]),
                      class = y, CancerType = x)
    return(res)
  })
  out <- do.call(rbind, out)
  return(out)
})
df12 <- do.call(rbind, df12)

df12$class <- factor(df12$class, levels = rev(c("Linear", "Complex-non-cyclic", "BFB", "ecDNA")))
df12$CancerType <- factor(df12$CancerType, levels = names(cancertype_colors))

p12.1 <- ggplot(df12, aes(x = class, y = CN)) + geom_jitter(height = 0, width = 0.2, alpha = 0.8, size = 4, aes(color = CancerType)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() + 
  scale_color_manual(values = cancertype_colors) + labs(x = "Amplicon Classification", y = "WGS-derived CN (Maximum segment)") + 
  scale_y_continuous(trans = "log2", breaks = c(0,2,4,8,16,32))
p12.2 <- ggplot(df12, aes(x = class, y = length)) + geom_jitter(height = 0, width = 0.2, alpha = 0.8, size = 4, aes(color = CancerType)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() + 
  scale_color_manual(values = cancertype_colors)  +
  scale_y_log10(breaks=10^(4:8), labels=trans_format("log10",math_format(10^.x))) +
  labs(x = "Amplicon Classification", y = "Amplicon length (bp)")
p12.3 <- ggplot(df12, aes(x = class, y = NumGenes)) + geom_jitter(height = 0, width = 0.2, alpha = 0.8, size = 4, aes(color = CancerType)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() + 
  scale_color_manual(values = cancertype_colors)  +
  labs(x = "Amplicon Classification", y = "# genes on amplicons")

pdf("output/Plots/03_Boxplot_Properties_CancerType_Mean.pdf", width = 4, height = 5)
p12.1
p12.2
p12.3
dev.off()
