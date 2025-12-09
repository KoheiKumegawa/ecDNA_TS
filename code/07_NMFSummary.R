#-------------------------------------------------------------------------------
# 07_NMFSummary.R
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

#------------ NMF genes ------------#
nmf_result <- lapply(c(2:20), function(i) readRDS(paste0("/mnt/host_mnt/Volumes/Shared2/Kume/Analysis/20251007_ecDNA_nmf/nmf_result_k", i, ".rds")))
names(nmf_result) <- paste0("k", c(2:20))
nrow(basis(nmf_result[[1]]))

#Num of genes :4842, top5% 242 genes
GEP_matrix <- sapply(paste0("k", c(9:18)), function(x){
  W <- basis(nmf_result[[x]])
  res <- apply(W, 2, function(n) rownames(W)[order(n, decreasing = T)[c(1:242)]])
  res <- apply(res, 2, function(x) rowData(se)[x,]$gene_name)
  colnames(res) <- paste0(x, "#", c(1:ncol(res)))
  return(res)
})
GEP_matrix <- do.call(cbind, GEP_matrix)
nrow(GEP_matrix)
ncol(GEP_matrix)

#jaccard
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
GEP_jaccard <- apply(GEP_matrix,2,function(x) apply(GEP_matrix,2,function(y) jaccard(x,y)))

col_fun1 <- colorRamp2(seq(0,1,0.1), c("white", rev(viridis(10, option = "A"))))
fh <- function(x) hclust(dist(x), method = "ward.D2")
ht1 <- Heatmap(GEP_jaccard, name = "Jaccard index", col = col_fun1, 
               cluster_rows = fh, cluster_columns = fh, 
               #bottom_annotation = ha1, right_annotation = ra1,
               row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6),
               use_raster = T)
p1 <- draw(ht1)
pdf("output/Plots/07_Heatmap_Jaccard_cNMF_prefilter.pdf", width = 10, height = 9)
p1
dev.off()
GEP_matrix_JIsort <- GEP_matrix[,column_order(p1)]
write.csv(GEP_matrix_JIsort, "output/Tables/07_GEP_JIsorted.csv")

validGEP <- apply(GEP_jaccard, 2, function(x) length(which(x > 0.66))) %>% sort
validGEP <- names(validGEP)[validGEP >= 5]

GEP_jaccard_valid <- GEP_jaccard[validGEP,validGEP]
ht2 <- Heatmap(GEP_jaccard_valid, name = "Jaccard index", col = col_fun1, 
               cluster_rows = fh, cluster_columns = fh, 
               #bottom_annotation = ha1, right_annotation = ra1,
               row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6),
               use_raster = T)
p2 <- draw(ht2)
pdf("output/Plots/07_Heatmap_Jaccard_cNMF_valid.pdf", width = 10, height = 9)
p2
dev.off()
GEP_matrix_JIsort2 <- GEP_matrix[,colnames(GEP_jaccard_valid)[column_order(p2)]]
write.csv(GEP_matrix_JIsort2, "output/Tables/07_GEP_JIsorted_valid.csv")

#GEP program
GEP_ls <- list(c(1:10),c(11:20),c(21:29),c(30:37),c(38:44),c(45:54),c(55:64),c(65:72),c(73:78),c(79:83),c(84:88),c(89:94),c(95:103),c(104:111))
GEP_genes <- lapply(GEP_ls, function(x){
  out <- as.character(GEP_matrix_JIsort2[,x])
  out <- unique(out) %>% sort
  return(out)
})
names(GEP_genes) <- paste0("GEP",c(1:14))
lapply(names(GEP_genes), function(x) write.table(GEP_genes[[x]], paste0("output/Tables/07_geneList_", x, ".txt"), row.names = F, col.names = F, quote = F))

#------------ ModuleScore ------------#
GEP_genes <- lapply(paste0("GEP",c(1:14)), function(x) read.table(paste0("output/Tables/07_geneList_",x,".txt"))[,1])
names(GEP_genes) <- paste0("GEP",c(1:14))

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

mat <-  read.csv("output/Tables/RNMF_InputMatrix.csv", row.names = 1)
rownames(mat) <- rowData(se)[rownames(mat),]$gene_name
GEP_moduleScore <- addGeneModuleScore(expr_matrix = mat, gene_list = GEP_genes)
GEP_moduleScore$CancerType <- se$CancerType2
GEP_moduleScore$ecDNA <- se$category1

##Difference from average
GEP_average <- sapply(GEP_genes, function(x) colMeans(mat[x,]))
p3 <- lapply(names(GEP_genes), function(x){
  df_tmp <- data.frame(Average = GEP_average[,x], ModuleScore = GEP_moduleScore[,x], CancerType = GEP_moduleScore$CancerType)
  out <- ggplot(df_tmp, aes(x = Average, y = ModuleScore, color = CancerType)) + 
    geom_point(alpha = 0.6, size = 2)  + 
    theme_cb() + scale_color_manual(values = cancertype_colors) + 
    geom_smooth(aes(x = Average, y = ModuleScore, color = NULL), method = "lm", se = F, color = "black") +
    theme(legend.key.size = unit(0.2, 'cm')) + ggtitle(x)
})
pdf("output/Plots/07_Scatterplot_Avrg_vs_Module.pdf", width = 4, height = 2.5)
p3
dev.off()

###Bubble plot###
GEP_SampleProp <- lapply(paste0("GEP", c(1:14)), function(x){
  tmp <- sapply(names(cancertype_colors), function(y) mean(GEP_moduleScore[,x][which(GEP_moduleScore$CancerType == y)])) %>% sort(.,decreasing = T)
  #cutoff <- median(GEP_moduleScore[,x])
  cutoff <- GEP_moduleScore[,x][which(GEP_moduleScore$CancerType == names(tmp)[1])] %>% quantile(., 0.25)

  tb <- table(GEP_moduleScore$CancerType[which(GEP_moduleScore[,x] > cutoff)])
  ntb <- names(tb)
  tb <- as.numeric(tb)
  names(tb) <- ntb
  out <- tb[names(cancertype_colors)]
  return(out)
})
GEP_SampleProp <- do.call(rbind, GEP_SampleProp)
GEP_SampleProp[is.na(GEP_SampleProp) == T] <- 0
rownames(GEP_SampleProp) <- paste0("GEP", c(1:14))
colnames(GEP_SampleProp) <- names(cancertype_colors)
GEP_SampleProp_Percent <- t(GEP_SampleProp) / as.numeric(table(se$CancerType2)[names(cancertype_colors)])

#calculate gini index
gini_index <- apply(t(GEP_SampleProp_Percent), 1, function(x) DescTools::Gini(x))
gini_index <- sort(gini_index)

#melt
GEP_SampleProp_Percent <- reshape2::melt(GEP_SampleProp_Percent)

#calculate mean usage per each GEP per each Sample
GEP_SampleMean <- lapply(names(cancertype_colors), function(x) colMeans(GEP_moduleScore[GEP_moduleScore$CancerType == x,c(1:14)]))
GEP_SampleMean <- do.call(rbind, GEP_SampleMean)
rownames(GEP_SampleMean) <- names(cancertype_colors)
colnames(GEP_SampleMean) <- paste0("GEP", c(1:14))
GEP_SampleMean <- scale(GEP_SampleMean)
GEP_SampleMean <- reshape2::melt(GEP_SampleMean)

#make visualize DF
GEP_UsageVis <- data.frame(CancerType = GEP_SampleMean$Var1, GEP = GEP_SampleMean$Var2, Mean = GEP_SampleMean$value, Prop = GEP_SampleProp_Percent$value)
GEP_UsageVis$GEP <- factor(GEP_UsageVis$GEP, levels = names(gini_index))

p4 <- ggplot(GEP_UsageVis, aes(x = CancerType, y = GEP, color = Mean, size = ifelse(Prop==0, NA, Prop))) + 
  scale_color_gradient2(low="#313695",mid="#FFFFBF",high="#A50026",midpoint = 0) +
  geom_point() + geom_point(aes(size=ifelse(Prop==0, NA, Prop)), shape = 21, colour="black", stroke=0.25) +
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color="black"), axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.key.size = unit(0.4, 'cm'))

col_fun1 <- colorRamp2(seq(0.2,1,0.1), c("white", rev(viridis(8, option = "A"))))
ht1 <- Heatmap(gini_index, name = "Gini Index", col = col_fun1, cluster_rows = F, cluster_columns = F) 
p5 <- draw(ht1)

pdf("output/Plots/07_BubblePlot_GEPScore.pdf", width = 7, height = 4)
p4
dev.off()
pdf("output/Plots/07_Heatmap_GiniIndex_GEPScore.pdf", width = 2, height = 3)
p5
dev.off()

GEP_moduleScore2 <- reshape2::melt(GEP_moduleScore)
GEP_moduleScore2$CancerType <- factor(GEP_moduleScore2$CancerType, levels = rev(names(cancertype_colors)))
GEP_moduleScore2$variable <- factor(GEP_moduleScore2$variable, levels = names(gini_index))

cutoff1 <- sapply(paste0("GEP", c(1:14)), function(x){
  tmp <- sapply(names(cancertype_colors), function(y) mean(GEP_moduleScore[,x][which(GEP_moduleScore$CancerType == y)])) %>% sort(.,decreasing = T)
  out <- GEP_moduleScore[,x][which(GEP_moduleScore$CancerType == names(tmp)[1])] %>% quantile(., 0.25)
  return(out)
})
df_25percentile <- data.frame(
  variable = factor(paste0("GEP", c(1:14)), levels = names(gini_index)),
  cutoff = as.numeric(cutoff1)
)

library(ggridges)
p6 <- ggplot(GEP_moduleScore2, aes(y=CancerType, x=value, fill = CancerType)) + 
  geom_density_ridges2() + 
  geom_vline(data = df_25percentile, aes(xintercept = cutoff, colour = "red")) + 
  scale_fill_manual(values = cancertype_colors) + facet_grid(~variable) + 
  theme(legend.position="none", axis.text.x = element_text(angle = 90, hjust = 1))
pdf("output/Plots/07_RidgePlot_GEPScore.pdf", width = 12, height = 3)
p6
dev.off()

#------------ MCPcounter ------------#
library(MCPcounter)
se_pc <- readRDS("rds/06_se_pc.rds")

norm_counts2 <- assays(se_pc)$log2tpm
rownames(norm_counts2) <- rowData(se_pc)$gene_name
MCPcounter_output <- MCPcounter.estimate(norm_counts2, featuresType = "HUGO_symbols")
MCPcounter_output <- as.data.frame(t(MCPcounter_output))

df_corTME <- cor(GEP_moduleScore[c(1:14)], MCPcounter_output, method = "pearson")
df_corTME <- reshape2::melt(df_corTME)

p7 <- ggplot(df_corTME, aes(x = factor(Var1, levels = names(gini_index)), y = value, color = Var2)) + geom_point(shape = 1, size = 2.5) + theme_cb() +
  labs(x = "GEP", y = "Pearson's correlation") + geom_hline(yintercept = 0, lty = "dotted") + ylim(-0.56,1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values = c("#FFB300","#803E75","#FF6800","#A6BDD7","#C10020","#CEA262","#817066","#007D34","#F6768E","#00538A"))
pdf("output/Plots/07_Dotplot_NMFscore_vs_TMEscore.pdf", width = 7, height = 3)
p7
dev.off()

p8 <- lapply(names(GEP_moduleScore), function(x){
  out <- lapply(colnames(MCPcounter_output), function(y){
    df7 <- data.frame(MCPscore = MCPcounter_output[,y], GEPscore = GEP_moduleScore[[x]], CancerType = GEP_moduleScore$CancerType)
    p <- ggplot(df7, aes(x = GEPscore, y = MCPscore, color = CancerType)) + geom_point(alpha = 0.6, size = 2)  + 
      theme_cb() + scale_color_manual(values = cancertype_colors) + 
      geom_smooth(aes(x = GEPscore, y = MCPscore, color = NULL), method = "lm", se = F, color = "black") +
      theme(legend.key.size = unit(0.2, 'cm')) + labs(x = paste0(x, " module score"), y = paste0(y, " score"))
    return(p)
  })
  return(out)
})
pdf("output/Plots/07_Scatterplot_NMFscore_TMEscore.pdf", width = 4, height = 2.5)
p8
dev.off()

#------------ GEP score visualization across cancer types ------------#
GEP_moduleScore
p9 <- lapply(paste0("GEP", c(1:14)), function(x){
  mtx_tmp <- data.frame(ModuleScore = GEP_moduleScore[,x])
  mtx_tmp$ecDNA <- factor(GEP_moduleScore$ecDNA, levels = c("ecDNA_n", "ecDNA_p"))
  out <- ggplot(mtx_tmp, aes(x=ecDNA, y=ModuleScore, fill=ecDNA)) + 
    geom_jitter_rast(height = 0, width = 0.2, alpha = 0.5, aes(color = ecDNA)) + 
    geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() +
    scale_color_manual(values = ecclass_colors) +
    labs(x = "", y = paste0("Module score")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    geom_signif(comparisons = list(c("ecDNA_n", "ecDNA_p")), textsize = 2) + ggtitle(x) + theme(legend.key.size = unit(0.3, 'cm'))
  return(out)
})
pdf("output/Plots/07_Boxplot_ModuleScore_ecDNA_ecDNAColored.pdf", width = 2.5, height = 3)
p9
dev.off()

p10 <- lapply(paste0("GEP", c(1:14)), function(x){
  mtx_tmp <- data.frame(ModuleScore = GEP_moduleScore[,x])
  mtx_tmp$ecDNA <- factor(GEP_moduleScore$ecDNA, levels = c("ecDNA_n", "ecDNA_p"))
  mtx_tmp$CancerType <- GEP_moduleScore$CancerType
  out <- ggplot(mtx_tmp, aes(x=ecDNA, y=ModuleScore, fill=ecDNA)) + 
    geom_jitter_rast(height = 0, width = 0.2, alpha = 0.5, aes(color = CancerType)) + 
    geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() +
    scale_color_manual(values = cancertype_colors) + guides(colour=F, fill = F) +
    labs(x = "", y = paste0("Module score")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    geom_signif(comparisons = list(c("ecDNA_n", "ecDNA_p")), textsize = 2) + ggtitle(x) + theme(legend.key.size = unit(0.3, 'cm'))
  return(out)
})
pdf("output/Plots/07_Boxplot_ModuleScore_ecDNA_CancertypeColored.pdf", width = 2, height = 3)
p10
dev.off()

p11 <- lapply(paste0("GEP", c(1:14)), function(x){
  mtx_tmp <- data.frame(ModuleScore = GEP_moduleScore[,x])
  mtx_tmp$ecDNA <- factor(GEP_moduleScore$ecDNA, levels = c("ecDNA_n", "ecDNA_p"))
  mtx_tmp$CancerType <- factor(GEP_moduleScore$CancerType, levels = names(cancertype_colors))
  out <- ggplot(mtx_tmp, aes(x=CancerType, y=ModuleScore, fill=ecDNA)) + 
    geom_jitter_rast(alpha = 0.5, size = 0.5, aes(color = ecDNA),
                     position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0)) + 
    geom_boxplot(outlier.shape = NA, alpha = 0) + theme_cb() +
    labs(x = "", y = paste0("Module score")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    scale_color_manual(values = ecclass_colors) + ggtitle(x) + theme(legend.key.size = unit(0.3, 'cm'))
  return(out)
})
pdf("output/Plots/07_Boxplot_ModuleScore_ecDNA_CancertypeFilled.pdf", width = 8, height = 3)
p11
dev.off()

GEP_moduleScore_corrected <- GEP_moduleScore[,c(1:14)]
GEP_moduleScore_corrected <- sapply(colnames(GEP_moduleScore_corrected), function(x){
  tmp <- GEP_moduleScore_corrected[,x]
  if(length(which(tmp < 0)) > 0){
    tmp <- tmp - min(tmp)
    return(tmp)
  } else {
    return(tmp)
  }
}) %>% as.data.frame
rownames(GEP_moduleScore_corrected) <- rownames(GEP_moduleScore)
GEP_moduleScore_corrected$CancerType <- GEP_moduleScore$CancerType
GEP_moduleScore_corrected$ecDNA <- GEP_moduleScore$ecDNA

GEP_VSecDNA <- lapply(paste0("GEP",c(1:14)), function(x){
  res <- lapply(names(cancertype_colors), function(y){
    value1 <- median(GEP_moduleScore_corrected[which(GEP_moduleScore_corrected$CancerType == y & GEP_moduleScore_corrected$ecDNA == "ecDNA_p"),x])
    value2 <- median(GEP_moduleScore_corrected[which(GEP_moduleScore_corrected$CancerType == y & GEP_moduleScore_corrected$ecDNA == "ecDNA_n"),x])
    out <- data.frame(ecDNA_p = value1, ecDNA_n = value2)
    return(out)
  }) %>% do.call(rbind,.)
  rownames(res) <- names(cancertype_colors)
  res <- res[-16,]
  res$FC <- res$ecDNA_p / res$ecDNA_n
  return(res$FC)
}) %>% do.call(cbind,.)
rownames(GEP_VSecDNA) <- names(cancertype_colors)[-16]
colnames(GEP_VSecDNA) <- paste0("GEP",c(1:14))

col_fun2 <- colorRamp2(c(0.5,1,2), c("#313695","#FFFFBF","#A50026"))
ht3 <- Heatmap(GEP_VSecDNA, name = "FC [median score ecDNA(+)/ecDNA(-)]",col = col_fun2, cluster_columns = fh, cluster_rows = F)
p12 <- draw(ht3)
pdf("output/Plots/07_Heatmap_FC_medianModuleScore_ecDNApn.pdf", width = 6.5, height = 4)
p12
dev.off()

GEP_VSecDNA_all <- lapply(paste0("GEP",c(1:14)), function(x){
  value1 <- median(GEP_moduleScore_corrected[which(GEP_moduleScore_corrected$ecDNA == "ecDNA_p"),x])
  value2 <- median(GEP_moduleScore_corrected[which(GEP_moduleScore_corrected$ecDNA == "ecDNA_n"),x])
  res <- value1 / value2
  return(res)
}) %>% do.call(cbind,.)
colnames(GEP_VSecDNA_all) <- paste0("GEP",c(1:14))

ht4 <- Heatmap(GEP_VSecDNA_all[,column_order(p12)], name = "FC [median score ecDNA(+)/ecDNA(-)]",col = col_fun2, cluster_columns = fh, cluster_rows = F)
p13 <- draw(ht4)
pdf("output/Plots/07_Heatmap_FC_medianModuleScore_ecDNApn_all.pdf", width = 4, height = 5)
p13
dev.off()

#wilcox test
GEP_ecDNA_wilcoxp <- sapply(paste0("GEP", c(1:14)), function(x){
  mtx_tmp <- data.frame(ModuleScore = GEP_moduleScore[,x])
  mtx_tmp$ecDNA <- factor(GEP_moduleScore$ecDNA, levels = c("ecDNA_n", "ecDNA_p"))
  res <- wilcox.test(mtx_tmp$ModuleScore~mtx_tmp$ecDNA)
  return(res$p.value)
})
GEP_ecDNA_wilcoxp <- -log10(GEP_ecDNA_wilcoxp)

col_fun3 <- colorRamp2(c(0,5,10,15,20), c(viridis(5, option = "D")))
ht5 <- Heatmap(GEP_ecDNA_wilcoxp[column_order(p12)], name = "-log10(P-value)",col = col_fun3, cluster_columns = F, cluster_rows = F)
p14 <- draw(ht5)
pdf("output/Plots/07_Heatmap_log2Pvalue_DiffecDNApn.pdf", width = 2, height = 4)
p14
dev.off()
