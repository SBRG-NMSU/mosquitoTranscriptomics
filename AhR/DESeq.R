############ Prereqs ############
options(stringsAsFactors = FALSE, scipen = 900)
oldPar <- par()
os <- Sys.info()["sysname"]
baseDir <- ifelse(os == "Windows", "C:/Users/ptrainor/gdrive/Mosquitos/AhR/", 
                  "~/gdrive/Mosquitos/AhR/")
setwd(baseDir)

library(tidyverse)
library(DESeq2)
library(ggrepel)

theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5))

############ Import and process raw counts ############
# Count data:
df1 <- read.csv("Data/AhR_RNAseq_RawCounts_20200802.csv", header = TRUE, check.names = FALSE)
names(df1) <- gsub(" - total_read_count", "", names(df1))
names(table(df1$Name))[table(df1$Name) > 1]
rmExtra1 <- which(df1$Name == "tRNA-Leu")[2]
df1 <- df1[-rmExtra1,]
# Get rid of extra part in name:
df1$Name <- gsub("AgaP_", "", df1$Name)
rownames(df1) <- df1$Name
df1$Name <- NULL

# Data with severe outlier removed:
df1b <- df1[, !(names(df1) ==  "dsAhR + Serratia Infection 2")]

# Column (phenotype) data:
colData1 <- data.frame(oldSampName = names(df1))
colData1$sampName <- gsub(" ", "_", gsub(" \\+ ", "_p_", colData1$oldSampName))
colData1$pheno <- substr(colData1$sampName, start = 1, stop = nchar(colData1$sampName) - 2)
colData1$rep <- substr(colData1$sampName, start = nchar(colData1$sampName), stop = nchar(colData1$sampName))
colData1$pheno <- as.factor(colData1$pheno)

# Phenotype data w/o outlier:
colData1b <- colData1 %>% filter(oldSampName != "dsAhR + Serratia Infection 2")

# Round fractional counts:
df1 <- round(df1) # 12,723 transcripts
df1b <- round(df1b) 

# Filter out transcripts w/ 0 counts:
filter1 <- apply(df1, 1, function(x) sum(x == 0) == 18)
filter1b <- apply(df1b, 1, function(x) sum(x == 0) == 17) 
df1 <- df1[!filter1, ] # 11,721
df1b <- df1b[!filter1b, ] # 11,719

############ Import annotation data ############
# Imports:
biomartDF <- read.csv("../Priming_Infection/biomart_agambiae_data.csv")
ensemblIDDF <- read.csv("../Priming_Infection/contrastDF2_ensembl_ids2.csv") %>% select(gene, ensembl_id) %>% unique()
# keggDF <- read.csv("./contrast_data_ensembl_kegg.csv") %>% select(gene, ensembl_id, kegg_id) %>% unique()
load("../Priming_Infection/kegg_pathway_data.Rdata")

# Join GO data to ensemble:
ensemblIDDF$GO_terms <- ""
for(i in 1:nrow(ensemblIDDF)){
  if(!is.na(ensemblIDDF$ensembl_id[i])){
    temp1 <- biomartDF[biomartDF$ensembl_gene_id == ensemblIDDF$ensembl_id[i],]
    temp2 <- temp1 %>% select(go_id, name_1006) %>% unique()
    ensemblIDDF$GO_terms[i] <- paste(paste(temp2$go_id, temp2$name_1006, sep = ":"), collapse = "; ")
    print(i)
  }
}
ensemblIDDF$GO_terms[ensemblIDDF$GO_terms == ":"] <- ""

# Process KEGG data:
pathway_list2 <- list()
for(i in 1:length(pathway_list)){
  pathway_list2[[i]] <- data.frame(gene = pathway_list[[i]], pathway = names(pathway_list)[i])
}
pathway_list2 <- do.call("rbind", pathway_list2)

# Join KEGG data:
ensemblIDDF$KEGGpaths1 <- ensemblIDDF$KEGGpaths2 <- ensemblIDDF$KEGGpaths <- ""
for(i in 1:nrow(ensemblIDDF)){
  # Match to ensembl
  if(!is.na(ensemblIDDF$ensembl_id[i])){
    temp1 <- pathway_list2[pathway_list2$gene == ensemblIDDF$ensembl_id[i],]
    temp1$pathway <- gsub("; ", ": ", temp1$pathway)
    ensemblIDDF$KEGGpaths1[i] <- paste(temp1$pathway, collapse = "; ")
  }

  # Match to gene name: 
  if(!is.na(ensemblIDDF$gene[i])){
    temp2 <- pathway_list2[pathway_list2$gene == ensemblIDDF$gene[i],]
    temp2$pathway <- gsub("; ", ": ", temp2$pathway)
    ensemblIDDF$KEGGpaths2[i] <- paste(temp2$pathway, collapse = "; ")
  }
  
  # Merge and check:
  if(ensemblIDDF$KEGGpaths1[i] == ensemblIDDF$KEGGpaths2[i]){
    ensemblIDDF$KEGGpaths[i] <- ensemblIDDF$KEGGpaths1[i]
  }else{
    if(ensemblIDDF$KEGGpaths1[i] == ""){
      ensemblIDDF$KEGGpaths[i] <- ensemblIDDF$KEGGpaths2[i]
    }else{
      ensemblIDDF$KEGGpaths[i] <- ensemblIDDF$KEGGpaths1[i]
    }
  }
    
  print(i)
}
table(ensemblIDDF$KEGGpaths1 != "", ensemblIDDF$KEGGpaths2 != "")
table(ensemblIDDF$KEGGpaths != "")

# Uniqueness:
keggGOAnno <- ensemblIDDF %>% select(gene, KEGGpaths, GO_terms) %>% unique()

############ Univariate analysis ############
# Set up DESeq2 analysis:
des1 <- DESeqDataSetFromMatrix(countData = df1, colData = colData1, design = ~ pheno)
des1 <- DESeq(des1)
des1b <- DESeqDataSetFromMatrix(countData = df1b, colData = colData1b, design = ~ pheno)
des1b <- DESeq(des1b) 

# Export:
rlog1b <- rlog(des1b, blind = TRUE)
rlog1b <- assay(rlog1b)
write.csv(as.data.frame(rlog1b), file = paste0("Results/rlogValues_", gsub("-", "", Sys.Date()), ".csv"))
rlog1b <- t(rlog1b)

# LRT:
des0 <- DESeq(des1b, test="LRT", reduced = ~ 1)
res0 <- as.data.frame(results(des0, cooksCutoff = Inf))

# Pre-filter low normalized counts:
lowNCounts <- rownames(res0)[res0$baseMean < 5]
des1 <- DESeqDataSetFromMatrix(countData = df1[!rownames(df1) %in% lowNCounts,], 
                               colData = colData1, design = ~ pheno) # 10,397
des1b <- DESeqDataSetFromMatrix(countData = df1b[!rownames(df1b) %in% lowNCounts,], 
                                colData = colData1b, design = ~ pheno) # 10,395
des1 <- DESeq(des1)
des1b <- DESeq(des1b) 

# Injury v. Serratia infection: infection altered genes
# AhR antagonist v. Serration infection: AhR dependent genes
# dsAhR v. dsGFP: AhR dependent genes
# dsTIEG v. dsGFP: TIEG dependent genes
# dsAhR v. dsTIEG: genes are affected by both AhR and TIEG.

# Make data.frame of comparisons between priming and Naive:
contrastDF1 <- data.frame(pheno1 = c("Serratia_Infection", "AhR_antagonist_p_Serratia_Infection", 
                                     "dsAhR_p_Serratia_Infection", "dsTieg_p_Serratia_Infection",
                                     "dsAhR_p_Serratia_Infection"),
                          pheno2 = c("Injury_Control", "Serratia_Infection", 
                                     "dsGFP_p_Serratia_Infection", "dsGFP_p_Serratia_Infection",
                                     "dsTieg_p_Serratia_Infection"))

# Make comparisons:
contrastDF2 <- list()
for(i in 1:nrow(contrastDF1)){
  temp1 <- as.data.frame(results(des1b, contrast = c("pheno", contrastDF1$pheno1[i], contrastDF1$pheno2[i]), cooksCutoff = Inf))
  temp1$lfcLower <- temp1$log2FoldChange - qnorm(0.975) * temp1$lfcSE
  temp1$lfcUpper <- temp1$log2FoldChange + qnorm(0.975) * temp1$lfcSE
  temp1$pheno1 <- contrastDF1$pheno1[i]
  temp1$pheno2 <- contrastDF1$pheno2[i]
  temp1$gene <- rownames(temp1)
  contrastDF2[[i]] <- temp1
}
contrastDF2 <- do.call("rbind", contrastDF2)

# Add FDR preserving q-values:
contrastDF2$qvalue <- qvalue::qvalue(contrastDF2$pvalue)$qvalues

# If p-adjusted is NA make it 1:
contrastDF2$padj[is.na(contrastDF2$padj)] <- 1.0
contrastDF2$qvalue[contrastDF2$padj == 1.0]

# Volcano Plots:
contrastDF2 <- contrastDF2 %>% mutate(Comparison = case_when(
  pheno2 == "Injury_Control" & pheno1 == "Serratia_Infection" ~ "Serratia infection vs. Injury control",
  pheno2 == "Serratia_Infection" & pheno1 == "AhR_antagonist_p_Serratia_Infection" ~ "AhR antagonist vs. Serration infection",
  pheno2 == "dsGFP_p_Serratia_Infection" & pheno1 == "dsAhR_p_Serratia_Infection" ~ "dsAhR vs. dsGFP",
  pheno2 == "dsGFP_p_Serratia_Infection" & pheno1 == "dsTieg_p_Serratia_Infection" ~ "dsTIEG vs. dsGFP",
  pheno2 == "dsTieg_p_Serratia_Infection" & pheno1 == "dsAhR_p_Serratia_Infection" ~ "dsAhR vs. dsTIEG"))

# Save it for LB:
# save("contrastDF2", "df1", file = paste0("Results/contrastDF2_", gsub("-", "", Sys.Date()), ".RData"))

# Add labels for volcano plot:
contrastDF2$volLabel <- ifelse(contrastDF2$qvalue < 0.05 & abs(contrastDF2$log2FoldChange) > 2.5, 
                               gsub("AgaP_", "", contrastDF2$gene), "")
comparisons <- unique(contrastDF2$Comparison)
for(i in 1:length(comparisons)){
  png(file = paste0("Plots/", gsub("\\.", "", gsub(" ", "_", comparisons[i])), gsub("-","", Sys.Date()),".png"),
      height = 4, width = 5, units = "in", res = 600)
  show(ggplot(contrastDF2 %>% filter(Comparison == comparisons[i]), 
              aes(x = log2FoldChange, y = -log10(qvalue), label = volLabel)) + 
    geom_point(shape = 21, color = "grey30", fill = "dodgerblue", alpha = .5, size = .65) +
    geom_text_repel(size = 1.05, segment.colour = "grey30", segment.alpha = .5, segment.size = .25) + 
    geom_vline(xintercept = -2.5, lwd = .25, lty = 2) + geom_vline(xintercept = 2.5, lwd = .25, lty = 2) +
    geom_hline(yintercept = 1.30103, lwd = .25, lty = 2) +
    labs(x = expression(paste(log[2], "(Fold Change)")), y = expression(paste(-log[10], "(q-value)")),
         title = comparisons[i]))
  dev.off()
}

# Wide contrast data.frame:
contrastDF2wFC <- contrastDF2 %>% select(gene, Comparison, log2FoldChange) %>% 
  spread(key = Comparison, value = log2FoldChange)
names(contrastDF2wFC)[names(contrastDF2wFC) != "gene"] <- 
  paste0("FC_", names(contrastDF2wFC)[names(contrastDF2wFC) != "gene"])
contrastDF2wQ <- contrastDF2 %>% select(gene, Comparison, qvalue) %>% spread(key = Comparison, value = qvalue)
names(contrastDF2wQ)[names(contrastDF2wQ) != "gene"] <- 
  paste0("Q_", names(contrastDF2wQ)[names(contrastDF2wQ) != "gene"])
contrastDF2w <- contrastDF2wFC %>% left_join(contrastDF2wQ)

# Join with KEGG and GO annotation data:
contrastDF2w <- contrastDF2w %>% left_join(keggGOAnno)

# Export:
# writexl::write_xlsx(contrastDF2w, path = paste0("Results/contrastDFw_", gsub("-", "", Sys.Date()), ".xlsx"))

rm(contrastDF1, contrastDF2w, contrastDF2wQ, contrastDF2wFC, lowNCounts, des0,
   res0, temp1, i, filter1, rmExtra1, filter1b)

############ Transformation & Entropy / Significance filtering ############
# Get regularized log expression:
rlog1 <- rlog(des1, blind = TRUE)
rlog1 <- assay(rlog1)
rlog1 <- t(rlog1)

rlog1b <- rlog(des1b, blind = TRUE)
rlog1b <- assay(rlog1b)
rlog1b <- t(rlog1b)

# Information gain:
rlog1Temp <- as.data.frame(rlog1b)

rlog1Temp$oldSampName <- rownames(rlog1Temp)

rlog1Temp <- colData1b %>% select(oldSampName, pheno) %>% left_join(rlog1Temp) %>% select(-oldSampName)

rlog1Entropy <- FSelectorRcpp::information_gain(pheno ~ ., data = rlog1Temp, type = "gainratio")

rlog1Entropy <- rlog1Entropy %>% select(gene = attributes, infoGain = importance)

# Entropy filter:
entropyFun <- function(x) entropy::entropy(entropy::discretize(x, 8), unit = "log2")
colEntropies <- apply(as.matrix(rlog1), 2, entropyFun)
rlog1Entropy <- rlog1Entropy %>% left_join(data.frame(gene = names(colEntropies), entropy = colEntropies))

rm(rlog1Temp, colEntropies)

# Histograms
hist(rlog1Entropy$entropy)
hist(rlog1Entropy$infoGain)
median(rlog1Entropy$entropy)
median(rlog1Entropy$infoGain)

############ Multivariate visualization ############
pca1 <- prcomp(scale(rlog1, scale = FALSE))
pca1Scores <- as.data.frame(pca1$x[,1:4])
pca1Scores$oldSampName <- rownames(pca1Scores)
pca1Scores <- colData1 %>% left_join(pca1Scores)

# png(filename = "Plots/PC1vsPC2_20200926.png", height = 4.5, width = 7.5, units = "in", res = 300)
set.seed(3)
ggplot(pca1Scores, aes(x = PC1, y = PC2, color = pheno, label = oldSampName)) + geom_point() +
  geom_text_repel(size = 2) + labs(x = "PC1 (60.4% of variance)", y = "PC2 (12.2% of variance)")
# dev.off()

# png(filename = "Plots/PC3vsPC4_20200926.png", height = 4.5, width = 7.5, units = "in", res = 300)
set.seed(3)
ggplot(pca1Scores, aes(x = PC3, y = PC4, color = pheno, label = oldSampName)) + geom_point() +
  geom_text_repel(size = 2) + labs(x = "PC3 (6.6% of variance)", y = "PC4 (5.2% of variance)")
# dev.off()

PCA3D <- plotly::plot_ly(pca1Scores, x = ~PC1, y = ~PC3, z = ~PC4, color = ~pheno)
PCA3D <- PCA3D %>% plotly::add_markers()
PCA3D <- PCA3D %>% plotly::layout(scene = list(xaxis = list(title = 'PC1'),
                                   yaxis = list(title = 'PC3'),
                                   zaxis = list(title = 'PC4')))
PCA3D

# png(filename = "Plots/hclust.png", height = 6.5, width = 6.5, units = "in", res = 300)
plot(hclust(dist(rlog1), method = "ward.D2"))
# dev.off()

############ Sensitivity analysis ############
pca2 <- prcomp(scale(rlog1b, scale = FALSE))
summary(pca2)
pca2Scores <- as.data.frame(pca2$x[,1:4])
pca2Scores$oldSampName <- rownames(pca2Scores)
pca2Scores <- colData1 %>% inner_join(pca2Scores)

# png(filename = "Plots/PC1vsPC2_sens_20200926.png", height = 4.5, width = 7.5, units = "in", res = 300)
set.seed(3)
ggplot(pca2Scores, aes(x = PC1, y = PC2, color = pheno, label = oldSampName)) + geom_point() +
  geom_text_repel(size = 2) + labs(x = "PC1 (60.1% of variance)", y = "PC2 (11.2% of variance)")
# dev.off()

# png(filename = "Plots/PC3vsPC4_sens_20200926.png", height = 4.5, width = 7.5, units = "in", res = 300)
set.seed(3)
ggplot(pca2Scores, aes(x = PC3, y = PC4, color = pheno, label = oldSampName)) + geom_point() +
  geom_text_repel(size = 2) + labs(x = "PC3 (7.3% of variance)", y = "PC4 (4.1% of variance)")
# dev.off()

# png(filename = "Plots/hclust_sens.png", height = 6.5, width = 6.5, units = "in", res = 300)
plot(hclust(dist(rlog1b), method = "ward.D2"))
# dev.off()

rm(pca1, pca2, pca1Scores, pca2Scores, PCA3D)

############ WGCNA "fitting" ############
# Entropy filters on both:
# rlog2b <- rlog1b[, (rlog1Entropyb$entropy > quantile(rlog1Entropyb$entropy)['50%'])] # 4,973
rlog2b <- rlog1b

powers <- 1:20
sft <- WGCNA::pickSoftThreshold(rlog2b, powerVector = powers, verbose = 5)
par(mfrow = c(1, 2))
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, col = "red")
# this line corresponds to using an R^2 cut-off of h
abline(h = 0.95, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = "red")
softPower <- 9
par(oldPar)

# Adjacency and distance matrices
adjacency <- WGCNA::adjacency(rlog2b, type = "signed", power = softPower) 
# Turn adjacency into topological overlap
TOM <- WGCNA::TOMsimilarity(adjacency)
# TOM2 <- WGCNA::TOMsimilarity(adjacency, TOMType = "signed") # 'twas the same
rownames(TOM) <- colnames(TOM) <- rownames(adjacency)
dissTOM <- 1 - TOM
tree <- hclust(as.dist(dissTOM), method = "average")

# Module identification using dynamic tree cut:
dynamicMods <- dynamicTreeCut::cutreeDynamic(dendro = tree, distM = dissTOM, deepSplit = 4, pamRespectsDendro = TRUE,
                           minClusterSize = 20, method = "hybrid", cutHeight = .925)
dynamicColors <- WGCNA::labels2colors(dynamicMods)

# Plot dendogram and module assignment:
png(filename = "Plots/WGCNA_Dendro1.png",height = 6, width = 10, units = "in", res = 600)
WGCNA::plotDendroAndColors(tree, cbind(dynamicColors), "Module Assignment", dendroLabels = FALSE, hang = 0.0, 
                           addGuide = TRUE, guideHang = 0.0, main = "Gene Modules")
dev.off()

# Plot distnace matrix heatmap w/ module assignment
png(filename = "Plots/WGCNA_Heatmap1.png", height = 6, width = 6, units = "in", res = 600)
WGCNA::TOMplot(dissTOM, tree, dynamicColors, main = "Module heatmap")
dev.off()

# Module-gene mapping:
modDF <- data.frame(gene = rownames(TOM), module = dynamicColors)

mEigen1 <- WGCNA::moduleEigengenes(rlog2b, dynamicColors, impute = FALSE, nPC = 1, align = "along average", 
                                   excludeGrey = TRUE, grey = if (is.numeric(colors)) 0 else "grey", 
                                   softPower = 9, scale = TRUE, verbose = 5, indent = 1)

writexl::write_xlsx(modDF, path = paste0("Results/WGCNA_Modules_", gsub("-", "", Sys.Date()), ".xlsx"))

############ FGSEA analysis ############
# Make a list of genes that are in each module
MElist <- list()
for(color in unique(modDF$module)){
  MElist[[color]] <- modDF$gene[modDF$module == color]
}

# GSEA analysis:
comps <- unique(contrastDF2$Comparison)
fgseaRes <- list()
for(i in 1:length(comps)){
  comp <- comps[i]
  comp2 <- contrastDF2[contrastDF2$Comparison == comp, c("gene", "stat")] 
  comp2$stat <- abs(comp2$stat)
  comp2 <- comp2 %>% deframe()
  fgsea1 <- fgsea::fgsea(pathways = MElist, stats = comp2, nperm = 100000)
  fgsea1$Comparison <- comp
  fgseaRes[[i]] <- fgsea1
  print(i)
}
fgseaRes0 <- fgseaRes
fgseaRes <- do.call("rbind", fgseaRes)
fgseaRes <- fgseaRes %>% group_by(pathway) %>% mutate(maxLogP = max(-log10(padj), na.rm = TRUE))
fgseaRes <- fgseaRes %>% filter(!pathway == "grey")

# Make one results dataset for export:
fgseaResP <- fgseaRes %>% select(pathway, Comparison, pval) %>% spread(key = "Comparison", value = "pval")
fgseaResNES <- fgseaRes %>% select(pathway, Comparison, NES) %>% spread(key = "Comparison", value = "NES")
fgseaRes2 <- fgseaResP %>% full_join(fgseaResNES, by = "pathway", suffix = c(".pValue", ".NES"))
writexl::write_xlsx(fgseaRes2, path = paste0("Results/WGCNA_GSEA_", gsub("-", "", Sys.Date()), ".xlsx"))

save.image(file = "Data/running_20200830.RData")
load(file = "Data/running_20200830.RData")

# Plot:
png(filename = "Plots/WGCNA_GSEARes.png",height = 6, width = 10, units = "in", res = 600)
ggplot(fgseaRes %>% filter(maxLogP > 1.30103), aes(x = Comparison, y = -log10(padj), fill = Comparison)) + 
  geom_bar(color = "black", position = "dodge", stat = "identity", width = .7) +
  facet_wrap(~pathway, ncol = 7) + theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = wesanderson::wes_palette("Royal1", 7, type = "continuous"))
dev.off()

# Make GSEA plots:
for(i in 1:length(comparisons)){
  comp2 <- contrastDF2[contrastDF2$Comparison == comparisons[i], c("gene", "stat")]
  comp2$stat <- abs(comp2$stat)
  comp2 <- comp2 %>% deframe()
  for(j in 1:length(names(MElist))){
    color <- names(MElist)[j]
    png(filename = paste0("Plots/moduleEnrichmentPlots", color, "_", gsub(" |\\.", "", comparisons[i]), ".png"), height = 5, width = 7, units = "in",
        res = 300)
    show(fgsea::plotEnrichment(MElist[[color]], comp2) + labs(title = comparisons[i], subtitle = paste0("Module: ", color)))
    dev.off()
  }
}

# Blue module:
blueTOMhClust <- hclust(as.dist(TOM[modDF$gene[modDF$module == "blue"], modDF$gene[modDF$module == "blue"]]),
                          method = "average")
blueExpression <- scale(rlog2b[match(colData1b$oldSampName, rownames(rlog2b)), modDF$gene[modDF$module == "blue"]])
png(filename = "Plots/modules/blue.png", height = 6, width = 24, units = "in", res = 600)
blueColors <- WGCNA::numbers2colors(blueExpression, signed = TRUE, commonLim = FALSE)
WGCNA::plotDendroAndColors(dendro = blueTOMhClust, colors = t(blueColors),
                           groupLabels = rownames(blueExpression),
                           cex.dendroLabels = 0.1, cex.colorLabels = .7, marAll = c(1, 6, 3, 1), main = "")
dev.off()

# magenta module:
magentaTOMhClust <- hclust(as.dist(TOM[modDF$gene[modDF$module == "magenta"], modDF$gene[modDF$module == "magenta"]]),
                        method = "average")
magentaExpression <- scale(rlog2b[match(colData1b$oldSampName, rownames(rlog2b)), modDF$gene[modDF$module == "magenta"]])
png(filename = "Plots/modules/magenta.png", height = 6, width = 9, units = "in", res = 600)
magentaColors <- WGCNA::numbers2colors(magentaExpression, signed = TRUE, commonLim = FALSE)
WGCNA::plotDendroAndColors(dendro = magentaTOMhClust, colors = t(magentaColors),
                           groupLabels = rownames(magentaExpression),
                           cex.dendroLabels = 0.1, cex.colorLabels = .55, marAll = c(1, 7, 3, 1), main = "")
dev.off()

# midnightblue module:
midnightblueTOMhClust <- hclust(as.dist(TOM[modDF$gene[modDF$module == "midnightblue"], modDF$gene[modDF$module == "midnightblue"]]),
                           method = "average")
midnightblueExpression <- scale(rlog2b[match(colData1b$oldSampName, rownames(rlog2b)), modDF$gene[modDF$module == "midnightblue"]])
png(filename = "Plots/modules/midnightblue.png", height = 6, width = 9, units = "in", res = 600)
midnightblueColors <- WGCNA::numbers2colors(midnightblueExpression, signed = TRUE, commonLim = FALSE)
WGCNA::plotDendroAndColors(dendro = midnightblueTOMhClust, colors = t(midnightblueColors),
                           groupLabels = rownames(midnightblueExpression),
                           cex.dendroLabels = 0.35, cex.colorLabels = .55, marAll = c(1, 7, 3, 1), main = "")
dev.off()
