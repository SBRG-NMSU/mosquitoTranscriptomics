############ Prereqs ############
options(stringsAsFactors = FALSE, scipen = 900)
oldPar <- par()
os <- Sys.info()["sysname"]
baseDir <- ifelse(os == "Windows", "C:/Users/ptrainor/gdrive/Mosquitos/", 
                  "~/gdrive/Mosquitos/")
setwd(baseDir)

library(tidyverse)
library(DESeq2)
library(ggrepel)

theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5))

############ Import and process raw counts ############
# Count data:
df1 <- read.csv("Data/Priming_RNAseq_RawCounts_5-31-20.csv", header = TRUE, check.names = FALSE)
names(df1) <- gsub(" - total_read_count", "", names(df1))
names(table(df1$Name))[table(df1$Name) > 1]
rmExtra1 <- which(df1$Name == "tRNA-Leu")[2]
df1 <- df1[-rmExtra1,]
rownames(df1) <- df1$Name
df1$Name <- NULL

# Data with severe outlier removed:
df1b <- df1[, !(names(df1) ==  "Ent Prim + Ser inf 2")]

# Column (phenotype) data:
colData1 <- data.frame(oldSampName = names(df1))
colData1$sampName <- gsub(" ", "_", gsub(" \\+ ", "_p_", colData1$oldSampName))
colData1$pheno <- substr(colData1$sampName, start = 1, stop = nchar(colData1$sampName) - 2)
colData1$rep <- substr(colData1$sampName, start = nchar(colData1$sampName), stop = nchar(colData1$sampName))
colData1$pheno <- as.factor(colData1$pheno)

# Phenotype data w/o outlier:
colData1b <- colData1 %>% filter(oldSampName != "Ent Prim + Ser inf 2")

# Round fractional counts:
df1 <- round(df1) # 12,723 transcripts
df1b <- round(df1b)

# Filter out transcripts w/ 0 counts:
filter1 <- apply(df1, 1, function(x) sum(x == 0) == 30)
filter1b <- apply(df1b, 1, function(x) sum(x == 0) == 29)
df1 <- df1[!filter1, ] # 12,104
df1b <- df1b[!filter1b, ] # 12,092

############ Outlier dis-agreement ############
# Data.frame for pairwise plots:
df1Disagree <- df1 %>% select("Ent Prim + Ser inf 1", "Ent Prim + Ser inf 2", "Ent Prim + Ser inf 3")

# Pairwise correlations:
cor(log10(df1Disagree$`Ent Prim + Ser inf 1` + 1), log10(df1Disagree$`Ent Prim + Ser inf 2` + 1))
cor(log10(df1Disagree$`Ent Prim + Ser inf 2` + 1), log10(df1Disagree$`Ent Prim + Ser inf 3` + 1))
cor(log10(df1Disagree$`Ent Prim + Ser inf 1` + 1), log10(df1Disagree$`Ent Prim + Ser inf 3` + 1))

# Pairwise plots
p1 <- ggplot(df1Disagree, aes(x = log10(`Ent Prim + Ser inf 2` + 1), y = log10(`Ent Prim + Ser inf 1` + 1))) + 
  geom_point(shape = 1, size = .75) + geom_smooth(method = "lm", se = FALSE) + xlim(0, 8) + ylim(0, 8) +
  annotate("label", x = 2, y = 6, label = expression(italic(r) == 0.780), size = 8) +
  ggtitle("1 versus 2")
p2 <- ggplot(df1Disagree, aes(x = log10(`Ent Prim + Ser inf 2`+ 1), y = log10(`Ent Prim + Ser inf 3`+ 1))) + 
  geom_point(shape = 1, size = .75) + geom_smooth(method = "lm", se = FALSE) + xlim(0, 8) + ylim(0, 8) +
  annotate("label", x = 2, y = 6, label = expression(italic(r) == 0.821), size = 8)+
  ggtitle("3 versus 2")
p3 <- ggplot(df1Disagree, aes(x = log10(`Ent Prim + Ser inf 1`+ 1), y = log10(`Ent Prim + Ser inf 3`+ 1))) + 
  geom_point(shape = 1, size = .75) + geom_smooth(method = "lm", se = FALSE) + xlim(0, 8) + ylim(0, 8) +
  annotate("label", x = 2, y = 6, label = expression(italic(r) == 0.978), size = 8)+
  ggtitle("3 versus 1")

# png(file = "Plots/Disagreement1.png", height = 4, width = 12, units = "in", res = 300)
gridExtra::grid.arrange(p3, p1, p2, ncol = 3)
# dev.off()

rm(df1Disagree, p1, p2, p3)

############ Univariate analysis ############
# Set up DESeq2 analysis:
des1 <- DESeqDataSetFromMatrix(countData = df1, colData = colData1, design = ~ pheno)
des1b <- DESeqDataSetFromMatrix(countData = df1b, colData = colData1b, design = ~ pheno)
des1 <- DESeq(des1)
des1b <- DESeq(des1b)

# Make data.frame of comparisons between priming and Naive:
contrastDF1 <- data.frame(pheno1 = c("Enterobacter_Priming", "Serratia_Priming", "Serratia_Priming",
                                     "Enterobacter_infection", "Serratia_Infection",
                                     "Ent_Prim_p_Ent_inf", "Ser_Prim_p_Ser_inf"),
                          pheno2 = c("Naive", "Naive", "Enterobacter_Priming", 
                                     "Injury_Control", "Injury_Control",
                                     "Enterobacter_infection", "Serratia_Infection"))

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
# save("contrastDF2", "df1b", file = "Results/contrastDF2_20200606.RData")

# If p-adjusted is NA make it 1:
contrastDF2$padj[is.na(contrastDF2$padj)] <- 1.0

# LRT:
des0 <- DESeq(des1b, test="LRT", reduced = ~ 1)
res0 <- as.data.frame(results(des0))

# Volcano Plots:
contrastDF2 <- contrastDF2 %>% mutate(Comparison = case_when(
  pheno2 == "Enterobacter_Priming" & pheno1 ==  "Serratia_Priming" ~ "Serratia vs. Enterobacter",
  pheno2 == "Naive" & pheno1 == "Serratia_Priming" ~ "Serratia vs. Naive",
  pheno2 == "Naive" & pheno1 == "Enterobacter_Priming" ~ "Enterobacter vs. Naive",
  pheno2 == "Injury_Control" & pheno1 == "Enterobacter_infection" ~ "Ent Inf vs Inj Ctrl",
  pheno2 == "Injury_Control" & pheno1 == "Serratia_Infection" ~ "Ser Inf vs Inj Ctrl",
  pheno2 == "Enterobacter_infection" & pheno1 == "Ent_Prim_p_Ent_inf" ~ "Ent Prim & Inf vs Ent Inf",
  pheno2 == "Serratia_Infection" & pheno1 == "Ser_Prim_p_Ser_inf" ~ "Ser Prim & Inf vs Ser Inf"))
contrastDF2$Comparison <- factor(contrastDF2$Comparison, levels = 
                                   c("Enterobacter vs. Naive", "Serratia vs. Naive", "Serratia vs. Enterobacter",
                                     "Ent Inf vs Inj Ctrl", "Ser Inf vs Inj Ctrl", "Ent Prim & Inf vs Ent Inf",
                                     "Ser Prim & Inf vs Ser Inf"))

# Add labels for volcano plot:
contrastDF2$volLabel <- ifelse(contrastDF2$padj < 0.05 & abs(contrastDF2$log2FoldChange) > 2.5, 
                               gsub("AgaP_", "", contrastDF2$gene), "")

# png(file = "Plots/PrimingVolcano1.png", height = 4, width = 12, units = "in", res = 600)
set.seed(3333)
ggplot(contrastDF2 %>% filter(Comparison %in% c("Enterobacter vs. Naive", "Serratia vs. Naive", "Serratia vs. Enterobacter") & padj < 1.0), 
       aes(x = log2FoldChange, y = -log10(padj), label = volLabel)) + 
  geom_point(shape = 21, color = "grey30", fill = "dodgerblue", alpha = .5, size = .65) + 
  geom_text_repel(size = 1.05, segment.colour = "grey30", segment.alpha = .5, segment.size = .25) +
  geom_vline(xintercept = -2.5, lwd = .25, lty = 2) + geom_vline(xintercept = 2.5, lwd = .25, lty = 2) +
  geom_hline(yintercept = 1.30103, lwd = .25, lty = 2) + 
  labs(x = expression(paste(log[2], "(Fold Change)")), y = expression(paste(-log[10], "(adjusted p-value)"))) +
  facet_grid(~Comparison)
# dev.off()

set.seed(3333)
ggplot(contrastDF2 %>% filter(Comparison %in% c("Ent Inf vs Inj Ctrl", "Ser Inf vs Inj Ctrl", "Ent Prim & Inf vs Ent Inf",
                                                "Ser Prim & Inf vs Ser Inf") & padj < 1.0), 
       aes(x = log2FoldChange, y = -log10(padj), label = volLabel)) + 
  geom_point(shape = 21, color = "grey30", fill = "dodgerblue", alpha = .5, size = .65) + 
  geom_text_repel(size = 1.05, segment.colour = "grey30", segment.alpha = .5, segment.size = .25) +
  geom_vline(xintercept = -2.5, lwd = .25, lty = 2) + geom_vline(xintercept = 2.5, lwd = .25, lty = 2) +
  geom_hline(yintercept = 1.30103, lwd = .25, lty = 2) + 
  labs(x = expression(paste(log[2], "(Fold Change)")), y = expression(paste(-log[10], "(adjusted p-value)"))) +
  facet_wrap(~Comparison)

# Concordant-discordant analysis:
contrastDF2wFC <- contrastDF2 %>% select(gene, Comparison, log2FoldChange) %>% spread(key = Comparison, value = log2FoldChange)
names(contrastDF2wFC)[names(contrastDF2wFC) != "gene"] <- paste0("FC_", names(contrastDF2wFC)[names(contrastDF2wFC) != "gene"])
contrastDF2wAdjP <- contrastDF2 %>% select(gene, Comparison, padj) %>% spread(key = Comparison, value = padj)
names(contrastDF2wAdjP)[names(contrastDF2wAdjP) != "gene"] <- 
  paste0("AdjP_", names(contrastDF2wAdjP)[names(contrastDF2wAdjP) != "gene"])
contrastDF2w <- contrastDF2wFC %>% left_join(contrastDF2wAdjP)

cor(x = contrastDF2w$`FC_Ent Inf vs Inj Ctrl`, y = contrastDF2w$`FC_Ser Inf vs Inj Ctrl`)
ggplot(contrastDF2w, aes(x = `FC_Ent Inf vs Inj Ctrl`, y = `FC_Ser Inf vs Inj Ctrl`)) + 
  geom_point(color = "dodgerblue", shape = 1, show.legend = FALSE) + stat_smooth(method = "lm")

# Priming concordant-discordant analysis:
contrastDF2a <- contrastDF2 %>% filter(!(pheno2 == "Enterobacter_Priming" & pheno1 ==  "Serratia_Priming"))  %>% select(-pheno2)
contrastDF2a1 <- contrastDF2a %>% filter(pheno1 == "Enterobacter_Priming")
contrastDF2a2 <- contrastDF2a %>% filter(pheno1 == "Serratia_Priming")
contrastDF2a <- contrastDF2a1 %>% full_join(contrastDF2a2 %>% select(-baseMean), by = "gene", suffix = c("_eP", "_sP"))
contrastDF2a <- contrastDF2a %>% filter(padj_eP < 0.05 | padj_sP < 0.05) # 1,114 genes
# Determine if concordant:
constrastDF2b <- contrastDF2a %>% mutate(catE = 
                         case_when(padj_eP < 0.05 & log2FoldChange_eP > 0 ~ "up05",
                                   padj_eP < 0.05 & log2FoldChange_eP < 0 ~ "down05",
                                   padj_eP >= 0.05 & log2FoldChange_eP > 0 ~ "upNS",
                                   padj_eP >= 0.05 & log2FoldChange_eP < 0 ~ "downNS"),
                         catS = 
                           case_when(padj_sP < 0.05 & log2FoldChange_sP > 0 ~ "up05",
                                     padj_sP < 0.05 & log2FoldChange_sP < 0 ~ "down05",
                                     padj_sP >= 0.05 & log2FoldChange_sP > 0 ~ "upNS",
                                     padj_sP >= 0.05 & log2FoldChange_sP < 0 ~ "downNS"))
xtabs(~catE + catS, data = constrastDF2b)

# Label some:
contrastDF2a$lab <- contrastDF2a$gene
contrastDF2a$lab <- gsub("AgaP_", "", contrastDF2a$lab)
contrastDF2a$lab[!(abs(contrastDF2a$log2FoldChange_eP) > 2.5 | abs(contrastDF2a$log2FoldChange_sP) > 2.5)] <- ""

# png(file = "Plots/PrimingAgreement1.png", height = 5, width = 6.25, units = "in", res = 300)
set.seed(3)
ggplot(contrastDF2a, aes(x = log2FoldChange_eP, y = log2FoldChange_sP, label = lab)) + 
  geom_point(color = "dodgerblue", shape = 1, show.legend = FALSE) + 
  geom_text_repel(size = 2, segment.colour = "grey30", segment.alpha = .5) +
  geom_vline(xintercept = 0, lty = 2, lwd = .25) + geom_hline(yintercept = 0, lty = 2, lwd = .25) +
  labs(x = "Enterobacter Priming / Naive (Log2 FC)", y = "Serratia Priming / Naive (Log2 FC)") 
# dev.off()

# Comparison of priming types:
contrastDF2c <- contrastDF2 %>% filter(pheno1 == "Enterobacter_Priming" & pheno2 == "Serratia_Priming" & padj < .05)

rm(constrastDF2b, contrastDF1, contrastDF2a, contrastDF2a1, contrastDF2a2, contrastDF2c, des0,
   res0, temp1, i, filter1, rmExtra1)

############ Transformation & Entropy filtering ############
# Get regularized log expression:
rlog1 <- rlog(des1, blind = TRUE)
rlog1b <- rlog(des1b, blind = TRUE)
rlog1 <- assay(rlog1)
rlog1b <- assay(rlog1b)
rlog1 <- t(rlog1)
rlog1b <- t(rlog1b)

# rlogPrime <- rlog1[, grepl("Serratia Priming", colnames(rlog1)) | grepl("Enterobacter Priming", colnames(rlog1)) |
#                                                       grepl("Naive", colnames(rlog1))]

# Information gain:
rlog1Temp <- as.data.frame(rlog1)
rlog1Tempb <- as.data.frame(rlog1b)

rlog1Temp$oldSampName <- rownames(rlog1Temp)
rlog1Tempb$oldSampName <- rownames(rlog1Tempb)

rlog1Temp <- colData1 %>% select(oldSampName, pheno) %>% left_join(rlog1Temp) %>% select(-oldSampName)
rlog1Tempb <- colData1b %>% select(oldSampName, pheno) %>% left_join(rlog1Tempb) %>% select(-oldSampName)

rlog1Entropy <- FSelectorRcpp::information_gain(pheno ~ ., data = rlog1Temp)
rlog1Entropyb <- FSelectorRcpp::information_gain(pheno ~ ., data = rlog1Tempb)

rlog1Entropy <- rlog1Entropy %>% select(gene = attributes, infoGain = importance)
rlog1Entropyb <- rlog1Entropyb %>% select(gene = attributes, infoGain = importance)

# Entropy filter:
entropyFun <- function(x) entropy::entropy(entropy::discretize(x, 8), unit = "log2")
colEntropies <- apply(as.matrix(rlog1), 2, entropyFun)
rlog1Entropy <- rlog1Entropy %>% left_join(data.frame(gene = names(colEntropies), entropy = colEntropies))

colEntropiesb <- apply(as.matrix(rlog1b), 2, entropyFun)
rlog1Entropyb <- rlog1Entropyb %>% left_join(data.frame(gene = names(colEntropiesb), entropy = colEntropiesb))

rm(rlog1Temp, rlog1Tempb, colEntropies, colEntropiesb)

# Histograms
hist(rlog1Entropy$entropy)
hist(rlog1Entropy$infoGain)
median(rlog1Entropy$entropy)
median(rlog1Entropy$infoGain)

############ Multivariate visualization ############
pca1 <- prcomp(scale(rlog2, scale = FALSE))
pca1Scores <- as.data.frame(pca1$x[,1:4])
pca1Scores$oldSampName <- rownames(pca1Scores)
pca1Scores <- colData1 %>% left_join(pca1Scores)

png(filename = "Plots/PC1vsPC2.png", height = 4.5, width = 6.5, units = "in", res = 300)
set.seed(3)
ggplot(pca1Scores, aes(x = PC1, y = PC2, color = pheno, label = oldSampName)) + geom_point() +
  geom_text_repel(size = 2)
dev.off()

png(filename = "Plots/PC3vsPC4.png", height = 4.5, width = 6.5, units = "in", res = 300)
set.seed(3)
ggplot(pca1Scores, aes(x = PC3, y = PC4, color = pheno, label = oldSampName)) + geom_point() +
  geom_text_repel(size = 2)
dev.off()

PCA3D <- plotly::plot_ly(pca1Scores, x = ~PC1, y = ~PC3, z = ~PC4, color = ~pheno)
PCA3D <- PCA3D %>% plotly::add_markers()
PCA3D <- PCA3D %>% plotly::layout(scene = list(xaxis = list(title = 'PC1'),
                                   yaxis = list(title = 'PC3'),
                                   zaxis = list(title = 'PC4')))
PCA3D

png(filename = "Plots/hclust.png", height = 6.5, width = 6.5, units = "in", res = 300)
plot(hclust(dist(rlog2), method = "ward.D2"))
dev.off()

############ Sensitivity analysis ############
pca2 <- prcomp(scale(rlog2[rownames(rlog2) != "Ent Prim + Ser inf 2",], scale = FALSE))
pca2Scores <- as.data.frame(pca2$x[,1:4])
pca2Scores$oldSampName <- rownames(pca2Scores)
pca2Scores <- colData1 %>% inner_join(pca2Scores)

png(filename = "Plots/PC1vsPC2_sens.png", height = 4.5, width = 6.5, units = "in", res = 300)
set.seed(3)
ggplot(pca2Scores, aes(x = PC1, y = PC2, color = pheno, label = oldSampName)) + geom_point() +
  geom_text_repel(size = 2)
dev.off()

png(filename = "Plots/PC3vsPC4_sens.png", height = 4.5, width = 6.5, units = "in", res = 300)
set.seed(3)
ggplot(pca2Scores, aes(x = PC3, y = PC4, color = pheno, label = oldSampName)) + geom_point() +
  geom_text_repel(size = 2)
dev.off()

png(filename = "Plots/hclust_sens.png", height = 6.5, width = 6.5, units = "in", res = 300)
plot(hclust(dist(rlog2[rownames(rlog2) != "Ent Prim + Ser inf 2",]), method = "ward.D2"))
dev.off()

rm(pca1, pca2, pca1Scores, pca2Scores, PCA3D, colEntropies, colEntropiesb)

############ WGCNA "fitting" ############

# Entropy filters on both:
rlog2 <- rlog1[, (rlog1Entropy$entropy > quantile(rlog1Entropy$entropy)['75%'] | 
                    rlog1Entropy$infoGain > quantile(rlog1Entropy$infoGain)['75%'])] # 4,973
rlog2b <- rlog1b[, (rlog1Entropyb$entropy > quantile(rlog1Entropyb$entropy)['75%'] | 
                      rlog1Entropyb$infoGain > quantile(rlog1Entropyb$infoGain)['75%'])] # 4,973

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
abline(h = 0.80, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = "red")
softPower <- 6
par(oldPar)

# Adjacency and distance matrices
adjacency <- WGCNA::adjacency(rlog2b, power = softPower) #corOptions="use = 'p', method = 'spearman'"
# Turn adjacency into topological overlap
TOM <- WGCNA::TOMsimilarity(adjacency)
rownames(TOM) <- colnames(TOM) <- rownames(adjacency)
dissTOM <- 1 - TOM
tree <- hclust(as.dist(dissTOM), method = "average")

# Module identification using dynamic tree cut:
dynamicMods <- dynamicTreeCut::cutreeDynamic(dendro = tree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = TRUE,
                           minClusterSize = 20, method = "hybrid")
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

############ Module eigengenes ############
mEigen1 <- WGCNA::moduleEigengenes(rlog2b, dynamicColors, impute = FALSE, nPC = 1, align = "along average", excludeGrey = TRUE, 
                 grey = if (is.numeric(colors)) 0 else "grey", softPower = 6, scale = TRUE,
                 verbose = 5, indent = 1)

save.image(file = "Data/running_20200606.RData")

# Make into long data.frame:
mEigen2 <- mEigen1$eigengenes
mEigen2$oldSampName <- rownames(mEigen2)
mEigen2 <- mEigen2 %>% gather(key = "ME", value = "value", -oldSampName)
mEigen2 <- colData1b %>% left_join(mEigen2)

# List of all module colors:
mEColors <- unique(mEigen2$ME)

# Make comparison of ME's:
comparisons <- data.frame(pheno1 = c("Naive", "Naive", "Enterobacter_Priming"), 
                          pheno2 = c("Enterobacter_Priming", "Serratia_Priming", "Serratia_Priming"))
dfTemp3 <- data.frame()
for(i in 1:nrow(comparisons)){
  for(j in 1:length(mEColors)){
    dfTemp1 <- mEigen2 %>% filter(pheno %in% c(comparisons$pheno1[i], comparisons$pheno2[i]) & ME == mEColors[j])
    dfTemp1$pheno <- factor(dfTemp1$pheno, levels = c(comparisons$pheno1[i], comparisons$pheno2[i]))
    glm1 <- glm(pheno ~ value, data = dfTemp1, family = "binomial")
    glm0 <- update(glm1, . ~ . -value)
    glm1LRT <- anova(glm0, glm1, test = "LRT")$`Pr(>Chi)`[2]
    dfTemp2 <- data.frame(pheno1 = comparisons$pheno1[i], pheno2 = comparisons$pheno2[i], ME = mEColors[j], lrtP = glm1LRT)
    dfTemp3 <- rbind(dfTemp3, dfTemp2)
    cat("i is ", i, "j is ", j, "\n")
  }
}

# MElightyellow heatmap
dfTemp4 <- rlog2b[grepl("Priming", rownames(rlog2b)) | grepl("Naive", rownames(rlog2b)), 
                  colnames(rlog2b) %in% 
                    modDF$gene[modDF$module %in% c("salmon", "grey60", "darkturquoise", "paleturquoise", "purple")]]
dfTemp4b <- rlog2b[, colnames(rlog2b) %in% 
                    modDF$gene[modDF$module %in% c("salmon")]]
heatmap(scale(dfTemp4))
heatmap(scale(dfTemp4b))
