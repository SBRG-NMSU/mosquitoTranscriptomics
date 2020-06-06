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

# Column (phenotype) data:
colData1 <- data.frame(oldSampName = names(df1))
colData1$sampName <- gsub(" ", "_", gsub(" \\+ ", "_p_", colData1$oldSampName))
colData1$pheno <- substr(colData1$sampName, start = 1, stop = nchar(colData1$sampName) - 2)
colData1$rep <- substr(colData1$sampName, start = nchar(colData1$sampName), stop = nchar(colData1$sampName))
colData1$pheno <- as.factor(colData1$pheno)

# Round fractional counts:
df1 <- round(df1) # 12,723 transcripts

# Filter out transcripts w/ 0 counts:
filter1 <- apply(df1, 1, function(x) sum(x == 0) == 30)
df1 <- df1[!filter1, ] 

############ Univariate analysis ############
des1 <- DESeqDataSetFromMatrix(countData = df1, colData = colData1, design = ~ pheno)
des1 <- DESeq(des1)
contrastDF1 <- data.frame(pheno1 = c("Enterobacter_Priming", "Serratia_Priming", "Serratia_Priming"),
                          pheno2 = c("Naive", "Naive", "Enterobacter_Priming"))

# Make comparisons:
contrastDF2 <- list()
for(i in 1:nrow(contrastDF1)){
  temp1 <- as.data.frame(results(des1, contrast = c("pheno", contrastDF1$pheno1[i], contrastDF1$pheno2[i]), cooksCutoff = Inf))
  temp1$lfcLower <- temp1$log2FoldChange - qnorm(0.975) * temp1$lfcSE
  temp1$lfcUpper <- temp1$log2FoldChange + qnorm(0.975) * temp1$lfcSE
  temp1$pheno1 <- contrastDF1$pheno1[i]
  temp1$pheno2 <- contrastDF1$pheno2[i]
  temp1$gene <- rownames(temp1)
  contrastDF2[[i]] <- temp1
}
contrastDF2 <- do.call("rbind", contrastDF2)
save("contrastDF2", "df1", file = "Results/contrastDF2_20200605.RData")

# If p-adjusted is NA make it 1:
contrastDF2$padj[is.na(contrastDF2$padj)] <- 1.0

# LRT:
des0 <- DESeq(des1, test="LRT", reduced = ~ 1)
res0 <- as.data.frame(results(des0))

# Volcano Plots:
contrastDF2 <- contrastDF2 %>% mutate(Comparison = case_when(
  pheno2 == "Enterobacter_Priming" & pheno1 ==  "Serratia_Priming" ~ "Serratia vs. Enterobacter",
  pheno2 == "Naive" & pheno1 == "Serratia_Priming" ~ "Serratia vs. Naive",
  pheno2 == "Naive" & pheno1 == "Enterobacter_Priming" ~ "Enterobacter vs. Naive"))
contrastDF2$Comparison <- factor(contrastDF2$Comparison, levels = 
                                   c("Enterobacter vs. Naive", "Serratia vs. Naive", "Serratia vs. Enterobacter"))

# Add labels for volcano plot:
contrastDF2$volLabel <- ifelse(contrastDF2$padj < 0.05 & abs(contrastDF2$log2FoldChange) > 2.5, 
                               gsub("AgaP_", "", contrastDF2$gene), "")

png(file = "Plots/PrimingVolcano1.png", height = 4, width = 12, units = "in", res = 600)
set.seed(3333)
ggplot(contrastDF2 %>% filter(padj < 1.0), aes(x = log2FoldChange, y = -log10(padj), label = volLabel)) + 
  geom_point(shape = 21, color = "grey30", fill = "dodgerblue", alpha = .5, size = .65) + 
  geom_text_repel(size = 1.05, segment.colour = "grey30", segment.alpha = .5, segment.size = .25) +
  geom_vline(xintercept = -2.5, lwd = .25, lty = 2) + geom_vline(xintercept = 2.5, lwd = .25, lty = 2) +
  geom_hline(yintercept = 1.30103, lwd = .25, lty = 2) + 
  labs(x = expression(paste(log[2], "(Fold Change)")), y = expression(paste(-log[10], "(adjusted p-value)"))) +
  facet_grid(~Comparison)
dev.off()

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

png(file = "Plots/PrimingAgreement1.png", height = 5, width = 6.25, units = "in", res = 300)
set.seed(3)
ggplot(contrastDF2a, aes(x = log2FoldChange_eP, y = log2FoldChange_sP, label = lab)) + 
  geom_point(color = "dodgerblue", shape = 1, show.legend = FALSE) + 
  geom_text_repel(size = 2, segment.colour = "grey30", segment.alpha = .5) +
  geom_vline(xintercept = 0, lty = 2, lwd = .25) + geom_hline(yintercept = 0, lty = 2, lwd = .25) +
  labs(x = "Enterobacter Priming / Naive (Log2 FC)", y = "Serratia Priming / Naive (Log2 FC)") 
dev.off()

# Comparison of priming types:
contrastDF2c <- contrastDF2 %>% filter(pheno1 == "Enterobacter_Priming" & pheno2 == "Serratia_Priming" & padj < .05)

rm(constrastDF2b, contrastDF1, contrastDF2a, contrastDF2a1, contrastDF2a2, contrastDF2c, des0,
   res0, temp1, i)

############ Transformation & Entropy filtering ############
rlog1 <- rlog(des1, blind = TRUE)
rlog1 <- assay(rlog1)
rlog1 <- t(rlog1)
# rlogPrime <- rlog1[, grepl("Serratia Priming", colnames(rlog1)) | grepl("Enterobacter Priming", colnames(rlog1)) |
#                                                       grepl("Naive", colnames(rlog1))]

entropyFun <- function(x) entropy::entropy(entropy::discretize(x, 10), unit = "log2")
colEntropies <- apply(as.matrix(rlog1), 2, entropyFun)
hist(colEntropies)
rlog2 <- rlog1[, colEntropies > 1] # 11,987

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
