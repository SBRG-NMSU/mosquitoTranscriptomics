############ Prereqs ############
options(stringsAsFactors = FALSE, scipen = 900)
oldPar <- par()
os <- Sys.info()["sysname"]
baseDir <- ifelse(os == "Windows", "C:/Users/ptrainor/gdrive/Mosquitos/Priming_Infection/", 
                  "~/gdrive/Mosquitos/Priming_Infection/")
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
# Get rid of extra part in name:
df1$Name <- gsub("AgaP_", "", df1$Name)
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

############ Import annotation data ############
# Imports:
biomartDF <- read.csv("./biomart_agambiae_data.csv")
ensemblIDDF <- read.csv("./contrastDF2_ensembl_ids2.csv") %>% select(gene, ensembl_id) %>% unique()
# keggDF <- read.csv("./contrast_data_ensembl_kegg.csv") %>% select(gene, ensembl_id, kegg_id) %>% unique()
load("kegg_pathway_data.Rdata")

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
                               colData = colData1, design = ~ pheno) # 10726
des1b <- DESeqDataSetFromMatrix(countData = df1b[!rownames(df1b) %in% lowNCounts,], 
                                colData = colData1b, design = ~ pheno) # 10714
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

# Add FDR preserving q-values:
contrastDF2$qvalue <- qvalue::qvalue(contrastDF2$pvalue)$qvalues

# If p-adjusted is NA make it 1:
contrastDF2$padj[is.na(contrastDF2$padj)] <- 1.0
contrastDF2$qvalue[contrastDF2$padj == 1.0]

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

# save("contrastDF2", "df1b", file = paste0("Results/contrastDF2_", gsub("-", "", Sys.Date()), ".RData"))

# Add labels for volcano plot:
contrastDF2$volLabel <- ifelse(contrastDF2$qvalue < 0.05 & abs(contrastDF2$log2FoldChange) > 2.5, 
                               gsub("AgaP_", "", contrastDF2$gene), "")

# png(file = paste0("Plots/PrimingVolcano1_", gsub("-","", Sys.Date()),".png"), height = 4, width = 8, units = "in", res = 600)
set.seed(3333)
ggplot(contrastDF2 %>% filter(Comparison %in% c("Enterobacter vs. Naive", "Serratia vs. Naive")), 
       aes(x = log2FoldChange, y = -log10(qvalue), label = volLabel)) + 
  geom_point(shape = 21, color = "grey30", fill = "dodgerblue", alpha = .5, size = .65) + 
  geom_text_repel(size = 1.05, segment.colour = "grey30", segment.alpha = .5, segment.size = .25) +
  geom_vline(xintercept = -2.5, lwd = .25, lty = 2) + geom_vline(xintercept = 2.5, lwd = .25, lty = 2) +
  geom_hline(yintercept = 1.30103, lwd = .25, lty = 2) + 
  labs(x = expression(paste(log[2], "(Fold Change)")), y = expression(paste(-log[10], "(q-value)"))) +
  facet_grid(~Comparison)
# dev.off()

# png(file = paste0("Plots/PrimingVolcano2_", gsub("-","", Sys.Date()),".png"), height = 4, width = 5, units = "in", res = 600)
set.seed(3333)
ggplot(contrastDF2 %>% filter(Comparison %in% c("Serratia vs. Enterobacter")), 
       aes(x = log2FoldChange, y = -log10(qvalue), label = volLabel)) + 
  geom_point(shape = 21, color = "grey30", fill = "dodgerblue", alpha = .5, size = .65) + 
  geom_text_repel(size = 2, segment.colour = "grey30", segment.alpha = .5, segment.size = .25) +
  geom_vline(xintercept = -2.5, lwd = .25, lty = 2) + geom_vline(xintercept = 2.5, lwd = .25, lty = 2) +
  geom_hline(yintercept = 1.30103, lwd = .25, lty = 2) + 
  labs(x = expression(paste(log[2], "(Fold Change)")), y = expression(paste(-log[10], "(q-value)")))
# dev.off()

# png(file = paste0("Plots/InfectionVolcano_", gsub("-","", Sys.Date()),".png"), height = 8, width = 9, units = "in", res = 600)
set.seed(3333)
prevLabs <- c("Ent Inf vs Inj Ctrl" = "Enterobacter Infection vs. Inj Ctrl", 
              "Ser Inf vs Inj Ctrl" = "Serratia Infection vs Inj Ctrl", 
              "Ent Prim & Inf vs Ent Inf" = "Enterobacter Priming & Infection vs. Infection Only",
              "Ser Prim & Inf vs Ser Inf" = "Serratia Priming & Infection vs. Infection Only")
ggplot(contrastDF2 %>% filter(Comparison %in% c("Ent Inf vs Inj Ctrl", "Ser Inf vs Inj Ctrl", "Ent Prim & Inf vs Ent Inf",
                                                "Ser Prim & Inf vs Ser Inf")), 
       aes(x = log2FoldChange, y = -log10(qvalue), label = volLabel)) + 
  geom_point(shape = 21, color = "grey30", fill = "dodgerblue", alpha = .5, size = .65) + 
  geom_text_repel(size = 1.05, segment.colour = "grey30", segment.alpha = .5, segment.size = .25) +
  geom_vline(xintercept = -2.5, lwd = .25, lty = 2) + geom_vline(xintercept = 2.5, lwd = .25, lty = 2) +
  geom_hline(yintercept = 1.30103, lwd = .25, lty = 2) + 
  labs(x = expression(paste(log[2], "(Fold Change)")), y = expression(paste(-log[10], "(q-value)"))) +
  facet_wrap(~Comparison, labeller = labeller(Comparison = prevLabs))
rm(prevLabs)
# dev.off()

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

# Infection concordance:
contrastDF2a <- contrastDF2w %>% filter(`Q_Ent Inf vs Inj Ctrl` < 0.05 | `Q_Ser Inf vs Inj Ctrl` < 0.05) # 3,961 genes
contrastDF2a$lab <- contrastDF2a$gene
contrastDF2a$lab[!(abs(contrastDF2a$`FC_Ent Inf vs Inj Ctrl`) > 2.5 | 
                     abs(contrastDF2a$`FC_Ser Inf vs Inj Ctrl`) > 2.5)] <- ""

png(file = "Plots/InfectionAgreement1.png", height = 5, width = 6, units = "in", res = 300)
ggplot(contrastDF2a , aes(x = `FC_Ent Inf vs Inj Ctrl`, y = `FC_Ser Inf vs Inj Ctrl`, label = lab)) + 
  geom_point(color = "dodgerblue", shape = 1, show.legend = FALSE) +
  geom_text_repel(size = 2, segment.colour = "grey30", segment.alpha = .5) +
  geom_vline(xintercept = 0, lty = 2, lwd = .25) + geom_hline(yintercept = 0, lty = 2, lwd = .25) +
  labs(x = "Enterobacter Infection / Injury Control (Log2 FC)", y = "Serratia Infection / Injury Control (Log2 FC)") 
dev.off()
rm(contrastDF2a)

# Priming concordant-discordant analysis:
contrastDF2a <- contrastDF2w %>% filter(`Q_Enterobacter vs. Naive` < 0.05 | `Q_Serratia vs. Naive` < 0.05) # 1,562 genes
# Determine if concordant:
constrastDF2b <- contrastDF2a %>% mutate(catE = 
                         case_when(`Q_Enterobacter vs. Naive` < 0.05 & `FC_Enterobacter vs. Naive` > 0 ~ "up05",
                                   `Q_Enterobacter vs. Naive` < 0.05 & `FC_Enterobacter vs. Naive` < 0 ~ "down05",
                                   `Q_Enterobacter vs. Naive` >= 0.05 & `FC_Enterobacter vs. Naive` > 0 ~ "upNS",
                                   `Q_Enterobacter vs. Naive` >= 0.05 & `FC_Enterobacter vs. Naive` < 0 ~ "downNS"),
                         catS = 
                           case_when(`Q_Serratia vs. Naive` < 0.05 & `FC_Serratia vs. Naive` > 0 ~ "up05",
                                     `Q_Serratia vs. Naive` < 0.05 & `FC_Serratia vs. Naive` < 0 ~ "down05",
                                     `Q_Serratia vs. Naive` >= 0.05 & `FC_Serratia vs. Naive`> 0 ~ "upNS",
                                     `Q_Serratia vs. Naive` >= 0.05 & `FC_Serratia vs. Naive` < 0 ~ "downNS"))
xtabs(~catE + catS, data = constrastDF2b)
writexl::write_xlsx(constrastDF2b, path = paste0("Results/contrastDF2b_Priming_", gsub("-", "", Sys.Date()), ".xlsx"))

# Label some:
contrastDF2a$lab <- contrastDF2a$gene
contrastDF2a$lab[!(abs(contrastDF2a$`FC_Enterobacter vs. Naive`) > 2.5 | 
                     abs(contrastDF2a$`FC_Serratia vs. Naive`) > 2.5)] <- ""

# png(file = "Plots/PrimingAgreement1.png", height = 5, width = 6, units = "in", res = 300)
set.seed(3)
ggplot(contrastDF2a, aes(x = `FC_Enterobacter vs. Naive`, y = `FC_Serratia vs. Naive`, label = lab)) + 
  geom_point(color = "dodgerblue", shape = 1, show.legend = FALSE) + 
  geom_text_repel(size = 2, segment.colour = "grey30", segment.alpha = .5) +
  geom_vline(xintercept = 0, lty = 2, lwd = .25) + geom_hline(yintercept = 0, lty = 2, lwd = .25) +
  labs(x = "Enterobacter Priming / Naive (Log2 FC)", y = "Serratia Priming / Naive (Log2 FC)") 
# dev.off()

# Comparison of priming types:
contrastDF2c <- contrastDF2 %>% filter(pheno1 == "Enterobacter_Priming" & pheno2 == "Serratia_Priming" & qvalue < .05)

rm(constrastDF2b, contrastDF1, contrastDF2a, contrastDF2w, contrastDF2wQ, contrastDF2wFC, lowNCounts, des0,
   res0, temp1, i, filter1, rmExtra1, filter1b, contrastDF2c)

############ Transformation & Entropy / Significance filtering ############
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

rlog1Entropy <- FSelectorRcpp::information_gain(pheno ~ ., data = rlog1Temp, type = "gainratio")
rlog1Entropyb <- FSelectorRcpp::information_gain(pheno ~ ., data = rlog1Tempb, type = "gainratio")

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
hist(rlog1Entropyb$entropy)
hist(rlog1Entropyb$infoGain)
median(rlog1Entropyb$entropy)
median(rlog1Entropyb$infoGain)

############ Multivariate visualization ############
pca1 <- prcomp(scale(rlog1, scale = FALSE))
pca1Scores <- as.data.frame(pca1$x[,1:4])
pca1Scores$oldSampName <- rownames(pca1Scores)
pca1Scores <- colData1 %>% left_join(pca1Scores)
summary(pca1)

# png(filename = "Plots/PC1vsPC2_20200926.png", height = 4.5, width = 6.5, units = "in", res = 300)
set.seed(3)
ggplot(pca1Scores, aes(x = PC1, y = PC2, color = pheno, label = oldSampName)) + geom_point() +
  geom_text_repel(size = 2) + labs(x = "PC1 (42.0% of variance)", y = "PC2 (15.4% of variance)")
# dev.off()

# png(filename = "Plots/PC3vsPC4_20200926.png", height = 4.5, width = 6.5, units = "in", res = 300)
set.seed(3)
ggplot(pca1Scores, aes(x = PC3, y = PC4, color = pheno, label = oldSampName)) + geom_point() +
  geom_text_repel(size = 2) + labs(x = "PC3 (11.4% of variance)", y = "PC4 (8.4% of variance)")
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

# png(filename = "Plots/PC1vsPC2_sens_20200926.png", height = 4.5, width = 6.5, units = "in", res = 300)
set.seed(3)
ggplot(pca2Scores, aes(x = PC1, y = PC2, color = pheno, label = oldSampName)) + geom_point() +
  geom_text_repel(size = 2) + labs(x = "PC1 (33.0% of variance)", y = "PC2 (22.1% of variance)")
# dev.off()

# png(filename = "Plots/PC3vsPC4_sens_20200926.png", height = 4.5, width = 6.5, units = "in", res = 300)
set.seed(3)
ggplot(pca2Scores, aes(x = PC3, y = PC4, color = pheno, label = oldSampName)) + geom_point() +
  geom_text_repel(size = 2) + labs(x = "PC3 (13.2% of variance)", y = "PC4 (6.2% of variance)")
# dev.off()

# png(filename = "Plots/hclust_sens.png", height = 6.5, width = 6.5, units = "in", res = 300)
plot(hclust(dist(rlog1b), method = "ward.D2"))
# dev.off()

rm(pca1, pca2, pca1Scores, pca2Scores, PCA3D, colEntropies, colEntropiesb)

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
dynamicMods <- dynamicTreeCut::cutreeDynamic(dendro = tree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = TRUE,
                           minClusterSize = 20, method = "hybrid", cutHeight = .9)
dynamicColors <- WGCNA::labels2colors(dynamicMods)

# Plot dendogram and module assignment:
# png(filename = "Plots/WGCNA_Dendro1.png",height = 6, width = 10, units = "in", res = 600)
WGCNA::plotDendroAndColors(tree, cbind(dynamicColors), "Module Assignment", dendroLabels = FALSE, hang = 0.0, 
                           addGuide = TRUE, guideHang = 0.0, main = "Gene Modules")
# dev.off()

# Plot distnace matrix heatmap w/ module assignment
# png(filename = "Plots/WGCNA_Heatmap1.png", height = 6, width = 6, units = "in", res = 600)
WGCNA::TOMplot(dissTOM, tree, dynamicColors, main = "Module heatmap")
# dev.off()

# Module-gene mapping:
modDF <- data.frame(gene = rownames(TOM), module = dynamicColors)

mEigen1 <- WGCNA::moduleEigengenes(rlog2b, dynamicColors, impute = FALSE, nPC = 1, align = "along average", 
                                   excludeGrey = TRUE, grey = if (is.numeric(colors)) 0 else "grey", 
                                   softPower = 9, scale = TRUE, verbose = 5, indent = 1)

writexl::write_xlsx(modDF, path = paste0("Results/WGCNA_Modules_", gsub("-", "", Sys.Date()), ".xlsx"))
save.image(file = "Data/running_20200612.RData")
load(file = "Data/running_20200612.RData")

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
  fgsea1 <- fgsea::fgsea(pathways = MElist, stats = comp2, nperm = 10000)
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

# Plot:
png(filename = "Plots/WGCNA_GSEARes.png",height = 6, width = 10, units = "in", res = 600)
ggplot(fgseaRes %>% filter(maxLogP > 1.30103), aes(x = Comparison, y = -log10(padj), fill = Comparison)) + 
  geom_bar(color = "black", position = "dodge", stat = "identity", width = .7) +
  facet_wrap(~pathway, ncol = 7) + theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = wesanderson::wes_palette("Royal1", 7, type = "continuous"))
dev.off()

comp2 <- contrastDF2[contrastDF2$Comparison == "Ser Inf vs Inj Ctrl", c("gene", "stat")]
comp2$stat <- abs(comp2$stat)
comp2 <- comp2 %>% deframe()

png(filename = "Plots/WGCNA_GSEA_GreenEnrich.png",height = 5, width = 7, units = "in", res = 600)
fgsea::plotEnrichment(MElist[["green"]], comp2) + labs(title = "Serratia Infection vs Injury Control: Green Module")
dev.off()

############ Module eigengenes ############
# Make into long data.frame:
mEigen2 <- mEigen1$eigengenes
mEigen2$oldSampName <- rownames(mEigen2)
mEigen2 <- mEigen2 %>% gather(key = "ME", value = "value", -oldSampName)
mEigen2 <- colData1b %>% left_join(mEigen2)

# List of all module colors:
mEColors <- unique(mEigen2$ME)

# ME factor variable:
mEigen2$pheno <- factor(mEigen2$pheno, levels = c("Naive","Enterobacter_Priming", "Serratia_Priming",
                  "Injury_Control", "Enterobacter_infection", "Serratia_Infection", "Ent_Prim_p_Ent_inf",
                  "Ser_Prim_p_Ser_inf", "Ent_Prim_p_Ser_inf", "Ser_Prim_p_Ent_inf"))
levels(mEigen2$pheno)[levels(mEigen2$pheno) == "Enterobacter_infection"] <- "Enterobacter_Infection"

# ANOVA of modules Eigengenes by group:
MElmDF <- data.frame(ME = mEColors, overallP = NA, primingVSnaive = NA, infectionVSinjury = NA)
for(i in 1:nrow(MElmDF)){
  lm0 <- lm(value ~ 1, data = mEigen2 %>% filter(ME == MElmDF$ME[i]))
  lm1 <- lm(value ~ pheno, data = mEigen2 %>% filter(ME == MElmDF$ME[i]))
  MElmDF$overallP[i] <- anova(lm0, lm1)$`Pr(>F)`[2]
  lsdTable <- as.data.frame(DescTools::PostHocTest(aov(lm1), method = "lsd")$pheno)
  lsdTable$compare <- rownames(lsdTable)
  lsdTable$pheno1 <- str_split(lsdTable$compare, "-", simplify = TRUE)[,1]
  lsdTable$pheno2 <- str_split(lsdTable$compare, "-", simplify = TRUE)[,2]
  MElmDF$primingVSnaive[i] <- min(lsdTable$pval[lsdTable$pheno2 == "Naive" & grepl("Priming", lsdTable$pheno1)])
  MElmDF$infectionVSinjury[i] <- min(lsdTable$pval[lsdTable$pheno2 == "Injury_Control" & grepl("Infection", lsdTable$pheno1)])
  
}
ggplot(mEigen2 %>% filter(ME == "MEsaddlebrown"), aes(x = pheno, y = value)) + geom_point()
ggplot(mEigen2 %>% filter(ME == "MEroyalblue"), aes(x = pheno, y = value)) + geom_point()
ggplot(mEigen2 %>% filter(ME == "MEgreen"), aes(x = pheno, y = value)) + geom_point()

colData1b$pheno <- factor(colData1b$pheno, levels = c("Naive","Enterobacter_Priming", "Serratia_Priming",
                                                      "Injury_Control", "Enterobacter_infection", "Serratia_Infection", "Ent_Prim_p_Ent_inf",
                                                      "Ser_Prim_p_Ser_inf", "Ent_Prim_p_Ser_inf", "Ser_Prim_p_Ent_inf"))
colData1b <- colData1b %>% arrange(pheno, rep)
colData1b$ord <- 1:nrow(colData1b)

# orange module:
orangeTOMhClust <- hclust(as.dist(TOM[modDF$gene[modDF$module == "orange"], modDF$gene[modDF$module == "orange"]]),
                         method = "average")
orangeExpression <- scale(rlog2b[match(colData1b$oldSampName, rownames(rlog2b)), 
                                modDF$gene[modDF$module == "orange"]])
orangeExpression <- orangeExpression[grepl("Injury|Infection|infection", rownames(orangeExpression)),]
png(filename = "Plots/WGCNA_Module_orange.png", height = 6, width = 12, units = "in", res = 600)
orangeColors <- WGCNA::numbers2colors(orangeExpression, signed = TRUE, commonLim = FALSE)
WGCNA::plotDendroAndColors(dendro = orangeTOMhClust, colors = t(orangeColors),
                           groupLabels = rownames(orangeExpression),
                           cex.dendroLabels = 0.75, main = "")
dev.off()

# Green module:
greenTOMhClust <- hclust(as.dist(TOM[modDF$gene[modDF$module == "green"], modDF$gene[modDF$module == "green"]]),
                         method = "average")
greenExpression <- scale(rlog2b[match(colData1b$oldSampName, rownames(rlog2b)), 
                                modDF$gene[modDF$module == "green"]])
greenExpression <- greenExpression[grepl("Injury|Infection|infection", rownames(greenExpression)),]
png(filename = "Plots/WGCNA_Module_Green.png", height = 6, width = 18, units = "in", res = 600)
greenColors <- WGCNA::numbers2colors(greenExpression, signed = TRUE, commonLim = FALSE)
WGCNA::plotDendroAndColors(dendro = greenTOMhClust, colors = t(greenColors),
                           groupLabels = rownames(greenExpression),
                           cex.dendroLabels = 0.25, main = "")
dev.off()

# cyan module:
cyanTOMhClust <- hclust(as.dist(TOM[modDF$gene[modDF$module == "cyan"], modDF$gene[modDF$module == "cyan"]]),
                         method = "average")
cyanExpression <- scale(rlog2b[match(colData1b$oldSampName, rownames(rlog2b)), 
                                modDF$gene[modDF$module == "cyan"]])
cyanExpression <- cyanExpression[grepl("Injury|Infection|infection", rownames(cyanExpression)),]
png(filename = "Plots/WGCNA_Module_cyan.png", height = 6, width = 12, units = "in", res = 600)
cyanColors <- WGCNA::numbers2colors(cyanExpression, signed = TRUE, commonLim = FALSE)
WGCNA::plotDendroAndColors(dendro = cyanTOMhClust, colors = t(cyanColors),
                           groupLabels = rownames(cyanExpression),
                           cex.dendroLabels = 0.4, main = "")
dev.off()

# Saddle brown module:
saddlebrownTOMhClust <- hclust(as.dist(TOM[modDF$gene[modDF$module == "saddlebrown"], modDF$gene[modDF$module == "saddlebrown"]]),
                         method = "average")
saddlebrownExpression <- scale(rlog2b[match(colData1b$oldSampName, rownames(rlog2b)), 
                                modDF$gene[modDF$module == "saddlebrown"]])
saddlebrownExpression <- saddlebrownExpression[grepl("Injury|Infection|infection", rownames(saddlebrownExpression)),]
png(filename = "Plots/WGCNA_Module_saddlebrown.png", height = 6, width = 12, units = "in", res = 600)
saddlebrownColors <- WGCNA::numbers2colors(saddlebrownExpression, signed = TRUE, commonLim = FALSE)
WGCNA::plotDendroAndColors(dendro = saddlebrownTOMhClust, colors = t(saddlebrownColors),
                           groupLabels = rownames(saddlebrownExpression),
                           cex.dendroLabels = 0.75, main = "")
dev.off()

# Royal blue module:
royalblueTOMhClust <- hclust(as.dist(TOM[modDF$gene[modDF$module == "royalblue"], modDF$gene[modDF$module == "royalblue"]]),
                         method = "average")
royalblueExpression <- scale(rlog2b[match(colData1b$oldSampName, rownames(rlog2b)), 
                                modDF$gene[modDF$module == "royalblue"]])
royalblueExpression <- royalblueExpression[grepl("Naive|Priming", rownames(royalblueExpression)),]
png(filename = "Plots/WGCNA_Module_RoyalBlue.png", height = 6, width = 18, units = "in", res = 600)
royalblueColors <- WGCNA::numbers2colors(royalblueExpression, signed = TRUE, commonLim = FALSE)
WGCNA::plotDendroAndColors(dendro = royalblueTOMhClust, colors = t(royalblueColors),
                           groupLabels = rownames(royalblueExpression),
                           cex.dendroLabels = 0.75, main = "")
dev.off()

# White module:
whiteTOMhClust <- hclust(as.dist(TOM[modDF$gene[modDF$module == "white"], modDF$gene[modDF$module == "white"]]),
                             method = "average")
whiteExpression <- scale(rlog2b[match(colData1b$oldSampName, rownames(rlog2b)), 
                                    modDF$gene[modDF$module == "white"]])
whiteExpression <- whiteExpression[grepl("Naive|Priming|Injury|Infection|infection", rownames(whiteExpression)),]
png(filename = "Plots/WGCNA_Module_white.png", height = 6, width = 12, units = "in", res = 600)
whiteColors <- WGCNA::numbers2colors(whiteExpression, signed = TRUE, commonLim = FALSE)
WGCNA::plotDendroAndColors(dendro = whiteTOMhClust, colors = t(whiteColors),
                           groupLabels = rownames(whiteExpression),
                           cex.dendroLabels = 0.75, main = "")
dev.off()

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
