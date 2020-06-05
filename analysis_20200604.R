############ Prereqs ############
options(stringsAsFactors = FALSE, scipen = 900)
oldPar <- par()
os <- Sys.info()["sysname"]
baseDir <- ifelse(os == "Windows", "C:/Users/ptrainor/gdrive/Mosquitos/", 
                  "~/gdrive/Mosquitos/")
setwd(baseDir)

library(tidyverse)

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
df1 <- round(df1)

############ Univariate analysis ############
des1 <- DESeq2::DESeqDataSetFromMatrix(countData = df1, colData = colData1, design = ~ pheno)
des1 <- DESeq2::DESeq(des1)
DESeq2::resultsNames(des1)
DESeq2::results(des1, contrast = c("pheno", "Naive", "Enterobacter_Priming"))

res1 <- as.data.frame(res1)

des0 <- DESeq2::DESeq(des1, test="LRT", reduced=~1)
res0 <- as.data.frame(DESeq2::results(des0))
