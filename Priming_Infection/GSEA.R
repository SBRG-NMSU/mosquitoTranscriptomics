##### GSEA Analysis for NMSU Mosquito project #####
  # Analysis by: Samantha Carlisle
  # Date: 6/6/20
  # Species: Anopheles gambiae

#### Loading libraries #####
library("biomaRt")
library("dplyr")
library("tidyr")
library("DESeq2")
library("stringr")
library("fgsea")
library("tibble")
library("ggplot2")
library("KEGGREST")

setwd("/Users/sambot/Desktop/Las_Cruces/Xu_mosquito_transcriptomics/")

# Note: To get the GO terms for GSEA we need to get the mappings per gene from biomaRt
# but biomaRt does not have the NCBI gene symbol annotations that were used for
# annotating the genes in the count and thus deseq2 results. We are using ensembl ids
# because they are the ids that have the most overlap between the deseq2 genes and
# the biomaRt genes. The NCBI mapping data didn't have many entrez gene ids that 
# mapped to the entrez gene ids in the biomaRt data. We used a combination of gene symbol,
# uniprot ids, and gene synonyms for mapping deseq2 data to biomart data. Used KEGGREST 
# for mapping KEGG IDs to gene ids-- there were 4975 genes mapped to KEGG ids.



##### Loading data and transforming #####
# Loading DESeq2 results from Patrick
load("contrastDF2_20200614.RData")
# Removing the leading AgaP_ from the gene names to get ensembl ids
# contrastDF2$gene <- gsub("AgaP_", "", contrastDF2$gene)

# Loading uniprot accession to id mapping file; retrieved from https://www.uniprot.org/uniprot/?query=taxonomy:7165 on 6/6/20
uniprot_taxonomy <- read.csv("uniprot_taxonomy.csv", stringsAsFactors=F)
# Isolating only the ensembl ids from the Gene.name column
uniprot_taxonomy$Gene.names <- str_extract(uniprot_taxonomy$Gene.names, "AGAP.*$")
# removing rows where ensembl id is NA
uniprot_taxonomy <- uniprot_taxonomy[!is.na(uniprot_taxonomy$Gene.names),]

# Loading NCBI mapping table for A. gambiae from https://www.ncbi.nlm.nih.gov/gene on 6/6/20;
# has entrez gene id and ncbi symbol
ncbi <- read.delim("ncbi_agambiae_gene_mappings.txt", stringsAsFactors = F)



###### Getting data (ensembl ids, entrez ids, uniprot ids, go terms) from BioMart #######
# First list available biomarts
listEnsemblGenomes()
# We need the metazoa mart
ensembl_metazoa <- useEnsemblGenomes(biomart = "metazoa_mart")
# Search the metazoa mart for the specific mosquito species used
searchDatasets(ensembl_metazoa, pattern = "Anopheles gambiae")
# Build mart for that species
ensembl_agambiae<- useEnsemblGenomes(biomart = "metazoa_mart", 
                                         dataset = "agambiae_eg_gene")

# Getting info about agambiae ensembl biomart (available data fields)
filters_agambiae <- listFilters(ensembl_agambiae)
attributes_agambiae <- listAttributes(ensembl_agambiae)

# Get data for agambiae
test <- getBM(attributes = c("uniprotsptrembl", "external_gene_name", "external_synonym", 'ensembl_gene_id', 'entrezgene_id', "description", "go_id",
                             "name_1006", "definition_1006", "namespace_1003"),
              #filters = 'external_gene_name',
              #values = unique(contrastDF2$gene), 
              mart = ensembl_agambiae)

# Only getting ensembl and entrez ids for merging with NCBI mapping file
test11 <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id'),
              #filters = 'external_gene_name',
              #values = unique(contrastDF2$gene), 
              mart = ensembl_agambiae)

test22 <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', "kegg", "kegg_enzyme"),
                                #filters = 'external_gene_name',
                                   #values = unique(contrastDF2$gene), 
                                 mart = ensembl_agambiae)



#### Merging ncbi and biomaRt mapping files ####
# merge entrez ids, ncbi symbols, ensembl ids
entrez_ensembl_ncbi_symbol <- full_join(dplyr::select(ncbi, GeneID, Symbol), 
                                             dplyr::select(test11, entrezgene_id, ensembl_gene_id),  
                                             by=c("GeneID"="entrezgene_id"))
# removing rows where Symbol is NA
entrez_ensembl_ncbi_symbol <- entrez_ensembl_ncbi_symbol[!is.na(entrez_ensembl_ncbi_symbol$Symbol) & !is.na(entrez_ensembl_ncbi_symbol$ensembl_gene_id),]

# merge the mapping file with the contrast file
contrast_ensembl <- merge(contrastDF2, dplyr::select(entrez_ensembl_ncbi_symbol, Symbol, ensembl_gene_id), by.x="gene", by.y="Symbol")
 #only 2532 genes remaining; this approach won't work



##### Converting gene symbols from DESeq2 results to ensembl ids (have greatest overlap bt our datasets) #####
  # first priming the new ensembl id column in the contrast file with "NA"
contrastDF2$ensembl_id <- "NA"

for (gene in contrastDF2$gene[grep("^AGAP", contrastDF2$gene, invert=T)]) {
  print(gene)
  ifelse(gene %in% toupper(test$external_gene_name), id <- test$ensembl_gene_id[toupper(test$external_gene_name)==gene], 
         ifelse(gene %in% test$external_synonym, id <- test$ensembl_gene_id[test$external_synonym==gene], 
                ifelse(gene %in% uniprot_taxonomy$Entry.name, id <- uniprot_taxonomy$Gene.names[uniprot_taxonomy$Entry.name==gene], 
                      ifelse(gene %in% entrez_ensembl_ncbi_symbol$Symbol, id <- entrez_ensembl_ncbi_symbol$ensembl_gene_id[entrez_ensembl_ncbi_symbol$Symbol==gene], 
                             id <- "NA"))))
  contrastDF2$ensembl_id[contrastDF2$gene==gene] <- unique(id)
}

for (gene in contrastDF2$gene[grep("^AGAP", contrastDF2$gene)]) {
  print(gene)
  contrastDF2$ensembl_id[contrastDF2$gene==gene] <- gene
}

#write.csv(contrastDF2, "contrastDF2_ensembl_ids2.csv", row.names=F)
contrast <- read.csv("contrastDF2_ensembl_ids2.csv", stringsAsFactors=F)

# Getting gene names from contrast file that didn't map to an ensembl id
test4 <- contrastDF2[contrastDF2$ensembl_id == "NA",]
    # 32 genes without ensembl ids
#write.csv(unique(test4$gene), "unmapped2.csv", row.names = F)

# Create a module list for GSEA based on GO terms from biomart
gene_list <- list()
for (id in unique(filter(test, go_id != "")$go_id)) {
  print(id)
  id2 <- paste(id, unique(test$name_1006[test$go_id == id]))
  df <- filter(dplyr::select(test, go_id, ensembl_gene_id), go_id == id)
  gene_list[[id2]] <- df$ensembl_gene_id
}
db_ensembl_ids_wo_go_id <- unique(filter(test, go_id == "")$ensembl_gene_id)
db_ensembl_ids_w_go_id <- unique(filter(test, go_id !="")$ensembl_gene_id)
contrast_ensembl_ids_wo_go_id <- sum(unique(contrastDF2$ensembl_id) %in% db_ensembl_ids_wo_go_id)
contrast_ensembl_ids_w_go_id <- sum(unique(contrastDF2$ensembl_id) %in% db_ensembl_ids_w_go_id)
genes_wo_ensembl_ids <- unique(filter(contrastDF2, ensembl_id=="NA")$gene)

# Creating sperate lists of GO terms based on category
GO_cellular_component <- data.frame("go_id"=dplyr::select(filter(test, namespace_1003 == "cellular_component"), go_id),
                               "pathname"=dplyr::select(filter(test, namespace_1003 == "cellular_component"), name_1006))
GO_cellular_component <- unique(paste(GO_cellular_component$go_id, GO_cellular_component$name_1006, sep=" "))

GO_molecular_function <- data.frame("go_id"=dplyr::select(filter(test, namespace_1003 == "molecular_function"), go_id),
                                    "pathname"=dplyr::select(filter(test, namespace_1003 == "molecular_function"), name_1006))
GO_molecular_function <- unique(paste(GO_molecular_function$go_id, GO_molecular_function$name_1006, sep=" "))

GO_biological_process <- data.frame("go_id"=dplyr::select(filter(test, namespace_1003 == "biological_process"), go_id),
                                    "pathname"=dplyr::select(filter(test, namespace_1003 == "biological_process"), name_1006))
GO_biological_process <- unique(paste(GO_biological_process$go_id, GO_biological_process$name_1006, sep=" "))

# Count of genes in each go_term list
#gene_count <- test %>% group_by(namespace_1003) %>% count(go_id)

#### Seperate contrast file into the different comparisons ####
enter_vs_naive_contrast <- filter(contrastDF2, Comparison == "Enterobacter vs. Naive")
serr_vs_naive_contrast <- filter(contrastDF2, Comparison == "Serratia vs. Naive")
serr_vs_enter_contrast <- filter(contrastDF2, Comparison == "Serratia vs. Enterobacter")
enter_inf_vs_inj_ctrl_contrast <- filter(contrastDF2, Comparison == "Ent Inf vs Inj Ctrl")
serr_inf_vs_inj_ctrl_contrast <- filter(contrastDF2, Comparison == "Ser Inf vs Inj Ctrl")
enter_prim_inf_vs_enter_inf_contrast <- filter(contrastDF2, Comparison == "Ent Prim & Inf vs Ent Inf")
serr_prim_inf_vs_serr_inf_contrast <- filter(contrastDF2, Comparison == "Ser Prim & Inf vs Ser Inf")

enter_vs_naive_contrast$stat <- abs(enter_vs_naive_contrast$stat)
serr_vs_naive_contrast$stat <- abs(serr_vs_naive_contrast$stat)
serr_vs_enter_contrast$stat <- abs(serr_vs_enter_contrast$stat)
enter_inf_vs_inj_ctrl_contrast$stat <- abs(enter_inf_vs_inj_ctrl_contrast$stat)
serr_inf_vs_inj_ctrl_contrast$stat <- abs(serr_inf_vs_inj_ctrl_contrast$stat)
enter_prim_inf_vs_enter_inf_contrast$stat <- abs(enter_prim_inf_vs_enter_inf_contrast$stat)
serr_prim_inf_vs_serr_inf_contrast$stat <- abs(serr_prim_inf_vs_serr_inf_contrast$stat)


# Create a preranked list of genes
ranked_gene_list1 <- enter_vs_naive_contrast %>% 
  dplyr::select(ensembl_id, stat) %>% 
  filter(ensembl_id != "NA") %>% 
  distinct() %>% 
  group_by(ensembl_id) %>% 
  summarize(stat=mean(stat)) %>%
  arrange(desc(stat))

ranked_gene_list1 <- deframe(ranked_gene_list1)

ranked_gene_list2 <- serr_vs_naive_contrast %>% 
  dplyr::select(ensembl_id, stat) %>% 
  filter(ensembl_id != "NA") %>% 
  distinct() %>% 
  group_by(ensembl_id) %>% 
  summarize(stat=mean(stat)) %>%
  arrange(desc(stat))

ranked_gene_list2 <- deframe(ranked_gene_list2)

ranked_gene_list3 <- serr_vs_enter_contrast %>% 
  dplyr::select(ensembl_id, stat) %>% 
  filter(ensembl_id != "NA") %>% 
  distinct() %>% 
  group_by(ensembl_id) %>% 
  summarize(stat=mean(stat)) %>%
  arrange(desc(stat))

ranked_gene_list3 <- deframe(ranked_gene_list3)

ranked_gene_list4 <- enter_inf_vs_inj_ctrl_contrast %>% 
  dplyr::select(ensembl_id, stat) %>% 
  filter(ensembl_id != "NA") %>% 
  distinct() %>% 
  group_by(ensembl_id) %>% 
  summarize(stat=mean(stat)) %>%
  arrange(desc(stat))

ranked_gene_list4 <- deframe(ranked_gene_list4)

ranked_gene_list5 <- serr_inf_vs_inj_ctrl_contrast %>% 
  dplyr::select(ensembl_id, stat) %>% 
  filter(ensembl_id != "NA") %>% 
  distinct() %>% 
  group_by(ensembl_id) %>% 
  summarize(stat=mean(stat)) %>%
  arrange(desc(stat))

ranked_gene_list5 <- deframe(ranked_gene_list5)

ranked_gene_list6 <- enter_prim_inf_vs_enter_inf_contrast %>% 
  dplyr::select(ensembl_id, stat) %>% 
  filter(ensembl_id != "NA") %>% 
  distinct() %>% 
  group_by(ensembl_id) %>% 
  summarize(stat=mean(stat)) %>%
  arrange(desc(stat))

ranked_gene_list6 <- deframe(ranked_gene_list6)

ranked_gene_list7 <- serr_prim_inf_vs_serr_inf_contrast %>% 
  dplyr::select(ensembl_id, stat) %>% 
  filter(ensembl_id != "NA") %>% 
  distinct() %>% 
  group_by(ensembl_id) %>% 
  summarize(stat=mean(stat)) %>%
  arrange(desc(stat))

ranked_gene_list7 <- deframe(ranked_gene_list7)


###### GSEA #######
fgseaRes_enter_vs_naive <- fgsea(pathways=gene_list, stats=ranked_gene_list1)
fgseaRes_serr_vs_naive <- fgsea(pathways=gene_list, stats=ranked_gene_list2)
fgseaRes_serr_vs_enter <- fgsea(pathways=gene_list, stats=ranked_gene_list3)
fgseaRes_enter_inf_vs_inj_ctrl <- fgsea(pathways=gene_list, stats=ranked_gene_list4)
fgseaRes_serr_inf_vs_inj_ctrl <- fgsea(pathways=gene_list, stats=ranked_gene_list5)
fgseaRes_enter_prim_inf_vs_enter_inf <- fgsea(pathways=gene_list, stats=ranked_gene_list6)
fgseaRes_serr_prim_inf_vs_serr_inf <- fgsea(pathways=gene_list, stats=ranked_gene_list7)

# write.csv(fgseaRes_enter_vs_naive[,1:7], "fgseaRes_enter_vs_naive_GO.csv", row.names=F)
# write.csv(fgseaRes_serr_vs_naive[,1:7], "fgseaRes_serr_vs_naive_GO.csv", row.names=F)
# write.csv(fgseaRes_serr_vs_enter[,1:7], "fgseaRes_serr_vs_enter_GO.csv", row.names=F)
# write.csv(fgseaRes_enter_inf_vs_inj_ctrl[,1:7], "fgseaRes_enter_inf_vs_inj_ctrl_GO.csv", row.names=F)
# write.csv(fgseaRes_serr_inf_vs_inj_ctrl[,1:7], "fgseaRes_serr_inf_vs_inj_ctrl_GO.csv", row.names=F)
# write.csv(fgseaRes_enter_prim_inf_vs_enter_inf[,1:7], "fgseaRes_enter_prim_inf_vs_enter_inf_GO.csv", row.names=F)
# write.csv(fgseaRes_serr_prim_inf_vs_serr_inf[,1:7], "fgseaRes_serr_prim_inf_vs_serr_inf_GO.csv", row.names=F)


###### Plots #######
#png("enter_vs_naive_GO_GSEA.png", width = 600, height = 400)
a <- ggplot(top_n(filter(fgseaRes_enter_vs_naive, padj<0.05 & pathway %in% GO_biological_process), -20, padj), aes(x=reorder(pathway, -padj), y=-log10(padj), color=NES)) +
  geom_point(aes(size=size)) + 
  scale_colour_gradient(low="green",high="red") +
  #scale_size_area(limits = c(0,225)) +
  scale_y_continuous(limits=c(0,7), expand=c(0,0)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=-log10(padj))) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) +
  labs(x="Pathway", y="-log10(q-value)",
       title="Biological Process GO terms Enrichment from GSEA \n(Enterobacter vs. Naive)") 
#dev.off()
b <- ggplot(top_n(filter(fgseaRes_enter_vs_naive, padj<0.05 & pathway %in% GO_cellular_component), -20, padj), aes(x=reorder(pathway, -padj), y=-log10(padj), color=NES)) +
  geom_point(aes(size=size)) + 
  scale_colour_gradient(low="green",high="red") +
  #scale_size_area(limits = c(0,225)) +
  scale_y_continuous(limits=c(0,7), expand=c(0,0)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=-log10(padj))) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) +
  labs(x="Pathway", y="-log10(q-value)",
       title="Cellular Component GO terms Enrichment from GSEA \n(Enterobacter vs. Naive)") 
c <- ggplot(top_n(filter(fgseaRes_enter_vs_naive, padj<0.05 & pathway %in% GO_molecular_function), -20, padj), aes(x=reorder(pathway, -padj), y=-log10(padj), color=NES)) +
  geom_point(aes(size=size)) + 
  scale_colour_gradient(low="green",high="red") +
  #scale_size_area(limits = c(0,225)) +
  scale_y_continuous(limits=c(0,7), expand=c(0,0)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=-log10(padj))) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) +
  labs(x="Pathway", y="-log10(q-value)",
       title="Molecular Function GO terms Enrichment from GSEA \n(Enterobacter vs. Naive)") 

ggplot(top_n(filter(fgseaRes_enter_inf_vs_inj_ctrl, padj<0.05), -20, padj), aes(x=reorder(pathway, -padj), y=-log10(padj), color=NES)) +
  geom_point(aes(size=size)) + 
  scale_colour_gradient(low="green",high="red") +
  #scale_size_area(limits = c(0,225)) +
  scale_y_continuous(limits=c(0,2), expand=c(0,0)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=-log10(padj))) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) +
  labs(x="Pathway", y="-log10(q-value)",
       title="GO terms Enrichment from GSEA \n(Ent Inf vs Inj Ctrl)") 

ggplot(top_n(filter(fgseaRes_serr_inf_vs_inj_ctrl, padj<0.05), -20, padj), aes(x=reorder(pathway, -padj), y=-log10(padj), color=NES)) +
  geom_point(aes(size=size)) + 
  scale_colour_gradient(low="green",high="red") +
  #scale_size_area(limits = c(0,225)) +
  scale_y_continuous(limits=c(0,6.7), expand=c(0,0)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=-log10(padj))) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) +
  labs(x="Pathway", y="-log10(q-value)",
       title="GO terms Enrichment from GSEA \n(Ser Inf vs Inj Ctrl)") 

ggplot(top_n(filter(fgseaRes_enter_prim_inf_vs_enter_inf, padj<0.05), -20, padj), aes(x=reorder(pathway, -padj), y=-log10(padj), color=NES)) +
  geom_point(aes(size=size)) + 
  scale_colour_gradient(low="green",high="red") +
  #scale_size_area(limits = c(0,225)) +
  scale_y_continuous(limits=c(0,3.7), expand=c(0,0)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=-log10(padj))) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) +
  labs(x="Pathway", y="-log10(q-value)",
       title="GO terms Enrichment from GSEA \n(Ent Prim & Inf vs Ent Inf)") 

ggplot(top_n(filter(fgseaRes_serr_prim_inf_vs_serr_inf, padj<0.05), -20, padj), aes(x=reorder(pathway, -padj), y=-log10(padj), color=NES)) +
  geom_point(aes(size=size)) + 
  scale_colour_gradient(low="green",high="red") +
  #scale_size_area(limits = c(0,225)) +
  scale_y_continuous(limits=c(0,1.6), expand=c(0,0)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=-log10(padj))) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) +
  labs(x="Pathway", y="-log10(q-value)",
       title="GO terms Enrichment from GSEA \n(Ser Prim & Inf vs Ser Inf)") 



#png("serr_vs_naive_GO_GSEA.png", width = 600, height = 400)
d <- ggplot(top_n(filter(fgseaRes_serr_vs_naive, padj<0.05 & pathway %in% GO_biological_process), -20, padj), aes(x=reorder(pathway, -padj), y=-log10(padj), color=NES)) +
  geom_point(aes(size=size)) + 
  scale_colour_gradient(low="green",high="red") +
  #scale_size_area(limits = c(0,225)) +
  scale_y_continuous(limits = c(0,6.5), expand=c(0,0)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=-log10(padj))) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) +
  labs(x="Pathway", y="-log10(q-value)",
       title="Biological Process GO terms Enrichment from GSEA \n(Serratia vs. Naive)") 
e <- ggplot(top_n(filter(fgseaRes_serr_vs_naive, padj<0.05 & pathway %in% GO_cellular_component), -20, padj), aes(x=reorder(pathway, -padj), y=-log10(padj), color=NES)) +
  geom_point(aes(size=size)) + 
  scale_colour_gradient(low="green",high="red") +
  #scale_size_area(limits = c(0,225)) +
  scale_y_continuous(limits = c(0,6.5), expand=c(0,0)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=-log10(padj))) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) +
  labs(x="Pathway", y="-log10(q-value)",
       title="Cellular Component GO terms Enrichment from GSEA \n(Serratia vs. Naive)") 
f <- ggplot(top_n(filter(fgseaRes_serr_vs_naive, padj<0.05 & pathway %in% GO_molecular_function), -20, padj), aes(x=reorder(pathway, -padj), y=-log10(padj), color=NES)) +
  geom_point(aes(size=size)) + 
  scale_colour_gradient(low="green",high="red") +
  #scale_size_area(limits = c(0,225)) +
  scale_y_continuous(limits = c(0,6.5), expand=c(0,0)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=-log10(padj))) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) +
  labs(x="Pathway", y="-log10(q-value)",
       title="Molecular Function GO terms Enrichment from GSEA \n(Serratia vs. Naive)") 
#dev.off()

png("serr_vs_enter_GO_GSEA.png", width = 600, height = 400)
g <- ggplot(top_n(filter(fgseaRes_serr_vs_enter, padj<0.05 & pathway %in% GO_biological_process), -20, padj), aes(x=reorder(pathway, -padj), y=-log10(padj), color=NES)) +
  geom_point(aes(size=size)) + 
  scale_colour_gradient(low="green",high="red") +
  #scale_size_area(limits = c(0,225)) +
  scale_y_continuous(limits = c(0,1.5), expand=c(0,0)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=-log10(padj))) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) +
  labs(x="Pathway", y="-log10(q-value)",
       title="GO terms Enrichment from GSEA \n(Serratia vs. Enterobacter)") 
h <- ggplot(top_n(filter(fgseaRes_serr_vs_enter, padj<0.05 & pathway %in% GO_cellular_component), -20, padj), aes(x=reorder(pathway, -padj), y=-log10(padj), color=NES)) +
  geom_point(aes(size=size)) + 
  scale_colour_gradient(low="green",high="red") +
  #scale_size_area(limits = c(0,225)) +
  scale_y_continuous(limits = c(0,1.5), expand=c(0,0)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=-log10(padj))) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) +
  labs(x="Pathway", y="-log10(q-value)",
       title="GO terms Enrichment from GSEA \n(Serratia vs. Enterobacter)") 
i <- ggplot(top_n(filter(fgseaRes_serr_vs_enter, padj<0.05 & pathway %in% GO_molecular_function), -20, padj), aes(x=reorder(pathway, -padj), y=-log10(padj), color=NES)) +
  geom_point(aes(size=size)) + 
  scale_colour_gradient(low="green",high="red") +
  #scale_size_area(limits = c(0,225)) +
  scale_y_continuous(limits = c(0,1.5), expand=c(0,0)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=-log10(padj))) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) +
  labs(x="Pathway", y="-log10(q-value)",
       title="GO terms Enrichment from GSEA \n(Serratia vs. Enterobacter)") 
dev.off()



# List of gene symbols in contrast that are not ensembl ids
gene_symbols <- data.frame("gene"=(contrastDF2$gene[grep("^[^AGAP]", contrastDF2$gene)]), stringsAsFactors=F)
# 2124 genes



# Using table downloaded from kegg to match kegg ids and gene ids
test33 <- read.csv("mosquito_kegg.csv", header=F)
test33$gene_id <- gsub("K.*|no| |aga:|AgaP_", "", test33$V1)
test33$kegg_id <- str_extract(test33$V1, "K.....")
test33$kegg_id <- gsub("KO ass", NA, test33$kegg_id)
test44 <- merge(contrast, dplyr::select(test33, kegg_id, gene_id), by.x="ensembl_id", by.y="gene_id")

### Using KeggRest
org <- keggList("organism")
org[grep("T01036",org),]

kegg_ids <- keggLink("T01036", "ko")
kegg_ids2 <- data.frame("kegg_id"=names(kegg_ids), "gene"=kegg_ids)
kegg_ids2$gene <- gsub("aga:|AgaP_", "", kegg_ids2$gene)
test55 <- merge(contrast, kegg_ids2, by.x="ensembl_id", by.y="gene")

length(unique(test55$kegg_id))
sum(unique(filter(test55, kegg_id != "")$kegg_id) %in% unique(contrastDF2$ensembl_id))

kegg_pathways <- keggLink("T01036", "pathway")
kegg_pathways2 <- data.frame("kegg_pathway"=names(kegg_pathways), "gene"=kegg_pathways)
kegg_pathways2$gene <- gsub("aga:|AgaP_", "", kegg_pathways2$gene)

kegg_pathway_names <- keggList("pathway", "aga")
kegg_pathway_names2 <- data.frame("kegg_pathway_id"=names(kegg_pathway_names), "pathway_name"=kegg_pathway_names)
kegg_pathway_names2$pathway_name <- gsub(" - Anopheles gambiae \\(mosquito\\)", "", kegg_pathway_names2$pathway_name)

# Create a pathway list for GSEA based on KEGG terms from KEGGREST
pathway_list <- list()
for (id in unique(kegg_pathways2$kegg_pathway)) {
  print(id) 
  id2 <- paste(id, kegg_pathway_names2$pathway_name[kegg_pathway_names2$kegg_pathway_id == id], sep="; ")
  df <- filter(kegg_pathways2, kegg_pathway == id)
  pathway_list[[id2]] <- df$gene
}



###### GSEA #######
fgseaRes_enter_vs_naive2 <- fgsea(pathways=pathway_list, stats=ranked_gene_list1)
fgseaRes_serr_vs_naive2 <- fgsea(pathways=pathway_list, stats=ranked_gene_list2)
fgseaRes_serr_vs_enter2 <- fgsea(pathways=pathway_list, stats=ranked_gene_list3)
fgseaRes_enter_inf_vs_inj_ctrl2 <- fgsea(pathways=pathway_list, stats=ranked_gene_list4)
fgseaRes_serr_inf_vs_inj_ctrl2 <- fgsea(pathways=pathway_list, stats=ranked_gene_list5)
fgseaRes_enter_prim_inf_vs_enter_inf2 <- fgsea(pathways=pathway_list, stats=ranked_gene_list6)
fgseaRes_serr_prim_inf_vs_serr_inf2 <- fgsea(pathways=pathway_list, stats=ranked_gene_list7)

#write.csv(fgseaRes_enter_vs_naive2[,1:7], "fgseaRes_enter_vs_naive_kegg.csv", row.names=F)
#write.csv(fgseaRes_serr_vs_naive2[,1:7], "fgseaRes_serr_vs_naive_kegg.csv", row.names=F)
#write.csv(fgseaRes_serr_vs_enter2[,1:7], "fgseaRes_serr_vs_enter_kegg.csv", row.names=F)
#write.csv(fgseaRes_enter_inf_vs_inj_ctrl2[,1:7], "fgseaRes_enter_inf_vs_inj_ctrl_kegg.csv", row.names=F)
#write.csv(fgseaRes_serr_inf_vs_inj_ctrl2[,1:7], "fgseaRes_serr_inf_vs_inj_ctrl_kegg.csv", row.names=F)
#write.csv(fgseaRes_enter_prim_inf_vs_enter_inf2[,1:7], "fgseaRes_enter_prim_inf_vs_enter_inf_kegg.csv", row.names=F)
#write.csv(fgseaRes_serr_prim_inf_vs_serr_inf2[,1:7], "fgseaRes_serr_prim_inf_vs_serr_inf_kegg.csv", row.names=F)


###### Plots #######
png("enter_vs_naive_KEGG_GSEA.png", width = 600, height = 400)
ggplot(top_n(filter(fgseaRes_enter_vs_naive2, padj<0.05 & pathway != "path:aga01100; Metabolic pathways"), -20, padj), aes(x=reorder(pathway, -padj), y=-log10(padj), color=NES)) +
  geom_point(aes(size=size)) + 
  scale_colour_gradient(low="green",high="red") +
  #scale_size_area(limits = c(0,225)) +
  scale_y_continuous(limits=c(0,3.4), expand=c(0,0)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=-log10(padj))) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) +
  labs(x="Pathway", y="-log10(q-value)",
       title="KEGG Pathway Enrichment from GSEA \n(Enterobacter vs. Naive)") 
dev.off()

png("serr_vs_naive_KEGG_GSEA.png", width = 600, height = 400)
ggplot(top_n(filter(fgseaRes_serr_vs_naive2, padj<0.05  & pathway != "path:aga01100; Metabolic pathways"), -20, padj), aes(x=reorder(pathway, -padj), y=-log10(padj), color=NES)) +
  geom_point(aes(size=size)) + 
  scale_colour_gradient(low="green",high="red") +
  #scale_size_area(limits = c(0,225)) +
  scale_y_continuous(limits=c(0,3), expand=c(0,0)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=-log10(padj))) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) +
  labs(x="Pathway", y="-log10(q-value)",
       title="KEGG Pathway Enrichment from GSEA \n(Serratia vs. Naive)") 
dev.off()

png("serr_vs_enter_KEGG_GSEA.png", width = 600, height = 400)
ggplot(top_n(filter(fgseaRes_serr_vs_enter2, padj<0.05 & pathway != "path:aga01100; Metabolic pathways"), -20, padj), aes(x=reorder(pathway, -padj), y=-log10(padj), color=NES)) +
  geom_point(aes(size=size)) + 
  scale_colour_gradient(low="green",high="red") +
  #scale_size_area(limits = c(0,225)) +
  scale_y_continuous(limits=c(0,2.2), expand=c(0,0)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=-log10(padj))) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) +
  labs(x="Pathway", y="-log10(q-value)",
       title="KEGG Pathway Enrichment from GSEA \n(Serratia vs. Enterobacter)") 
dev.off()

ggplot(top_n(filter(fgseaRes_enter_inf_vs_inj_ctrl2, padj<0.05 & pathway != "path:aga01100; Metabolic pathways"), -20, padj), aes(x=reorder(pathway, -padj), y=-log10(padj), color=NES)) +
  geom_point(aes(size=size)) + 
  scale_colour_gradient(low="green",high="red") +
  #scale_size_area(limits = c(0,225)) +
  scale_y_continuous(limits=c(0,6.2), expand=c(0,0)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=-log10(padj))) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) +
  labs(x="Pathway", y="-log10(q-value)",
       title="KEGG Pathway Enrichment from GSEA \n(Ent Inf vs Inj Ctrl)") 

ggplot(top_n(filter(fgseaRes_serr_inf_vs_inj_ctrl2, padj<0.05 & pathway != "path:aga01100; Metabolic pathways"), -20, padj), aes(x=reorder(pathway, -padj), y=-log10(padj), color=NES)) +
  geom_point(aes(size=size)) + 
  scale_colour_gradient(low="green",high="red") +
  #scale_size_area(limits = c(0,225)) +
  scale_y_continuous(limits=c(0,5), expand=c(0,0)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=-log10(padj))) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) +
  labs(x="Pathway", y="-log10(q-value)",
       title="KEGG Pathway Enrichment from GSEA \n(Ser Inf vs Inj Ctrl)") 

ggplot(top_n(filter(fgseaRes_enter_prim_inf_vs_enter_inf2, padj<0.05 & pathway != "path:aga01100; Metabolic pathways"), -20, padj), aes(x=reorder(pathway, -padj), y=-log10(padj), color=NES)) +
  geom_point(aes(size=size)) + 
  scale_colour_gradient(low="green",high="red") +
  #scale_size_area(limits = c(0,225)) +
  scale_y_continuous(limits=c(0,4.2), expand=c(0,0)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=-log10(padj))) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) +
  labs(x="Pathway", y="-log10(q-value)",
       title="KEGG Pathway Enrichment from GSEA \n(Ent Prim & Inf vs Ent Inf)") 

ggplot(top_n(filter(fgseaRes_serr_prim_inf_vs_serr_inf2, padj<0.05 & pathway != "path:aga01100; Metabolic pathways"), -20, padj), aes(x=reorder(pathway, -padj), y=-log10(padj), color=NES)) +
  geom_point(aes(size=size)) + 
  scale_colour_gradient(low="green",high="red") +
  #scale_size_area(limits = c(0,225)) +
  scale_y_continuous(limits=c(0,5.9), expand=c(0,0)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=-log10(padj))) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) +
  labs(x="Pathway", y="-log10(q-value)",
       title="KEGG Pathway Enrichment from GSEA \n(Ser Prim & Inf vs Ser Inf)") 

library(nord)
library(purrr)

par(mfrow=c(8, 2), lheight = 2, mar=rep(1, 4), adj = 0)

walk(names(nord_palettes), nord_show_palette)





