
suppressPackageStartupMessages({
library(tidyverse)
library(ggrepel)
library(BiocParallel)
library(ggpubr)
library(magrittr)
library(broom)
library(data.table)
library(cowplot)
library(BiocSingular)
library(clusterProfiler)
library(enrichR)
})
rm(list=ls())

dir.create("dge/functional_enrichment/")

# ad10 vs wt10
dge <- read.table("dge/ad10vswt10_Dge_HumanID.txt",header=T)
l <- split(dge, dge$Direction)
l <- l[sapply(l, nrow)>0] #remove objects with less than one gene

GOI <- list()
GeneOnto <- list()

for(i in 1:length(l)){
GOI[[i]] <- bitr(as.character(l[[i]]$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Hs.eg.db::org.Hs.eg.db)

GeneOnto[[i]] <- enrichGO(gene = unique(GOI[[i]]$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 0.2, 
                     qvalueCutoff = 0.2, 
                     readable = TRUE)

openxlsx::write.xlsx(as.data.frame(GeneOnto[[i]]), 
                     file = sprintf("dge/functional_enrichment/%s_ad10vswt10_GO.xlsx", names(l)[[i]]), 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)


PLOT <- dotplot(GeneOnto[[i]])
print(PLOT)
ggsave(sprintf("dge/functional_enrichment/%s_ad10vswt10_dotplot.pdf", names(l)[[i]]))


PLOT_2 <- barplot(GeneOnto[[i]])
print(PLOT_2)
ggsave(sprintf("dge/functional_enrichment/%s_ad10vswt10_dotplot.pdf", names(l)[[i]]))

}


# ENRICHR
dge <- read.table("dge/ad10vswt10_Dge_HumanID.txt",header=T)

dbs <- c('GO_Biological_Process_2018','GO_Cellular_Component_2018',
         'GO_Molecular_Function_2018')


collapsed_output <- data.frame()
for(cur in as.character(unique(dge$Direction))){
  print(cur)
  # select genes
  cur_genes <- dge %>%
    subset(Direction == cur) %>% .$Gene

  # run enrichR on different gene sets:
  cur_result <- enrichr(cur_genes, dbs)

  # collapse results into one dataframe
  for(db in dbs){
    cur_result[[db]]$cluster <- cur
    cur_result[[db]]$db <- db
    cur_result[[db]]$Diagnosis <- 'HET'
    collapsed_output <- rbind(collapsed_output, cur_result[[db]])
  }
}

collapsed_output %>%
  write.csv(file='dge/functional_enrichment/ENRICHR_ad10vswt10_GO_terms.csv')


collapsed_output <- read_csv('dge/functional_enrichment/ENRICHR_ad10vswt10_GO_terms.csv')

input_bub_up <- collapsed_output %>% 
    #filter(str_detect(Term, "cytokine-mediated signaling pathway|inflammatory response|negative regulation of apoptotic process|positive regulation of chemokine|response to interleukin-1")) %>% 
    filter(cluster == "Upreg",str_detect(Term, "GO:0001959|GO:0050727|GO:0070555")) %>%
    filter(db == "GO_Biological_Process_2018", cluster %in% c("Downreg","Upreg")) %>% 
    mutate(log = -log10(P.value)) %>%
    group_by(cluster) %>%
    #top_n(3,log) %>%
    mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Term, perl=TRUE)) %>%    
    mutate(Term2 = as.factor(Term2), cluster = as.factor(cluster))


input_bub_down <- collapsed_output %>% 
    #filter(str_detect(Term, "cytokine-mediated signaling pathway|inflammatory response|negative regulation of apoptotic process|positive regulation of chemokine|response to interleukin-1")) %>% 
    filter(cluster == "Downreg",str_detect(Term, "GO:0006304|GO:0006285|GO:0042981")) %>%
    filter(db == "GO_Biological_Process_2018", cluster %in% c("Downreg","Upreg")) %>% 
    mutate(log = -log10(P.value)) %>%
    group_by(cluster) %>%
    #top_n(3,log) %>%
    mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Term, perl=TRUE)) %>%    
    mutate(Term2 = as.factor(Term2), cluster = as.factor(cluster))

input_bub <- rbind(input_bub_up,input_bub_down) %>% mutate(cluster = fct_relevel(cluster,c("Downreg","Upreg")),Term2 = fct_relevel(Term2,c("regulation of cytokine-mediated signaling pathway","regulation of inflammatory response","response to interleukin-1","DNA modification","base-excision repair","regulation of apoptotic process")))

colors <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF","#E16462FF", "#FCA636FF", "#F0F921FF")

pdf("dge/functional_enrichment/ENRICHR_ad10vswt10_GO_bubblechart.pdf", width = 4, height = 5)
ggballoonplot(input_bub %>% as.data.frame(), x = "cluster", y = "Term2",
              size = "log", fill = "log") +
   scale_fill_gradientn(colors = colors) +
     scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  guides(size = FALSE)
dev.off()



# ad3 vs wt3
dge <- read.table("dge/ad3vswt3_Dge_HumanID.txt",header=T)
l <- split(dge, dge$Direction)
l <- l[sapply(l, nrow)>0] #remove objects with less than one gene

GOI <- list()
GeneOnto <- list()

for(i in 1:length(l)){
GOI[[i]] <- bitr(as.character(l[[i]]$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Hs.eg.db::org.Hs.eg.db)

GeneOnto[[i]] <- enrichGO(gene = unique(GOI[[i]]$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 0.2, 
                     qvalueCutoff = 0.2, 
                     readable = TRUE)

openxlsx::write.xlsx(as.data.frame(GeneOnto[[i]]), 
                     file = sprintf("dge/functional_enrichment/%s_ad3vswt3_GO.xlsx", names(l)[[i]]), 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)


PLOT <- dotplot(GeneOnto[[i]])
print(PLOT)
ggsave(sprintf("dge/functional_enrichment/%s_ad3vswt3_dotplot.pdf", names(l)[[i]]))


PLOT_2 <- barplot(GeneOnto[[i]])
print(PLOT_2)
ggsave(sprintf("dge/functional_enrichment/%s_ad3vswt3_dotplot.pdf", names(l)[[i]]))

}


# ENRICHR
dge <- read.table("dge/ad3vswt3_Dge_HumanID.txt",header=T)

dbs <- c('GO_Biological_Process_2018','GO_Cellular_Component_2018',
         'GO_Molecular_Function_2018')


collapsed_output <- data.frame()
for(cur in as.character(unique(dge$Direction))){
  print(cur)
  # select genes
  cur_genes <- dge %>%
    subset(Direction == cur) %>% .$Gene

  # run enrichR on different gene sets:
  cur_result <- enrichr(cur_genes, dbs)

  # collapse results into one dataframe
  for(db in dbs){
    cur_result[[db]]$cluster <- cur
    cur_result[[db]]$db <- db
    cur_result[[db]]$Diagnosis <- 'HET'
    collapsed_output <- rbind(collapsed_output, cur_result[[db]])
  }
}

collapsed_output %>%
  write.csv(file='dge/functional_enrichment/ENRICHR_ad3vswt3_GO_terms.csv')


input_bub <- collapsed_output %>% 
    filter(db == "GO_Biological_Process_2018") %>% 
    mutate(log = -log10(P.value)) %>%
    group_by(cluster) %>%
    top_n(3,log) %>%
    mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Term, perl=TRUE)) %>%    
    mutate(Term2 = as.factor(Term2), cluster = as.factor(cluster))


colors <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF","#E16462FF", "#FCA636FF", "#F0F921FF")

pdf("dge/functional_enrichment/ENRICHR_ad3vswt3_GO_bubblechart.pdf", width = 4, height = 5)
ggballoonplot(input_bub, x = "cluster", y = "Term2",
              size = "log", fill = "log") +
   scale_fill_gradientn(colors = colors) +
     scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  guides(size = FALSE)
dev.off()


# wt10 vs wt3
dge <- read.table("dge/wt10vswt3_Dge_HumanID.txt",header=T)
l <- split(dge, dge$Direction)
l <- l[sapply(l, nrow)>0] #remove objects with less than one gene

GOI <- list()
GeneOnto <- list()

for(i in 1:length(l)){
GOI[[i]] <- bitr(as.character(l[[i]]$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Hs.eg.db::org.Hs.eg.db)

GeneOnto[[i]] <- enrichGO(gene = unique(GOI[[i]]$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 0.2, 
                     qvalueCutoff = 0.2, 
                     readable = TRUE)

openxlsx::write.xlsx(as.data.frame(GeneOnto[[i]]), 
                     file = sprintf("dge/functional_enrichment/%s_wt10vswt3_GO.xlsx", names(l)[[i]]), 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)


PLOT <- dotplot(GeneOnto[[i]])
print(PLOT)
ggsave(sprintf("dge/functional_enrichment/%s_wt10vswt3_dotplot.pdf", names(l)[[i]]))


PLOT_2 <- barplot(GeneOnto[[i]])
print(PLOT_2)
ggsave(sprintf("dge/functional_enrichment/%s_wt10vswt3_dotplot.pdf", names(l)[[i]]))

}


# ENRICHR
dge <- read.table("dge/wt10vswt3_Dge_HumanID.txt",header=T)

dbs <- c('GO_Biological_Process_2018','GO_Cellular_Component_2018',
         'GO_Molecular_Function_2018')


collapsed_output <- data.frame()
for(cur in as.character(unique(dge$Direction))){
  print(cur)
  # select genes
  cur_genes <- dge %>%
    subset(Direction == cur) %>% .$Gene

  # run enrichR on different gene sets:
  cur_result <- enrichr(cur_genes, dbs)

  # collapse results into one dataframe
  for(db in dbs){
    cur_result[[db]]$cluster <- cur
    cur_result[[db]]$db <- db
    cur_result[[db]]$Diagnosis <- 'HET'
    collapsed_output <- rbind(collapsed_output, cur_result[[db]])
  }
}

collapsed_output %>%
  write.csv(file='dge/functional_enrichment/ENRICHR_wt10vswt3_GO_terms.csv')


input_bub <- collapsed_output %>% 
    filter(db == "GO_Biological_Process_2018") %>% 
    mutate(log = -log10(P.value)) %>%
    group_by(cluster) %>%
    top_n(3,log) %>%
    mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Term, perl=TRUE)) %>%    
    mutate(Term2 = as.factor(Term2), cluster = as.factor(cluster))


colors <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF","#E16462FF", "#FCA636FF", "#F0F921FF")

pdf("dge/functional_enrichment/ENRICHR_wt10vswt3_GO_bubblechart.pdf", width = 4, height = 5)
ggballoonplot(input_bub, x = "cluster", y = "Term2",
              size = "log", fill = "log") +
   scale_fill_gradientn(colors = colors) +
     scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  guides(size = FALSE)
dev.off()
