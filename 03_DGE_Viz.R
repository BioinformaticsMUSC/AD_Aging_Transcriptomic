# Load libraries
suppressPackageStartupMessages({
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(data.table)
library(RColorBrewer)
library(tidyverse)
library(preprocessCore)
library(future.apply)
library(DESeq2)
library(pheatmap)
library(sva)
library(viridis)
library(limma)
library(janitor)
library(UpSetR)
})


# Load data
load("dge/AD_Dge_Defined.RData")
load("futcounts/Expression_Input.RData")
load("utils/geneset/GeneSets_Senescence.RData")

# Convert senescence in mgi

senescence <- do.call(rbind,GeneSets)

human = biomaRt::useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
mouse = biomaRt::useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")

MGI = biomaRt::getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = senescence$Gene,mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)

comb <- merge(senescence,MGI,by.x="Gene",by.y="HGNC.symbol",all=F)


# Create input expression for viz
pd <- data.frame(row.names=colnames(exp), 
                 Genotype = as.factor(c(rep("WT_3m",3),rep("WT_10m",3),rep("AD_3m",3),rep("AD_10m",3))),
                 Time = as.factor(c(rep("3M",3),rep("10M",3),rep("3M",3),rep("10M",3))),
                 Condition = as.factor(c(rep("WT",3),rep("WT",3),rep("AD",3),rep("AD",3))))

# Filter the expression by condition
filter=apply(rpkm, 1, function(x) (all(x[1:3] >= 0.5) | all(x[4:6] >= 0.5) | all(x[7:9] >= 0.5) | all(x[10:12] >= 0.5)))
count_filt <- exp[filter,]
rpkm_filt <- rpkm[filter,]

logCPM <- log2(rpkm_filt+1)
p <- normalize.quantiles(as.matrix(logCPM))
rownames(p) <- rownames(logCPM)
colnames(p) <- colnames(logCPM)



# Input for viz Plot for ad10 vs wt10
df <- ad10vswt10 %>% 
        mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>% 
        mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg", logFC < -0.3 & FDR < 0.05 ~ "DownReg"))



top_labelled <- df %>% 
                  filter(Gene %in% c("Cdkn2a","Icam1","Ccl3","Cxcl10","Cd68","Gfap")) %>%
                  group_by(Direction) %>% 
                  na.omit() 
                  #%>%
                  #arrange(ABS) %>%
                  #top_n(n = 5, wt = ABS)

#  boxplots
pd_filt <- pd %>% 
            filter(Genotype %in% c("WT_10m","AD_10m"))

mat <- p[rownames(p)%in% top_labelled$Gene,colnames(p) %in% rownames(pd_filt)] %>%
        t() %>%
        as.data.frame() %>%
        mutate(Genotype = pd_filt$Genotype) %>%
        pivot_longer(!Genotype, names_to = "Gene", values_to="Exp")

pdf("dge/Boxplots_TopGenes_ad10vswt10.pdf",width=8,height=2,useDingbats=FALSE)
ggboxplot(mat, "Genotype", "Exp", color = "Genotype",
 palette = c("red", "black")) +
      xlab("")+ 
      ylab("log2(Expression Adjusted)")+
theme_classic() + 
facet_wrap(.~Gene,scales="free",ncol=6,nrow=1) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
dev.off()

pdf("dge/Vulcano_Plot_ad10vswt10.pdf",width=6,height=6,useDingbats=FALSE)
ggscatter(df, 
            x = "logFC", 
            y = "LOG",
            color = "Threshold",
            palette=c("grey","red"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
      #geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      #geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                      mapping = aes(label = Gene), 
                      size = 5,
                      #box.padding = unit(0.4, "lines"),
                      point.padding = unit(0.4, "lines"), max.overlaps = Inf,
                      nudge_x = .15,
                      box.padding = 0.5,
                      nudge_y = 1,
                      segment.curvature = -0.1,
                      segment.ncp = 3,
                      segment.angle = 20)+
      theme(legend.position="none")+
      ylim(0,12) + xlim(-5,+5)
dev.off()

# heatmap
mat <- p[rownames(p)%in% ad10vswt10_Sign$Gene,colnames(p) %in% rownames(pd_filt)]
anno <- pd_filt
Condition        <- c("red", "black")
names(Condition) <- c("AD", "WT")
anno_colors <- list(Condition = Condition)
pdf("dge/Heatmap_ad10vswt10.pdf",width=4,height=6)
pheatmap(mat,scale="row",show_rownames = F,annotation=anno,annotation_colors = anno_colors)
dev.off()


# Input for viz Plot for ad3 vs wt3
df <- ad3vswt3 %>% 
        mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>% 
        mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg", logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

top_labelled <- df %>% 
                  group_by(Direction) %>% 
                  na.omit() %>%
                  arrange(FDR) %>%
                  top_n(n = 5, wt = LOG)

#  boxplots
pd_filt <- pd %>% 
            filter(Genotype %in% c("WT_3m","AD_3m"))

mat <- p[rownames(p)%in% top_labelled$Gene,colnames(p) %in% rownames(pd_filt)] %>%
        t() %>%
        as.data.frame() %>%
        mutate(Genotype = pd_filt$Genotype) %>%
        pivot_longer(!Genotype, names_to = "Gene", values_to="Exp")

pdf("dge/Boxplots_TopGenes_ad3vswt3.pdf",width=6,height=5,useDingbats=FALSE)
ggboxplot(mat, "Genotype", "Exp", color = "Genotype",
 palette = c("black", "red")) +
      xlab("")+ 
      ylab("log2(Expression Adjusted)")+
theme_classic() + 
facet_wrap(.~Gene,scales="free",ncol=4,nrow=4) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
dev.off()

pdf("dge/Vulcano_Plot_ad3vswt3.pdf",width=6,height=6,useDingbats=FALSE)
ggscatter(df, 
            x = "logFC", 
            y = "LOG",
            color = "Threshold",
            palette=c("grey","red"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
      #geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      #geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                      mapping = aes(label = Gene), 
                      size = 5,
                      box.padding = unit(0.4, "lines"),
                      point.padding = unit(0.4, "lines"), max.overlaps = Inf)+
      theme(legend.position="none")+
      ylim(0,10) + xlim(-5,+5)
dev.off()

# heatmap
mat <- p[rownames(p)%in% ad3vswt3_Sign$Gene,colnames(p) %in% rownames(pd_filt)]
anno <- pd_filt
Condition        <- c("red", "black")
names(Condition) <- c("AD", "WT")
anno_colors <- list(Condition = Condition)
pdf("dge/Heatmap_ad3vswt3.pdf",width=4,height=6)
pheatmap(mat,scale="row",show_rownames = F,annotation=anno,annotation_colors = anno_colors)
dev.off()


# Input for viz Plot for wt10 vs wt3
df <-  wt10vswt3 %>% 
        mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>% 
        mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg", logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

top_labelled <- df %>% 
                  group_by(Direction) %>% 
                  na.omit() %>%
                  arrange(FDR) %>%
                  top_n(n = 5, wt = LOG)

#  boxplots
pd_filt <- pd %>% 
            filter(Genotype %in% c("WT_10m","WT_3m"))

mat <- p[rownames(p)%in% top_labelled$Gene,colnames(p) %in% rownames(pd_filt)] %>%
        t() %>%
        as.data.frame() %>%
        mutate(Genotype = pd_filt$Genotype) %>%
        pivot_longer(!Genotype, names_to = "Gene", values_to="Exp")

pdf("dge/Boxplots_TopGenes_wt10vswt3.pdf",width=6,height=5,useDingbats=FALSE)
ggboxplot(mat, "Genotype", "Exp", color = "Genotype",
 palette = c("black", "red")) +
      xlab("")+ 
      ylab("log2(Expression Adjusted)")+
theme_classic() + 
facet_wrap(.~Gene,scales="free",ncol=4,nrow=4) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
dev.off()

pdf("dge/Vulcano_Plot_wt10vswt3.pdf",width=6,height=6,useDingbats=FALSE)
ggscatter(df, 
            x = "logFC", 
            y = "LOG",
            color = "Threshold",
            palette=c("grey","red"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
      #geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      #geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                      mapping = aes(label = Gene), 
                      size = 5,
                      box.padding = unit(0.4, "lines"),
                      point.padding = unit(0.4, "lines"), max.overlaps = Inf)+
      theme(legend.position="none")+
      ylim(0,10) + xlim(-5,+5)
dev.off()


# UpSet
m3 <- read.table("dge/ad3vswt3_Dge_MouseID.txt",header=T)
m10 <- read.table("dge/ad10vswt10_Dge_MouseID.txt",header=T)

m3$Class <- "3m"
m10$Class <- "10m"

m3 <- m3 %>%
      unite("Def", Direction:Class, na.rm = TRUE, remove = T) %>%
      filter(!(Def == "All_3m"))

m10 <- m10 %>%
      unite("Def", Direction:Class, na.rm = TRUE, remove = T) %>%
      filter(!(Def == "All_10m"))

tmp_upset <- rbind(m3,m10)

l <- split(as.character(tmp_upset$Gene),tmp_upset$Def)
Class <- names(l)
ToTGene <- as.numeric(sapply(l, length))
metadata <- as.data.frame(cbind(Class, ToTGene))
names(metadata) <- c("Class", "ToTGene")
metadata$ToTGene <- as.numeric(as.character(metadata$ToTGene))

pdf("dge/Upset_Plot_Intersection.pdf", width = 6, height = 4)
upset(fromList(l),,nsets = 4, set.metadata = list(data = metadata, plots = list(list(type = "hist", 
    column = "ToTGene", assign = 20), 
    list(type = "matrix_rows", column = "sets", colors = c(Upreg_3m = "#89E651", Downreg_3m = "#DB8BE4",Upreg_10m = "#DDA187", Downreg_10m="#DDD3DD"), 
    alpha = 0.5))))
dev.off()














# Boxplot
posthoc_p <- read.table("dge/AD_Dge_Interaction_Posthoc_Pvalues.txt") 

sign_posthoc_p <- posthoc_p %>%
                   filter(if_any(is.numeric, ~ .x < 0.05)) %>%
                   rownames_to_column("Gene")

mat <- p[rownames(p) %in% sign_posthoc_p$Gene,] %>% 
        as.data.frame() %>%
        rownames_to_column("Gene") %>%
        melt() %>%
        as.data.frame() 

tmp <- merge(mat,pd,by.x="variable",by.y="row.names") %>%
        mutate(Gene == as.factor(Gene)) %>%
        mutate(Condition = fct_relevel(Condition, c("WT","AD")), 
               Time = fct_relevel(Time, c("3M","10M")))

dir.create("dge/Plot")

cl <- colors(distinct = TRUE)
set.seed(15887) # to set random generator seed
cols <- c("#454b87","#bd3106")
doPlot = function(sel_name) 
{
    df = subset(tmp, Gene == sel_name)
    PLOT= ggboxplot(df, "Time", "value", color = "Condition",
              palette = cols,
              outlier.shape = NA) +
      xlab("")+ 
      ylab("Gene Expression")+
        theme_classic() +
        rotate_x_text(angle = 45)

    print(PLOT)
    ggsave(sprintf("dge/Plot/%s.pdf", sel_name),width=4,height=4)
 }

lapply(unique(tmp$Gene), doPlot)


















