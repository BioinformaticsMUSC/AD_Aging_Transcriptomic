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
})

dir.create("dge")

load("futcounts/Expression_Input.RData")


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


pdf("dge/PCA.pdf",width=6,height=6,useDingbats=FALSE)
pca.Sample<-prcomp(t(p))
PCi<-data.frame(pca.Sample$x,Genotype=pd$Genotype,ID = rownames(pd))
eig <- (pca.Sample$sdev)^2
variance <- eig*100/sum(eig)
ggscatter(PCi, x = "PC1", y = "PC2",
          color = "Genotype",palette="paired", 
          shape = "Genotype", size = 3,label = "ID")+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )"))+
theme_classic()
dev.off()

# Surrogate variable
mod <- model.matrix(~., pd)
mod0 <- model.matrix(~ 1, pd)

svaobj <- sva(as.matrix(p),mod,mod0,n.sv=NULL,B=100)

svaobj$sv <- data.frame(svaobj$sv)
colnames(svaobj$sv) = c(paste0('SV',seq(svaobj$n.sv)))
pdSv <- cbind(pd,svaobj$sv)


# Regression Expression
pd_sva <- pdSv %>%
            dplyr::select(-Genotype) %>% #Removing ID, Region, Diagnosis from the model
            droplevels()

betas<-future_lapply(1:nrow(p), function(x)
            {
                lm(unlist(p[x,])~., data = pd_sva)
            })

residuals<-future_lapply(betas, function(x)residuals(summary(x)))
residuals<-do.call(rbind, residuals)
p_regressed <- residuals+matrix(future_apply(p, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
rownames(p_regressed)<-rownames(p)

write.table(p_regressed,"dge/expression_regressed.txt",sep="\t",quote=F)

pdf("dge/PCA_AdjustedForConfound.pdf",width=8,height=8,useDingbats=FALSE)
pca.Sample<-prcomp(t(p_regressed))
PCi<-data.frame(pca.Sample$x,Genotype=pd$Genotype,ID = rownames(pd))
eig <- (pca.Sample$sdev)^2
variance <- eig*100/sum(eig)
ggscatter(PCi, x = "PC1", y = "PC2",
          color = "Genotype",palette="paired", 
          shape = "Genotype", size = 3,label = "ID")+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )"))+
theme_classic()  +
theme(legend.position= "none")
dev.off()



# Emmeans post-hoc
# Modeling
model1 <- 'geneExpr ~ Condition * Time' # Null model

# Function for fitting the two models with/without the "diagnosis".
emmeans_lm_p <- function(vectorizedExpression) {
  tmpMetaData <- cbind(pd, data.frame(geneExpr = unname(vectorizedExpression)))
  residuals1 <- lm(model1,data=tmpMetaData)
  pval <- emmeans::emmeans(residuals1,pairwise ~ Condition * Time,adjust="tukey")$contrasts %>%
          broom::tidy() %>%
          select(contrast,adj.p.value) %>% 
          pivot_wider(names_from = contrast, values_from = adj.p.value) %>%
          as.data.frame()
       
}

emmeans_lm_p_fun <- function(vectorizedExpression) {
  tryCatch(emmeans_lm_p(vectorizedExpression))
}


posthoc_p <- future_apply(p, 1, emmeans_lm_p_fun) %>% 
                      bind_rows() %>%
                      clean_names() %>%
                      as.data.frame() 

rownames(posthoc_p) <- rownames(p)

# Eff Size
emmeans_lm_eff <- function(vectorizedExpression) {
  tmpMetaData <- cbind(pd, data.frame(geneExpr = unname(vectorizedExpression)))
  residuals1 <- lm(model1,data=tmpMetaData)
  
  effect_size <- emmeans::emmeans(residuals1, adjust="tukey",pairwise ~ Condition * Time)$contrasts %>%
          broom::tidy() %>%
          select(contrast,estimate) %>% 
          pivot_wider(names_from = contrast, values_from = estimate) %>%
          as.data.frame()
       
}

emmeans_lm_eff_fun <- function(vectorizedExpression) {
  tryCatch(emmeans_lm_eff(vectorizedExpression))
}


posthoc_eff <- future_apply(p, 1, emmeans_lm_eff_fun) %>% 
                      bind_rows() %>%
                      clean_names() %>%
                      as.data.frame() 

rownames(posthoc_eff) <- rownames(p)


write.table(posthoc_p,"dge/AD_Dge_Interaction_Posthoc_Pvalues.txt",sep="\t",quote=F)
write.table(posthoc_eff,"dge/AD_Dge_Interaction_Posthoc_EffSize.txt",sep="\t",quote=F)


# Filter for significant
sign_posthoc_p <- posthoc_p %>%
                   filter(if_any(is.numeric, ~ .x < 0.05)) %>%
                   rownames_to_column("Gene")

sign_posthoc_eff <- posthoc_eff[rownames(posthoc_eff) %in% sign_posthoc_p$Gene,] %>%
                    clean_names() %>%
                   rownames_to_column("Gene")


# Save
save(posthoc_eff,posthoc_p,sign_posthoc_p,sign_posthoc_eff,  file = "dge/AD_Dge_Interaction_Posthoc.RData")


openxlsx::write.xlsx(sign_posthoc_p, 
                     file = "dge/AD_Dge_Interaction_Posthoc_Pvalues_Sign.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")


openxlsx::write.xlsx(sign_posthoc_eff, 
                     file = "dge/AD_Dge_Interaction_Posthoc_EffSize_Sign.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")

# Combine

ad10vswt10 <- data.frame(Gene = rownames(p), logFC = fc$ad_10m_wt_10m, FDR = p$ad_10m_wt_10m)
ad3vswt3 <- data.frame(Gene = rownames(p), logFC = fc$ad_3m_wt_3m, FDR = p$ad_3m_wt_3m)
wt10vswt3 <- data.frame(Gene = rownames(p), logFC = fc$wt_10m_wt_3m, FDR = p$wt_10m_wt_3m)

ad10vswt10_Sign <- ad10vswt10 %>% filter(FDR < 0.05, abs(logFC) > 0.3)
ad3vswt3_Sign <- ad3vswt3 %>% filter(FDR < 0.05, abs(logFC) > 0.3)
wt10vswt3_Sign <- wt10vswt3 %>% filter(FDR < 0.05, abs(logFC) > 0.3)

save(ad10vswt10,ad3vswt3,wt10vswt3,ad10vswt10_Sign,ad3vswt3_Sign,wt10vswt3_Sign, file = "dge/AD_Dge_Defined.RData")


