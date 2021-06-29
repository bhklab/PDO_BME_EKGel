library(ggrepel)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(snowfall)
library(GSA)
library(piano)

setwd("PDO_BME_EKGel") # set to path of cloned project in your machine


out_PDO <- readRDS("Data/UHN_BR_PDOs_hg38_kallisto_mixed.rda") # RNAseq data

PDOs <- out_PDO$rnaseq@assayData$exprs  ## RNAseq data [normalization: log2(TPM+0.001)]

features_gene <- out_PDO$rnaseq@featureData@data

proteinCoding <- rownames(features_gene)[features_gene$gene_type=="protein_coding"]



xx <- cor(PDOs,method = "s")
col_fun = colorRamp2(c(0.7,1), c("white", "red"))

pdf("Plots/Fig4A.pdf",width = 10,height = 8)
Heatmap(xx,show_row_names = T,show_column_names = F,column_title = "Pairwise correlation [All genes]"
        ,col = col_fun,name = "Spearman",)
dev.off()

############################
############################

data_counts <- out_PDO$rnaseq.counts ## RNAseq data [normalization: log2(counts+1)]



pData(data_counts) <- data.frame("sampleid"=colnames(exprs(data_counts)),stringsAsFactors = F)
pData(data_counts)$patientID <- factor(c(rep("BPTO_51",3), rep("BPTO_95",3),rep("BXTO_64",5)))
pData(data_counts)$batch <- factor(c("BME","CNC","Tumour_Tissue","BME","CNC","Tumour_Tissue"
                                     ,"BME","BME","CNC","CNC","Tumour_Tissue"))
pData(data_counts)$type <- factor(c("ORG","ORG","Tumour_Tissue",
                                    "ORG","ORG","Tumour_Tissue",
                                    "ORG","PDX","ORG","PDX","Tumour_Tissue"))

data_counts_sub <- data_counts[,data_counts$type%in%c("ORG")]
data_counts_sub <- data_counts_sub[,data_counts_sub$batch%in%c("BME","CNC")]

data_counts_sub$batch <- factor(pData(data_counts_sub)$batch)
data_counts_sub$patientID <- factor(pData(data_counts_sub)$patientID)
dds <- DESeqDataSetFromMatrix(countData = round(2^(exprs(data_counts_sub))-1),
                              colData = pData(data_counts_sub),
                              design= ~ patientID + batch)


keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
dds <- DESeq(dds)


resultsNames(dds) # lists the coefficients
res <- results(dds, name="batch_CNC_vs_BME")

# or to shrink log fold changes association with condition:
resLF <- lfcShrink(dds, coef="batch_CNC_vs_BME", type="apeglm")

res <- res[order(abs(res$stat),decreasing = T),]
res <- res[order(abs(res$padj),decreasing = F),]

plotMA(resLF)

res[1:10,]
res$gene <- features_gene[rownames(res),"gene_name"]


res_df <- as.data.frame(res,stringsAsFactors=F)
res_df[,"FDR < 0.05"] <- factor(ifelse(!is.na(res_df$padj) &res_df$padj<0.05,"Yes","No"))


pdf("Plots/Fig4B.pdf",width = 8,height = 6)
ggplot(res_df,aes(log2FoldChange,-log10(pvalue),col=`FDR < 0.05`)) +
  geom_point() +
  scale_colour_manual(values = c("black","red")) +
  geom_text_repel(data = res_df[!is.na(res_df$padj) & res_df$padj<0.05,]
                  ,aes(log2FoldChange,-log10(pvalue),label=gene)
  )  +
  ggtitle("BME vs CNC Differential Expression Analysis [ORG]") +
  theme_bw() + theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank())


dev.off()




###################################
# Pathway analysis
set.seed(3425)


nPerm <- 10000
gsc1 <- loadGSC("Data/c2.cp.reactome.v7.4.symbols.gmt")


sigGenes <- res[!is.na(res$pvalue) & res$pvalue < 0.05,]
finalSetOfGenes_tmp <- sigGenes$gene

gsea_out <- piano::runGSAhyper(genes = finalSetOfGenes_tmp, gsc=gsc1,universe = unique(features_gene[proteinCoding,"gene_name"]))
gsea_out_final <- as.data.frame(gsea_out$resTab)
gsea_out_final$fdr <- p.adjust(gsea_out_final$`p-value`,method = "fdr")
gsea_out_final <- gsea_out_final[order(gsea_out_final$`p-value`),]

gsea_out_final1 <- gsea_out_final[gsea_out_final$fdr<0.05,]

gsea_out_final1 <- gsea_out_final1[order(gsea_out_final1$fdr,decreasing = T),]



write.table(gsea_out_final,"Hypergeometric_gene_set_analysis.csv",sep = ",",row.names = T,col.names = NA)

