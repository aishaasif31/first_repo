
library(DESeq2)

FC <- 2 # Fold-change cutoff
FDR <- 0.05 # FDR cutoff
alpha <- 0.1 # independent filtering

raw_counts = read.csv("converted_counts_data.csv")
row.names(raw_counts) <- raw_counts$User_ID
raw_counts <- raw_counts[, -(1:3)] 
str(raw_counts)

col_data <- data.frame(
  "sample" = c("NF10", "NF11", "NF12", "NF13", "NF14", "NF1", "NF3", "NF5", "NF7", "NF15", "NF2", "NF4", "NF6", "NF8", "DCM8", "DCM14", "DCM12", "DCM13", "DCM11", "DCM17", "DCM16", "DCM18", "DCM19", "DCM21", "DCM20", "DCM22", "DCM24", "DCM26", "DCM28", "DCM29", "DCM2", "DCM31", "DCM32", "DCM33", "DCM34", "DCM35", "DCM38", "DCM39", "DCM3", "DCM43", "DCM44", "DCM46", "ICM47", "ICM49", "DCM4", "ICM50", "ICM52", "ICM53", "ICM54", "ICM55", "ICM56", "ICM58", "ICM59", "ICM61", "ICM63", "ICM64", "DCM6", "DCM72", "DCM73", "DCM74", "DCM75", "DCM76", "DCM77", "DCM78"),
  "groups" = c("NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "ICM", "ICM", "DCM", "ICM", "ICM", "ICM", "ICM", "ICM", "ICM", "ICM", "ICM", "ICM", "ICM", "ICM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM")
)

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = col_data,
  design = ~groups
)
dds <- DESeq2::DESeq(dds)

# Comparison 1 of 3:  DCM-ICM
res <- DESeq2::results(dds, 
                       contrast = c("groups", "DCM", "ICM"),
                       independentFiltering = TRUE,
                       alpha = alpha
)
summary(res)
plotMA(res)
plotCounts(dds, gene = which.min(res$padj), intgroup = colnames(col_data)[2])
res <- subset(res, padj < FDR & abs(log2FoldChange) > log2(FC)) # Select
table(sign(res$log2FoldChange)) # N. of genes Down, Up
res <- res[order(-res$log2FoldChange), ] #sort
head(res) #top upregulated
tail(res) #top downregulated

res <- as.data.frame(res) %>% mutate(Expression = case_when(log2FoldChange >= 1.5 & pvalue <= 0.05 ~ "Up-regulated",log2FoldChange <= -1.5 & pvalue <= 0.05 ~ "Down-regulated",TRUE ~ "Unchanged"))

write.csv(res, file = "GSE116250_DEG.csv", row.names = TRUE)

# Comparison 2 of 3:  DCM-NF
ress <- DESeq2::results(dds, 
                       contrast = c("groups", "DCM", "NF"),
                       independentFiltering = TRUE,
                       alpha = alpha
)
summary(ress)
plotMA(ress)
plotCounts(dds, gene = which.min(res$padj), intgroup = colnames(col_data)[2])
ress <- subset(ress, padj < FDR & abs(log2FoldChange) > log2(FC)) # Select
table(sign(ress$log2FoldChange)) # N. of genes Down, Up
ress <- ress[order(-ress$log2FoldChange), ] #sort
head(ress) #top upregulated
tail(ress) #top downregulated
ress <- as.data.frame(ress) %>% mutate(Expression = case_when(log2FoldChange >= 1.5 & pvalue <= 0.05 ~ "Up-regulated",log2FoldChange <= -1.5 & pvalue <= 0.05 ~ "Down-regulated",TRUE ~ "Unchanged"))

write.csv(ress, file = "GSE116250_DEG2.csv", row.names = TRUE)

# Comparison 3 of 3:  ICM-NF
resss <- DESeq2::results(dds, 
                       contrast = c("groups", "ICM", "NF"),
                       independentFiltering = TRUE,
                       alpha = alpha
)
summary(resss)
plotMA(resss)
plotCounts(dds, gene = which.min(resss$padj), intgroup = colnames(col_data)[2])
resss <- subset(resss, padj < FDR & abs(log2FoldChange) > log2(FC)) # Select
table(sign(resss$log2FoldChange)) # N. of genes Down, Up
resss <- resss[order(-resss$log2FoldChange), ] #sort
head(resss) #top upregulated
tail(resss) #top downregulated

resss <- as.data.frame(resss) %>% mutate(Expression = case_when(log2FoldChange >= 1.5 & pvalue <= 0.05 ~ "Up-regulated",log2FoldChange <= -1.5 & pvalue <= 0.05 ~ "Down-regulated",TRUE ~ "Unchanged"))

write.csv(resss, file = "GSE116250_DEG3.csv", row.names = TRUE)


#----------------------

tT <- as.data.frame(tT) %>% mutate(Expression = case_when(logFC >= 1.5 & P.Value <= 0.05 ~ "Up-regulated",logFC <= -1.5 & P.Value <= 0.05 ~ "Down-regulated",TRUE ~ "Unchanged"))

write.csv(resss, file = "GSE51472_DEG1.csv", row.names = TRUE)


