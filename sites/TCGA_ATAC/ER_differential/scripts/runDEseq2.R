library(DESeq2)

## count tables
counts <- read.table('../all_sites/BRCA_peak_counts_for_DEseq2.txt', header = TRUE, row.names = 1)
status <- factor(read.table('../all_sites/labels_for_DEseq2.txt')$V1)
## object construction
dds <- DESeqDataSetFromMatrix(counts, DataFrame(status), ~ status)
## standard analysis
dds <- DESeq(dds)
res <- results(dds)
## moderated log2 fold changes
resLFC <- lfcShrink(dds, coef=2, type="apeglm")

pdf('../plots/ER_differential_MA_plot.pdf')
plotMA(resLFC)
dev.off() 

write.table(as.data.frame(resLFC), file="../all_sites/ER_diff_DEseq2_results.txt", quote=FALSE, sep="\t")
