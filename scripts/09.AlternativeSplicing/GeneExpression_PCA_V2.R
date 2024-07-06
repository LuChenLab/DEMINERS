library(data.table)
library(GenomicFeatures)
library(GenomicAlignments)
library(BiocParallel)
library(ggplot2)

TxDb <- rtracklayer::readGFFAsGRanges("/mnt/raid61/Personal_data/songjunwei/reference/plasmoDB/PlasmoDB-53_PbergheiANKA.gff")
TxDb <- TxDb[mapply(length, with(TxDb, Parent)) == 0]
table(TxDb$type)

ol <- findOverlaps(TxDb, TxDb)
ol <- ol[queryHits(ol) != subjectHits(ol)]
TxDb <- TxDb[-queryHits(ol)]
names(TxDb) <- TxDb$ID



bamFile <- list.files("/mnt/raid61/Personal_data/tangchao/Temp/20211025/02.SplitFastq", "bam$", full.names = T)
bamFile <- grep("RTA", bamFile, value = T)

flag <- scanBamFlag(isSecondaryAlignment = FALSE,
                    isNotPassingQualityControls = FALSE,
                    isUnmappedQuery = FALSE,
                    isDuplicate = FALSE)
sbp <- ScanBamParam(flag = flag, mapqFilter = 1)

bamLst <- BamFileList(bamFile, yieldSize = 2000000)
options(srapply_fapply = "parallel", mc.cores = 6)
se_count1 <- summarizeOverlaps(features = TxDb, 
                               reads = bamLst, 
                               mode = "Union",
                               ignore.strand = FALSE, 
                               inter.feature = FALSE, 
                               singleEnd = TRUE,
                               fragments = FALSE, 
                               param = sbp, 
                               preprocess.reads = NULL)
colnames(se_count1) <- gsub(".bam", "", colnames(se_count1))
colSums(assay(se_count1))


meta <- data.frame(Sample = colnames(se_count1), Stage = rep(c("Trophozoite", "Schizont"), each = 3), row.names = colnames(se_count1))

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = assay(se_count1), colData = DataFrame(meta), design = ~ Stage)

hist(log2(1 + rowSums(counts(dds))), breaks = 100, main = "Histogram of log2 gene counts")
abline(v = log2(1 + 10), col = 2)

sum(rowSums(counts(dds)) > 10)
length(dds)

dds <- dds[rowSums(counts(dds)) > 10, ]

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

ntd <- normTransform(dds)
rld <- rlog(dds, blind = FALSE)




library("RColorBrewer")
library(pheatmap)
sampleDists <- dist(t(assay(ntd)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, annotation_row = meta,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, main = "ntd")


sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, annotation_row = meta,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, main = "rld")


plotPCA(ntd, intgroup = c("Stage"))
plotPCA(ntd, intgroup = c("Stage"), ntop = 100)
plotPCA(ntd, intgroup = c("Stage"), ntop = length(ntd))

plotPCA(rld, intgroup = c("Stage"))
plotPCA(rld, intgroup = c("Stage"), ntop = 100)
plotPCA(rld, intgroup = c("Stage"), ntop = length(rld))

plotPCA(rld, intgroup = c("Stage"), ntop = length(rld)) + 
  scale_colour_manual(values = RColorBrewer::brewer.pal(n = 3, "Dark2")[1:2]) + 
  theme_classic(base_size = 16) + 
  guides(colour = guide_legend(title = "Stage")) + 
  theme(legend.position = "top")



CorMat <- cor(assay(rld), method = "spearman")
CorMat <- 1 - CorMat
loc <- cmdscale(CorMat)
colnames(loc) <- c("Dim1", "Dim2")
loc <- merge(as.data.table(loc, keep.rownames = "Sample"), meta, by = "Sample")

library(ggrepel)
ggplot(loc, aes(x = Dim1, y = Dim2, color = Stage))+
  geom_point(size = 2) + #Size and alpha just for fun
  geom_text_repel(aes(label = Sample)) + 
  # scale_color_manual(values = TissueCols) +
  theme_bw() +
  xlab("Dim1") +
  ylab("Dim2")


ggplot(loc, aes(x = Dim1, y = Dim2, color = Stage)) + 
  geom_point() + 
  geom_text_repel(aes(label = Sample)) +
  scale_colour_manual(values = RColorBrewer::brewer.pal(n = 3, "Dark2")[1:2]) + 
  theme_bw(base_size = 15) + 
  guides(colour = guide_legend(title = "Stage")) + 
  theme(legend.position = "top")
ggsave("/mnt/raid61/Personal_data/tangchao/Temp/20211025/02.SplitFastq/GeneExpression_MDS.pdf", width = 4.2, height = 4.2)




tsne_res <- Rtsne::Rtsne(t(assay(rld)), perplexity = 1)
tsne_res <- tsne_res$Y
row.names(tsne_res) <- colnames(assay(rld))
colnames(tsne_res) <- c("t-SNE 1", "t-SNE 2")

tsne_res <- merge(as.data.table(tsne_res, keep.rownames = "Sample"), meta, by = "Sample")

ggplot(tsne_res, aes(x = `t-SNE 1`, y = `t-SNE 2`, color = Stage))+
  geom_point(size = 2) + #Size and alpha just for fun
  geom_text_repel(aes(label = Sample)) + 
  # scale_color_manual(values = TissueCols) +
  theme_bw()

ggplot(tsne_res, aes(x = `t-SNE 1`, y = `t-SNE 2`, color = Stage)) + 
  geom_point() + 
  geom_text_repel(aes(label = Sample)) +
  scale_colour_manual(values = RColorBrewer::brewer.pal(n = 3, "Dark2")[1:2]) + 
  theme_bw(base_size = 15) + 
  guides(colour = guide_legend(title = "Stage")) + 
  theme(legend.position = "top")
ggsave("/mnt/raid61/Personal_data/tangchao/Temp/20211025/02.SplitFastq/GeneExpression_tSNE.pdf", width = 4.2, height = 4.2)


