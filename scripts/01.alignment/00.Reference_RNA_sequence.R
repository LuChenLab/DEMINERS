library(openxlsx)
Mat <- openxlsx::read.xlsx("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/Meta/Barcode_RNA-更正终版.xlsx", sheet = 2)
Mat <- Mat[, c(1, 2, 5, 7, 8)]
Mat1 <- subset(Mat, !grepl("\\/", Gene.Name))
Mat2 <- subset(Mat, grepl("\\/", Gene.Name))
Mat1 <- Mat1[, 1:4]

library(stringr)
Mat1$RNA.Sequnce <- str_remove_all(Mat1$RNA.Sequnce, "\r\n")

Mat2$RNA.Sequnce <- str_remove_all(Mat2$RNA.Sequnce, "\r\n")
Mat2$X8 <- str_remove_all(Mat2$X8, "\r\n")

library(Biostrings)
Seq1 <- DNAStringSet(Mat1$RNA.Sequnce)
names(Seq1) <- Mat1$Gene.Name

Seq2 <- DNAStringSet(c(Mat2$RNA.Sequnce, Mat2$X8))
names(Seq2) <- unlist(strsplit(Mat2$Gene.Name, "/"))

Seqs <- c(Seq1, Seq2)

writeXStringSet(Seqs, filepath = "/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/Reference/RNAsequence.fa")
