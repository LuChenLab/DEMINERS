library(openxlsx)
library(data.table)

openxlsx::getSheetNames("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/Meta/Barcode_and_RNA_information.xlsx")

# 20200605 ----

Mat <- as.data.table(openxlsx::read.xlsx("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/Meta/Barcode_and_RNA_information.xlsx", sheet = 1))
Mat <- Mat[, c(1, 2, 5, 7)]
Mat1 <- Mat[, 1:4]

library(stringr)
Mat1$RNA.Sequnce <- str_remove_all(Mat1$RNA.Sequnce, "\r\n")

library(Biostrings)
Seq1 <- DNAStringSet(Mat1$RNA.Sequnce)
names(Seq1) <- gsub("\\(与后期不一样\\)", "", Mat1$Gene.Name)

writeXStringSet(Seq1, filepath = "/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/Reference/RNAsequence_20200605.fa")


# 20200620 ----


Mat <- as.data.table(openxlsx::read.xlsx("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/Meta/Barcode_and_RNA_information.xlsx", sheet = 3))
Mat <- Mat[, c(1, 2, 5, 7)]
Mat1 <- Mat[, 1:4]

Mat1$RNA.Sequnce <- str_remove_all(Mat1$RNA.Sequnce, "\r\n")

Seq1 <- DNAStringSet(Mat1$RNA.Sequnce)
names(Seq1) <- gsub("\\(与后期不一样\\)", "", Mat1$Gene.Name)

writeXStringSet(Seq1, filepath = "/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/Reference/RNAsequence_20200620.fa")



# 20200902 ----


Mat <- as.data.table(openxlsx::read.xlsx("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/Meta/Barcode_and_RNA_information.xlsx", sheet = 3))
Mat <- Mat[, c(1, 2, 5, 7)]
Mat1 <- Mat[, 1:4]

Mat1$RNA.Sequnce <- str_remove_all(Mat1$RNA.Sequnce, "\r\n")

Seq1 <- DNAStringSet(Mat1$RNA.Sequnce)
names(Seq1) <- gsub("\\(与后期不一样\\)", "", Mat1$Gene.Name)

writeXStringSet(Seq1, filepath = "/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/Reference/RNAsequence_20200902.fa")





# 20210703 ----

Mat <- openxlsx::read.xlsx("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/Rawdata/20210703/Barcode RNA-20210703测序对应关系.xlsx", sheet = 2)

library(stringr)
Mat$RNA.Sequnce <- str_remove_all(Mat$RNA.Sequnce, "\r\n")

library(Biostrings)
Seqs <- DNAStringSet(Mat$RNA.Sequnce)
names(Seqs) <- Mat$Gene.Name

writeXStringSet(Seqs, filepath = "/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/Reference/RNAsequence_20210703.fa")


