library(data.table)
library(parallel)
library(caret)
library(doParallel)
library(pbapply)

setwd("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex")

SigTab <- readRDS("./analysis/07.VirusClassification/03.Merge/00.BarcodeSignal/BarcodeSignal.Rds")

for(i in seq(2, 24, 2)) {
  print(i)
  
  read2gene <- readRDS("./analysis/07.VirusClassification/03.Merge/00.BarcodeSignal/Reads2Barcode.Rds")
  setkey(read2gene, read)
  
  rbt2u <- read2gene[, .N, Barcode][order(N, decreasing = T)][seq_len(i), Barcode]
  read2gene <- read2gene[Barcode %in% rbt2u]
  
  set.seed(9560)
  R2B <- as.data.table(downSample(x = read2gene, y = factor(read2gene$Barcode)))

  set.seed(3456)
  trainIndex <- createDataPartition(R2B$Class, p = .7, list = TRUE, times = 1)[[1]]
  TrainingReads <- R2B[trainIndex, ]
  TestReads <- read2gene[!read %in% TrainingReads$read, ]
  TestReads[, Class := factor(Barcode)]
  
  TrainingData <- SigTab[TrainingReads[, read], ]
  stopifnot(identical(row.names(TrainingData), TrainingReads$read))
  TrainingData$Class <- TrainingReads$Class
  TestData <- SigTab[TestReads[, read], ]
  stopifnot(identical(row.names(TestData), TestReads$read))
  TestData$Class <- TestReads$Class
  
  out <- paste0("./analysis/11.Classifiers/01.ModelingData/Barcode_", i, ".RData")
  save(TrainingReads, TestReads, TrainingData, TestData, file = out)
}



setwd("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex")

N_Reads <- lapply(seq(2, 24, 2), function(i) {
  load(paste0("./analysis/11.Classifiers/01.ModelingData/Barcode_", i, ".RData"))
  reads <- merge(TrainingReads[, .(train = .N), Barcode], TestReads[, .(test = .N), Barcode], by = "Barcode")
  data.table(N = i, reads)
})
N_Reads <- do.call(rbind, N_Reads)

openxlsx::write.xlsx(N_Reads, "./analysis/11.Classifiers/01.ModelingData/Classifier_Reads.xlsx")
