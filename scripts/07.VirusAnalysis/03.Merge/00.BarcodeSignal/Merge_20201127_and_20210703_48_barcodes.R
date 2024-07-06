setwd("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex")
library(data.table)
library(parallel)

R2B_1 <- fread("analysis/02.BamReadsSplit/read2gene_20201127.txt")
R2B_2 <- fread("analysis/02.BamReadsSplit/read2gene_20210703.txt")

R2B <- rbind(data.table(R2B_1, batch = "2020-11-27"), data.table(R2B_2, batch = "2021-07-03"))
R2B[, read := paste("read", read, sep = "_")]

# 20201127
dir_bs <- "./analysis/03.BarcodeProcessing/02.BarcodeSignalBinExp/20201127"
files <- list.files(path = dir_bs, pattern = ".signal", full.names = TRUE)

sigs1 <- mclapply(FUN = function(x) {
  load(x)
  barcode
}, files, mc.cores = 4)
sigs1 <- do.call("c", sigs1)

# 20210703
dir_bs <- "./analysis/03.BarcodeProcessing/02.BarcodeSignalBinExp/20210703"
files <- list.files(path = dir_bs, pattern = ".signal", full.names = TRUE)

sigs2 <- mclapply(FUN = function(x) {
  load(x)
  barcode
}, files, mc.cores = 4)
sigs2 <- do.call("c", sigs2)

sigs <- c(sigs1, sigs2)

SigTab <- as.data.frame(do.call(rbind, sigs))
colnames(SigTab) <- paste0("BIN", sprintf("%03d", seq_len(ncol(SigTab))))

SigTab <- subset.data.frame(SigTab, row.names(SigTab) %in% R2B[, read])
R2B <- R2B[read %in% row.names(SigTab)]

saveRDS(SigTab, file = "./analysis/07.VirusClassification/03.Merge/00.BarcodeSignal/BarcodeSignal.Rds")
saveRDS(R2B, file = "./analysis/07.VirusClassification/03.Merge/00.BarcodeSignal/Reads2Barcode.Rds")
