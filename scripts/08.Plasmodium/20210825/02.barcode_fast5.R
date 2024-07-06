setwd("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex")

aligns <- readRDS("./analysis/07.VirusClassification/02.Prediction/20210825/AlignmentResult.Rds")
aligns[, qname := paste0("read_", qname)]
aligns <- aligns[!flag %in% c(2048, 2064)]
aligns <- aligns[Species == "PbergheiANKA"]
aligns <- aligns[Length >= 200]

Preds <- readRDS("./analysis/07.VirusClassification/02.Prediction/20210825/BarcodesPrediction.Rds")
colnames(Preds) <- gsub("\\.", "-", colnames(Preds))
Preds <- as.data.table(Preds, keep.rownames = "qname")
Preds <- Preds[, .(qname, `RTA-08`, `RTA-10`, `RTA-27`, `RTA-33`, `RTA-37`)]
Preds <- merge(Preds, melt.data.table(Preds)[,.(pred = variable[which.max(value)]), qname], by = "qname")
Preds$PP <- apply(Preds[, grepl("RTA", colnames(Preds)), with = F], 1, max)
Preds <- Preds[pred %in% c("RTA-27", "RTA-08")]
Preds[PP > 0.1, table(pred)]
Preds <- Preds[PP > 0.1, ]

Mat <- merge(aligns, Preds, by = "qname")
Mat[, table(pred)]
Mat[mapq == 60, table(pred)]
Mat[mapq == 60, prop.table(table(pred))]
Mat[, prop.table(table(pred))]
Mat[, qname := gsub("^read_", "", qname)]

single_fast5 <- list.files("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/08.m6A/01.Plasmodium/20210825/01.single_fast5/fast5_pass", "fast5$", recursive = T, full.names = T)
names(single_fast5) <- gsub(".fast5$", "", basename(single_fast5))

Mat <- Mat[qname %in% names(single_fast5)]
Mat[, table(pred)]
Mat[, prop.table(table(pred))]

Mat[, flag := NULL]
Mat[, mapq := NULL]
Mat <- unique(Mat)
Mat[, table(pred)]
Mat[, prop.table(table(pred))]


sh1 <- paste0("cp ", single_fast5[Mat[pred == "RTA-08", qname]], " /mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/08.m6A/01.Plasmodium/20210825/02.barcode_fast5/RTA08")
write.table(sh1, "/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/08.m6A/01.Plasmodium/20210825/02.barcode_fast5/RTA08/README", row.names = F, col.names = F, quote = F)

sh2 <- paste0("cp ", single_fast5[Mat[pred == "RTA-27", qname]], " /mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/08.m6A/01.Plasmodium/20210825/02.barcode_fast5/RTA27")
write.table(sh2, "/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/08.m6A/01.Plasmodium/20210825/02.barcode_fast5/RTA27/README", row.names = F, col.names = F, quote = F)



