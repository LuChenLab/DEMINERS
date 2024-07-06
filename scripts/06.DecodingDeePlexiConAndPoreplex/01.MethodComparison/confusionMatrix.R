library(caret)
library(data.table)

read2gene <- readRDS(file = "analysis/06.DecodingDeePlexiConAndPoreplex/00.ModelTraining2/00.RData/read2gene.Rds")

# DecodeR-deeplexicon

DecodeR <- fread("/mnt/raid61/Personal_data/tangchao/DecodeR/BarcodeDecomplex/analysis/10.MethodComparison/05.Accuracy/DecodeR/DecodeR_20200605_output.tsv")
DecodeR[, mean(Read %in% read2gene$read)]
DecodeR <- merge(DecodeR, read2gene[, .(read, Name)], by.x = "Read", by.y = "read")
DecodeR[, mean(Barcode == Name)]

DecodeR[, pred := factor(Barcode, levels = c("D-BC1", "D-BC2", "D-BC3", "D-BC4"))]
DecodeR[, truth := factor(Name, levels = c("D-BC1", "D-BC2", "D-BC3", "D-BC4"))]

DecodeR[, confusionMatrix(pred, truth)]


# deeplexicon

deeplexicon <- fread("/mnt/raid62/BetaCoV/Person/tangchao/analysis/DrictRNA/deeplexicon/conda/output1.tsv")
deeplexicon[, ReadID := paste0("read_", ReadID)]
PPs <- melt.data.table(deeplexicon[, .(ReadID, P_bc_1, P_bc_2, P_bc_3, P_bc_4)], id.vars = "ReadID", variable.name = "Barcode", value.name = "Probability")
PPs <- PPs[, .SD[which.max(Probability), ], ReadID]
deeplexicon <- merge(deeplexicon[, .(ReadID, P_bc_1, P_bc_2, P_bc_3, P_bc_4)], PPs, by = "ReadID")

deeplexicon <- merge(deeplexicon, read2gene[, .(read, Name)], by.x = "ReadID", by.y = "read")
deeplexicon[, Name := plyr::mapvalues(Name, c("D-BC1", "D-BC2", "D-BC3", "D-BC4"), c("P_bc_1", "P_bc_2", "P_bc_3", "P_bc_4"))]
deeplexicon[, mean(Barcode == Name)]

deeplexicon[, pred := factor(Barcode, levels = c("P_bc_1", "P_bc_2", "P_bc_3", "P_bc_4"))]
deeplexicon[, truth := factor(Name, levels = c("P_bc_1", "P_bc_2", "P_bc_3", "P_bc_4"))]

deeplexicon[, confusionMatrix(pred, truth)]


# Poreplex

bams <- list.files("/mnt/raid62/BetaCoV/Person/tangchao/analysis/DrictRNA/Poreplex/guppy_basecaller/output", "sorted.bam", full.names = TRUE)

lapply(bams, function(x) {
  bam <- readGAlignments(x, param = Rsamtools::ScanBamParam(what = c("qname", "flag", "mapq")), use.names = T)
  bam <- bam[with(bam, mapq == 60 & flag == 0)]
  data.table(ReadID = names(bam), Gene = as.character(seqnames(bam)), pred = gsub(".sorted.bam", "", basename(x)))
}) -> poreplex_output
poreplex_output <- do.call(rbind, poreplex_output)
poreplex_output[, ReadID := paste0("read_", ReadID)]

R2BC <- readRDS("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/06.DecodingDeePlexiConAndPoreplex/00.ModelTraining/00.RData/Poreplex_TestSet_Meta.Rds")
poreplex_output <- merge(poreplex_output, R2BC[, .(read, Name, gene)], by.x = "ReadID", by.y = "read")

setnames(poreplex_output, "Name", "obs")
poreplex_output[, pred := plyr::mapvalues(x = pred, from = c("BC1", "BC2", "BC3", "BC4"), to = c("P-BC1", "P-BC2", "P-BC3", "P-BC4"))]

table(poreplex_output$pred)

mean(poreplex_output$pred != "undetermined")

poreplex_output[, pred := factor(pred, levels = c("P-BC1", "P-BC2", "P-BC3", "P-BC4", "undetermined"))]
poreplex_output[, obs := factor(obs, levels = c("P-BC1", "P-BC2", "P-BC3", "P-BC4", "undetermined"))]

postResample(pred = poreplex_output$pred, obs = poreplex_output$obs)
mean(poreplex_output$pred == poreplex_output$obs)

mean(poreplex_output[pred != "undetermined", ]$pred == poreplex_output[pred != "undetermined", ]$obs)

confusionMatrix(data = poreplex_output$pred, reference = poreplex_output$obs)


# DecodeR-Poreplex

load(file = "/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/06.DecodingDeePlexiConAndPoreplex/00.ModelTraining2/01.Model/Poreplex_Classifier.RData")
TestData <- readRDS("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/06.DecodingDeePlexiConAndPoreplex/00.ModelTraining2/00.RData/Poreplex_TestSet.Rds")

ROC_Tab <- data.frame(obs = TestData$Class, 
                      predict(Fit1, TestData, type = "prob"), 
                      pred = predict(Fit1, newdata = TestData))
ROC_Tab <- as.data.table(ROC_Tab, keep.rownames = "read")





F1 <- function(cM) {
  TestAccuracy <- cM$overall[1]
  Sensitivity <- mean(na.omit(cM$byClass)[,1])
  Specificity <- mean(na.omit(cM$byClass)[,2])
  Precision <- mean(na.omit(cM$byClass)[,5])
  Recall <- mean(na.omit(cM$byClass)[,6])
  F1 <- mean(na.omit(cM$byClass)[,7])
  data.table(TestAccuracy, Sensitivity, Specificity, Precision, Recall, F1)
}


F1(DecodeR[, confusionMatrix(pred, truth)])
F1(deeplexicon[, confusionMatrix(pred, truth)])
F1(poreplex_output[, confusionMatrix(pred, obs)])
F1(ROC_Tab[, confusionMatrix(pred, obs)])

confM <- rbind(F1(DecodeR[, confusionMatrix(pred, truth)]), 
               F1(deeplexicon[, confusionMatrix(pred, truth)]),
               F1(ROC_Tab[, confusionMatrix(pred, obs)]),
               F1(poreplex_output[, confusionMatrix(pred, obs)]))
confM <- data.table(Method = c("DecodeR(DeePlexiCon)", "DeePlexiCon", "DecodeR(Poreplex)", "Poreplex"), confM)
openxlsx::write.xlsx(confM, "./analysis/06.DecodingDeePlexiConAndPoreplex/01.Comparison/confusionMatrix.xlsx", overwrite = T)

