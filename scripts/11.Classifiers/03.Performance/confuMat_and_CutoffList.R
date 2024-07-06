setwd("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex")
list.files("./analysis/11.Classifiers/02.Classifier", ".Rds")
library(caret)
library(data.table)

Perform <- mclapply(seq(2, 22, 2), FUN = function(i) {
  load(paste0("./analysis/11.Classifiers/01.ModelingData/Barcode_", i, ".RData"))
  Fit <- readRDS(paste0("./analysis/11.Classifiers/02.Classifier/Barcode_", i, ".Rds"))
  
  ROC_Tab <- data.frame(predict(Fit, TestData, type = "prob"), 
                        pred = predict(Fit, newdata = TestData))
  ROC_Tab <- as.data.table(ROC_Tab, keep.rownames = "read")
  ROC_Tab$PP <- apply(ROC_Tab[, grepl("RTA", colnames(ROC_Tab)), with = F], 1, max)
  ROC_Tab <- merge(ROC_Tab, TestReads, by = "read")
  cM <- ROC_Tab[, confusionMatrix(pred, Class)]
  
  if(i == 2) {
    TrainAccuracy <- Fit$results$Accuracy
    TestAccuracy <- cM$overall[1]
    Sensitivity <- cM$byClass[1]
    Specificity <- cM$byClass[2]
    Precision <- cM$byClass[5]
    Recall <- cM$byClass[6]
    F1 <- cM$byClass[7]
  } else {
    TrainAccuracy <- Fit$results$Accuracy
    TestAccuracy <- cM$overall[1]
    Sensitivity <- mean(cM$byClass[, 1])
    Specificity <- mean(cM$byClass[2])
    Precision <- mean(cM$byClass[5])
    Recall <- mean(cM$byClass[6])
    F1 <- mean(cM$byClass[7])
  }
  confuMat <- data.table(Barcodes = i, TrainAccuracy, TestAccuracy, Sensitivity, Specificity, Precision, Recall, F1)
  
  CutoffList <- do.call(rbind, lapply(1:10/10, FUN = function(ct) {
    data.table(Cutoff = ct, Unclassified = ROC_Tab[, mean(PP < ct) * 100], Accuracy = ROC_Tab[PP >= ct, mean(pred == Class) * 100])
  }))
  CutoffList <- data.table(Barcodes = i, CutoffList)
  
  list(confuMat, CutoffList)
}, mc.cores = 4)


confuMat <- do.call(rbind, lapply(Perform, function(x) x[[1]]))
CutoffList <- do.call(rbind, lapply(Perform, function(x) x[[2]]))

openxlsx::write.xlsx(confuMat1, "./analysis/11.Classifiers/03.Performance/confuMat1.xlsx")
openxlsx::write.xlsx(CutoffList1, "./analysis/11.Classifiers/03.Performance/CutoffList1.xlsx")







# 24

load("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/04.RandomForest/01.ClassifierTraining/RData/BIN_100_24barcodes_TestData_V2.RData")
load(file = "/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/04.RandomForest/01.ClassifierTraining/BIN_100_24barcodes_Classifier_V2.RData")

ROC_Tab <- data.frame(predict(Fit1, TestData, type = "prob"), 
                      pred = predict(Fit1, newdata = TestData))
ROC_Tab <- as.data.table(ROC_Tab, keep.rownames = "read")
ROC_Tab$PP <- apply(ROC_Tab[, grepl("RTA", colnames(ROC_Tab)), with = F], 1, max)
ROC_Tab <- merge(ROC_Tab, metaInfo, by = "read")
ROC_Tab[, pred := as.character(pred)]
ROC_Tab[, Barcode := as.character(Barcode)]
od <- sort(ROC_Tab[, unique(Barcode)])
ROC_Tab[, pred := factor(pred, levels = od)]
ROC_Tab[, Barcode := factor(Barcode, levels = od)]
cM <- ROC_Tab[, confusionMatrix(pred, Barcode)]


TrainAccuracy <- Fit1$results$Accuracy
TestAccuracy <- cM$overall[1]
Sensitivity <- mean(cM$byClass[, 1])
Specificity <- mean(cM$byClass[2])
Precision <- mean(cM$byClass[5])
Recall <- mean(cM$byClass[6])
F1 <- mean(cM$byClass[7])

confuMat <- data.table(Barcodes = 24, TrainAccuracy, TestAccuracy, Sensitivity, Specificity, Precision, Recall, F1)

CutoffList <- do.call(rbind, lapply(1:10/10, FUN = function(ct) {
  data.table(Cutoff = ct, Unclassified = ROC_Tab[, mean(PP < ct) * 100], Accuracy = ROC_Tab[PP >= ct, mean(pred == Barcode) * 100])
}))
CutoffList <- data.table(Barcodes = 24, CutoffList)

list(confuMat, CutoffList)

openxlsx::write.xlsx(confuMat, "./analysis/11.Classifiers/03.Performance/confuMat_24.xlsx")
openxlsx::write.xlsx(CutoffList, "./analysis/11.Classifiers/03.Performance/CutoffList_24.xlsx")






