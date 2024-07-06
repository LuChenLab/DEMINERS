library(data.table)
library(parallel)
library(caret)
library(doParallel)
library(pbapply)

setwd("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex")

for(i in seq(14, 24, 2)) {
  print(i)
  load(paste0("./analysis/11.Classifiers/01.ModelingData/Barcode_", i, ".RData"))
  rm(list = c("TrainingReads", "TestReads", "TestData")); gc()
  
  cl <- makePSOCKcluster(6)
  registerDoParallel(cl)
  fitControl <- trainControl(method = "repeatedcv",
                             number = 10,
                             repeats = 10)
  set.seed(825)
  Fit1 <- train(Class ~ ., data = TrainingData,
                preProc = c("center", "scale", "YeoJohnson", "nzv"),
                method = "rf",
                trControl = fitControl,
                verbose = FALSE,
                # to evaluate:
                # tuneGrid = expand.grid(mtry = 2),
                tuneLength = 1,
                metric = "Accuracy",
                allowParallel = TRUE)
  saveRDS(Fit1, file = paste0("./analysis/11.Classifiers/02.Classifier/Barcode_", i, ".Rds"))
  stopCluster(cl)
  rm(list = ls()); gc()
}


# 3102

for(i in rev(seq(18, 22, 2))) {
  print(i)
  load(paste0("/mnt/raid61/Personal_data/tangchao/Temp/01.ModelingData/Barcode_", i, ".RData"))
  rm(list = c("TrainingReads", "TestReads", "TestData")); gc()
  
  set.seed(3456)
  trainIndex <- createDataPartition(TrainingData$Class, 
                                    p = 10000 / min(table(TrainingData$Class)), 
                                    list = TRUE, 
                                    times = 1)
  TrainingData <- TrainingData[trainIndex[[1]], ]
  gc()
  
  cl <- makePSOCKcluster(10)
  registerDoParallel(cl)
  fitControl <- trainControl(method = "repeatedcv",
                             number = 10,
                             repeats = 5)
  set.seed(825)
  Fit1 <- train(Class ~ ., data = TrainingData,
                preProc = c("center", "scale", "YeoJohnson", "nzv"),
                method = "rf",
                trControl = fitControl,
                verbose = FALSE,
                # to evaluate:
                # tuneGrid = expand.grid(mtry = 2),
                tuneLength = 1,
                metric = "Accuracy",
                allowParallel = TRUE)
  saveRDS(Fit1, file = paste0("/mnt/raid61/Personal_data/tangchao/Temp/01.ModelingData/Barcode_", i, ".Rds"))
  stopCluster(cl)
  rm(list = ls()); gc()
}


setwd("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/11.Classifiers/02.Classifier")

files <- list.files(".", "Rds")
files <- paste0("Barcode_", seq(14, 22, 2), ".Rds")
for(i in files) {
  Fit1 <- readRDS(i)
  Fit1$trainingData <- NULL
  gdata::mv("Fit1", gsub(".Rds", "barcodes", gsub("^Barcode", "Model", i)))
}

Fit1 <- butcher::butcher(Fit1) # better


save(Model_2barcodes, Model_4barcodes, Model_6barcodes, Model_8barcodes, Model_10barcodes, Model_12barcodes, file = "Models.RData")
save(Model_2barcodes, file = "Model_2barcodes.RData")
save(Model_4barcodes, file = "Model_4barcodes.RData")
save(Model_6barcodes, file = "Model_6barcodes.RData")
save(Model_8barcodes, file = "Model_8barcodes.RData")
save(Model_10barcodes, file = "Model_10barcodes.RData")
save(Model_12barcodes, file = "Model_12barcodes.RData")
save(Model_14barcodes, file = "Model_14barcodes.RData")
save(Model_16barcodes, file = "Model_16barcodes.RData")
save(Model_18barcodes, file = "Model_18barcodes.RData")
save(Model_20barcodes, file = "Model_20barcodes.RData")
save(Model_22barcodes, file = "Model_22barcodes.RData")












load(file = "/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/04.RandomForest/01.ClassifierTraining/BIN_100_24barcodes_Classifier_V2.RData")
Fit1$trainingData <- NULL
gdata::mv("Fit1", "Model_24barcodes")
save(Model_24barcodes, file = "Model_24barcodes.RData")
