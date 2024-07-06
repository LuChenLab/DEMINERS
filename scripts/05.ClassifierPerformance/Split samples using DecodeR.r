---
title: "Split samples using DecodeR "
author: "junwei"
date: "2022/9/2"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(changepoint)
library(data.table)
library(randomForest)
library(smoother)
library(caret)
library(rhdf5)
library(DecodeR)
library(data.table)
```

#split RTA
```{r}
Model4 <- readRDS("/mnt/raid61/Personal_data/tangchao/Temp/01.ModelingData/Barcode_4.Rds")
files0 <- list.files("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/01Rawdata/Tumor4/raw_f5/20211128_0852_MN31508_FAQ93968_a0437aad/fast5_pass", full.names = T)
files1 <- list.files("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/01Rawdata/Tumor4/raw_f5/20211129_0337_MN29097_FAQ89846_fe3787be/fast5_pass", full.names = T)
files2 <- list.files("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/01Rawdata/Tumor4/raw_f5/20211129_1201_MN29097_FAQ89846_2cd5a1e3/fast5_pass", full.names = T)
files <- c(files0,files1, files2)

lapply(files, function(x) {
  print(paste0(which(files == x), " of ", length(files)))
  pred <- DecodeR(fast5 = x, model = Model4, NT = 10)
  write.table(pred, paste0("/mnt/raid61/Personal_data/tangchao/Temp/SJW/", basename(x), ".txt"), sep = "\t", quote = F, row.names = F)
})


pred <- lapply(list.files("/mnt/raid61/Personal_data/tangchao/Temp/SJW", full.names = T), fread)
pred <- do.call(rbind, pred)

hist(pred$Probability, breaks = 100)
pred[, .N, Barcode]
pred[Probability > 0.6, .N, Barcode]

```

```{r}
library(PorexploreR)

# get example file from package
fast5file <- system.file("extdata/demo2_0.fast5", package = "DecodeR")

# load in the model, limited by file size only the 2 barcodes model were built into the package
data("Model_2barcodes")

# predict the barcode of example fast5 file
pred <- DecodeR(fast5 = fast5file, model = Model_2barcodes)
write.csv(pred,"/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/pred_decodeR.csv",row.names = F)
# histogram of predicted probability
hist(pred$Probability, xlab = "Probability", main = "Histogram of Probability")

# number of each barcode
table(pred$Barcode)

# set cutoff for unclassified read
pred2 <- DecodeR(fast5 = fast5file, model = Model_2barcodes, cutoff = 0.8)
table(pred2$Barcode)
```


#SARS
## 24sam
```{r}
Model24 <- readRDS("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/11.Classifiers/02.Classifier/Barcode_24.Rds")

files <- list.files("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/01Rawdata/batch5_SARS/230808_24sam/second_lib/",pattern = ".fast5",full.names = T,recursive = T)
files <- list.files("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/01Rawdata/batch5_SARS/230808_24sam/fast5/",pattern = ".fast5",full.names = T,recursive = T)
lapply(files, function(x) {
  print(paste0(which(files == x), " of ", length(files)))
  pred <- DecodeR(fast5 = x, model = Model24, NT = 20)
  write.table(pred, paste0("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/batch5_SARS/230808_24sam/dec_res", basename(x), ".txt"), sep = "\t", quote = F, row.names = F)
})
```


#0X Batch1234 profiling
```{r}
res123_path <- list.files("/mnt/raid5/Personal/minion/DRS_mul/other/nanom6a/",pattern = "*resFin",full.names = T)
filename <- gsub("_resFin","",basename(res123_path))
m6a_res123 <- m6a_process(res123_path,filename = filename)

res4_path <- list.files("/mnt/raid5/Personal/minion/DRS_mul/batch4/nanom6a/",pattern = "*result_final",full.names = T)
pb_filename <-c("batch4_Pb_tro1","batch4_Pb_tro2","batch4_Pb_tro3","batch4_Pb_sch1","batch4_Pb_sch2","batch4_Pb_sch3")
RTA <- c("RTA03","RTA10","RTA16","RTA17","RTA24","RTA32")
m6a_res4_n <- m6a_process(res4_path[1:6],filename = RTA)

test <- data.table(m6a_res4$ratio)[sample=="batch4_Pb_tro1"]
test[ID=="PbANKA_12_v3:1735993"]
test[which(duplicated(test$gID))]
length(unique(test$gID))

data.table(m6a_res4$geno_abun_m6a)[sample=="batch4_Pb_tro1"][,.N,by=gene][gene=="PBANKA_1101300"]
test1 <- data.table(m6a_res4$geno_abun_m6a)[sample=="batch4_Pb_tro1"][,.N,by=gID]
test2 <- merge(test,test1)
table(test2$m6aN==test2$N)

write.table(m6a_res4_n$ratio,"/mnt/raid61/Personal_data/songjunwei/DRS_RTA/01analysis/rep4_pb_6sample_m6atable.txt",quote = F,col.names = T,row.names = F)
```

##merge ratio
```{r }
res_all <- rbind(m6a_res123$ratio,m6a_res4$ratio)
res_all <-res_all[sample %in% res_all[,.N,by=sample][N>5]$sample]
res_all$sample <- gsub("batch","Rep",res_all$sample)
table(res_all$sample)
res_all[sample %in% c("Rep1_SVV","Rep2_PRRSV","Rep2_SVV"), species := "Virus"]
res_all[sample %in% c("Rep3_Ecoli","Rep3_S_enter"), species := "Bacteria"]
res_all[sample %in% c("Rep3_S_cere"), species := "Fungi"]
res_all[sample %in% c("Rep2_Pb_sch4","Rep4_Pb_sch1","Rep4_Pb_sch2","Rep4_Pb_sch3","Rep4_Pb_tro1","Rep4_Pb_tro2","Rep4_Pb_tro3"), species := "Plasmodium"]
table(is.na(res_all$species))
res_all[, ratio := as.numeric(ratio)]
res_all[, m6aN := as.numeric(m6aN)]
res_all$species <- factor(res_all$species)

order_name <- unique(res_all$sample)
order_name[c(1,4,3,5,7,6,8,9,10,11,12,13,2)]
res_all$sample <- factor(res_all$sample,levels = order_name[c(1,4,3,5,7,6,8,9,10,11,12,13,2)])
res_all$species <- factor(res_all$species,levels =c("Virus","Bacteria","Fungi","Plasmodium") )
write.table(res_all,"./m6a_all_table.txt",col.names = T,row.names = F)
```

```{r fig.height=5,fig.width=10}
pdf("/mnt/raid5/Personal/minion/DRS_mul/other/analysis/m6a_ratio_box.pdf",height = 5,width = 10)
compaired <- list(c("Rep2_SVV","Rep2_PRRSV"),c("Rep1_SVV","Rep2_PRRSV"))
ggplot(data = res_all,aes(x = sample, y = ratio)) +
  geom_boxplot(aes(fill = species))+
  geom_jitter(size=0.1)+
  theme_bw(base_size = 10) + 
    scale_colour_viridis_c(direction = -1) + 
    theme_bw() + 
    theme(axis.title = element_blank(), 
          axis.text.y = element_text(size = 12), 
          axis.text.x = element_text(size = 10, angle = 30, hjust = 1),
          legend.key.size = unit(15,"pt"))+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = T,test = "wilcox.test")
dev.off()
```

```{r}
m6a_yeast <- fread("/mnt/raid5/Personal/minion/Visium_DRS/DownData/m6A_Yeast_Basic_Site_Information.txt",fill = T)
table(m6a_yeast$V16)
```

## merge genomo
```{r fig.height=8,fig.width=10}
geno_all <- rbind(m6a_res123$geno_abun,m6a_res4$geno_abun)
geno_all$sample <- gsub("batch","Rep",geno_all$sample)

filename <- unique(res_all$sample)[c(1,4,3,5,7,6,8,9,10,11,12,13,2)]
kmer_fig <- lapply(filename, function(file){
  kmer <- geno_all[sample==file]$kmer
  kfig <- ggseqlogo(kmer, method = 'bits' )+labs(title = file)
  return(kfig)
})
pdf("/mnt/raid5/Personal/minion/DRS_mul/other/analysis/m6a_motif.pdf",height = 8,width = 10)
ggarrange(kmer_fig[[1]],kmer_fig[[2]],kmer_fig[[3]],kmer_fig[[4]],kmer_fig[[5]],
          kmer_fig[[6]],kmer_fig[[7]],kmer_fig[[8]],kmer_fig[[9]],kmer_fig[[10]],
          kmer_fig[[11]],kmer_fig[[12]],kmer_fig[[13]])
dev.off()
```

