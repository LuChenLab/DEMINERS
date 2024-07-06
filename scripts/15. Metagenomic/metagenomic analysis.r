---
title: "Performance of classifiers"
author: "Junwei Song (crChao Tang)"
date: 'Report created: `r Sys.Date()`'
output: 
  html_document: 
    code_folding: "hide"
    toc: true
    toc_depth: 4
    toc_float: 
      collapsed: false
    number_sections: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, 
                      warning = FALSE, 
                      message = FALSE, 
                      fig.align = 'center')
knitr::opts_knit$set(root.dir = "/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/")

library(data.table)
library(parallel)
library(cowplot)
library(patchwork)
library(tibble)
library(Biostrings)
library(GenomicRanges)
library(GenomicAlignments)
library(ggpubr)
library(DecodeR)
library(ShortRead)

library(multiROC)
library(ggvenn)
library(ggplot2)
library(ggalluvial)
source("/mnt/raid61/Personal_data/songjunwei/AllScript/00Function/DIyplot.R")
setwd("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/")
```


# Karenk2
```{linux}
#Kraken2
## build index
/mnt/raid61/Personal_data/songjunwei/software/kraken2-master/install_path/kraken2-build --standard --threads 24 --db $DBNAME
#db can download https://benlangmead.github.io/aws-indexes/k2

cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/batch123/karenk2
kraken2="/mnt/raid61/Personal_data/songjunwei/software/kraken2-master/install_path/kraken2"
index="/mnt/raid61/Personal_data/songjunwei/reference/Kraken2_index/standard_all"
fq_b1="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/01Rawdata/batch1/guppy/DRS_hac.fastq.gz"
fq_b2="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/01Rawdata/batch2_fq/DRS_virus2.fastq"
fq_b3="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/01Rawdata/batch3_fq/DRS_multi3_6sample.fastq"

seqkit seq $fq_b1 --rna2dna > b1_dna.fq
seqkit seq $fq_b2 --rna2dna > b2_dna.fq
seqkit seq $fq_b3 --rna2dna > b3_dna.fq

for i in b1 b2 b3;
do
$kraken2 --db $index --threads 30 --classified-out ${i}_st_class.fq  --unclassified-out ${i}_st_unclass.fq --output ${i}_st_reads.txt --report ${i}_st_report.txt  \
  --use-names --confidence 0 --minimum-hit-groups 2 --report-zero-counts ${i}_dna.fq
done
#b1  9008 sequences classified (87.55%) 1281 sequences unclassified (12.45%)
#b2  48312 sequences classified (51.49%) 45518 sequences unclassified (48.51%)
#b3  464316 sequences classified (90.52%) 48634 sequences unclassified (9.48%)
```

```{r}
report_files <- list.files("batch123/karenk2/",pattern = "st_report.txt",full.names = T)
report_b2 <- kraken2_read_report(report_files[[2]])

lapply(report_files, function(x) {
  tmp <- data.table(kraken2_read_report(x))
  tmp[tax_count>1][order(tax_count,decreasing = T)]
})

report_list <- lapply(report_files, function(x){
  report_tmp <- data.table(kraken2_read_report(x))
  report_tmp <- report_tmp[clade_count!=0]
  tmp_G <- report_tmp[rank %in% grep("G",report_tmp$rank,value = T)][order(percent,decreasing = T)]
  tmp_S <- report_tmp[rank %in% grep("S",report_tmp$rank,value = T)][order(percent,decreasing = T)]
  tmp_S <- tmp_S[percent>0.1]
  tmp_S$sample <- gsub("b","Exp",gsub("_st_report.txt","",basename(x)))
  tmp_S
})
View(report_list[[3]])
report_tab <- rbind(report_list[[1]],report_list[[2]][percent>0.5],report_list[[3]][percent>1])
p1 <- ggplot(report_tab,aes(x=sample,y=percent,fill=name)) + geom_bar(stat = "identity") + mytheme+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "right") + xlab(NULL) + ylab("Percentage%")
p1
write.csv(rbindlist(report_list),"batch123/karenk2/exp123_kraken2_report.csv",row.names = F)

reads_files <- list.files("batch123/karenk2/",pattern = "st_reads.txt",full.names = T)
```

#BugSeq
##report
```{r fig.width=12,fig.height=6}
bugseq_report_files <- c(
    "batch123/bugseq/batch1/metagenomic_classification/kreports/DRS_hac.kreport",
    "batch123/bugseq/batch2/metagenomic_classification/kreports/DRS_virus2.kreport",
    "batch123/bugseq/batch3_5/metagenomic_classification/kreports/DRS_multi3_6sample.kreport")

bugseq_report_list <- lapply(bugseq_report_files, function(x){
  report_tmp <- data.table(kraken2_read_report(x))[clade_count>1][order(percent,decreasing = T)]
  tmp <- report_tmp[rank %in% grep("S",report_tmp$rank,value = T)][order(percent,decreasing = T)]  #S or G
  #tmp <- report_tmp[rank %in% grep("G",report_tmp$rank,value = T)][order(percent,decreasing = T)] 
  tmp
})
bugseq_report_list[[1]]$sample <- "Exp1"
bugseq_report_list[[2]]$sample <- "Exp2"
bugseq_report_list[[3]]$sample <- "Exp3"
bugseq_report_tab <- rbind(bugseq_report_list[[1]][clade_count>=10],bugseq_report_list[[2]][clade_count>=10],bugseq_report_list[[3]][clade_count>=10])
write.csv(bugseq_report_tab,"batch123/bugseq/bathc123_bugseq.csv",row.names = F)
table(bugseq_report_tab$sample)

bugseq_report_tab <- read.csv("batch123/bugseq/bathc123_bugseq.csv") #reads大于10的种S
ggplot(bugseq_report_tab,aes(x=sample,y=percent,fill=name)) + 
    geom_bar(stat = "identity",width = 0.5) + mytheme+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "right") + xlab(NULL) + ylab("Percentage%")
```

## Reads and Report
```{r}
bugseq_reads_files <- c(
    "batch123/bugseq/batch1/metagenomic_classification/read-based/details/DRS_hac.csv",
    "batch123/bugseq/batch2/metagenomic_classification/read-based/details/DRS_virus2.csv",
    "batch123/bugseq/batch3_5/metagenomic_classification/read-based/details/DRS_multi3_6sample.csv")
bugseq_reads_list <- lapply(bugseq_reads_files, function(x){
  tab <- read.csv(x)
  tab[,c("read_id","tax_id")] })
bugseq_reads_list[[1]]$sample <- "Exp1"
bugseq_reads_list[[2]]$sample <- "Exp2"
bugseq_reads_list[[3]]$sample <- "Exp3"
bugseq_reads_tab <- rbindlist(bugseq_reads_list)
table(bugseq_reads_tab$sample)

#######
bugseq_report_files <- c(
    "batch123/bugseq/batch1/metagenomic_classification/kreports/DRS_hac.kreport",
    "batch123/bugseq/batch2/metagenomic_classification/kreports/DRS_virus2.kreport",
    "batch123/bugseq/batch3_5/metagenomic_classification/kreports/DRS_multi3_6sample.kreport")

bugseq_report_list <- lapply(bugseq_report_files, function(x){
  tmp <- data.table(kraken2_read_report(x))[clade_count>1][order(percent,decreasing = T)]
  #tmp <- tmp[rank %in% grep("S",tmp$rank,value = T)][order(percent,decreasing = T)]  #S or G
  #tmp <- tmp[rank %in% grep("G",tmp$rank,value = T)][order(percent,decreasing = T)] 
  tmp
})
bugseq_report_list[[1]]$sample <- "Exp1"
bugseq_report_list[[2]]$sample <- "Exp2"
bugseq_report_list[[3]]$sample <- "Exp3"
bugseq_report_tab <- bugseq_report_list %>% rbindlist()

tax_id_name <- unique(bugseq_report_tab[,c("rank","tax_id","name")])
colnames(tax_id_name)[3] <- "tax_name"
table(bugseq_reads_tab$tax_id %in% tax_id_name$tax_id)

bugseq_reads_tab <- merge(bugseq_reads_tab,tax_id_name,by="tax_id")
write.table(bugseq_reads_tab,"batch123/Exp123_bugseq_res.txt")
```


```{r}
files <- list.files("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/02align/batch123_new/",full.names = T)
param <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isSupplementaryAlignment = FALSE,
                                                               isPaired = FALSE,
                                                               isSecondaryAlignment = FALSE),
                                 mapqFilter = 0,what=c("qname", "flag", "mapq"))

Mat0 <- mclapply(files, function(x) {
  map0 <- readGAlignments(x, param = param,use.names = T)
  data.frame(read= names(map0),seqname=map0@seqnames,Batch = gsub(".bam","",basename(x)),mapq= mcols(map0)$mapq, Length = qwidth(map0), Width = width(map0),flag = mcols(map0)$flag)
},mc.cores = 10)
minimap_tab <- Mat0 %>% rbindlist()
table(minimap_tab$flag)
hist(minimap_tab$mapq)

minimap_tab <- data.table(minimap_tab)
filter_highest_mapq <- function(data_table) {
  # 按 read 分组，然后获取每组中 mapq 最大的行
  result <- data_table[, .SD[which.max(mapq)], by = read]
  return(result)
}

filtered_minimap_tab <- filter_highest_mapq(minimap_tab)
table(duplicated(filtered_minimap_tab$read))
table(filtered_minimap_tab$Batch)
colnames(filtered_minimap_tab)[3] <- "minimap_name"

#add sample id
fq_files <- list.files("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/batch123/karenk2/",pattern = "dna.fq",full.names = T)
fq_read <- lapply(fq_files, function(x) {
  b1_fq <- readFastq(x)
  readid <- lapply(strsplit(b1_fq@id," "), function(x) x[1] %>% as.character()) %>% unlist()
  data.frame(read=readid,sample=basename(x))
  })
fq_read_tab <- fq_read %>% rbindlist()
table(fq_read_tab$sample)
fq_read_tab$sample <- gsub("_dna.fq","",fq_read_tab$sample)
#b1_dna.fq b2_dna.fq b3_dna.fq 
#    10289     93830    512950 
filtered_minimap_tab <- merge(filtered_minimap_tab,fq_read_tab,by="read",all.x=T)
filtered_minimap_tab[sample=="b1"]$sample <- "Exp1"
filtered_minimap_tab[sample=="b2"]$sample <- "Exp2"
filtered_minimap_tab[sample=="b3"]$sample <- "Exp3"
table(filtered_minimap_tab$sample)
#   b1     b2     b3 
#  8775  87920 447389 
filtered_minimap_tab$minimap_name[grep("S_enter",filtered_minimap_tab$minimap_name)] <- "Salmonella enteritidis"
filtered_minimap_tab$minimap_name[grep("S_cere",filtered_minimap_tab$minimap_name)] <- "Saccharomyces cerevisiae"
filtered_minimap_tab$minimap_name[grep("Ecoli",filtered_minimap_tab$minimap_name)] <- "Escherichia coli"
fwrite(filtered_minimap_tab,"batch123/Exp123_minimap2.txt",row.names = F)
```

# add DecodeR
```{r}
Pred <- readRDS("batch123/Exp123_DecodeR_new.rds")
Pred_exp1 <- Pred[[1]]

Pred_exp2 <- readRDS("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/07.VirusClassification/03.Merge/02.Prediction/Prediction_Table_20210825.Rds")
Pred_exp3 <- readRDS("/mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/analysis/07.VirusClassification/03.Merge/02.Prediction/Prediction_Table_20211008.Rds")


colnames(Pred_exp2)[1] <- c("Read")
colnames(Pred_exp2)[5] <- c("Barcode")
colnames(Pred_exp2)[6] <- c("Probability")
Pred_exp2$Read <- gsub("read_","",Pred_exp2$Read)
Pred_exp2$sample <- "Exp2"
colnames(Pred_exp3)[1] <- c("Read")
colnames(Pred_exp3)[8] <- c("Barcode")
colnames(Pred_exp3)[9] <- c("Probability")
Pred_exp3$Read <- gsub("read_","",Pred_exp3$Read)
Pred_exp3$sample <- "Exp3"

Pred_list <- list(Exp1=Pred_exp1,Exp2=Pred_exp2,Exp3=Pred_exp3)
saveRDS(Pred_list,"batch123/Exp123_DecodeR_n1o23.rds")
```

## map & pred
```{r}
species_merge2 <- fread("batch123/Exp123_bugseq_minimap2_unique.txt")
Pred <- readRDS("batch123/Exp123_DecodeR_new.rds") #new123
Pred_list <- readRDS("batch123/Exp123_DecodeR_n1o23.rds") # new1 old 123 use

mapp_pred_list <- lapply(c("Exp1","Exp2","Exp3"), function(x){
  merge(species_merge2[sample.y==x & mapq==60],Pred_list[[x]],by.x="read_id",by.y="Read")
})
names(mapp_pred_list) <- c("Exp1","Exp2","Exp3")
saveRDS(mapp_pred_list,"batch123/Exp123_mapped_DecodeR_n1o23.rds")


pred_sta_list <- lapply(mapp_pred_list, function(tab){
  tab <- data.table(tab)
  sta1 <- tab[, .N, by = c("sample.y", "minimap_name")][order(N, decreasing = T)][order(sample.y)]
  colnames(sta1)[3] <- "sum"
  sta2 <- tab[, .N, by = c("sample.y", "Barcode", "minimap_name")][order(Barcode, N, decreasing = T)][order(sample.y)]
  sta3 <- merge(sta1, sta2, by = c("sample.y", "minimap_name"))
  sta3$precent <- sta3$N / sta3$sum
  return(list(sta2=sta2,sta3=sta3))
})

target <- rbind(data.frame(exp="Exp1",Barcode=c("RTA-08","RTA-27","RTA-33"),Speceis=c("GETA","PEDV","SVV")),
            data.frame(exp="Exp2",Barcode=c("RTA-08","RTA-27","RTA-10","RTA-33"),Speceis=c("PbergheiANKA","PbergheiANKA","PRRSV","SVV")),
            data.frame(exp="Exp3",Barcode=c("RTA-10","RTA-16","RTA-17","RTA-24","RTA-32"), 
                       Speceis= c("PRRSV","Salmonella enteritidis","Escherichia coli",
                     "PbergheiANKA","Saccharomyces cerevisiae"))) %>% data.table()
target$id <- paste0(target$Barcode,"_",target$Speceis)

pred_sta <- lapply(c("Exp1","Exp2","Exp3"), function(x) {
  tab <- data.table(pred_sta_list[[x]]$sta3)
  tab$id <- paste0(tab$Barcode,"_",tab$minimap_name)
  target_sub <- target[exp == x]
  tab[id %in% target_sub$id]
}) %>% rbindlist()
mean(pred_sta$precent)
sum(pred_sta$N)
write.csv(pred_sta,"batch123/exp123_pred_map_res_sheet1.csv",row.names = F)
```
