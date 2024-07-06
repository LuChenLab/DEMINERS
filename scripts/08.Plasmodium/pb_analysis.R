---
title: "DecodeR_pb_m6a_analysis.Rmd"
author: "junwei"
date: "2021/11/12"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/rep4_pb_analysis/")
source("/mnt/raid61/Personal_data/songjunwei/AllScript/00Function/m6a_process.R")
source("/mnt/raid61/Personal_data/songjunwei/AllScript/00Function/Isoform_psi_func.R")
source("/mnt/raid61/Personal_data/songjunwei/AllScript/00Function/Ontology_hostgene.R")
source("/mnt/raid61/Personal_data/songjunwei/AllScript/00Function/DIyplot.R")
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)
library(GenomicRanges)
library(DESeq2)
library(sva)
library(data.table)
library(stringr)
library(tibble)
library(dplyr)

library(ggpubr)
library(ggplot2)
library(ggbio)
library(ggthemes)
library(seqLogo)
library(ggseqlogo)
library(pheatmap)
library(RColorBrewer)
library(patchwork)
library(FactoMineR)


library(plyr)
library(tidyr)
library("AnnotationDbi")
library("rtracklayer")
library("isoband")
library("infotheo")
# library("org.Pf.plasmo.db")
# library("org.PbANKA.eg.db")
# library("nlcor")

setwd("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/rep4_pb_analysis/")
```

# 0ref
```{r}
refpath="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/"
#gtf import
gtf_pb <- rtracklayer::import(paste0(refpath,"PbergheiANKA/PlasmoDB-53_PbergheiANKA.gff"))
gtf_pb_df <- as.data.frame(gtf_pb)
gtftmp <- unique(na.omit(data.frame(strand=gtf_pb_df$strand,gene=gtf_pb_df$gene_id)))


table(gtf_pb$type)
gtf_pb[,c("ID","Name")]
gtf_pb_gene <- gtf_pb[which(gtf_pb$type %in% c("mRNA","pseudogene","ncRNA_gene","ncRNA","tRNA","rRNA","snRNA","snoRNA" )), ][,c("type","ID","Parent","gene_id","Name")]
```

gene2transcrit
```{r}
library(data.table)
library(dplyr)
# grep ">" PlasmoDB-53_PbergheiANKA_AnnotatedTranscripts.fasta|awk -F '\\|' '{print $1,$2}'>ref_transcripts_name.tx
cdna_name <- fread("/mnt/raid5/Personal/minion/DRS_mul/ref/ref_transcripts_name.txt",sep = " ",fill = T,header = F)
names(cdna_name) <- c("transcript","gene")
cdna_name[,transcript:=gsub(">","",transcript)]
cdna_name[,gene:=gsub("gene=","",gene)]
g2t <- cdna_name %>% group_by(gene) %>% summarise(new_transcript = paste(transcript, collapse = " "))
write.table(g2t,"/mnt/raid5/Personal/minion/DRS_mul/ref/PlasmoDB-53_PbergheiANKA_gene2transcript1.txt",row.names = F,col.names = F,quote = F)
```

# 0data
```{r}
meta <- data.frame(sample=c("RTA03", "RTA10", "RTA16", "RTA17", "RTA24", "RTA32"),
           type=c("tro", "tro", "tro", "sch", "sch", "sch"),
           type1=c("Trophozoite", "Trophozoite", "Trophozoite", "Schizont", "Schizont", "Schizont"))
rownames(meta) <- paste0(meta$type,"_",meta$sample)
#meta <- meta[c(2:6),] #####2 vs 3 样本！
annotation_col <- meta
annotation_col$sample <- NULL
annotation_col$type <- NULL

setwd("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/rep4_pb_analysis/")
gene_tab <- readRDS("01gene_rawtab.rds")
table(gene_tab$sample)
isoform_tab <- readRDS("01isoform_rawtab.rds")
table(isoform_tab$sample)
m6A_Tab <- data.table(read.table("01m6a_res_rawtab.txt",header = T))
table(m6A_Tab$sample)
polya_tab <- read.table("01polyA_rawtab.txt")
table(polya_tab$sample)
Mat_List <- readRDS("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/rep4_pb_analysis/01Expresion_mat.rds")

gene_mat <- readRDS("01gene_mat.rds")
iso_mat <- readRDS("01isoform_mat.rds")
polya_mat <- readRDS("01polyA_mat.rds")
m6aN_mat <- readRDS("01m6a_N_mat.rds")
m6aR_mat <- readRDS("01m6a_Ratio_mat.rds")
ASSJ_mat <- readRDS("01SJ_PSI_mat.rds")
Mat_List <- list(gene=gene_mat,isoform=iso_mat,polya=polya_mat,m6aN=m6aN_mat,m6aR=m6aR_mat,ASSJ_psi=ASSJ_mat)
saveRDS(Mat_List,"./01Expresion_mat.rds")

Mat_List[[5]]
```

#01 make Mat
## isoform
```{r}
isoform_files <- list.files("/mnt/raid62/Personal_data/zhangdan/song_DRS/flair", pattern = ".isoform.read.map.txt", full.names = TRUE)
isoform <- lapply(isoform_files, function(x){
  tmp = read.table(x, stringsAsFactors = FALSE, header = FALSE)
  readid = strsplit(tmp$V2, ",")
  name_length = unlist(lapply(readid, function(x){length(x)}))
  name = rep(tmp$V1, name_length)
  res = data.table(readname = unlist(readid), isoform_gene = name)
  res2 = res[!grep(":", res$isoform_gene), ]
  tmp2 = strsplit(res2$isoform_gene, "_P")
  tmp3 = data.table(isoform = unlist(lapply(tmp2, function(x){x[[1]]})), gene = paste0("P", unlist(lapply(tmp2, function(x){x[[2]]}))))
  res2$isoform = tmp3$isoform
  res2$gene <- tmp3$gene
  res2$sample <- gsub(".isoform.read.map.txt","",gsub("flair.quantify.","",basename(x)))
  res2
})
iso_rb <- rbindlist(isoform)
saveRDS(iso_rb,"/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/rep4_pb_analysis//01isoform_rawlist.rds")

gene_tab <- data.table(rbindlist(isoform))[,.N,by=c("sample","gene")]
saveRDS(gene_tab,"./01gene_rawtab.rds")

isoform_tab <- data.table(rbindlist(isoform))[,.N,by=c("sample","isoform")]
isoform_tab$pbtyep <- gsub("RTA-[0-9][0-9].","",isoform_tab$sample)
table(isoform_tab$pbtyep)
isoform_tab$iso_type <- ifelse(isoform_tab$isoform %in% grep("PBANKA*",isoform_tab$isoform,value = T),"Anno","Novel")
table(isoform_tab$iso_type,isoform_tab$pbtyep)
saveRDS(isoform_tab,"./01isoform_rawtab.rds")




iso_mat <- dcast(data = isoform_tab,formula = isoform ~ sample,value.var = "N")
iso_mat[is.na(iso_mat)] <- 0
iso_mat %>% as.data.frame() %>% tibble::column_to_rownames("isoform") -> iso_mat
colnames(iso_mat) <- c("tro_RTA03" ,"tro_RTA10" ,"tro_RTA16", "sch_RTA17" ,"sch_RTA24" ,"sch_RTA32")
dim(iso_mat)
saveRDS(iso_mat,"./01isoform_mat.rds")
```

## gene 
```{r}
gene_tab <- readRDS("./01gene_rawtab.rds")
table(gene_tab$sample)
gene_mat <- dcast(data = gene_tab,formula = gene ~ sample,value.var = "N")
colnames(iso_mat) <- c("tro_RTA03" ,"tro_RTA10" ,"tro_RTA16", "sch_RTA17" ,"sch_RTA24" ,"sch_RTA32")
gene_mat[is.na(gene_mat)] <- 0
gene_mat %>% as.data.frame() %>% tibble::column_to_rownames("gene") -> gene_mat
colnames(gene_mat) <- rownames(meta)
dim(gene_mat)
saveRDS(gene_mat,"./01gene_mat.rds")
```

## polya
```{r}
files = list.files("/mnt/raid62/Personal_data/zhangdan/song_DRS/polya", pattern = "_results.tsv", full.names = TRUE)
polya_length <- lapply(files, function(x){
  tmp = read.table(x, stringsAsFactors = FALSE, header = TRUE)
  tmp = tmp[which(tmp$qc_tag == "PASS"), c("readname", "polya_length")]
  tmp$sample = gsub("_polya_results.tsv", "", basename(x))
  tmp
})
polya_tab = Reduce(rbind, polya_length)
table(polya_tab$sample)
polya_tab$type = plyr::mapvalues(from = metainfo$RTA.ID, metainfo$group,x = polya_length2$sample)

iso_rb <- readRDS("01isoform_rawlist.rds")
iso_rb$sample <- str_split_fixed(iso_rb$sample,"\\.",2)[,1]
polya_tab <- data.table(polya_tab)
iso_rb <- data.table(iso_rb)

new_polyA <- lapply(unique(iso_rb$sample), function(i){
  polyA = merge(polya_tab[sample==i], iso_rb[sample==i], by = c("readname","sample"))
})
new_polyA2 <- Reduce(rbind, new_polyA)
new_polyA2$group = paste0(new_polyA2$sample, ":", new_polyA2$gene)
write.table(new_polyA2,"01polyA_rawtab.txt")
```

```{r}
new_polyA2_median = aggregate(new_polyA2["polya_length"],by = new_polyA2[c("group")],FUN = median) #取miadn
new_polyA2_median$sample = unlist(lapply(strsplit(new_polyA2_median$group, ":"), function(x){x[[1]]}))
new_polyA2_median$gene = unlist(lapply(strsplit(new_polyA2_median$group, ":"), function(x){x[[2]]}))
new_polyA2_median$group <- NULL
new_polyA2_median <- new_polyA2_median[, c(2,3,1)]
new_polyA2_median_matrix <-dcast(data = new_polyA2_median, gene~sample )
rownames(new_polyA2_median_matrix) <- new_polyA2_median_matrix$gene
new_polyA2_median_matrix$gene <- NULL
new_polyA2_median_matrix[is.na(new_polyA2_median_matrix)] = 0
new_polyA2_median_matrix2 <- lapply(new_polyA2_median_matrix, function(x){tmp = as.integer(x)})
new_polyA2_median_matrix2 <- Reduce(cbind, new_polyA2_median_matrix2)
rownames(new_polyA2_median_matrix2) <- rownames(new_polyA2_median_matrix)
colnames(new_polyA2_median_matrix2) <- colnames(new_polyA2_median_matrix)
new_polyA2_median_matrix2 <- as.data.frame(new_polyA2_median_matrix2)
polya_mat <- new_polyA2_median_matrix2[rowSums(new_polyA2_median_matrix2[, 1:3]) != 0 & rowSums(new_polyA2_median_matrix2[, 4:6]) != 0, ]
colnames(polya_mat) <- rownames(meta)
polya_mat <- mat_filt(polya_mat,3,10)
saveRDS(polya_mat,"./01polyA_mat.rds")
```

## m6A 
```{r}
result_path <- list.files("/mnt/raid5/Personal/minion/DRS_mul/batch4/nanom6a/",pattern = "*result_final",full.names = T)
m6a_res <- m6a_process(resFin_path = result_path[1:6],cores = 20)
saveRDS(m6a_res,"./01m6a_ratio_geno.rds")

#input result
result_path <- list.files("/mnt/raid5/Personal/minion/DRS_mul/batch4/nanom6a//",pattern = "*result_final",full.names = T)
result_path <- paste0(result_path,"/ratio.0.5.tsv")
file.exists(result_path)
tsv <- lapply(result_path[1:6], function(x) fread(x, fill = T,sep = "\t",header = F))

tsv_split <- lapply(tsv,function(data){
  library(tidyr)
  tmp1 <- unite(data[,-1],col=pos,sep = "\t")
  tmp2 <- data.frame(gene=data$V1,pos=tmp1)
  tmp3 <- separate_rows(tmp2,pos,sep= "\t")
  tmp3 <- tmp3[which(tmp3$pos!=""),]
  library(stringr)
  tmp <-  cbind(data.frame(str_split_fixed(tmp3$gene, "\\|", 2)),data.frame(str_split_fixed(tmp3$pos, "\\|", 4)))
  names(tmp) <- c("gene","chr","pos","m6aN","readN","ratio")
  tmp$strand <- plyr::mapvalues(from = gtftmp$gene,to = gtftmp$strand,x=tmp$gene)
  tmp$strand <- gsub("1","+",tmp$strand)
  tmp$strand <- gsub("2","-",tmp$strand)
  tmp$name <- paste0(tmp$chr,":",tmp$pos,":",tmp$strand)
  return(tmp)
})

RTA <- c("RTA03","RTA10","RTA16","RTA17","RTA24","RTA32")
for (i in 1:6){
  tsv_split[[i]]$barcode <- RTA[i]
}
m6A_Tab <- rbindlist(tsv_split)
table(m6A_Tab$barcode)
m6A_Tab[barcode %in% c("RTA03","RTA10","RTA16"), type := "Trophozoite"]
m6A_Tab[barcode %in% c("RTA17","RTA24","RTA32"), type := "Schizont"]
m6A_Tab[, ratio := as.numeric(ratio)]
m6A_Tab[, m6aN := as.numeric(m6aN)]
m6A_Tab[, readN := as.numeric(readN)]
m6A_Tab[, pos := as.numeric(pos)]

### check gene name!!
m6A_gr <- GRanges(seqnames = m6A_Tab$chr,ranges = IRanges(start = m6A_Tab$pos,end = m6A_Tab$pos),strand = m6A_Tab$strand)
find <- findOverlaps(m6A_gr,gtf_pb_gene)
m6A_Tab$gene_gtf <- NA
m6A_Tab$gene_gtf[queryHits(find)] <- gtf_pb_gene[subjectHits(find)]$Parent
idx_na <- is.na(m6A_Tab$gene_gtf)
m6A_Tab$gene_gtf[idx_na] <- m6A_Tab[idx_na]$gene
m6A_Tab$gene_gtf <- as.character(m6A_Tab$gene_gtf)
table(m6A_Tab$gene==m6A_Tab$gene_gtf)
m6A_Tab$sample <- paste0(m6A_Tab$type,"_",m6A_Tab$barcode)
m6A_Tab$id <- paste0(m6A_Tab$gene_gtf,"..",gsub("PbANKA_","",m6A_Tab$name))
m6A_Tab$id_s <- paste0(m6A_Tab$gene_gtf,"..",gsub("PbANKA_","",m6A_Tab$name),m6A_Tab$sample)
m6A_Tab <- m6A_Tab[!duplicated(m6A_Tab$id_s)]
m6A_Tab$id_s <- NULL
write.table(m6A_Tab,"01m6a_res_rawtab.txt",row.names = F,sep = "\t",quote = F)

m6A_Tab[,.(sum(readN)),by=barcode]
```



### N
```{r}
m6A_Tab <- data.table(read.table("01m6a_res_rawtab.txt",header = T))
table(m6A_Tab$sample)
dim(m6A_Tab)
m6aN_mat <- dcast(data = m6A_Tab, formula = id ~ barcode, value.var = "m6aN")
m6aN_mat[is.na(m6aN_mat)] <- 0
m6aN_mat %>% as.data.frame() %>% tibble::column_to_rownames("id") -> m6aN_mat
colnames(m6aN_mat) <- c("tro_RTA03" ,"tro_RTA10" ,"tro_RTA16", "sch_RTA17" ,"sch_RTA24" ,"sch_RTA32")
dim(m6aN_mat)
m6aN_mat <- mat_filt(m6aN_mat,n=2)
saveRDS(m6aN_mat,"./01m6a_N_mat.rds")
```

### Ratio
```{r}
m6A_Tab[, ratio := as.numeric(ratio)]
m6aR_mat <- dcast(data = m6A_Tab, formula = id ~ barcode, value.var = "ratio")
m6aR_mat %>% as.data.frame() %>% tibble::column_to_rownames("id") -> m6aR_mat
# hv <- data.table(ID = rownames(m6aR_mat), sd = apply(m6aR_mat, 1, sd, na.rm = T))
# hv <- hv[!is.na(sd)][order(sd, decreasing = T)]
# m6aR_mat[hv$ID,]
m6aR_mat[is.na(m6aR_mat)] <- 0
m6aR_mat <- mat_filt(m6aR_mat,n=2)
colnames(m6aR_mat) <- c("tro_RTA03" ,"tro_RTA10" ,"tro_RTA16", "sch_RTA17" ,"sch_RTA24" ,"sch_RTA32")
saveRDS(m6aR_mat,"./01m6a_Ratio_mat.rds")
```

##SJ
```{r}
SJs <- list.files("/mnt/raid61/Personal_data/tangchao/Temp/20211025/05.TranscriptClean", "SJ.out.tab$", full.names = T)
meta <- data.frame(Stage = c("Trophozoite", "Trophozoite", "Trophozoite", "Schizont", "Schizont", "Schizont"),
                   row.names = c("RTA-03", "RTA-10", "RTA-16", "RTA-17", "RTA-24", "RTA-32"))
SJs_list <- lapply(SJs, function(x) {
  tmp <- read.table(x)
  colnames(tmp) <- c("chr","start","end","strand","intron_motif","annotation","uniq_map_num","multi_map_num","overhang")
  tmp
})
load(file = "/mnt/raid61/Personal_data/tangchao/Temp/PSI_Tab.RData")
dim(PSI_Tab)
ASSJ_mat <- data.frame(PSI_Tab)
ASSJ_mat[is.na(ASSJ_mat)] <- 0
colnames(ASSJ_mat) <- c("tro_RTA03" ,"tro_RTA10" ,"tro_RTA16", "sch_RTA17" ,"sch_RTA24" ,"sch_RTA32")
saveRDS(ASSJ_mat,"./01SJ_PSI_mat.rds")
```

#02 statistics
##pca
直接跑
```{r fig.width=10, fig.height=3}
# gene_mat <- mat_filt(gene_mat, 5, 40) 
pca_list <- lapply(names(Mat_List),function(x){
  mat <- Mat_List[[x]]
  mat <- mat[,rownames(meta)]
  pca_res <- PCA(t(mat), ncp = 3, graph = F)
  pca_result <- data.frame(pca_res$svd$U, Name = colnames(mat))
  pca_result$Stage <- plyr::mapvalues(pca_result$Name, c("tro_RTA03", "tro_RTA10", "tro_RTA16", "sch_RTA17", "sch_RTA24", "sch_RTA32"), 
                                    c("Trophozoite", "Trophozoite", "Trophozoite", "Schizont", "Schizont", "Schizont"))
  p1 <- ggplot(pca_result, aes(x = X1, y = X2, color = Stage)) +
    geom_point(size = 2) + #Size and alpha just for fun
    ggrepel::geom_text_repel(aes(label = Name)) +
    theme_bw() +
    xlab(paste("PC1(", round(pca_res$eig[, 2][1], 2), "%)", sep = "")) +
    ylab(paste("PC2(", round(pca_res$eig[, 2][2], 2), "%)", sep = "")) + ggtitle(x)
  p1
})
ggarrange(plotlist = pca_list)

pca_list <- lapply(names(Mat_List)[c(1,3,4)],function(x){
  mat <- Mat_List[[x]]
  mat <- mat[,rownames(meta)]
  pca_res <- PCA(t(mat), ncp = 3, graph = F)
  pca_result <- data.frame(pca_res$svd$U, Name = colnames(mat))
  pca_result$Stage <- plyr::mapvalues(pca_result$Name, c("tro_RTA03", "tro_RTA10", "tro_RTA16", "sch_RTA17", "sch_RTA24", "sch_RTA32"), 
                                    c("Trophozoite", "Trophozoite", "Trophozoite", "Schizont", "Schizont", "Schizont"))
  p1 <- ggplot(pca_result, aes(x = X1, y = X2, color = Stage)) +
    geom_point(size = 2) + #Size and alpha just for fun
    #ggrepel::geom_text_repel(aes(label = Stage)) +
    theme_bw() + scale_color_manual(values =c("#d05e15","#239873"))+
    xlab(paste("PC1(", round(pca_res$eig[, 2][1], 2), "%)", sep = "")) +
    ylab(paste("PC2(", round(pca_res$eig[, 2][2], 2), "%)", sep = "")) + ggtitle(x)
  p1
})

plot <- ggarrange(plotlist = pca_list,nrow = 1,common.legend = T)
ggsave(filename = "./02pb_PCA_3.pdf",plot = plot,width = 10,height = 3)
```
唐老师代码跑
```{r}
Meth <-  data.table(read.table("01m6a_res_rawtab.txt",header = T))
Meth[, ratio := as.numeric(ratio)]
# hv <- Meth[, sd(ratio, na.rm = T), name][!is.na(V1)][order(V1, decreasing = T)][1:50, name]
MethTab <- dcast(data = Meth, formula = name + gene ~ barcode, value.var = "ratio")
# MethTab <- MethTab[name %in% hv]
MethTab <- data.frame(MethTab[, -c(1:2)], row.names = paste0(MethTab[[1]], "_", MethTab[[2]]))

m6aR_mat[m6aR_mat==0]=NA
hv <- data.table(ID = rownames(m6aR_mat), sd = apply(m6aR_mat, 1, sd, na.rm = T))
hv <- hv[!is.na(sd)][order(sd, decreasing = T)]
length(hv)

m6aR_mat[is.na(m6aR_mat)] <- 0

library(FactoMineR)
pca_res <- PCA(t(m6aR_mat[hv$ID, ]), ncp = 3, graph = F)
library(FactoMineR)
pca_result <- data.frame(pca_res$svd$U, Name = colnames(MethTab))
pca_result$name <- gsub("RTA", "RTA-", pca_result$Name)
pca_result$Stage <- plyr::mapvalues(pca_result$Name, c("RTA03", "RTA10", "RTA16", "RTA17", "RTA24", "RTA32"), 
                                    c("Trophozoite", "Trophozoite", "Trophozoite", "Schizont", "Schizont", "Schizont"))
library(ggrepel)
ggplot(pca_result, aes(x = X1, y = X2, color = Stage))+
  geom_jitter(size = 2) + #Size and alpha just for fun
  geom_text_repel(aes(label = name)) + 
  theme_bw(base_size = 15) +
  xlab(paste("PC1(", round(pca_res$eig[,2][1]), "%)", sep = "")) +
  ylab(paste("PC2(", round(pca_res$eig[,2][2]), "%)", sep = "")) + 
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 3, name = "Dark2")[1:2]) + 
  theme(legend.position = "top")
ggsave("/mnt/raid61/Personal_data/tangchao/Temp/20211025/05.TranscriptClean/m6A_AS/m6A_PCA.pdf", width = 4.2, height = 4.2)
```

##isoform
```{r fig.width=6,fig.height=5}
isoform_tab_uniq <- unique(isoform_tab[,c("isoform","pbtyep","iso_type")])
table(isoform_tab_uniq$pbtyep)
table(isoform_tab_uniq$pbtyep,isoform_tab_uniq$iso_type)

iso_list1 <- list(Trophozoite=isoform_tab[pbtyep=="Trophozoite"]$isoform,
                 Schizont=isoform_tab[pbtyep=="Schizont"]$isoform)
p1 <- ggvenn::ggvenn(iso_list1,fill_color = c("#D05E15","#239873")) +labs(title = "All isoform")

novel_iso <- data.table(isoform_tab[isoform_tab$iso_type=="Novel",])
iso_list2 <- list(Trophozoite=novel_iso[pbtyep=="Trophozoite"]$isoform,
                 Schizont=novel_iso[pbtyep=="Schizont"]$isoform)
p2 <- ggvenn::ggvenn(iso_list2,fill_color = c("#D05E15","#239873")) +labs(title = "Novel isoform")

pdf("./02isoform_venn.pdf",6,5)
p1+p2
dev.off()
```

##m6a
```{r}
m6A_Tab <- data.table(read.table("01m6a_res_rawtab.txt",header = T))
m6A_Tab_uniq <- data.table(unique(m6A_Tab[,c("name","type")]))
table(m6A_Tab_uniq$type)
ggvenn::ggvenn(list(tro=m6A_Tab_uniq[type=="tro"]$name,sch=m6A_Tab_uniq[type=="sch"]$name))
```




#03 Analysis
## isoform
```{r}
iso_mat <- iso_mat[,rownames(meta)]
# iso_mat <- mat_filt(iso_mat, n = 5, s = 40) # 40 improve p-value
# rownames(iso_mat) <- plyr::mapvalues(x=rownames(iso_mat),from = gtf_pb_gene$ID,to = gtf_pb_gene$Name)
print(dim(iso_mat))
ddsMat <- DESeqDataSetFromMatrix(countData = iso_mat,
                                 colData = meta,
                                 design = ~ type)
dds <- DESeq(ddsMat)
res <- results(dds)
print(res[grep("^PBANKA_0406650",rownames(iso_mat),value=T), ])


iso_res <- subset(res,pvalue<0.05)
iso_res$type <- ifelse(iso_res$log2FoldChange>0,"Trophozoite","Schizont")
table(iso_res$type)
iso_res <- data.frame(iso_res) %>% rownames_to_column("isoform") 
novel_l <- unique(iso_rb[,c("isoform","gene")])
iso_res$gene <- plyr::mapvalues(x=iso_res$isoform,from = novel_l$isoform,to = novel_l$gene)
write.csv(iso_res,"./03diff_iso_5PBsam_p005.csv",row.names = F)

iso_vst <- rlog(ddsMat)
assay(iso_vst)["PBANKA_0406650",]
pheatmap(
  t(assay(iso_vst)[iso_res$isoform,]),
  scale = "column",
  show_rownames = T,
  cluster_cols = T,
  border_color = NA,
  color = colorRampPalette(c("#188AD7", "white", "#DA222C"))(100)
  )

```

## gene
```{r fig.width=6,fig.height=3}
gene_mat <- gene_mat[,rownames(meta)]
gene_mat <- mat_filt(gene_mat, n = 5, s = 40) # 40 improve p-value
#rownames(gene_mat) <- plyr::mapvalues(x=rownames(gene_mat),from = gtf_pb_gene$ID,to = gtf_pb_gene$Name)
print(dim(gene_mat))
ddsMat <- DESeqDataSetFromMatrix(countData = gene_mat,
                                 colData = meta,
                                 design = ~ type)
gene_vst <- rlog(ddsMat)
dds <- DESeq(ddsMat)
res <- results(dds)
print(res["PBANKA_0406650", ])


gene_res <- subset(res,pvalue<0.05)
gene_res$type <- ifelse(gene_res$log2FoldChange>0,"Trophozoite","Schizont")
table(gene_res$type)
gene_res <- data.frame(gene_res) %>% rownames_to_column("gene") 
write.csv(gene_res,"./03diff_gene_5PBsam_p005.csv",row.names = F)


assay(gene_vst)["PBANKA_0406650",]

pdf("./05n5s40_heatmap_pb_genemv03_res_p005_t.pdf",6,3)
pheatmap(
  t(assay(gene_vst)[gene_res$gene,]),
  scale = "column",
  show_rownames = T,
  cluster_cols = T,
  border_color = NA,
  color = colorRampPalette(c("#188AD7", "white", "#DA222C"))(100)
  )
# pheatmap(
#   assay(gene_vst)[gene_res_sub$gene,],
#   scale = "row",
#   show_rownames = T,
#   cluster_cols = T,
#   annotation_col = annotation_col,border_color = NA,
#   # annotation_colors = ann_colors,
#   color = colorRampPalette(c("#188AD7", "white", "#DA222C"))(100)
#   )
dev.off()

write.csv(gene_res,"./05n5s40_pb_genemv03_res_p005.csv",row.names = F)
```


```{r fig.width=6,fig.height=3}
gene_feat <- c("PBANKA_0305300","PBANKA_1338600","PBANKA_1354400","PBANKA_1406800")

exp <- assay(gene_vst)[gene_feat,]
p <- pheatmap(
  exp,
  scale = "row",
  show_rownames = T,
  cluster_cols = T,
  border_color = NA,
  color = colorRampPalette(c("#188AD7", "white", "#DA222C"))(100)
  )
saveRDS(list(exp=exp,plot=p),"/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/rep4_pb_analysis/gene_heatmap.rds")
```

## polya 
```{r fig.width=6,fig.height=3}
polya_mat <- readRDS("01polyA_mat.rds")
polya_mat <- polya_mat[,rownames(meta)]
dim(polya_mat)

#rownames(polya_mat) <- plyr::mapvalues(x=rownames(polya_mat),from = gtf_pb_gene$ID,to = gtf_pb_gene$Name)
polya_mat <- mat_filt(polya_mat,3,0)
polya_mat["PBANKA_0406650", ]
dim(polya_mat)
table(rownames(polya_mat) %in% rownames(gene_mat))
ddsMat <- DESeqDataSetFromMatrix(countData = polya_mat,
                         colData = meta,
                         design = ~ type)
dds <- DESeq(ddsMat)
res <- results(dds)
print(res["PBANKA_0406650", ])
polya_res <- subset(res,pvalue<0.05)
polya_res$type <- ifelse(polya_res$log2FoldChange>0,"Trophozoite","Schizont")
table(polya_res$type)
polya_res <- data.frame(polya_res) %>% rownames_to_column("gene") 
write.csv(polya_res,"./03diff_polya_5PBsam_p005.csv",row.names = F)



pdf("./05n5s40_5heatmap_pb_polyamv03_res_p005_t.pdf",6,3)
polya_vst <- rlog(ddsMat)
pheatmap(
  t(assay(polya_vst)[polya_res$gene,]),
  scale = "column",
  show_rownames = T,
  cluster_cols = T,
  #annotation_col = annotation_col,
  border_color = NA,
  # annotation_colors = ann_colors,
  color = colorRampPalette(c("#188AD7", "white", "#DA222C"))(100)
  )
dev.off()

write.csv(polya_res,"./05n5s40_pb_polyamv03_res_p005.csv",row.names = F)
```

## m6A_N
```{r fig.width=5,fig.height=5}
m6aN_mat <- readRDS("01m6a_N_mat.rds")
m6aN_mat <- m6aN_mat[,rownames(meta)]
dim(m6aN_mat)
m6aN_mat[grep("^PBANKA_0406650*",rownames(m6aN_mat),value = T), ]
#m6aN_mat <- mat_filt(m6aN_mat,3,0)
ddsMat <- DESeqDataSetFromMatrix(countData = m6aN_mat,
                         colData = meta,
                         design = ~ type)
dds <- DESeq(ddsMat)
res <- results(dds)
print(res[grep("^PBANKA_1245821*",rownames(m6aN_mat),value = T), ])
res_sub <- subset(res,pvalue<0.05)

res_sub$type <- ifelse(res_sub$log2FoldChange>0,"Trophozoite","Schizont")
table(res_sub$type)
m6aN_res_sub <- data.frame(res_sub) %>% rownames_to_column("id") 
str_tmp <- str_split_fixed(m6aN_res_sub$id,"\\.\\.",2)
m6aN_res_sub$gene <- str_tmp[,1]
write.csv(m6aN_res_sub,"./03diff_m6aN_6PBsam_p005.csv",row.names = F)

write.csv(m6aN_res_sub,"./05n5s40_pb_m6aNmv03_res_p005.csv",row.names = F)


m6aN_vst <- rlog(ddsMat)
assay(m6aN_vst)["PBANKA_0406650_04_v3:245948:+", ]
pdf("./05heatmap_pb_m6aNmv03_res_p005.pdf",5,5)
pheatmap(
  assay(m6aN_vst)[m6aN_res_sub$id,],
  scale = "row",
  show_rownames = F,
  cluster_cols = T,
  annotation_col = annotation_col,border_color = NA,
  # annotation_colors = ann_colors,
  color = colorRampPalette(c("#188AD7", "white", "#DA222C"))(100)
  )
dev.off()
```



## m6aRatio
```{r fig.width=5,fig.height=5}
m6aR_mat <- readRDS("./01m6a_Ratio_mat.rds")
m6aR_mat <- m6aR_mat[,rownames(meta)]
RawTab <- m6A_Tab[,c("sample","id","m6aN","readN")]
RawTab$sample <- gsub("Trophozoite","tro",RawTab$sample)
RawTab$sample <- gsub("Schizont","sch",RawTab$sample)
fm_iso <- FindMarkers_PII(
  PIIratio = m6aR_mat,
  meta = meta,
  design = "type",
  ident.1 = "tro",
  ident.2 = "sch",
  min.pct = 0,
  delta.threshold = 0,
  only.pos = F,
  min.cells.group = 0,
  test.use = "prop",
  RawTab =  RawTab,
  NT = 10
)
fm_iso <- FindMarkers_PII(
  PIIratio = m6aR_mat,
  meta = meta,
  design = "type",
  ident.1 = "tro",
  ident.2 = "sch",
  min.pct = -Inf,
  delta.threshold = 0,
  only.pos = F,
  min.cells.group = 0,
  test.use = "t.test",
  NT = 10
)
res_sub<-subset(fm_iso,p_val < 0.05)
hist(res_sub$p_val)
res_sub$type <- ifelse(res_sub$deltaPSI<0,"Schizont","Trophozoite")
table(res_sub$type)
m6aR_res_sub <- data.frame(res_sub) %>% rownames_to_column("id") 
str_tmp <- str_split_fixed(m6aR_res_sub$id,"\\.\\.",2)
m6aR_res_sub$gene <- str_tmp[,1]
write.csv(m6aR_res_sub,"./03diff_m6a_Ratio_6PBsam_p005.csv",row.names = F)


pdf("./05heatmap_pb_m6aR_res_p005.pdf",5,5)
pheatmap(
  m6aR_mat[m6aR_res_sub$id, ],
  scale = "row",
  show_rownames = F,
  cluster_cols = T,
  annotation_col = annotation_col,
  border_color = NA,
  # annotation_colors = ann_colors,
  color = colorRampPalette(c("#188AD7", "white", "#DA222C"))(100)
  )
dev.off()
write.csv(m6aR_res_sub,"./05pb_m6aR_res_p005.csv",row.names = F)
```

#04 example box
```{r}
resfile <- list.files("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/rep4_pb_analysis/",pattern = "03diff_*",full.names = T)
res1 <- lapply(resfile, function(x){
  tmp <- data.table(read.csv(x))
  tmp1 <- tmp[,.N,by="type"]
  tmp2 <- t(tmp1) %>% data.frame() 
  colnames(tmp2) <- tmp2[1,]
  tmp3 <- data.frame(tmp2[-1,])
  colnames(tmp3) <- tmp2[1,]
  rownames(tmp3) <- gsub(".csv","",gsub("03diff_","",basename(x)))
  tmp3
})
res2 <- rbindlist(res1,fill = T)
res2[is.na(res2)]=0
res2$Diff_type <- gsub(".csv","",gsub("03diff_","",basename(resfile)))
```

## 3 vs 3
```{r}
res_tab <- lapply(resfile[c(2,4,5,7,9)], function(x) data.table(read.csv(x)))
names(res_tab) <- gsub("_p005.csv","",gsub("03diff_","",basename(resfile[c(2,4,5,7,9)])))
res_gene <- res_tab$gene_6PBsam
res_m6a <- res_tab$m6aN_6PBsam
res_m6aR <- res_tab$m6a_Ratio_6PBsam
table(res_m6a$type)

pheatmap(m6aR_mat[res_m6aR$id,])

## m6a counts
# gene_feat <- intersect(res_gene[type=="Trophozoite"]$gene,res_m6a$gene)
# m6a_feat <- res_m6a[gene %in% gene_feat]$id
# 
# assay(gene_vst)[gene_feat,]
# pheatmap(assay(gene_vst)[gene_feat,],scale = "row",cluster_cols = F,cluster_rows = F)
# assay(m6aN_vst)[m6a_feat,]
# pheatmap(assay(m6aN_vst)[m6a_feat,],scale = "row",cluster_cols = F,cluster_rows = F)
# m6aR_mat[m6a_feat,]
# pheatmap(m6aR_mat[m6a_feat,],cluster_cols = F,cluster_rows = F)


## m6a ratio 
gene_feat <-
  c(
    intersect(res_gene[type == "Trophozoite"]$gene, res_m6aR[type == "Schizont"]$gene),
    intersect(res_gene[type == "Schizont"]$gene, res_m6aR[type == "Trophozoite"]$gene)
  )
## input pathway
# GOinfo <- readRDS("/mnt/raid61/Personal_data/yangqingxin/data/Gene_ontology_info/PBANKA_GO_info_from_PlasmoDBv56/PBANK_GO_info_from_PlasmoDB.Rds")
# gene_feat <- intersect(gene_feat,GOinfo$PlasmoDB_name)
# lapply(gene_feat, function(x) gtf_pb[gtf_pb$gene_id %in% x,])
m6a_feat <- res_m6aR[gene %in% gene_feat]$id
#都有多个exon
m6aR_mat[grep("PBANKA_1245821",rownames(m6aR_mat)),]
m6A_Tab[gene_gtf=="PBANKA_1245821"]
m6A_Tab[id=="PBANKA_1245821..12_v3:1735101:+"]
```


```{r}
#box plot
ddsMat <- DESeqDataSetFromMatrix(countData = gene_mat,
                                 colData = meta,
                                 design = ~ type)
gene_vst <- rlog(ddsMat)
# exp <- rbind(assay(gene_vst)[gene_feat, ],
#              assay(m6aN_vst)[m6a_feat,])
exp <- rbind(assay(gene_vst)[gene_feat, ],
             m6aR_mat[m6a_feat,])
data <- data.frame(exp) %>% gather()
data$Name <- rep(c(gene_feat,m6a_feat),6)
data$type <- ifelse(1:nrow(data) %in% grep("tro",data$key),"Trophozoite","Schizont")
data$type <- factor(data$type,levels = c("Trophozoite","Schizont"))
# data$compare <- NA
# for (i in 1:length(gene_feat)) {
#   data$compare[grep(gene_feat[i],data$Name)] <- i
# }
data$tech <- ifelse(1:nrow(data) %in% grep(":",data$Name),"m6A","gene")
data <- data.table(data)
tmp <- data

# lapply(split(data,f = data$compare), function(tmp){
  pgene <- ggplot(data = tmp[tech=="gene"], aes(x = type, y = value, fill = type)) +
    geom_boxplot() + ylim(0,10.5)+
    geom_point(position = position_jitterdodge())+
    theme(legend.position = 'none')+xlab(label = NULL)+
    labs(
      y = "Normalized gene expression",
      title = unique(tmp[tech=="gene"]$Name),
      subtitle =  paste0("Wald.test=", round(res_gene[res_gene$gene == unique(tmp[tech=="gene"]$Name), "pvalue"], 6))
    )+ scale_fill_manual(values = c("#239873","#D05E15"))+mytheme 
  
  pm6a <- ggplot(data = tmp[tech=="m6A"], aes(x = type, y = value, fill = type)) +
    geom_boxplot() +
    geom_point(position = position_jitterdodge()) +
    theme(legend.position = 'none')+xlab(label = NULL)+
    mytheme +  labs(
      y = "m6A Ratio",
      title = unique(tmp[tech=="m6A"]$Name),
      subtitle =  paste0("t.test=", round(res_m6aR[res_m6aR$id == unique(tmp[tech=="m6A"]$Name), "p_val"], 6))
    )+ scale_fill_manual(values = c("#239873","#D05E15"))
  p <- pgene+pm6a
# })

ggsave(filename = "./04Pb1234821_box3.pdf",plot = p,width = 7.5,height = 5)
```




## motif
```{r}
m6a_res <- readRDS("./01m6a_ratio_geno.rds")
m6a_res$geno_abun$IDs <- paste0(m6a_res$geno_abun$chr,"_",m6a_res$geno_abun$pos)
geno_tab <- data.table(m6a_res$geno_abun)[gene==gene_feat][IDs=="PbANKA_12_v3_1735101"]
table(geno_tab$kmer)
table(geno_tab$sample)
pmotif <- ggseqlogo(geno_tab$kmer, method = 'prob') 
ggsave(plot = pmotif,filename = "./04Pb1234821_motif.pdf",width = 4,height = 3)
```


##IGV
```{r}
gene_feat="PBANKA_1245821" #PBANKA_0406650(2vs3)
iso_table <- data.table(readRDS("./01isoform_rawlist.rds"))
iso_table$sample <- str_split_fixed(iso_table$sample,"\\.",2)[,1]
iso_table <- iso_table[gene==gene_feat]
table(iso_table$isoform)
table(iso_table$sample)
iso_table <- iso_table[,c(1,3,5,4)]
colnames(iso_table) <- c("Reads","trans","sample","symbol")
work = "./IGV/"
dir = paste0(work, gene_feat)
if (file.exists(dir) == F) {
  icesTAF::mkdir(dir)
}
sam_name <- unique(iso_table$sample)
rtracklayer::export(
  subset(gtf_pb, gene_id %in% gene_feat),
  paste0(work, gene_feat, "/", gene_feat, "_anno.gtf")
)
lapply(sam_name, function(sam) {
          write.table(
            data.table(iso_table)[sample == sam][, "Reads"],
            paste0(work, "/", gene_feat,"/", gene_feat, "@", sam, "@readID.txt"),
            quote = F,
            sep = '\t',
            col.names = F,
            row.names = F
          )
  })

sashimi_gene <- as.data.table(gtf_pb)[ID %in% gene_feat]
sashimi_gene <- tibble::column_to_rownames(sashimi_gene, var = "ID")
cor <- paste0(rownames(sashimi_gene),"_",sashimi_gene$seqnames,":",sashimi_gene$start,"-",sashimi_gene$end,":",sashimi_gene$strand)
write.table(cor,paste0("./IGV/",gene_feat,"/_gene_cor.txt"),col.names = F,quote = F,row.names = F)
```

```{linux}
conda active
IGVwork="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/rep4_pb_analysis/IGV/"
bamPath="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/batch4_6pb"
picard="/mnt/raid61/Personal_data/songjunwei/software/picard.jar"

#Split sample bam
for gene in $(ls $IGVwork); do
cd $IGVwork/$gene
for samp in $(ls *readID.txt |awk -F '@' '{print $2}') ; do

java -jar $picard FilterSamReads \
I=$bamPath/${samp}.bam    O=${samp}_s.bam \
READ_LIST_FILE=${gene}@${samp}@readID.txt \
FILTER=includeReadList

bedtools bamtobed -bed12 -split -i ${samp}_s.bam |bedToGenePred /dev/stdin /dev/stdout| \
genePredToGtf "file" stdin ${samp}_plot.gtf
done
done
```

```{r fig.width=6,fig.height=4}
work = "./IGV/"
gene = "PBANKA_1245821" #PbANKA_12_v3:1735034-1737075:+
file <- rev(list.files(paste0(work, gene), pattern = "plot.gtf", full.names = T))
tx_list <- GRangesList(tx_list_fun(file))
names(tx_list) <- as.character(seq(1:length(names(tx_list))))

tx_tab <- data.table(data.frame(unlist(tx_list)))
tx_tab[Bam %in% paste0(c("RTA-03","RTA-10","RTA-16"),"_plot.gtf")]$method="#D05E15" #D05E15,红，tro
tx_tab[Bam %in% paste0(c("RTA-17","RTA-24","RTA-32"),"_plot.gtf")]$method="#239873" #239873,绿，sch

igv_color <- unique(tx_tab[,c("Bam","method")])$method
igv_lab <- unique(tx_tab[,c("Bam","method")])$Bam

coord <- read.table(paste0(work, gene, "/_gene_cor.txt"))$V1
p_raw <-
  ggbio::autoplot(
    tx_list,
    chr = "PbANKA_12_v3",
    gap.geom = "chevron",
    label.color = "black",
    aes(fill = Bam),
    show.coverage = T,
    indName = "grl_name",
    group.selfish = T
  )

p0 <- p_raw@ggplot +  theme_clear() +
  theme_classic2() +  xlim(as.numeric(1735034),as.numeric(1737075)) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) + scale_fill_manual(values = igv_color)


p1 <- ggbio:: autoplot(gtf_pb[gtf_pb$ID %in% gene])
p2 <- p1@ggplot +
  theme_clear() +
  scale_fill_manual(values = c("black")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) +
  labs(y = "GTF",x=gene)

p_all <- (p0 / p2) + plot_layout(ncol = 1, heights  = c(8, 1))
ggsave(
  p_all,
  filename = paste0("./IGV/", gene, "_igv_gtf.pdf"),
  width = 10,
  height = 5.5
)
```

```{r}
m6a_mat <- readRDS("./01m6a_Ratio_mat.rds")
m6a_mat["PBANKA_0406650",]
assay(gene_vst)["PBANKA_0406650",]
```



#07 m6a analysis
##venn
```{r fig.width=6,fig.height=5}
id_list <- list(Trophozoite=m6A_Tab[type=="Trophozoite"]$name,
                Schizont=m6A_Tab[type=="Schizont"]$name)
p1 <- ggvenn::ggvenn(id_list,fill_color = c("#D05E15","#239873"))+labs(title = "All m6A")

id_list2 <- list(Trophozoite=m6A_Tab[barcode=="RTA03"]$name,
                Schizont=m6A_Tab[barcode=="RTA17"]$name)
p2 <- ggvenn::ggvenn(id_list2,fill_color =c("#D05E15","#239873"))+labs(title = "One sample m6A")

pdf("./07m6a_venn.pdf",6,5)
p1+p2
dev.off()
```

##metagene
```{r fig.width=5,fig.height=4}
#gtf
gtf_pb <- rtracklayer::import(paste0(refpath,"PbergheiANKA/PlasmoDB-53_PbergheiANKA.gff"))
naidx <- is.na(gtf_pb$gene_id)
gtf_pb$gene_id[naidx] <- gtf_pb$ID[naidx]
gtf_pb_clean <- unique(gtf_pb[,c("type","gene_id","ID")])
View(data.frame(gtf_pb))
## over geneome 
gtf_pb_gene <- gtf_pb_clean[gtf_pb_clean$type %in% c("protein_coding_gene","pseudogene","ncRNA_gene"),]
table(m6A_Tab$gene %in% gtf_pb_gene$gene_id)
gtf_pb_gene <- data.table(data.frame(gtf_pb_gene))

m6A_gene_spl <- split(m6A_Tab,f = m6A_Tab$gene)
m6A_gene_rb <- rbindlist(lapply(m6A_gene_spl, function(g){
  if (length(unique(g$strand))!=1) {
    stop("strand?")
  }
  gend_tmp <- gtf_pb_gene[gene_id==unique(g$gene)]
  if (unique(g$strand)=="+") {
   g$gtftype <-  ifelse(g$pos>gend_tmp$end,"over_3end",
                        ifelse(g$pos<gend_tmp$start,"over_5end",
                               ifelse(g$pos>=gend_tmp$start & g$pos<=gend_tmp$end,"inside","ques")))
  }else{
     g$gtftype <-  ifelse(g$pos<gend_tmp$start,"over_3end",
           ifelse(g$pos>gend_tmp$end,"over_5end",
                               ifelse(g$pos>=gend_tmp$start & g$pos<=gend_tmp$end,"inside","ques")))
  }
  return(g)
}))
table(m6A_gene_rb$gtftype)

#m6A
m6A_pos <- unique(m6A_gene_rb[,c("chr","pos","strand","type","gtftype","gene")])
m6A_pos$pos2 <- m6A_pos$pos
m6A_gr <- GRanges(seqnames = m6A_pos$chr,ranges = IRanges(start = m6A_pos$pos,end = m6A_pos$pos2),strand = m6A_pos$strand)
mcols(m6A_gr) <- m6A_pos[,c("type","gtftype","gene")]

## cds 3UTR
gtf_pb_nei <- gtf_pb_clean[gtf_pb_clean$type %in% c("three_prime_UTR","CDS","exon"),]
table(gtf_pb_nei$type)
hts <- findOverlaps(m6A_gr,gtf_pb_nei)

find <- data.frame(m6A_gr[queryHits(hts)],gtf_pb_nei[subjectHits(hts)])
find$id <- paste0(find$gene,"_",gsub("PbANKA_","",find$seqnames),":",find$start,":",find$strand)
find <- find[,c("id","type.1")] %>% unique() 
table(find$gtftype)
table(find$type.1)
table(duplicated(find$id))
m6A_gene_new <- merge(m6A_gene_rb,find,by="id",all=T)
table(m6A_gene_new$type.1)
m6A_gene_new$type.1 <- as.character(m6A_gene_new$type.1)
naidx <- is.na(m6A_gene_new$type.1)
m6A_gene_new$type.1[naidx] <- m6A_gene_new$gtftype[naidx]
table(m6A_gene_new$type.1)
table(m6A_gene_new$type.1,m6A_gene_new$type)

m6A_gene_new <- m6A_gene_new[type.1!="three_prime_UTR"]
m6A_gene_new_tro <- m6A_gene_new[type == "tro"]
data <- m6A_gene_new_tro[,.N,by="type.1"]
data[5,] <- c("three_prime_UTR",0)
data$precent = round(data$N/sum(data$N),2)*100
pie1 <- ggplot(data, aes(x="", y=precent, fill=type.1)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  geom_text(aes(label = paste0(precent, "%")), position = position_stack(vjust=0.5)) +
  labs(x = NULL, y = NULL, fill = NULL)+ theme_bw() + theme(legend.text = element_text(size = 14))+
  ggtitle(label = "Trophozoite")
m6A_gene_new_sch <- m6A_gene_new[type == "sch"]
data_sch <- m6A_gene_new_sch[,.N,by="type.1"]
data_sch$precent = round(data_sch$N/sum(data_sch$N),2)*100
pie2 <- ggplot(data_sch, aes(x="", y=precent, fill=type.1)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  geom_text(aes(label = paste0(precent, "%")), position = position_stack(vjust=0.5)) +
  labs(x = NULL, y = NULL, fill = NULL)+ theme_bw() + theme(legend.text = element_text(size = 14))+
  ggtitle(label = "Schizont ")

ggpubr::ggarrange(pie1,pie2,common.legend = T)
```

#08 model trainning
```{r}
pb_align <-
  readGAlignments("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/02align/batch3_9Pb/PbergheiANKA.bam",
    param = ScanBamParam(
      flag = scanBamFlag(
        isSupplementaryAlignment = FALSE,
        isPaired = FALSE,
        isSecondaryAlignment = FALSE
      ),
      mapqFilter = 60
    ),
    use.names = T
  )

write.table(names(pb_align)[1:12000],"/mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/01rawdata/Train_pb_readid.txt",row.names = F,col.names = F,quote = F)
write.table(names(pb_align)[12001:15000],"/mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/01rawdata/Valid_pb_readid.txt",row.names = F,col.names = F,quote = F)
```


```{linux}
conda activate base #2101
f5="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/01Rawdata/batch3_fq/DRS_multiplex3/Pb_9sample/20210918_2113_MN27836_FAQ60861_83a343d5"
pb_bam="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/02align/batch3_9Pb/PbergheiANKA.bam"
pb_ref="/mnt/raid61/Personal_data/songjunwei/reference/plasmoDB/PlasmoDB-53_PbergheiANKA_Genome.fasta"
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model
for i in Train  Valid; 
do
fast5_subset -i $f5 \
            -s 01rawdata/${i}_pb_fast5 \
            --batch_size  2000 \
            --read_id_list 01rawdata/${i}_pb_readid.txt \
            --recursive  --ignore_symlinks \
            -t 30
seqkit grep -f 01rawdata/${i}_pb_readid.txt /mnt/raid61/Personal_data/songjunwei/DRS_RTA/01Rawdata/batch3_fq/DRS_multi3_9Pb.fastq > 01rawdata/${i}_pb_guppy.fq
minimap2 -ax map-ont -t 30  $pb_ref 01rawdata/${i}_pb_guppy.fq |  samtools sort -@ 30 | samtools view -b > 01rawdata/${i}_pb_guppy.bam
samtools index 01rawdata/${i}_pb_guppy.bam
done

#### tayaki ####
conda activate base #2048
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model
tayaki="/mnt/raid5/Personal/minion/software/taiyaki/bin/"
pb_ref="/mnt/raid61/Personal_data/songjunwei/reference/plasmoDB/PlasmoDB-53_PbergheiANKA_Genome.fasta"

for i in Train Valid ;
do
python $tayaki/generate_per_read_params.py 01rawdata/${i}_pb_fast5 --output 02tayaki/${i}_pb_readid.tsv
python $tayaki/get_refs_from_sam.py $pb_ref  01rawdata/${i}_pb_guppy.bam  --output 02tayaki/${i}_pb_tayaki.fa --reverse
python $tayaki/prepare_mapped_reads.py 01rawdata/${i}_pb_fast5 02tayaki/${i}_pb_readid.tsv /mnt/raid5/Personal/minion/DRS_mul/pb_train_mod/${i}_tayaki.h5 \
        /mnt/raid5/Personal/minion/software/taiyaki/models/r941_rna_minion_upgrae.checkpoint     02tayaki/${i}_pb_tayaki.fa --jobs 30 --overwrite 
cp /mnt/raid5/Personal/minion/DRS_mul/pb_train_mod/${i}_tayaki.h5 02tayaki/${i}_tayaki.h5
done

cd /mnt/raid5/Personal/minion/DRS_mul/pb_train_mod/
rodan="/mnt/raid61/Personal_data/songjunwei/software/RODAN/"
for i in Train Valid ;
do
mkdir rodan_${i}
python  $rodan/gendata.py -i ${i}_tayaki.h5 --outdir rodan_${i}
mv rodan_${i}/train.hdf5 rna-${i}.hdf5
rm -rf rodan_${i}
done
mkdir runs
#vim rna.config
python  $rodan/model.py -c ./rna.config -a $rodan/rnaarch -n Pb --arch $rodan/rnaarch --savedir ./runs --rna -w 30 -l
```


```{linux}
#### basecalling ####

#rodan
conda activate base #2101
multi_to_single_fast5 --input_path /mnt/raid61/Personal_data/songjunwei/DRS_RTA/01Rawdata/batch3_fq/DRS_multiplex3/Pb_9sample  --save_path /mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/04f5sg --recursive -t 5

f5="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/04f5sg"
rodan="/mnt/raid61/Personal_data/songjunwei/software/RODAN/"
pb_ref="/mnt/raid61/Personal_data/songjunwei/reference/plasmoDB/PlasmoDB-53_PbergheiANKA_Genome.fasta"
cd /mnt/raid5/Personal/minion/DRS_mul/pb_train_mod/03accuracy/
for model in /mnt/raid5/Personal/minion/DRS_mul/pb_train_mod/runs/Pb-epoch*.torch; do
    epoch=$(basename "$model" | sed 's/Pb-epoch\([0-9]*\).torch/\1/')
    output_fa="fafiles/Pb_rodon_ep${epoch}.fa"
    output_sam="samfiles/Pb_rodon_ep${epoch}.sam"
    accuracy_result_file="accufiles/accuracy_epoch${epoch}.txt"
    python $rodan/basecall_mod.py -m $model $f5 --arch $rodan/rnaarch -b 250 > $output_fa
    minimap2 --secondary=no -ax map-ont -t 32 --cs $pb_ref $output_fa > $output_sam
    python $rodan/accuracy.py $output_sam $pb_ref > $accuracy_result_file
done
cp -r /mnt/raid5/Personal/minion/DRS_mul/pb_train_mod/03accuracy/ /mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/03accuracy

## rodan raw
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/03accuracy
rodan="/mnt/raid61/Personal_data/songjunwei/software/RODAN/"
pb_ref="/mnt/raid61/Personal_data/songjunwei/reference/plasmoDB/PlasmoDB-53_PbergheiANKA_Genome.fasta"
python $rodan/basecall_mod.py -m $rodan/rna.torch $f5 --arch $rodan/rnaarch -b 200 > /mnt/raid5/Personal/minion/DRS_mul/pb_train_mod/03accuracy/fafiles/Pb_rodan_raw.fa
minimap2 --secondary=no -ax map-ont -t 32 --cs $pb_ref fafiles/Pb_rodan_raw.fa > samfiles/Pb_rodan_raw.sam
python $rodan/accuracy.py samfiles/Pb_rodan_raw.sam $pb_ref > accufiles/accuracy_rodan.txt
# Total: 20892 Median accuracy: 0.891363662295334 Average accuracy: 0.8628132824601701 std: 0.08594840858519276
# Median  - Mismatch: 0.022126047577183276 Deletions: 0.056451612903225805 Insertions: 0.02063886678400006
# Average - Mismatch: 0.02532956776709016 Deletions: 0.08870056421248074 Insertions: 0.02315658556025895

#guppy 2048
guppy_basecaller \
-i /mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/04f5sg \
--save_path /mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/01rawdata/guppy_fq/ --device cuda:0 \
-c rna_r9.4.1_70bps_hac_prom.cfg \
--min_qscore 7  --verbose_logs  \
--chunks_per_runner 256  --chunk_size 100 \
--num_callers 4 --gpu_runners_per_device 8 -r
cat /mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/01rawdata/guppy_fq/pass/* > /mnt/raid5/Personal/minion/DRS_mul/pb_train_mod/03accuracy/fafiles/Pb_guppy.fastq

cd /mnt/raid5/Personal/minion/DRS_mul/pb_train_mod/03accuracy/
pb_ref="/mnt/raid61/Personal_data/songjunwei/reference/plasmoDB/PlasmoDB-53_PbergheiANKA_Genome.fasta"
minimap2 --secondary=no -ax map-ont -t 30 --cs $pb_ref fafiles/Pb_guppy.fastq > samfiles/Pb_guppy.sam
python $rodan/accuracy.py samfiles/Pb_guppy.sam $pb_ref > accufiles/accuracy_guppy.txt

#guppy 6.5.1
# Total: 18290 Median accuracy: 0.8504049180327868 Average accuracy: 0.8284826892161672 std: 0.08036716002367697
# Median  - Mismatch: 0.038365042147282064 Deletions: 0.06972702661165428 Insertions: 0.027944111776447105
# Average - Mismatch: 0.04176736477915147 Deletions: 0.09872129812093006 Insertions: 0.031028647883751332

# Total: 18581 Median accuracy: 0.8994974874371859 Average accuracy: 0.8680436405131445 std: 0.09119572458812626
# Median  - Mismatch: 0.018666666666666668 Deletions: 0.05086848635235732 Insertions: 0.022290545734050732
# Average - Mismatch: 0.020604501145293615 Deletions: 0.08782336398897211 Insertions: 0.023528494352589798

```


```{r}
files <- list.files(path = "/mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/03accuracy/accufiles/",pattern = "accuracy_.*\\.txt$",full.names = T)

results_df <- data.frame()

# 遍历每个文件，读取内容，并提取所需信息
for (file in files) {
  epoch <- gsub("accuracy_epoch([0-9]+)\\.txt", "\\1", basename(file))  # 提取 epoch 编号
  content <- readLines(file)
  
  # 提取 Total, Median accuracy, Average accuracy, std
  total <- as.numeric(gsub("Total: ([0-9]+).*", "\\1", content[1]))
  median_accuracy <- as.numeric(gsub(".*Median accuracy: ([0-9.]+).*", "\\1", content[1]))
  average_accuracy <- as.numeric(gsub(".*Average accuracy: ([0-9.]+).*", "\\1", content[1]))
  std <- as.numeric(gsub(".*std: ([0-9.]+)", "\\1", content[1]))
  
  # 提取 Mismatch, Deletions, Insertions 的 Median 和 Average 值
  mismatch_median <- as.numeric(gsub(".*Mismatch: ([0-9.]+).*", "\\1", content[2]))
  deletions_median <- as.numeric(gsub(".*Deletions: ([0-9.]+).*", "\\1", content[2]))
  insertions_median <- as.numeric(gsub(".*Insertions: ([0-9.]+)", "\\1", content[2]))
  
  mismatch_average <- as.numeric(gsub(".*Mismatch: ([0-9.]+).*", "\\1", content[3]))
  deletions_average <- as.numeric(gsub(".*Deletions: ([0-9.]+).*", "\\1", content[3]))
  insertions_average <- as.numeric(gsub(".*Insertions: ([0-9.]+)", "\\1", content[3]))
  
  # 将提取的数据添加到数据帧
  results_df <-
    rbind(
      results_df,
      data.frame(
        Epoch = epoch,
        Total = total,
        MedianAccuracy = median_accuracy,
        AverageAccuracy = average_accuracy,
        Std = std,
        MismatchMedian = mismatch_median,
        DeletionsMedian = deletions_median,
        InsertionsMedian = insertions_median,
        MismatchAverage = mismatch_average,
        DeletionsAverage = deletions_average,
        InsertionsAverage = insertions_average
      )
    )
}

rownames(results_df) <- results_df$Epoch
results_df$Epoch <- as.numeric(as.character(results_df$Epoch))
results_df <- results_df[order(results_df$Epoch),]
write.csv(results_df,"/mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/03accuracy/Epochs_res.csv",row.names = T)
```

```{r}
results_com <- results_df[c("34","accuracy_guppy.txt","accuracy_rodan.txt"),]
results_com <- data.table(results_com[,-1])
results_com$type <- c("Rodan_v1.0_Train","Guppy_v6.5.7","Rodan_v1.0_Raw")

pdata_com <- melt(results_com[,-c("Total","Std")])
pdata_com$type <- factor(pdata_com$type,levels = c("Rodan_v1.0_Train","Rodan_v1.0_Raw","Guppy_v6.5.7"))
pdata_com$value <- round(pdata_com$value*100,2)
pdata_com_aver <- pdata_com[grep("Average",pdata_com$variable),]
pdata_com_aver$variable <- gsub("Average","",pdata_com_aver$variable)
pdata_com_medi <- pdata_com[grep("Median",pdata_com$variable),]
pdata_com_medi$variable <- gsub("Median","",pdata_com_medi$variable)

p1 <- ggplot(pdata_com_medi, aes(x = variable, y = value, fill = type)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_text(
    aes(label = value),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    size = 3.5
  ) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)) +  
  xlab(NULL) + ylab("Percentage%") + ggtitle("Comparison of median")+
  mytheme
p2 <- ggplot(pdata_com_aver, aes(x = variable, y = value, fill = type)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_text(
    aes(label = value),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    size = 3.5
  ) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)) +  
  xlab(NULL) + ylab("Percentage%") + ggtitle("Comparison of average")+
  mytheme

ggpubr::ggarrange(p1,p2,nrow  = 1,common.legend = T)
ggsave("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/03accuracy/Comparison_res.pdf",width = 9,height = 5)
```


```{r}
results_df <- na.omit(results_df)
# 1. Total 随 Epoch 变化的图
p1 <- ggplot(data = results_df, aes(x = Epoch, y = Total)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Total Over Epochs", x = "Epoch", y = "Total")

# 2. Median 和 Average Accuracy 随 Epoch 变化的图
accuracy_df <- melt(results_df[, c("Epoch", "MedianAccuracy", "AverageAccuracy")], id.vars = "Epoch")
p2 <- ggplot(data = accuracy_df, aes(x = Epoch, y = value, color = variable)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Accuracy Over Epochs", x = "Epoch", y = "Accuracy") +
  scale_color_brewer(palette = "Set1")

# 3. Median 和 Average Mismatch 随 Epoch 变化的图
mismatch_df <- melt(results_df[, c("Epoch", "MismatchMedian", "MismatchAverage")], id.vars = "Epoch")
p3 <- ggplot(data = mismatch_df, aes(x = Epoch, y = value, color = variable)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Mismatch Over Epochs", x = "Epoch", y = "Mismatch") +
  scale_color_brewer(palette = "Set1")

# 4. Median 和 Average Deletions 随 Epoch 变化的图
deletions_df <- melt(results_df[, c("Epoch", "DeletionsMedian", "DeletionsAverage")], id.vars = "Epoch")
p4 <- ggplot(data = deletions_df, aes(x = Epoch, y = value, color = variable)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Deletions Over Epochs", x = "Epoch", y = "Deletions") +
  scale_color_brewer(palette = "Set1")

# 5. Median 和 Average Insertions 随 Epoch 变化的图
insertions_df <- melt(results_df[, c("Epoch", "InsertionsMedian", "InsertionsAverage")], id.vars = "Epoch")
p5 <- ggplot(data = insertions_df, aes(x = Epoch, y = value, color = variable)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Insertions Over Epochs", x = "Epoch", y = "Insertions") +
  scale_color_brewer(palette = "Set1")

pall <- ggarrange(p2,p3,p4,p5)
ggsave("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/03accuracy/Epochs_total.pdf",p1,width = 4,height = 4)
ggsave("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/03accuracy/Epochs_accuracy_res.pdf",pall,width = 8,height = 6)
```
