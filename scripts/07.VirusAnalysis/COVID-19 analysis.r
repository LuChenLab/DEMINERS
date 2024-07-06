---
title: "COVID-19 analysis"
author: "junwei"
date: "2023-08-01"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/")
library(data.table)
library(Rsamtools)
library(ShortRead)
library(GenomicRanges)
library(GenomicAlignments)
library(ggplot2)
library(ggvenn)
library(vcfR)
library(dplyr)
library(ggseqlogo)
setwd("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/")
source("/mnt/raid61/Personal_data/songjunwei/AllScript/00Function/DIyplot.R")
source("/mnt/raid61/Personal_data/songjunwei/AllScript/00Function/Ontology_hostgene.R")
source("/mnt/raid61/Personal_data/songjunwei/AllScript/00Function/m6a_process.R")
source("/mnt/raid61/Personal_data/songjunwei/AllScript/00Function/Isoform_psi_func.R")
```

```{r}
sars_gtf <-  rtracklayer::import('/mnt/raid61/Personal_data/songjunwei/reference/SARS_CoV2/Sars_cov_2.ASM985889v3.101.gtf')
sampleinfo <- read.csv("24sample_info.csv")
```

#00 Barcode
```{r }
#230808_24sam 
pred <- lapply(list.files("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/batch5_SARS/230808_24sam", pattern = "fast5",full.names = T), fread) %>% rbindlist()
pred$Read <- gsub("read_","",pred$Read)

hist(pred$Probability, breaks = 100)
pred_tab <- pred[Probability > 0.3, .N, Barcode]
pred_tab$pro <- round(pred_tab$N/sum(pred_tab$N),4)
pred_tab$Barcode <- gsub("-","_",pred_tab$Barcode)
View(pred_tab)
```

#01 fastq/align
```{r}
fastq <- readFastq(c("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/01Rawdata/batch5_SARS/230808_24sam_all_dna.fq.gz","/mnt/raid61/Personal_data/songjunwei/DRS_RTA/01Rawdata/batch5_SARS/230808_24sam2_all_dna.fq.gz"))
fastq_qc <- fread("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/QC/24sam_NanoComp-data.tsv.gz",header = T)

sars_24sam_fq <- data.table(data.frame(readid=ShortRead::id(fastq),seq=sread(fastq),width=width(fastq),fastq_qc[dataset %in% c("1st_all","2nd_all")]))
sars_24sam_fq$readid <- stringr::str_split_fixed(sars_24sam_fq$readid," ",3)[,1]
table(sars_24sam_fq$width==sars_24sam_fq$lengths)

#merge SARS
param <- ScanBamParam(flag = scanBamFlag(isSupplementaryAlignment = FALSE,isPaired = FALSE,isSecondaryAlignment = FALSE),mapqFilter = 0,what = "seq")
map1_sars <- readGAlignments("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/02align/batch5/230808_24sam_SARS.bam", param = param,use.names = T)
map2_sars <- readGAlignments("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/02align/batch5/230808_24sam2_SARS.bam", param = param,use.names = T)
map_sars <- c(map1_sars,map2_sars)
sars_24sam_fq$align <- ifelse(sars_24sam_fq$readid %in% names(map_sars),"SARS",NA)
table(sars_24sam_fq$align)

#merge pred
sars_24sam_fq_RTA <- merge(sars_24sam_fq,pred[,c("Read","Barcode","Probability")],by.x="readid",by.y="Read",all=T) 
pred <- pred[,c("ReadID","RTA","Probability")]
colnames(pred) <- c("ReadID","Barcode_lla","Probability")
sars_24sam_fq_RTA1 <- merge(sars_24sam_fq_RTA,pred,by.x="readid",by.y="ReadID",all.x=T) 
sars_24sam_fq_RTA1$width <- NULL
sars_24sam_fq_RTA1$length_filter <- NULL
saveRDS(sars_24sam_fq_RTA1,"01filter_rawdata/Combin_QC_align_raw.rds")

#filter
hist(sars_24sam_fq_RTA1$quals)
sars_24sam_fq_RTA1[quals>5,.N]/nrow(sars_24sam_fq_RTA1)
sars_24sam_fq_RTA1[align=="SARS"][quals>5,.N]
sars_24sam_fq_RTA1[lengths>100,.N]/nrow(sars_24sam_fq_RTA1)
sars_24sam_fq_RTA1[align=="SARS"][lengths>100,.N]

sars_24sam_sub <- sars_24sam_fq_RTA1[lengths>100 & quals>5]
write.table(sars_24sam_sub$readid,"01filter_rawdata/readid_filter.txt",quote = F,col.names = F,row.names = F)
saveRDS(sars_24sam_sub,"01filter_rawdata/Combin_QC_align_filter.rds")
```



```{r}
###### 过滤保存
sars_24sam_sub2 <- sars_24sam_sub[Probability > 0.3]
table(sars_24sam_sub2$Barcode_lla)
lapply(unique(sars_24sam_sub2$Barcode_lla), function(x){
  write.table(sars_24sam_sub2[Barcode_lla==x]$readid,paste0("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata/00split_barcode/",x,"_readid_filter.txt"),quote = F,col.names = F,row.names = F)
})

sars_24sam_sub3 <- sars_24sam_sub[Probability > 0.3][align=="SARS"]
table(sars_24sam_sub3$Barcode_lla)
lapply(unique(sars_24sam_sub3$Barcode_lla), function(x){
  write.table(sars_24sam_sub3[Barcode_lla==x]$readid,paste0("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/05m6A/01align_SARS_f5/",x,"_SARS_readid.txt"),quote = F,col.names = F,row.names = F)
})
```

## ⭐️STAR
```{linux}
## mk fastq
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata
cat /mnt/raid61/Personal_data/songjunwei/DRS_RTA/01Rawdata/batch5_SARS/230808_24sam*_all.fq.gz > 230808_24sam_comb.fq.gz
seqkit grep -f readid_filter.txt 30808_24sam_comb.fq.gz > 230808_24sam_comb_filter.fq
seqkit seq 230808_24sam_comb_filter.fq --rna2dna > 230808_24sam_comb_filter_dna.fq
fast5_subset -i /mnt/raid61/Personal_data/songjunwei/DRS_RTA/01Rawdata/batch5_SARS/230808_24sam \
            -s 230808_24sam_comb_filter_fast5 \
            --batch_size  5000 \
            --read_id_list readid_filter.txt \
            --recursive  --ignore_symlinks \
            -t 10

## align
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata
minimap2 -ax map-ont -t 30  /mnt/raid62/BetaCoV/Person/zhangyiming/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta \
       230808_24sam_comb_filter.fq |  samtools sort -@ 30 | samtools view -b > 230808_24sam_comb_filter_sars.bam
minimap2 -ax splice -uf -k14 -t 30  /mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/HomSap.fasta \
       230808_24sam_comb_filter.fq |  samtools sort -@ 30 | samtools view -b > 230808_24sam_comb_filter_homsap.bam
minimap2 -ax splice -uf -k14 -t 30 /mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/S_cere/GCF_000146045.2_R64_genomic.fa \
       230808_24sam_comb_filter.fq |  samtools sort -@ 30 | samtools view -b > 230808_24sam_comb_filter_Sere.bam

## mut
bcftools mpileup -f /mnt/raid62/BetaCoV/Person/zhangyiming/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta \
   $output/01filter_rawdata/230808_24sam_comb_filter_sars.bam  -Q 0 | bcftools call -mv -o $output/01filter_rawdata/230808_24sam_SARS_mut.vcf

## coverage
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata/02align_filter/
samtools depth 230808_24sam_comb_filter_sars.bam > 24sam_ft_sars.depth
bedtools genomecov -bga -pc -ibam 230808_24sam_comb_filter_sars.bam > 24sam_ft_sars.frag.cov

## split
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata/00split_barcode
for i in $( ls *_readid_filter.txt | awk -F '_read' '{print $1}');
do 
seqkit grep -f ${i}_readid_filter.txt /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata/230808_24sam_comb.fq.gz > ${i}_ft.fq
done

##QC
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata/00split_barcode
NanoComp --fastq  $(ls *_ft.fq) -o /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/06microorganism/00QC -t 30  -p Comb_filt_fq_ --raw

cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata/
NanoComp --fastq 230808_24sam_comb.fq.gz 230808_24sam_comb_filter.fq -n comb_raw comb_filt -o QC -t 30 -p comb.fq_ --raw 

cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/03NGS_res
NanoComp --bam  $(ls 20*bam) -o /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/00result -t 30  -p NGS_bam_ --raw
```

##QC/Gene/Coverage
QC
```{r }
Nano_res <- fread("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/06microorganism/00QC/Comb_filt_fq_NanoComp-data.tsv.gz",header = T) 
#reads number
Nano_res$dataset <- gsub("_ft.fq","",Nano_res$dataset)
table(Nano_res$length_filter)
Reads_number <- Nano_res[,.N,by=dataset]
Reads_number$dataset <- factor(Reads_number$dataset,levels =Reads_number[order(N,decreasing = T)]$dataset )
pnum <- ggplot(Reads_number,aes(x=dataset,y=N,fill=dataset))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")+
  geom_bar(stat = 'identity',position = 'stack')+ mytheme+xlab("Samples")+
  ylab("Reads number")
ggsave("00result/QC_reads_number.pdf",plot = pnum,width = 6,height = 4)

#quality
pqual <- ggplot(data=Nano_res,mapping = aes(x=dataset,y=quals,fill=dataset))+
       geom_violin(trim = FALSE,color="white") +
       geom_boxplot(outlier.colour = NA,position = position_dodge(0.9))+
      xlab("Samples")+ ylab("Reads Quality")+ mytheme +
      theme(
            panel.grid.minor = element_blank(),
            text = element_text(size = 18),
            axis.text.x= element_text(size=12,family="ArialMT",angle =  45,hjust = 1),
            legend.position = "none")
ggsave("00result/QC_quality.pdf",plot = pqual,width = 6,height = 4)

#length
plen <- ggplot(data=Nano_res,mapping = aes(x=lengths)) + 
  geom_density(aes(fill = dataset),alpha=0.3)+ scale_x_log10()+
  geom_vline(data = Nano_res, aes(xintercept = mean(lengths), color=dataset),linetype="dashed")+
  mytheme+ 
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank(),
        text = element_text(size = 12),
        axis.text.x= element_text(size=12,family="ArialMT"))+
    xlab("Length")+ ylab("Denisty")
ggsave("00result/QC_length.pdf",plot = plen,width = 8,height = 4)
```

gene
```{r fig.width=4.5,fig.height=3}
sars_align <-
  readGAlignments(
    "/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata/02align_filter/230808_24sam_comb_filter_sars.bam",
    param = ScanBamParam(
      flag = scanBamFlag(
        isSupplementaryAlignment = FALSE,
        isPaired = FALSE,
        isSecondaryAlignment = FALSE
      ),
      mapqFilter = 0
    ),
    use.names = T
  )
sars_gtf <-  rtracklayer::import('/mnt/raid61/Personal_data/songjunwei/reference/SARS_CoV2/Sars_cov_2.ASM985889v3.101.gtf')
sars_gtf_gene <- sars_gtf[sars_gtf$type=="gene",]
sars_gtf_gene <- sars_gtf_gene[-2]
overlaps <- findOverlaps(sars_align,sars_gtf_gene)
genes <- sars_gtf_gene[subjectHits(overlaps)]$gene_name

gene_df <- data.frame(table(genes))
gene_df <- gene_df[order(gene_df$Freq,decreasing = T),]
gene_df$genes <- factor(gene_df$genes,levels = gene_df$genes)
pgene <- ggplot(gene_df,aes(genes,Freq,fill=genes))+geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")+ xlab(NULL)+ylab("Gene frequency")+
  mytheme
ggsave("00result/SARS-CoV2-gene_freq.pdf",pgene,width = 4.5,height = 3)
```

coverage
```{r }
gff <- ape::read.gff("/mnt/raid62/BetaCoV/Person/tangchao/data/reference/MN908947.genome.gff")
gff <- gff[gff$type == "gene", ]
gff$Gene <- sapply(gff$attributes, function(x) {
  x = str_split(x, ";")[[1]]
  for (i in x) {
    if (str_detect(i, "Name")) {
      i = str_split(i, "=")[[1]]
      return(i[2])
    }
  }
})
gff$Gene = as.character(gff$Gene)
gff$molecule = "Ref"
gff$direction = sapply(gff$strand, function(x) { ifelse(x == "+", 1, -1) })
gff$strand = sapply(gff$strand, function(x) { ifelse(x == "+", "forward", "reverse") })
gff$Gene <- factor(gff$Gene, levels = gff$Gene)

library(gggenes)
p0 <- ggplot(gff, aes(xmin = start, xmax = end, y = molecule, fill = Gene, label = Gene)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), 
                  arrowhead_width = unit(1, "mm")) +
  geom_gene_label(align = "centre") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 11, name = "Paired")) +
  theme_genes() +
  theme(
    legend.position = "top", 
    legend.text = element_text(size = 8), 
    legend.title = element_text(size = 8), 
    axis.line.x = element_line(), 
    axis.title.y = element_blank(),
    axis.text = element_text(size = 8), 
    axis.title = element_text(size = 8)
  ) +  scale_x_continuous(breaks = c(seq(0, 30000, 5000)), position = "bottom") 
  guides(fill = guide_legend(nrow = 1))

cov_depth <- read.table("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata/02align_filter/24sam_ft_sars.depth")
cov_depth$V4 <-  log10(cov_depth$V3)
cov_depth$V4 <- ifelse(cov_depth$V4 == -Inf,NA,cov_depth$V4)

p1 <- ggplot(cov_depth, aes(x = V2, y = V4)) + 
  geom_step(aes(colour = V1)) +   theme_classic() +
  geom_hline(yintercept = 0, colour = "white", lwd = 1) + 
  theme(legend.position = "none")+ ylab("Log10(Coverage)")+
  mytheme
  # annotate("rect", xmin = 13342, xmax = 13460, ymin = 0, ymax = Inf, alpha = .2) + 
  # annotate("rect", xmin = 28287, xmax = 28358, ymin = 0, ymax = Inf, alpha = .2) + 
  # annotate("rect", xmin = 29164, xmax = 29230, ymin = 0, ymax = Inf, alpha = .2) + 
  # annotate("rect", xmin = 28881, xmax = 28979, ymin = 0, ymax = Inf, alpha = .2) 
pcov <- p0 / p1 + patchwork::plot_layout(heights = c(1, 8)) 

ggsave("00result/SARS-CoV2-seqcov_log.pdf",pcov,width = 8,height = 4)
```

Reads分布
```{r}
param <- ScanBamParam(flag = scanBamFlag(isSupplementaryAlignment = FALSE,isPaired = FALSE,isSecondaryAlignment = FALSE),mapqFilter = 0,what = "seq")
bampath <- list.files("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata",pattern = "bam$",full.names = T)
Map_list <- lapply(bampath, function(x) readGAlignments(x, param = param,use.names = T) )
names(Map_list) <- gsub("230808_24sam_comb_filter_","",basename(bampath))
sars_24sam_sub$align <- NA
sars_24sam_sub$align[sars_24sam_sub$readid %in% names(Map_list$sars.bam)] <- "SARS"
sars_24sam_sub$align[sars_24sam_sub$readid %in% names(Map_list$Sere.bam)] <- "Sere"
sars_24sam_sub$align[sars_24sam_sub$readid %in% names(Map_list$homsap.bam)] <- "homsap"
sars_24sam_sub$align[is.na(sars_24sam_sub$align)] <- "unknown"
table(sars_24sam_sub$align)

p_rta_align_lla <- sars_24sam_sub[Probability > 0.3, .N, c("Barcode_lla","align")]
p1 <- ggplot(p_rta_align_lla, aes(Barcode_lla, N, fill = align)) + geom_bar(stat = "identity") +
  mytheme + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- ggplot(p_rta_align_lla, aes(Barcode_lla, N, fill = align)) + geom_bar(position = "fill",stat = "identity") +
  mytheme + theme(axis.text.x = element_text(angle = 45, hjust = 1)) # +scale_fill_manual(c('darkgreen','red','darkblue','black'))
p1/p2
```


#02 SARS_algin_mut
```{linux}
for i in $( ls *_readid_filter.txt | awk -F '_read' '{print $1}');
do 
output="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/"
#align
minimap2 -ax map-ont -t 30  \
       /mnt/raid62/BetaCoV/Person/zhangyiming/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta \
       ${i}_ft.fq |  samtools sort -@ 30 | samtools view -b > $output/02split_SARS_align/${i}_sars.bam
samtools index $output/02split_SARS_align/${i}_sars.bam
## mutation
bcftools mpileup -f /mnt/raid62/BetaCoV/Person/zhangyiming/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta \
   $output/02split_SARS_align/${i}_sars.bam  -Q 0 | bcftools call -mv -o $output/02split_SARS_align/${i}_sars.vcf
done

cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata/02align_filter/
/mnt/raid61/Personal_data/songjunwei/software/qualimap_v2.3/qualimap bamqc \
                    -bam  230808_24sam_comb_filter_sars.bam -outdir ./bam_QC 

```


## mut DecodeR
```{r}
#ALL sample
DecodeR_snp_merge <- read.vcfR("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata/230808_24sam_SARS_mut.vcf")
DecSnp <- data.table(data.frame(DecodeR_snp_merge@fix))
DecSnp$POS <- as.numeric(DecSnp$POS)
DecSnp$QUAL <- as.numeric(DecSnp$QUAL)
DecSnp[QUAL>20]
DecSnp$ID <- paste0(DecSnp$POS,":",DecSnp$REF,">",DecSnp$ALT)

#split sample
snp_list <- list.files("02split_SARS_align/",pattern = ".vcf",full.names = T)
DecSnp_list <- rbindlist(lapply(snp_list, function(x){
  tmp1 <- read.vcfR(x)
  tmp2 <- data.table(data.frame(tmp1@fix))
  tmp2$sample <- gsub("_sars.vcf","",basename(x))
  return(tmp2)
}))
table(DecSnp_list$sample)
DecSnp_list$POS <- as.numeric(DecSnp_list$POS)
DecSnp_list$QUAL <- as.numeric(DecSnp_list$QUAL)
DecSnp_list[QUAL>20]
table(DecSnp_list[QUAL>20]$sample)
DecSnp_list$ID <- paste0(DecSnp_list$POS,":",DecSnp_list$REF,">",DecSnp_list$ALT)
DecSnp_list$FILTER <- NULL
saveRDS(DecSnp_list,"04SNPs_analy/DecodeR_SNPs_raw.rds")
ggvenn(list(merge=DecSnp$ID,split=DecSnp_list$ID))

#BA5.2 XBB.1 sample
snp_list <- list.files("02split_SARS_align/",pattern = "class.vcf",full.names = T)
DecSnp_line <- rbindlist(lapply(snp_list, function(x){
  tmp1 <- read.vcfR(x)
  tmp2 <- data.table(data.frame(tmp1@fix))
  tmp2$class <- gsub("_class.vcf","",basename(x))
  return(tmp2)
}))
table(DecSnp_line$class)
DecSnp_line$POS <- as.numeric(DecSnp_line$POS)
DecSnp_line$QUAL <- as.numeric(DecSnp_line$QUAL)
DecSnp_line[QUAL>20]
```

#03 NGS
```{r}
sampleinfo <- read.csv("24sample_info.csv")
sta_read <- sars_24sam_sub[align=="SARS"][Probability > 0.3, .N, c("Barcode_lla")]
sampleinfo <- merge(sampleinfo,sta_read,by.x="RTA_ID",by.y="Barcode_lla",all=T)
write.csv(sampleinfo,"24sample_info_read.csv")

table(sampleinfo$class1)
sampleinfo$fa_path <- NA
sampleinfo$fa_path[1:11] <- paste0("/mnt/raid64/COVID19/Personal/lhkp/HX_NGS/20230426/04results/04fasta/",sampleinfo$sample_id[1:11],".fa")
sampleinfo$fa_path[12:24] <- paste0("/mnt/raid64/COVID19/Personal/lhkp/HX_NGS/20230616/04results/04fasta/",sampleinfo$sample_id[12:24],".fa")
table(file.exists(sampleinfo$fa_path))
NGS_fa <-lapply(sampleinfo$fa_path, readDNAStringSet)
NGS_fa = do.call(c,NGS_fa)
writeFasta(NGS_fa,"03NGS_res//NGS_strainGV.fa")

sampleinfo$bam_path <- NA
sampleinfo$bam_path[1:11] <- paste0("ln -s ","/mnt/raid64/COVID19/Personal/lhkp/HX_NGS/20230426/02align/Merge_",sampleinfo$Index[1:11],".bam ./",sampleinfo$sample_id[1:11],".bam")
sampleinfo$bam_path[12:24] <- paste0("ln -s ","/mnt/raid64/COVID19/Personal/lhkp/HX_NGS/20230616/02align/Merge_",sampleinfo$Index[12:24],".bam ./",sampleinfo$sample_id[12:24],".bam")
write.table(sampleinfo$bam_path,"./03NGS_res/NGS_bam.txt",quote = F,col.names = F,row.names = F)


sampleinfo$bai_path <- NA
sampleinfo$bai_path[1:11] <- paste0("ln -s ","/mnt/raid64/COVID19/Personal/lhkp/HX_NGS/20230426/02align/Merge_",sampleinfo$Index[1:11],".bam.bai ./",sampleinfo$sample_id[1:11],".bam.bai")
sampleinfo$bai_path[12:24] <- paste0("ln -s ","/mnt/raid64/COVID19/Personal/lhkp/HX_NGS/20230616/02align/Merge_",sampleinfo$Index[12:24],".bam.bai ./",sampleinfo$sample_id[12:24],".bam.bai")
write.table(sampleinfo$bai_path,"./03NGS_res/NGS_bai.txt",quote = F,col.names = F,row.names = F)


sampleinfo$fq_path <- NA
sampleinfo$fq_path[1:11] <- paste0("gzip -c ","/mnt/raid64/COVID19/Personal/lhkp/HX_NGS/20230426/01trim/Merge_",sampleinfo$Index[1:11],".trim.fq  > ./NGS_COVID19_",sampleinfo$RTA_ID[1:11],".fq.gz")
sampleinfo$fq_path[12:24] <- paste0("gzip -c ","/mnt/raid64/COVID19/Personal/lhkp/HX_NGS/20230616/01trim/Merge_",sampleinfo$Index[12:24],".trim.fq > ./NGS_COVID19_",sampleinfo$RTA_ID[12:24],".fq.gz")
write.table(sampleinfo$fq_path,"/mnt/raid61/Personal_data/songjunwei/DRS_RTA/Upload/tmp3/NGS_fastq.txt",quote = F,col.names = F,row.names = F)


sampleinfo$bam_path1 <- NA
sampleinfo$bam_path1[1:11] <- paste0("/mnt/raid64/COVID19/Personal/lhkp/HX_NGS/20230426/02align/Merge_",sampleinfo$Index[1:11],".bam")
sampleinfo$bam_path1[12:24] <- paste0("/mnt/raid64/COVID19/Personal/lhkp/HX_NGS/20230616/02align/Merge_",sampleinfo$Index[12:24],".bam")
sampleinfo[,"sample_type"] <- gsub(" swab","",sampleinfo[,"sample_type"])
write.table(sampleinfo[,c("RTA_ID","bam_path1","sample_type")],"./03NGS_res/NGS_name_bam.txt",quote = F,col.names = F,row.names = F,sep = "\t")
```

```{linux}
#NGS bam QC
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/03NGS_res
/mnt/raid61/Personal_data/songjunwei/software/qualimap_v2.3/qualimap multi-bamqc \
                    -d NGS_name_bam.txt  -outdir ./bam_QC 
 
  
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/03NGS_res
for i in $( ls *bam | awk -F '.' '{print $1}');
do 
bcftools mpileup -f /mnt/raid62/BetaCoV/Person/zhangyiming/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta \
   ${i}.bam  -Q 0 | bcftools call -mv -o ${i}_sars.vcf
done
```

#04 SNPs
```{r}
ngs_snp_list <- list.files("03NGS_res",pattern = ".vcf",full.names = T)
NGSSnp_list <- rbindlist(lapply(ngs_snp_list, function(x){
  tmp1 <- read.vcfR(x)
  tmp2 <- data.table(data.frame(tmp1@fix))
  tmp2$sample <- gsub("_sars.vcf","",basename(x))
  return(tmp2)
}))
table(NGSSnp_list$sample)
NGSSnp_list$POS <- as.numeric(NGSSnp_list$POS)
NGSSnp_list$QUAL <- as.numeric(NGSSnp_list$QUAL)
NGSSnp_list[QUAL>20]
table(NGSSnp_list[QUAL>20]$sample)
NGSSnp_list$ID <- paste0(NGSSnp_list$POS,":",NGSSnp_list$REF,">",NGSSnp_list$ALT)
saveRDS(NGSSnp_list,"04SNPs_analy/NGS_SNPs_raw.rds")


```

```{r fig.width=6}
sampleinfo <- read.csv("24sample_info.csv")
NGSSnp_list <- readRDS("04SNPs_analy/NGS_SNPs_raw.rds")
DecSnp_list <- readRDS("04SNPs_analy/DecodeR_SNPs_raw.rds")

DecSnp_list$mtype <- NA
DecSnp_list$mtype[grep("INDEL",DecSnp_list$INFO)] <- "INDEL"
DecSnp_list$mtype[is.na(DecSnp_list$mtype)] <- "SNP"
DecSnp <- DecSnp_list[mtype=="SNP"][QUAL>20]

NGSSnp_list$mtype <- NA
NGSSnp_list$mtype[grep("INDEL",NGSSnp_list$INFO)] <- "INDEL"
NGSSnp_list$mtype[is.na(NGSSnp_list$mtype)] <- "SNP"
NGSSnp <- NGSSnp_list[mtype=="SNP"][QUAL>20]

# vari_gff <- readGFF("04SNPs_analy/2019-nCoV_EPI_ISL_17953380_variants.gff3")
# vari_gff$ID <- paste0(vari_gff$start,":",vari_gff$REF,">",vari_gff$ALT)
# ggvenn(list(DecodeR=DecSnp$ID,nCoV=vari_gff$ID))


pSNP <- lapply(unique(DecSnp$sample), function(x){
  ggvenn(list(DecodeR=DecSnp[sample==x]$ID,NGS=NGSSnp[RTA_ID==x]$ID))+ggtitle(x)
})
ggarrange(plotlist = pSNP,nrow = 1)

SNPs <- c(
    split.data.frame(DecSnp, DecSnp$sample),
    split.data.frame(NGSSnp[RTA_ID %in% c("RTA_08", "RTA_09", "RTA_27")], NGSSnp[RTA_ID %in% c("RTA_08", "RTA_09", "RTA_27")]$RTA_ID))
names(SNPs) <- c(paste0("Decode",names(SNPs)[1:3]),paste0("NGS",names(SNPs)[1:3]))
library(UpSetR)
upset(fromList(lapply(SNPs,function(x) x$ID)), nsets = 6,order.by = "freq")

ggvenn(list(DecodeR=DecSnp$ID,NGS=NGSSnp[RTA_ID %in% c("RTA_08", "RTA_09", "RTA_27")]$ID))


data <- m6a_ratio_DecodeR[,.N,by="type"]
pieSNP <- ggplot(data, aes(x="", y=precent, fill=type)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  geom_text(aes(label = paste0(precent, "%")), position = position_stack(vjust=0.5)) +
  labs(x = NULL, y = NULL, fill = NULL)+ theme_bw() + theme(legend.text = element_text(size = 14))
ggsave("00result/SARS_SNPs_vail.pdf",pieSNP,width = 6,height = 4)

DecSnp <- DecSnp_list[QUAL>20]
DecSnp %>% group_by(ID) %>% summarise(
  genome=CHROM,
  site=POS,
  REF=REF,
  ALT=ALT,
  sample=paste0(unique(sample),collapse = " ")
)-> DecSnp1

NGSSnp <- NGSSnp_list[QUAL>20]
NGSSnp %>% group_by(ID) %>% summarise(
  genome=CHROM,
  site=POS,
  REF=REF,
  ALT=ALT,
  sample=length(unique(RTA_ID))
)-> NGSSnp1

Snp_sta <- unique(merge(DecSnp1,NGSSnp1[,c("ID","sample")],by="ID"))
colnames(Snp_sta)[6] <- "DecodeR Identi Samples"
colnames(Snp_sta)[7] <- "NGS Identi Samples number"
write.csv(Snp_sta,"00result/04SNPs_DEL_sta.csv")
```


#05 m6A
```{linux}
#2101
conda create -n tom python=3.6
conda install ont-tombo

cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/05m6A/01align_SARS_f5
for i in RTA_08 RTA_09 RTA_10 RTA_15 RTA_19 RTA_26 RTA_27 RTA_29 RTA_32 RTA_36 RTA_40;
do 
conda activate base #2101
fast5_subset -i /mnt/raid61/Personal_data/songjunwei/DRS_RTA/01Rawdata/batch5_SARS/230808_24sam \
            -s ${i}_fast5 \
            --batch_size  1000 \
            --read_id_list ${i}_SARS_readid.txt \
            --recursive  --ignore_symlinks \
            -t 10
seqkit grep -f ${i}_SARS_readid.txt /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata/230808_24sam_comb.fq.gz > ${i}_ft.fq
multi_to_single_fast5 --input_path ${i}_fast5  --save_path ${i}_f5Sg --recursive -t 10 #线程不能太大 卡住跑不过去
ll ${i}_f5Sg/0|wc -l
tombo preprocess annotate_raw_with_fastqs --fast5-basedir ${i}_f5Sg --fastq-filenames ${i}_ft.fq --processes 10 --overwrite
done

conda activate tom #2101
for i in RTA_08 RTA_09 RTA_10 RTA_15 RTA_19 RTA_26 RTA_27 RTA_29 RTA_32 RTA_36 RTA_40;
do 
echo $i
tombo resquiggle ${i}_f5Sg /mnt/raid61/Personal_data/songjunwei/reference/SARS_CoV2/nCoV-2019.reference.fasta \
 --q-score 0  --signal-matching-score 4  \
 --rna --overwrite --processes 20  --ignore-read-locks --num-most-common-errors 5
done

# --signal-matching-score 类似于q值，该值越高，获得的电信号越多
```

## nanom6A process
```{r}
ref_cds <- readDNAStringSet("/mnt/raid61/Personal_data/songjunwei/reference/SARS_CoV2/SARS-CoV-2_NC_045512_CDS.fasta")
cds_n <- names(ref_cds)
cds_n <- gsub(" \\[Severe acute respiratory syndrome coronavirus 2\\]","",cds_n)
write.table(cds_n,"/mnt/raid61/Personal_data/songjunwei/reference/SARS_CoV2/cds_transcript.txt",quote = F,col.names = F,row.names = F)

cds_tran <- read.table("/mnt/raid61/Personal_data/songjunwei/reference/SARS_CoV2/cds_transcript.txt")
cds_tran <- cds_tran %>% group_by(V2) %>% summarise(new_transcript = paste(V1, collapse = " "))
write.table(cds_tran,"/mnt/raid61/Personal_data/songjunwei/reference/SARS_CoV2/cds2transcript.txt",quote = F,col.names = F,row.names = F)
```

```{linux}
conda activate nanom6a  #2048
refpath="/mnt/raid61/Personal_data/songjunwei/reference/SARS_CoV2/"
Nanom6ASof="/mnt/raid5/Personal/minion/software/nanom6A_2021_10_22" #2048
nanom6apath="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/05m6A/02nanom6a"
#genome.fa 建立index
# java -jar /mnt/raid61/Personal_data/songjunwei/software/picard.jar CreateSequenceDictionary \
#       R=$refpath/nCoV-2019.reference.fasta \
#       O=$refpath/nCoV-2019.reference.dict

for i in RTA_08 RTA_09 RTA_10 RTA_15 RTA_19 RTA_26 RTA_27 RTA_29 RTA_32 RTA_36 RTA_40;
do 
#find fast5
find /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/05m6A/01align_SARS_f5/${i}_f5Sg -name "*.fast5" > $nanom6apath/${i}_f5name.txt
echo $i
# extracting signals
python $Nanom6ASof/extract_raw_and_feature_fast.py --cpu=20 --fl=$nanom6apath/${i}_f5name.txt -o $nanom6apath/${i}_result --clip=10 
# predicting m6A site cds2transcritpt 见下面的代码
python $Nanom6ASof/predict_sites.py --cpu 20 --model $Nanom6ASof/bin/model \
        -r $refpath/SARS-CoV-2_NC_045512_CDS.fasta \
        -g $refpath/nCoV-2019.reference.fasta -b $refpath/cds2transcript.txt \
        -i $nanom6apath/${i}_result \
        -o $nanom6apath/${i}_result_final 
# Visualization of m6A sites 
python $Nanom6ASof/nanoplot.py --input $nanom6apath/${i}_result_final \
   -o $nanom6apath/${i}_plot
done
```


## tandemMod
```{linux}
# git clone https://github.com/yulab2021/TandemMod
# pip install torch==1.9.1
# usage https://yulab2021.github.io/TandemMod_document/usage.html

conda activate base # 2101
input="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/05m6A/01align_SARS_f5"
TandemMod="/mnt/raid61/Personal_data/songjunwei/software/TandemMod/scripts/"
reffa="/mnt/raid62/BetaCoV/Person/zhangyiming/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta"
model_file="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/05m6A/04tandemMod/RNA_modifi_motif.txt"
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/05m6A/04tandemMod
for i in RTA_08 RTA_09 RTA_10 RTA_15 RTA_19 RTA_26 RTA_27 RTA_29 RTA_32 RTA_36 RTA_40;
do 
# Extract signal files and features from resquiggled fast5 files using the following python scripts.
minimap2 -ax map-ont $reffa -t 30   $input/${i}_ft.fq  > ${i}_ft.sam
python $TandemMod/extract_signal_from_fast5.py -p 40 --fast5 $input/${i}_f5Sg \
      --reference $reffa --sam ${i}_ft.sam --output ${i}_signal.tsv --clip 10

#Predict modifications using pretrained model.
cat $model_file | while read rows; #
do
type=$(echo $rows | awk -F ' ' '{print $1}')
motif=$(echo $rows | awk -F ' ' '{print $2}')
model=$(echo $rows | awk -F ' ' '{print $3}')
echo $type
echo $motif
python $TandemMod/extract_feature_from_signal.py  --signal_file ${i}_signal.tsv \
      --clip 10 --output ${i}_${type}_feature.tsv --motif $motif
python $TandemMod/TandemMod.py --run_mode predict \
    --pretrained_model $TandemMod/../models/${model} \
    --feature_file ${i}_${type}_feature.tsv \
    --predict_result ${i}_${type}_prediction.tsv
done
done


## analysis
```{r}
## nanom6A 
nanom6A_8 <- m6a_process(resFin_path = "/mnt/raid61/Personal_data/hezhifeng/nanom6A/result/RTA_08_result_final/")
nanom6A_9 <- m6a_process(resFin_path = "/mnt/raid61/Personal_data/hezhifeng/nanom6A/result/RTA_09_result_final/")
nanom6a_res <- rbind(nanom6A_8$ratio,nanom6A_9$ratio)
nanom6a_res$IDs <- paste0(nanom6a_res$gene,":",nanom6a_res$pos,":+")
nanom6a_res <- nanom6a_res[readN>10 & gene != "ORF10"]
nanom6a_res$IDs <- paste0(nanom6a_res$chr,":",nanom6a_res$pos,":+")
nanom6a_tab <- dcast(unique(nanom6a_res[,c("sample","IDs","ratio")]),formula = IDs ~ sample,value.var = "ratio")
write.csv(nanom6a_tab,"05m6A/nanom6a_tab_ratio.csv")
# test <- merge(m6a_tab,nanom6a_tab,by="IDs",all=T)

## tandemod
tande_res <- tande_process("05m6A/04tandemMod/")
tande_res$IDs <- paste0("MN908947.3",":",tande_res$site,":+")
write.table(tande_res,"05m6A/tande_res.txt",quote = F,row.names = F)

tande_res_sub <- tande_res[sum>=20]
table(tande_res_sub$sample,tande_res_sub$mod_type)
tande_tab_list <- lapply(split.data.frame(tande_res_sub,tande_res_sub$mod_type), function(tande_res_sub_sp){
  tande_tab <- dcast(unique(tande_res_sub_sp[,c("sample","IDs","ratio")]),formula = IDs ~ sample,value.var = "ratio")
  write.csv(tande_tab,paste0("05m6A/tande_",unique(tande_res_sub_sp$mod_type),"_tab_ratio.csv"))
  tande_tab
})
View(tande_tab_list$m7G)

#m6A-net
m6a_net_sars <- read.csv("/mnt/raid61/Personal_data/liudefu/m6Anet_m6A_tab_ratio.csv")

# cell publish
mmc5 <- data.table::data.table(openxlsx::read.xlsx("cell_2020_SARS_mmc5.xlsx"))
mmc5$IDs <- paste0("MN908947.3:",mmc5$`pos+1`,":+")

tande_tab_list
```

```{r }
## plot
ggvenn(list(
  nano08 = nanom6a_res[sample == "RTA_08"]$IDs,
  nano09 = nanom6a_res[sample == "RTA_09"]$IDs,
  tand08 = tande_res_sub[sample == "RTA_08" & mod_type == "m6A"]$IDs,
  tand09 = tande_res_sub[sample == "RTA_09" & mod_type == "m6A"]$IDs
            ))
#tande
ggvenn(lapply(split.data.frame(tande_res_sub,f = tande_res_sub$mod_type), function(x) x$IDs))
ggvenn(lapply(split.data.frame(tande_res_sub[mod_type=="m1A"],f = tande_res_sub[mod_type=="m1A"]$sample), function(x) x$IDs))+ggtitle("m1A_tandem") +
  ggvenn(lapply(split.data.frame(tande_res_sub[mod_type=="m6A"],f = tande_res_sub[mod_type=="m6A"]$sample), function(x) x$IDs))+ggtitle("m6A_tandem") +
  ggvenn(lapply(split.data.frame(tande_res_sub[mod_type=="m5C"],f = tande_res_sub[mod_type=="m5C"]$sample), function(x) x$IDs))+ggtitle("m5C_tandem") +
  ggvenn(lapply(split.data.frame(tande_res_sub[mod_type=="m7G"],f = tande_res_sub[mod_type=="m7G"]$sample), function(x) x$IDs))+ggtitle("m7G_tandem")

#publish
ggvenn(list(
    #nanom6a = m6a_res$pos,
    m6A_tand = tande_res_sub[mod_type == "m6A"]$site,
    m1A_tand = tande_res_sub[mod_type == "m1A"]$site,
    Cell2020 = mmc5[mod_base == "A"]$`pos+1`
  )) 
ggvenn(list(
    m5C_tand=tande_res_sub[ mod_type == "m5C"]$site,
    Cell2020 = mmc5[mod_base == "C"]$`pos+1`
  ))+
ggvenn(list(
   m7G_tand=tande_res_sub[ mod_type == "m7G"]$site,
    Cell2020 = mmc5[mod_base == "G"]$`pos+1`
  ))

#motif
pmoti <- lapply(unique(tande_res_sub$mod_type), function(x){
  ggseqlogo(tande_res_sub[mod_type==x]$motif, method = 'prob') + ggtitle(x)
})
p_moti <- ggarrange(plotlist = pmoti)
ggsave("00result/SARS-CoV2-RNA_modtif.pdf",plot = p_moti,width = 5,height = 4)
```

```{r}
cell_Res <-paste0("MN908947.3:",c("29428","29450","29522"),":+")## Cell Research Published
vail_exp <- c(intersect(mmc5$IDs,tande_res_sub$IDs),cell_Res)
dcast(tande_res_sub[IDs %in% vail_exp ],IDs + mod_type ~ sample ,value.var = "ratio")
dcast(m6a_res[IDs %in% vail_exp ],IDs  ~ sample ,value.var = "ratio")
write.table(vail_exp,"05m6A/05tombo_plot/vail_siteFrom_Publish.txt",quote = F,col.names = F,row.names = F) #两个样本都画

m6A_corrd <- setdiff(c(tande_res_sub[ratio >= 0.5 & mod_type == "m1A"]$IDs,
            tande_res_sub[ratio >= 0.8 & mod_type == "m7G"]$IDs,
            tande_res_sub[ratio >= 0.6 & mod_type == "m5C"]$IDs),
          vail_exp)
write.table(m6A_corrd,"05m6A/05tombo_plot/Top_m1A_m5C_m7G.txt",quote = F,col.names = F,row.names = F)

ggvenn(list(tand=tande_tab_list$m6A$IDs,
            nanom6A=nanom6a_tab$IDs))
comm <- intersect(tande_tab_list$m6A$IDs,nanom6a_tab$IDs)
all <- unique(c(comm,vail_exp,m6A_corrd))
write.table(all,"05m6A/05tombo_plot/all_site.txt",quote = F,col.names = F,row.names = F)
```



## SNPs & m6A
```{r}
m6a_net_sars <- read.csv("/mnt/raid61/Personal_data/liudefu/m6Anet_m6A_tab_ratio.csv")
m6a_net_sars$site <- as.numeric(str_split_fixed(m6a_net_sars$IDs,":",3)[,2])
m6a_net_sars$chr <- str_split_fixed(m6a_net_sars$IDs,":",3)[,1]
m6a_net_sars$IDs <- paste0(m6a_net_sars$chr,":",m6a_net_sars$site+1,":+")

sars_net <- m6a_net_sars[rowSums(m6a_net_sars[,3:4],na.rm = T)>0,]
sars_net$N <- apply(sars_net[,3:4],1,function(x){sum(x>0,na.rm = T)})
sars_tande <- tande_res_sub[mod_type=="m6A"][,.N,by=IDs]
sars_nano <- nanom6a_res[,.N,by=IDs]
sars_m6a_3soft <- merge(merge(sars_net[,c("IDs","N")],sars_nano,by="IDs",all=T),sars_tande,by="IDs",all=T)
colnames(sars_m6a_3soft) <- c("m6A IDs","m6A-net","nanom6A","TandeMod")
sars_m6a_3soft$Software_number <- apply(sars_m6a_3soft[,2:4],1,function(x){sum(x>0,na.rm = T)})
sars_m6a_3soft <- sars_m6a_3soft[order(sars_m6a_3soft$`m6A IDs`),]
write.csv(sars_m6a_3soft,"00result/SARS-CoV2-m6A_list.csv")

m6a_data <- data.frame(IDs=data.table(sars_m6a_3soft)[Software_number==3 &TandeMod==2 ]$`m6A IDs`)
m6a_data <- data.frame(IDs=intersect(intersect(sars_nano[N>1]$IDs,sars_tande[N>1]$IDs),sars_net[N>1]$IDs))
m6a_data$chr <- str_split_fixed(m6a_data$IDs,":",3)[,1]
m6a_data$pos <- str_split_fixed(m6a_data$IDs,":",3)[,2]
m6a_data$pos <- as.numeric(m6a_data$pos)
```


```{r}
# 以m6A为中心
snp_dist<- sapply(snp_data$POS, function(site) {
  tmp0 <- site - m6a_data$pos
  tmp0[tmp0>=-500 & tmp0<=500]
}) %>% unlist()
hist(snp_dist)
snp_dist <- data.frame(distance_to_m6a=snp_dist)
set.seed(275)
random_site <- sample(0:30000, 300)
random_site <- setdiff(random_site,snp_data$POS)
random_distances <- sapply(random_site, function(site) {
  tmp0 <- site - m6a_data$pos
  tmp0[tmp0>=-500 & tmp0<=500]
}) %>% unlist()
# hist(random_distances)
random_background <- data.frame(distance_to_m6a=random_distances)
p <- wilcox.test(snp_dist$distance_to_m6a,random_background$distance_to_m6a)
plot_m6a6<- ggplot() +
  geom_density(data = snp_dist, aes(x = distance_to_m6a, color = "m6A"), alpha = 0.8) +
  geom_density(data = random_background, aes(x = distance_to_m6a, color = "Random"), alpha = 0.5) +
  xlab("Distance to m6A Site") +
  ylab("Density") +
  ggtitle(i) + mytheme+
  labs(subtitle = (paste0("Wilcoxon, P = ",round(p$p.value,4))))+
  scale_color_manual(values = c("m6A" = "red", "Random" = "gray"))
print(plot_m6a6)
ggsave("00result/SARS-CoV2-SNPs_m6A_corr111.pdf",plot_m6a6,width = 6,height = 4)
```


#06 microorganism
## kraken2
```{linux }
#Kraken2
## build index
/mnt/raid61/Personal_data/songjunwei/software/kraken2-master/install_path/kraken2-build --standard --threads 24 --db $DBNAME
#db can download https://benlangmead.github.io/aws-indexes/k2

## search
kraken2="/mnt/raid61/Personal_data/songjunwei/software/kraken2-master/install_path/kraken2"
fq="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata/230808_24sam_comb_filter_dna.fq"
### viral
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/06microorganism/01kraken2_res
index="/mnt/raid61/Personal_data/songjunwei/reference/Kraken2_index/viral"
$kraken2 --db $index --threads 10 --classified-out viral_class.fq  --output viral_reads.txt --report viral_report.txt  \
  --use-names --confidence 0 --minimum-hit-groups 2 $fq
### standard
index="/mnt/raid61/Personal_data/songjunwei/reference/Kraken2_index/standard_all"
$kraken2 --db $index --threads 10 --classified-out standard_class.fq  --unclassified-out standard_unclass.fq --output standard_reads.txt --report standard_report.txt  \
  --use-names --confidence 0 --minimum-hit-groups 2 --report-zero-counts $fq
### PlusPF_16
index="/mnt/raid61/Personal_data/songjunwei/reference/Kraken2_index/PlusPF_16/"
$kraken2 --db $index --threads 10 --classified-out PlusPF_class.fq --unclassified-out PlusPF_unclass.fq  --output PlusPF_reads.txt --report PlusPF_report.txt  \
  --use-names --confidence 0 --minimum-hit-groups 2 --report-zero-counts $fq
### standard过滤后
$kraken2 --db $index --threads 10  --output standard_unclass_PlusPF_reads.txt --report standard_unclass_PlusPF_report.txt  \
  --use-names --confidence 0 --minimum-hit-groups 2 standard_unclass.fq
```

report
```{r fig.width=4,fig.height=4}
kraken2_read_report <- function(report.file){
  fread(file = normalizePath(report.file),
        sep = "\t",
        data.table = F,
        col.names = c("percent", "clade_count", "tax_count", "rank", "tax_id", "name")) -> tbl
  return(tbl)
}
# (U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G) enus，或（S）pecies
# U未分类，R根，D域，k界，P门，C纲，O目，F科，G属，S种

stand_report <- data.table(kraken2_read_report("06microorganism/01kraken2_res/standard_report.txt"))
write.csv(stand_report,"06microorganism/00result/01standard_sta.csv",row.names = F)
stand_report[rank %in% grep("R",stand_report$rank,value = T)]
stand_report[rank %in% grep("D",stand_report$rank,value = T)]
stand_report[rank %in% grep("K",stand_report$rank,value = T)]

# PlusPF_report <- data.table(kraken2_read_report("06microorganism/01kraken2_res/PlusPF_report.txt"))
# PlusPF_report[rank %in% grep("R",PlusPF_report$rank,value = T)]
# PlusPF_report[rank %in% grep("D",PlusPF_report$rank,value = T)]
# PlusPF_report[rank %in% grep("K",PlusPF_report$rank,value = T)]
# write.csv(PlusPF_report,"06microorganism/00result/01PlusPF_sta.csv",row.names = F)

class_sta <- stand_report[rank %in% c("U","D")]
stand_report[name=="Archaea"]
class_sta$category <- "category"
class_sta$percent <- class_sta$clade_count/sum(class_sta$clade_count)*100
class_sta$name <- paste0(class_sta$name,"(",round(class_sta$percent,2),"%)")
class_sta$name <- factor(class_sta$name,levels = c("Bacteria(59.66%)","Eukaryota(15.77%)","Viruses(0.1%)","Archaea(0.01%)","unclassified(24.47%)"))
sum(class_sta$percent)
sum(class_sta$clade_count)

p_all_sta <- ggplot(class_sta,aes(x=category,y=percent,fill=name)) + geom_bar(stat = "identity") +
  mytheme + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab(NULL) + ylab("Percentage%")
ggsave("06microorganism/00result/01standard_catology.pdf",width = 4,height = 4)
```

reads
```{r}
sars_24sam_sub <- readRDS("01filter_rawdata/Combin_QC_align_filter.rds")
stand_blast <- fread("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/06microorganism/01kraken2_res/standard_reads.txt")
colnames(stand_blast) <- c("Tag","readid","species","align_region","LCA")

stand_blast$species[stand_blast$species=="unclassified (taxid 0)"] <- "unclassified"
# stand_blast$species <- ifelse(stand_blast$species %in% stand_blast[,.N,by="species"][N<=10]$species,"Others",stand_blast$species)
sars_24sam_kraken <- merge(sars_24sam_sub,stand_blast[,c("readid","species")])
ggvenn(list(minimap = sars_24sam_sub[align == "SARS"]$readid, kraken =
                      stand_blast[grep("Severe acute respiratory syndrome",stand_blast$species),]$readid))+ggtitle("SARS-CoV-2")
ggvenn(list(minimap = sars_24sam_sub[align == "homsap"]$readid, kraken =
                      stand_blast[species == "Homo sapiens (taxid 9606)"]$readid))+ggtitle("Human")
length(unique(sars_24sam_kraken$species))

allsample_sta <- sars_24sam_kraken[,.N,by=species][order(N,decreasing = T)][N>10]
write.csv(allsample_sta,"00result/micro_species_sta.csv")
View(sars_24sam_kraken[,.N,by=species][order(N,decreasing = T)])
```


```{r fig.width=12,fig.height=4}
#合并某些属
N=12
sars_24sam_kraken$class <- NA
sars_24sam_kraken$class[grep("Veillonella",sars_24sam_kraken$species)] <- "Veillonella"
sars_24sam_kraken$class[grep("Streptococcus",sars_24sam_kraken$species)] <- "Streptococcus"
sars_24sam_kraken$class[grep("Prevotella",sars_24sam_kraken$species)] <- "Prevotella"
sars_24sam_kraken$class[grep("Neisseria gonorrhoeae",sars_24sam_kraken$species)] <- "Neisseria gonorrhoeae"
sars_24sam_kraken$class[sars_24sam_kraken$align=="SARS"] <- "SARS-CoV-2"

# other 未过滤
sars_24sam_kraken$class[sars_24sam_kraken$species=="Bacteria (taxid 2)"] <- "Other bacteria"
index <- is.na(sars_24sam_kraken$class)
sars_24sam_kraken$class[index] <- sars_24sam_kraken$species[index] 
major <- c(head(sars_24sam_kraken[,.N,by=class][order(N,decreasing = T)],N)$class,"SARS-CoV-2")
sars_24sam_kraken$class <- ifelse(sars_24sam_kraken$class %in% major,sars_24sam_kraken$class,"Other bacteria")
major <- c(head(sars_24sam_kraken[,.N,by=class][order(N,decreasing = T)],N)$class,"SARS-CoV-2")

sars_24sam_kraken[class=="Other bacteria"][,.N,by=c("species")][N>10]
#Other bacteria(n=706)

sta_tab_class <- sars_24sam_kraken[,.N,by=c("class")][order(N,decreasing = T)]
sta_tab_class$taxid <- gsub(")","",str_split_fixed(sta_tab_class$class,"\\(taxid ",2)[,2])
sta_tab_class$class1 <- str_split_fixed(sta_tab_class$class,"\\(taxid ",2)[,1]
# lapply(c("Veillonella","Streptococcus","Prevotella","Neisseria gonorrhoeae","Severe acute"), function(x){
#   id <- unique(sars_24sam_kraken$species[grep(x,sars_24sam_kraken$species)])
#   id
# })
sta_tab_class$taxid[1] <- "2"
sta_tab_class$taxid[2] <- "0"
sta_tab_class$taxid[4] <- "838"
sta_tab_class$taxid[5] <- "29465"
sta_tab_class$taxid[6] <- "1301"
sta_tab_class$taxid[9] <- "485"
sta_tab_class$taxid[13] <- "2697049"
sta_tab_class$class1[1] <- "Other bacteria(n=705)"
sta_tab_class$class1[13] <- "Severe acute respiratory syndrome coronavirus 2"
sta_tab_class$`Precentage\\%` <- round(sta_tab_class$N/sum(sta_tab_class$N)*100,2)
write.csv(sta_tab_class,"06microorganism/00result/01classified_sta_5type.csv")
#add Archaea

## plot
class_sta <-data.frame(
  name = c("Bacteria(59.79%)","Eukaryota(15.71%)","Viruses(0.1%)","Archaea(0.01%)","unclassified(24.38%)"),
  count = c(223744, 58778, 406, 37, 91239))
class_sta$category <- "category"
class_sta$percent <- class_sta$count/sum(class_sta$count)*100
class_sta$name <- factor(class_sta$name,levels = c("Bacteria(59.79%)","Eukaryota(15.71%)","Viruses(0.1%)","Archaea(0.01%)","unclassified(24.38%)"))
sum(class_sta$percent)
sum(class_sta$count)
p1 <- ggplot(class_sta,aes(x=category,y=percent,fill=name)) + 
  geom_bar(stat = "identity") +
  mytheme + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab(NULL) + ylab("Percentage%")


sta_tab_class$class1 <- paste0(sta_tab_class$class1,"(taxid ",sta_tab_class$taxid,")")
p_kraken <- sars_24sam_kraken[Probability>0.3][class != "unclassified"][,.N, c("Barcode_lla","class")]
p_kraken$class1 <- plyr::mapvalues(p_kraken$class,from = sta_tab_class$class,to = sta_tab_class$class1)
p_kraken$class1 <- factor(p_kraken$class1,levels = rev(sta_tab_class$class1))
p2 <- ggplot(p_kraken, aes(Barcode_lla, N, fill = class1)) + 
  geom_bar(position = "fill", stat = "identity") +
  mytheme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = scales::brewer_pal(palette = "BrBG")(12)) + ylab("Percentage%")

p <- p1+p2+plot_layout(widths = c(1,10))
ggsave("06microorganism/00result/02standard_top10_plot.pdf",plot = p,device = "pdf",width = 17,height = 4)
# sta_tab <- dcast(p_kraken,formula = Barcode_lla ~ class)
# sta_tab[is.na(sta_tab)] <- 0
# write.csv(sta_tab,"06microorganism/00result/02standard_top10_RTA_martix.csv")
```

```{r fig.width=10,fig.height=5}
data_kraken <- dcast(p_kraken,formula = Barcode_lla ~ class1,value.var = "N")
data_kraken[is.na(data_kraken)] <- 0
tmp <- data.frame(data_kraken)[,2:13] %>% as.matrix()
mx_karken <- prop.table(tmp, margin = 1)
rownames(mx_karken) <- data_kraken$Barcode_lla

meta_1 <- sample_de
meta_1 <- meta_1[rownames(mx_karken),]
anno_col <- data.frame(type=meta_1[,"sample_type"])
rownames(anno_col) <- rownames(meta_1)
pheatmap::pheatmap(t(data.frame(mx_karken)),cluster_cols = T,cluster_rows = T,annotation_col = anno_col)
```


## 微生态差异
```{r fig.width=9,fig.height=4}
sars_24sam_kraken$class2 <- ifelse(sars_24sam_kraken$species %in%  sars_24sam_kraken[,.N,by=species][N>100]$species,sars_24sam_kraken$species,"Other bacteria")
sars_24sam_kraken[class2=="Bacteria (taxid 2)"]$class2 <- "Other bacteria"
sta_RTA <- sars_24sam_kraken[Probability > 0.3, .N, by = c("class2", "Barcode_lla")][order(N, decreasing = T)]

sta_RTA_tab <- dcast(sta_RTA,formula = Barcode_lla ~ class2)
sta_RTA_tab[is.na(sta_RTA_tab)] <- 0
sta_RTA_tab <- sta_RTA_tab %>% tibble::column_to_rownames("Barcode_lla")
sta_RTA_tab <- data.frame(t(sta_RTA_tab))
saveRDS(sta_RTA_tab,"06microorganism/00result/03RTA_species_counts.rds")


library(DESeq2)
sample_de <- read.csv("24sample_info.csv")
rownames(sample_de) <- sample_de$RTA_ID
dds <- DESeqDataSetFromMatrix(countData = sta_RTA_tab,
                              colData = sample_de,
                              design = ~ sample_type )
dds <- DESeq(dds)
ddsv <- varianceStabilizingTransformation(dds, blind = FALSE)
res <- results(dds)
res_sub <- res[res$pvalue<0.05,]
sta_RTA[class2 %in% rownames(res_sub)][,sum(N),by=class2]
dcast(sta_RTA[class2 %in% rownames(res_sub)],formula = Barcode_lla ~ class2)

## heatmap
mat <- sta_RTA_tab[rownames(res_sub),]
mat <- mat[,colSums(mat)>=10]
anno_col <- data.frame(type=sample_de[,"sample_type"])
rownames(anno_col) <- rownames(sample_de)
pde_heat <- pheatmap::pheatmap(plot_mat,scale = "column",cluster_cols = T,cluster_rows = T,annotation_col = anno_col)
ggsave("06microorganism/00result/03DE_species_col0_vst.pdf",pde_heat,width = 10,height = 4)


#PCA
counts <- assay(ddsv)
sample_de$sample_type <- gsub(" swab","",sample_de$sample_type)
meta <- sample_de
mat <- mat[,colSums(mat)>=10]
meta <- sample_de[colnames(mat),]
counts <- counts[,colnames(mat)]
meta$sample_type <- factor(meta$sample_type)
p_pca <- DEGreport::degPCA(counts,metadata=meta,condition="sample_type") + coord_fixed() +theme_bw()+
  scale_fill_manual(labels =c(nasopharyngeal="#00BFC4",Oropharyngeal="#F8766D"))+
                ggrepel::geom_text_repel(aes(label = colnames(counts))) #+theme(legend.position="none")
ggsave("06microorganism/00result/03species_PCA.pdf",p_pca,width = 6,height = 4)
```

## RNA修饰
```{r}
View(sars_24sam_kraken[,.N,by=species][order(N,decreasing = T)])
sars_24sam_kraken[,.N,by=class][order(N,decreasing = T)]
# max reads: "Homo sapiens (taxid 9606)" "Veillonella (taxid 29465)" ##? 不知道参考"Streptococcus (taxid 1301)"

kraken_24sam <- sars_24sam_kraken[Probability > 0.3]
homo <- kraken_24sam[species=="Homo sapiens (taxid 9606)"]
table(homo$Barcode_lla,homo$species)
write.table(unique(homo$Barcode_lla),"24sample_name.txt",quote = F,row.names = F,col.names = F)


Veill <- kraken_24sam[class=="Veillonella"]
Strep <- kraken_24sam[class=="Streptococcus"]
Prevo <- kraken_24sam[class=="Prevotella"]
Bacill <- kraken_24sam[species %in% c("Bacillus thuringiensis (taxid 1428)","Bacillus (taxid 1386)")]
Fusob <- kraken_24sam[species %in% "Fusobacterium pseudoperiodonticum (taxid 2663009)"]

x1 <- cbind(table(Prevo$Barcode_lla, Prevo$class),
      table(Strep$Barcode_lla, Strep$class),
      table(Veill$Barcode_lla, Veill$class))
x1 <- data.frame(x1) %>% tibble::rownames_to_column("Barcode_lla")
x2 <- Bacill[,.N,by=Barcode_lla]
x3 <- Fusob[,.N,by=Barcode_lla]
x4 <- merge(x2,x3,by="Barcode_lla",all=T)
colnames(x4) <- c("Barcode_lla","Bacill","Fusob")
micro_tab <- merge(x1,x4,by="Barcode_lla",all=T)


lapply(unique(homo$Barcode_lla), function(x){
  write.table(homo[Barcode_lla==x]$readid,paste0("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/06microorganism/04fast5/",x,"_homo_readid.txt"),quote = F,col.names = F,row.names = F)
  write.table(Veill[Barcode_lla==x]$readid,paste0("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/06microorganism/04fast5/",x,"_Veill_readid.txt"),quote = F,col.names = F,row.names = F)
  write.table(Strep[Barcode_lla==x]$readid,paste0("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/06microorganism/04fast5/",x,"_Strep_readid.txt"),quote = F,col.names = F,row.names = F)
  write.table(Prevotella[Barcode_lla==x]$readid,paste0("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/06microorganism/04fast5/",x,"_Prevo_readid.txt"),quote = F,col.names = F,row.names = F)
  write.table(Bacill[Barcode_lla==x]$readid,paste0("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/06microorganism/04fast5/",x,"_Bacill_readid.txt"),quote = F,col.names = F,row.names = F)
  write.table(Fusob[Barcode_lla==x]$readid,paste0("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/06microorganism/04fast5/",x,"_Fusob_readid.txt"),quote = F,col.names = F,row.names = F)
})

table(kraken_24sam)

```

###tombo
```{linux}
conda activate base #2101
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/06microorganism/04fast5/
name="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/24sample_name.txt"
for species in Prevo Bacill Fusob;do #homo Veill Strep
for i in  $(cat $name);
do 
fast5_subset -i /mnt/raid61/Personal_data/songjunwei/DRS_RTA/01Rawdata/batch5_SARS/230808_24sam \
            -s ${i}_${species}_fast5 \
            --batch_size  1000 \
            --read_id_list ${i}_${species}_readid.txt \
            --recursive  --ignore_symlinks \
            -t 10
seqkit grep -f ${i}_${species}_readid.txt /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata/230808_24sam_comb.fq.gz > ${i}_${species}_ft.fq
multi_to_single_fast5 --input_path ${i}_${species}_fast5  --save_path ${i}_${species}_f5Sg --recursive -t 5 #线程不能太大 卡住跑不过去
ll ${i}_${species}_f5Sg/0|wc -l
tombo preprocess annotate_raw_with_fastqs --fast5-basedir ${i}_${species}_f5Sg --fastq-filenames ${i}_${species}_ft.fq --processes 10 --overwrite
done
done

conda activate tom #2101
ref_file="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/ref_path1"
name="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/24sample_name.txt"
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/06microorganism/
cat $ref_file | while read rows; 
do
species=$(echo $rows | awk -F ' ' '{print $1}')
ref=$(echo $rows | awk -F ' ' '{print $2}')
echo $species
echo $ref

mkdir -p 05align/${species}/bam/
for i in $(cat $name);
do
## Alignment ## 
# base 2101
minimap2 -ax map-ont -t 30 $ref/genomic.fa 04fast5/${i}_${species}_ft.fq  > 05align/${species}/${i}_${species}.sam ## 跑tandemMod必须是-ax map-ont
samtools sort -@ 30 05align/${species}/${i}_${species}.sam | samtools view -b > 05align/${species}/bam/${i}_${species}.bam
samtools index 05align/${species}/bam/${i}_${species}.bam

## Tombo re-squiggle ## 
#tom 2101
echo $i
tombo resquiggle 04fast5/${i}_${species}_f5Sg $ref/genomic.fa --q-score 0  --signal-matching-score 4 \
 --rna --overwrite --processes 20  --ignore-read-locks --num-most-common-errors 5
done
done

##QC
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/06microorganism/05align
NanoComp --bam  $(ls */bam/*bam) -o /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/06microorganism/00QC -t 30  -p miciro_bam_ 
for i in $(ls */bam/*bam);do 
  echo $i
  samtools flagstat $i | sed -n '5p'
done
```

###tandemMod 
```{linux}
conda activate base # 2101
TandemMod="/mnt/raid61/Personal_data/songjunwei/software/TandemMod/scripts/"
name="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/24sample_name.txt"
model_file="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/05m6A/04tandemMod/RNA_modifi_motif.txt"
ref_file="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/ref_path1"
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/06microorganism/
cat $ref_file | while read rows; 
do
species=$(echo $rows | awk -F ' ' '{print $1}')
ref=$(echo $rows | awk -F ' ' '{print $2}')
echo $species
echo $ref
for i in $(cat $name);
do
python $TandemMod/extract_signal_from_fast5.py -p 30 --fast5 04fast5/${i}_${species}_f5Sg \
      --reference $ref/genomic.fa --sam 05align/${species}/${i}_${species}.sam \
      --output 06tandMod/${i}_${species}_signal.tsv --clip 10
cat $model_file | while read rows; 
do
type=$(echo $rows | awk -F ' ' '{print $1}')
motif=$(echo $rows | awk -F ' ' '{print $2}')
model=$(echo $rows | awk -F ' ' '{print $3}')
echo $type
echo $motif
python $TandemMod/extract_feature_from_signal.py  --signal_file 06tandMod/${i}_${species}_signal.tsv \
      --clip 10 --output 06tandMod/${i}_${type}_${species}_feature.tsv --motif $motif
python $TandemMod/TandemMod.py --run_mode predict \
    --pretrained_model $TandemMod/../models/${model} \
    --feature_file 06tandMod/${i}_${type}_${species}_feature.tsv \
    --predict_result 06tandMod/${i}_${type}_${species}_prediction.tsv
done
done
done


###nanom6A 
```{r}
#sed 's/lcl|//g' cds_from_genomic.fa > cds_from_genomic_rename.fa
cds_path <- list.files("/mnt/raid61/Personal_data/songjunwei/reference/ncbi_genome/",recursive=T,pattern = "cds_from_genomic_rename.fa",full.names = T)
lapply(cds_path, function(x){
  cds <- read.fasta(x)
  cds_n <- strsplit(names(cds)," \\[")
  head(cds_n)
  cds_df0 <- lapply(cds_n, function(x){
    x <- grep("gene=",x,invert = T,value = T)
    x <- gsub("]","",x)
    x <- gsub("db_xref=","",x)
    data.frame(t(x))[,1:4]
  })
  cds_df <- rbindlist(cds_df0,fill = T)
  colnames(cds_df) <- c("CDS_name","locus_tag","db_xref","protein") #,"protein_id","location","gbkey"
  table(duplicated(cds_df$db_xref))
  write.table(cds_df,gsub("cds_from_genomic_rename.fa","cds_table.txt",x),quote = F,row.names = F)
  write.table(cds_df[,c("db_xref","CDS_name")],gsub("cds_from_genomic_rename.fa","cds2transcritp.txt",x),quote = F,row.names = F,col.names = F)
})
# gtf <- rtracklayer::import('/mnt/raid61/Personal_data/songjunwei/reference/ncbi_genome/Veillonella_parvula_NCTC11810/GCF_900186885.1/genomic.gtf')
# gtf_tmp <- unique(gtf[,c("gene_id","gene","db_xref")])
# length(unique(gtf_tmp$db_xref))
# table(cds_df$db_xref %in% gtf_tmp$db_xref )
```

```{linux}
conda activate nanom6a  #2048
name="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/24sample_name.txt"
Nanom6ASof="/mnt/raid5/Personal/minion/software/nanom6A_2021_10_22" #2048
nanom6apath="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/06microorganism/07nanom6A/"
ref_file="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/ref_path1"
cat $ref_file | while read rows; 
do
species=$(echo $rows | awk -F ' ' '{print $1}')
ref=$(echo $rows | awk -F ' ' '{print $2}')
echo $species
echo $ref
# java -jar /mnt/raid61/Personal_data/songjunwei/software/picard.jar CreateSequenceDictionary \
#       R=$ref/genomic.fa \
#       O=$ref/genomic.dict
# samtools faidx $ref/genomic.fa
for i in $(cat $name);
do 
#find fast5
find /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/06microorganism/04fast5/${i}_${species}_f5Sg -name "*.fast5" > $nanom6apath/${i}_${species}_f5name.txt
echo $i
# extracting signals
python $Nanom6ASof/extract_raw_and_feature_fast.py --cpu=20 --fl=$nanom6apath/${i}_${species}_f5name.txt -o $nanom6apath/${i}_${species}_res --clip=10 
# predicting m6A site cds2transcritpt
python $Nanom6ASof/predict_sites.py --cpu 20 --model $Nanom6ASof/bin/model \
        -r $ref/cds_from_genomic_rename.fa \
        -g $ref/genomic.fa \
        -b $ref/cds2transcritp.txt \
        -i $nanom6apath/${i}_${species}_res \
        -o $nanom6apath/${i}_${species}_res_final
# Visualization of m6A sites 
# python $Nanom6ASof/nanoplot.py --input $nanom6apath/${i}_${species}_res_final \
#    -o $nanom6apath/${i}_${species}_plot
done
done
```

###analysis
```{r}
## tandMod
micro_tande_res <- tande_process(tande_res_path = "06microorganism/06tandMod/")
micro_tande_res$modif <- str_split_fixed(micro_tande_res$mod_type,"_",2)[,1]
micro_tande_res$species <- str_split_fixed(micro_tande_res$mod_type,"_",2)[,2]
saveRDS(micro_tande_res,"06microorganism/00result/06tande_tab.rds")

micro_tande_res <- readRDS("06microorganism/00result/06tande_tab.rds")
micro_tande_res$IDs <- paste0(micro_tande_res$transcript_id,":",micro_tande_res$site,":+")
micro_tande_res <- micro_tande_res[sum>10]
micro_tande_m6a <- micro_tande_res[modif=="m6A"]
table(micro_tande_res$modif,micro_tande_res$sample,micro_tande_res$species)
# lapply(split.data.frame(micro_tande_res,micro_tande_res$species), function(spe_table){
#   dcast(spe_table,transcript_id+site  ~ sample,value.var = "ratio")
# })

## m6a-net
micro_m6anet_list <- lapply(list.files("/mnt/raid61/Personal_data/liudefu/bacteria_m6anet_inference",pattern = "csv",full.names = T),read.csv)
micro_m6anet_res <- rbindlist(micro_m6anet_list)
micro_m6anet_res$sample <- paste0("RTA_",str_split_fixed(micro_m6anet_res$Type,"_",3)[,2])
micro_m6anet_res$species <- str_split_fixed(micro_m6anet_res$Type,"_",3)[,3]
micro_m6anet_res$site <- as.numeric(str_split_fixed(micro_m6anet_res$IDs,":",3)[,2])
micro_m6anet_res$chr <- str_split_fixed(micro_m6anet_res$IDs,":",3)[,1]
micro_m6anet_res$site <- micro_m6anet_res$site +1
micro_m6anet_res$IDs <- paste0(micro_m6anet_res$chr,":",micro_m6anet_res$site,":+")
table(micro_m6anet_res$sample,micro_m6anet_res$species)

## nanom6A
sb_path <- list.files("/mnt/raid61/Personal_data/hezhifeng/nanom6A/result_bact",pattern = "res_final",full.names = T)
micro_m6a_res_list <- list()
for (i in 1:length(sb_path)) {
  micro_m6a_res_list[[i]] <- m6a_process(resFin_path =   sb_path[i])
}
saveRDS(micro_nanom6a_ratio,"06microorganism/00result/06nanom6a_tab.rds")

micro_nanom6a_ratio <- rbindlist(lapply(micro_m6a_res_list, function(x) { x$ratio}))
micro_nanom6a_abun <- rbindlist(lapply(micro_m6a_res_list, function(x) { x$geno_abun}))
micro_nanom6a_ratio <- micro_nanom6a_ratio[readN>10]
micro_nanom6a_ratio$sample <- gsub("_res_final","",micro_nanom6a_ratio$sample)
micro_nanom6a_ratio$samples <- paste0("RTA_",str_split_fixed(micro_nanom6a_ratio$sample,"_",3)[,2])
micro_nanom6a_ratio$species <- str_split_fixed(micro_nanom6a_ratio$sample,"_",3)[,3]
micro_nanom6a_ratio$IDs <- paste0(micro_nanom6a_ratio$chr,":",micro_nanom6a_ratio$pos,":+")
table(micro_nanom6a_ratio$samples,micro_nanom6a_ratio$species)
```


```{r fig.width=12}
# 5细菌软件 overlap
spec_sub_list <- list(micro_tande_m6a,micro_m6anet_res,micro_nanom6a_ratio)
names(spec_sub_list) <- c("tande","m6anet","nanom6a")
spec_sub <-rbindlist(lapply(c("tande","m6anet","nanom6a"),function(x){
  spec_sub <- spec_sub_list[[x]][,c("IDs","species")]
  spec_sub$soft <- x
  spec_sub
}))
spec_sub_ID <- lapply(unique(micro_tande_m6a$species), function(x){
  tmp <- split.data.frame(spec_sub[species==x],spec_sub[species==x]$soft)
  ggvenn(lapply(tmp,function(x) x$ID))+ggtitle(x)
})
ggarrange(plotlist = spec_sub_ID)

hist(micro_m6anet_res$mod_ratio) 
hist(micro_tande_m6a$ratio)
```
### Strep
#### 统计
```{r}
strep_tand <- micro_tande_m6a[mod_type=="m6A_Strep"]
plog <- ggseqlogo(strep_tand$motif,method = 'prob' )
ggsave("06microorganism/00result/Strep_m6a_seqlogo.pdf",plog,width = 4,height = 3)

strep_net <- micro_m6anet_res[species=="Strep"]
pven <- ggvenn(list(TandemMod=strep_tand$IDs,m6anet=strep_net$IDs),)
ggsave("06microorganism/00result/Strep_m6a_venn.pdf",pven,width = 4,height = 4)

library(eulerr)
## intersect m6A 统计每个样本数目
micor_int <- strep_tand[IDs %in% intersect(strep_tand$IDs,strep_net$IDs)]
length(unique(micor_int$IDs))
table(micor_int$sample)
p_int <- micor_int[,.N,by=sample][order(N,decreasing = T)]
p_int$sample <- factor(p_int$sample,levels =p_int$sample )
p_stre_m6a <- ggplot(p_int,aes(sample,N,fill=sample))+geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")+ 
  xlab(NULL)+ylab("Streptococcus m6A frequency")+
  mytheme
ggsave("06microorganism/00result/Strep_m6a_sam.pdf",p_stre_m6a,width = 6,height = 4)
```

####  diff m6a
```{r}
ggvenn::ggvenn(list(tande=micro_tande_m6a[species=="Strep"]$IDs,net=micro_m6anet_res[species=="Strep"]$IDs))
Strep_inter <- intersect(micro_tande_m6a[species=="Strep"]$IDs,micro_m6anet_res[species=="Strep"]$IDs)

#mod N filter
m6a_isoTab <- micro_tande_m6a[species=="Strep"]
m6A_matN <- reshape2::dcast(m6a_isoTab,IDs ~ sample, value.var = "mod")
m6A_matN[is.na(m6A_matN)]=0
m6A_matN <- data.frame(m6A_matN) %>% tibble::column_to_rownames("IDs") 
#filter reads>10 samples>3
dim(m6A_matN)
m6A_matN <- m6A_matN[apply(m6A_matN, 1, function(x){
  sum(x>0) >= 3 & sum(x)>=10
}),]
dim(m6A_matN)
m6A_matN <- m6A_matN[,colSums(m6A_matN)>100]
hv <- data.table(ID = rownames(m6A_matN), sd = apply(m6A_matN, 1, sd, na.rm = T))
hist(hv$sd)
m6A_matN <- m6A_matN[hv[!is.na(sd) & sd< 300 ][,ID], ]
dim(m6A_matN) 
# ggvenn(list(inter=Strep_inter,filter=rownames(m6A_matN)))

#ratio_mat
Strep_m6a <- dcast(micro_tande_m6a[species=="Strep"],IDs ~ sample, value.var = "ratio")
Strep_m6a[is.na(Strep_m6a)] <- 0
Strep_m6a <- Strep_m6a %>% tibble::column_to_rownames("IDs")
Strep_m6a <- Strep_m6a[rownames(m6A_matN),colnames(m6A_matN)]
```


```{r}
## ratio
#PCA
metap <- sample_de[colnames(Strep_m6a),]
metap$sample_type <- factor(metap$sample_type)
table(metap$sample_type)
DEGreport::degPCA(Strep_m6a,metadata=metap,condition="sample_type") + coord_fixed() +theme_bw()+
  scale_fill_manual(labels =c(nasopharyngeal="#00BFC4",Oropharyngeal="#F8766D"))+
                ggrepel::geom_text_repel(aes(label = colnames(Strep_m6a))) #+theme(legend.position="none")

# diff
fm_m6a <- FindMarkers_PII(
  PIIratio = Strep_m6a,
  meta = metap,
  design = "sample_type",
  ident.1 = "nasopharyngeal",
  ident.2 = "Oropharyngeal",
  min.pct = 0,
  delta.threshold = 0,
  only.pos = F,
  min.cells.group=0,
  test.use = "prop",
  RawTab =  m6a_isoTab[,c("sample","IDs","mod","sum")]
)
hist(fm_m6a$p_val)
fm_m6a_sub<-subset(fm_m6a,p_val < 0.05)
fm_m6a_sub$type <- ifelse(fm_m6a_sub$deltaPSI<0,"Oropharyngeal","nasopharyngeal")
table(fm_m6a_sub$type)

anno_col <- data.frame(type=metap[,"sample_type"])
rownames(anno_col) <- rownames(metap)
pheatmap::pheatmap(Strep_m6a[rownames(fm_m6a_sub),],cluster_cols = T,cluster_rows = T,annotation_col = anno_col)
```

```{r}
## m6A counts
metap <- sample_de[colnames(m6A_matN),]
metap$sample_type <- factor(metap$sample_type)
dds <- DESeqDataSetFromMatrix(countData = m6A_matN,
                              colData = metap,
                              design = ~ sample_type )
dds <- DESeq(dds)
ddsv <- varianceStabilizingTransformation(dds, blind = FALSE)
res <- results(dds)
res_sub <- res[res$pvalue<0.05,]

#PCA
table(metap$sample_type)
DEGreport::degPCA(assay(ddsv),metadata=metap,condition="sample_type") + coord_fixed() +theme_bw()+
  scale_fill_manual(labels =c(nasopharyngeal="#00BFC4",Oropharyngeal="#F8766D"))+
                ggrepel::geom_text_repel(aes(label = colnames(m6A_matN))) #+theme(legend.position="none")
anno_col <- data.frame(type=metap[,"sample_type"])
rownames(anno_col) <- rownames(metap)
pheatmap::pheatmap(assay(ddsv)[rownames(res_sub),],cluster_cols = T,cluster_rows = T,annotation_col = anno_col)
```



#07humanCoV
```{linux}
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/07humanCoV
for i in $(ls *fa | awk -F '.fa' '{print $1}')
do
minimap2 -ax map-ont -t 30  ${i}.fa \
       /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata/230808_24sam_comb_filter.fq |  samtools sort -@ 30 | samtools view -b > ${i}.bam
samtools index ${i}.bam
done
```

```{r}
CovPath <- list.files("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/07humanCoV/",pattern = "bam$",full.names = T)

CovBam <- lapply(CovPath, function(x){
  align <-
  readGAlignments(x,
    param = ScanBamParam(
      flag = scanBamFlag(
        isSupplementaryAlignment = FALSE,
        isPaired = FALSE,
        isSecondaryAlignment = FALSE
      ),
      mapqFilter = 0
    ),
    use.names = T
  )
  align
})
names(CovBam) <- gsub(".bam","",basename(CovPath))

CovBam_sub <-CovBam[c(1,6,7)] 
ggvenn(lapply(CovBam_sub, function(c) names(c)))
```
