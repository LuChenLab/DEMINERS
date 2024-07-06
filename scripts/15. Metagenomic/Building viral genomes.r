---
title: "Building viral genomes"
author: "junwei"
date: "2023-08-01"
output: html_document
---


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, 
                      warning = FALSE, 
                      message = FALSE, 
                      fig.align = 'center')
knitr::opts_knit$set(root.dir = "/mnt/raid61/Personal_data/songjunwei/DRS_RTA")

library(data.table)
library(ggplot2)
library(parallel)
library(cowplot)
library(patchwork)
library(Biostrings)
library(GenomicRanges)
library(GenomicAlignments)
library(ggpubr)
source("/mnt/raid61/Personal_data/songjunwei/AllScript/00Function/DIyplot.R")
```

#1 Length
2021-08-11
```{r}
SVV_fa <- readDNAStringSet("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/SVV.fasta")
PRRSV_fa <- readDNAStringSet("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/PRRSV.fasta")
Ecoli_fa <- readDNAStringSet("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/Ecoli/GCF_000005845.2_ASM584v2_genomic.fa")
S_enter_fa <- readDNAStringSet("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/S_enter/GCF_000006945.2_ASM694v2_genomic.fa")
S_cere_fa <- readDNAStringSet("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/S_cere/GCF_000146045.2_R64_genomic.fa")
```


```{r}
SVV_gr <- GenomicRanges::GRanges(seqnames = mapply(function(x) x[1], strsplit(names(SVV_fa), " ")), 
                                 ranges = IRanges(start = 1, width = width(SVV_fa)))

PRRSV_gr <- GenomicRanges::GRanges(seqnames = mapply(function(x) x[1], strsplit(names(PRRSV_fa), " ")), 
                                   ranges = IRanges(start = 1, width = width(PRRSV_fa)))

Ecoli_gr <- GenomicRanges::GRanges(seqnames = mapply(function(x) x[1], strsplit(names(Ecoli_fa), " ")), 
                                   ranges = IRanges(start = 1, width = width(Ecoli_fa)))
```



```{r}
param <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isSupplementaryAlignment = FALSE,
                                                               isPaired = FALSE,
                                                               isSecondaryAlignment = FALSE),
                                 mapqFilter = 0)
file1 <- list.files("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align", "SVV.bam$", full.names = TRUE)
file2 <- list.files("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align", "PRRSV.bam$", full.names = TRUE)
files <- c(file1,file2)

Mat0 <- lapply(files, function(x) {
  map0 <- readGAlignments(x, param = param)
  data.table(Batch = gsub(".bam", "", basename(x)), Length = qwidth(map0), Width = width(map0))
})
Mat0 <- do.call(rbind, Mat0)
table(Mat0$Batch)
Mat0$Batch <- factor(Mat0$Batch,levels = c("batch1_SVV","batch2_SVV","batch2_PRRSV"))
```


```{r fig.width=5, fig.height=4}
plength <- ggplot(Mat0, aes(x = Batch, y = Width)) +
  geom_violin() +
  geom_hline(yintercept = width(SVV_fa),
             lty = 2,
             color = "grey50") +
  geom_hline(yintercept = width(PRRSV_fa),
             lty = 2,
             color = "grey50") +
  stat_summary(fun.data = "mean_sd") +
  theme_bw(base_size = 15) +
  labs(y = "Read length (nt)") +
  scale_y_continuous(breaks = c(0, 2000, 4000, 6000, width(SVV_fa), width(PRRSV_fa)))+
  theme( axis.text.x = element_text(size=12,angle=45,hjust=1))

pdf("./SVV_PRRSV_length.pdf",5,4)
plength
dev.off()
```

```{r}
lapply(split(Mat0,f = Mat0$Batch)[1:2], function(tmp){
  tmp1 <- tmp[,.N,Width>=7000]
  tmp1$Percentage <- tmp1$N/sum(tmp1$N)*100
  print(max(tmp$Width)-width(SVV_fa))
  tmp1
})
lapply(split(Mat0,f = Mat0$Batch)[3], function(tmp){
  tmp1 <- tmp[,.N,Width>=15000]
  tmp1$Percentage <- tmp1$N/sum(tmp1$N)*100
  print(max(tmp$Width)-width(PRRSV_fa))
  tmp1
})

```

#2 DELs
```{r}
Map0 <- lapply(files, function(x) readGAlignments(x, param = param))
Map0 <- do.call(c, Map0)

DELs <- unlist(GenomicAlignments::cigarRangesAlongReferenceSpace(cigar = cigar(Map0), ops = "D", pos = start(Map0)))
DELs <- DELs[width(DELs) >= 3]
DELs <- data.table(pos = seq_along(as.numeric(coverage(DELs))), DEL = as.numeric(coverage(DELs)))
```

```{r}
cov <- coverage(Map0)[[1]]
cov <- data.table(pos = seq_along(as.numeric(cov)), D = as.numeric(cov))

DELs <- merge(DELs, cov, by = "pos")
DELss <- reduce(DELs[DEL / D > 0.3, IRanges(pos, pos)])
```


```{r}
Ms <- unlist(GenomicAlignments::cigarRangesAlongReferenceSpace(cigar = cigar(Map0), ops = "M", pos = start(Map0)))
Ms <- data.table(pos = seq_along(as.numeric(coverage(Ms))), Ms = as.numeric(coverage(Ms)))
```

```{r}
i <- 1
ggplot(Ms[pos >= start(DELss[i]) - 20 & pos <= end(DELss[i]) + 20]) +
  geom_step(aes(x = pos, y = Ms)) +
  theme_light(base_size = 15) + 
  scale_y_continuous(limits = c(0, Ms[pos >= start(DELss[i]) - 20 & pos <= end(DELss[i]) + 20, max(Ms)] + 100))
```

#3 SNVs
```{liunx}
conda activate redtools2 #2048
reditools="/mnt/raid61/Personal_data/songjunwei/software/reditools2.0/"
SOURCE_BAM_PATH="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/batch2_PRRSV.bam" #must have index
REFERENCE="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/PRRSV.fasta"
SIZE_FILE="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/PRRSV.fasta.fai"
OUTPUT_PATH="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/SNPs"

i="batch2_PRRSV"
mkdir $OUTPUT_PATH/${i}
bash $reditools/extract_coverage.sh $SOURCE_BAM_PATH ${OUTPUT_PATH}/${i}/coverage/ $SIZE_FILE
mpirun -np 20   $reditools/src/cineca/parallel_reditools.py -S \
-f $SOURCE_BAM_PATH \
-o ${OUTPUT_PATH}/${i}/parallel_table.txt.gz -r $REFERENCE \
-t $OUTPUT_PATH/${i}/tmp -Z $SIZE_FILE \
-G ${OUTPUT_PATH}/${i}/coverage/${i}.cov -D ${OUTPUT_PATH}/${i}/coverage/ \
--min-read-length 10 \
-q 10 --min-base-quality 10 \
--min-base-position 5 \
--max-base-position 5 \
--min-edits-per-nucleotide 1 -V

```


# PRRSV (DecodeR vs Sanger)
##call snps
```{linux}
#Sanger seq res
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/PRRSV_sanger
minimap2 -t 4 -ax splice -k14 --secondary=no --eqx /mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/PRRSV.fasta PRRSV_Sanger_res.fa \
 | samtools view -Sb | samtools sort  -o PRRSV_Sanger.bam
samtools index PRRSV_Sanger.bam
bcftools mpileup -f /mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/PRRSV.fasta PRRSV_Sanger.bam | bcftools call -mv -o PRRSV_Sanger_SNP.vcf

#DecodeR
bcftools mpileup -f /mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/PRRSV.fasta /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/batch2_PRRSV.bam | bcftools call -mv -o PRRSV_DecodeR_SNP.vcf 
bcftools mpileup -f  /mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/PRRSV.fasta /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/batch2_PRRSV.bam  -Q 0 | bcftools call -mv -o PRRSV_DecodeR_SNP.vcf11

#
# sambamba depth base -F '' -L 10:100021916-100021916 test1.bam
# samtools mpileup -go samtools.bcf -f /mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/PRRSV.fasta /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/batch2_PRRSV.bam
# bcftools call -vmO z -o bcftools_raw.vcf.gz samtools.bcf
```

##snps valid
```{r}
DecodeR_snp <- read.vcfR("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/batch123/PRRSV_sanger/PRRSV_DecodeR_SNP.vcf")
Sanger_snp <- read.vcfR("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/batch123/PRRSV_sanger/PRRSV_Sanger_SNP.vcf")
DecodeR_snp1 <- read.vcfR("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/batch123/PRRSV_sanger/PRRSV_DecodeR_SNP.vcf1")
```

```{r}
#Sanger
SanSnp <- data.table(data.frame(Sanger_snp@fix))
SanSnp$SNP <- paste0(SanSnp$POS,":",SanSnp$REF,">",SanSnp$ALT)

#Dec0
DecSnp <- data.table(data.frame(DecodeR_snp@fix))
DecSnp$POS <- as.numeric(DecSnp$POS)
DecSnp$QUAL <- as.numeric(DecSnp$QUAL)
DecSnp_sub <- DecSnp[POS>13941][POS<15340]
DecSnp_sub$SNP <- paste0(DecSnp_sub$POS,":",DecSnp_sub$REF,">",DecSnp_sub$ALT)
DecSnp_sub$vail <- ifelse(DecSnp_sub$SNP %in% SanSnp$SNP,"Sanger_vail","Sanger_novail")
DecSnp_sub[vail=="Sanger_novail"]

#Dec1
DecSnp1 <- data.table(data.frame(DecodeR_snp1@fix))
DecSnp1$POS <- as.numeric(DecSnp1$POS)
DecSnp1$QUAL <- as.numeric(DecSnp1$QUAL)
DecSnp_sub1 <- DecSnp1[POS>13941][POS<15340]
DecSnp_sub1$SNP <- paste0(DecSnp_sub1$POS,":",DecSnp_sub1$REF,">",DecSnp_sub1$ALT)
DecSnp_sub1$vail <- ifelse(DecSnp_sub1$SNP %in% SanSnp$SNP,"Sanger_vail","Sanger_novail")
DecSnp_sub1[vail=="Sanger_novail"]
table(DecSnp_sub1$vail)

ggvenn(list(DecodeR=DecSnp_sub$SNP,DecodeR1=DecSnp_sub1$SNP,Sanger=SanSnp$SNP))+ggtitle("SNPs Validated")
```

```{r fig.width=6,fig.height=4}
DecSnp_sub1$vail <- factor(DecSnp_sub1$vail,levels = c("Sanger_vail","Sanger_novail"))
pvio <- ggplot(DecSnp_sub1, aes(vail, QUAL)) +
  geom_violin(aes(fill = vail), cex = 1.2) +
  #scale_fill_manual(values = c('#FB5554','#868B31','#42F203','#579ABB','#B978AE'))+
  geom_boxplot(width = 0.1, cex = 1.2) +
  theme_classic(base_size = 14)  +geom_signif(comparisons =  list(c("Sanger_vail","Sanger_novail")))+
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'none') + xlab(NULL) + ylab("Quality of variant")
ggsave(filename = "/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/PRRSV_sanger/DecodeRvail_vio.pdf",plot = pvio,width = 6,height = 4)
```

```{r}
SanSnp$vail <- ifelse(SanSnp$SNP %in% DecSnp_sub$SNP,"DRSvail","DRS_novail")
SanSnp[vail=="DRS_novail"]

data <- DecSnp_sub1[,.N,by="vail"]
data$precent = round(data$N/sum(data$N),2)*100
pie1 <- ggplot(data, aes(x="", y=precent, fill=vail)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  geom_text(aes(label = paste0(precent, "%")), position = position_stack(vjust=0.5)) +
  labs(x = NULL, y = NULL, fill = NULL)+ theme_bw() + theme(legend.text = element_text(size = 14),legend.position = 'none')+
  ggtitle(label = "PRRSV SNPs validated")

ggsave(filename = "/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/PRRSV_sanger/DecodeRvail_pie.pdf",plot = pie1,width = 3,height = 3)
```

## 测序深度
```{linux}
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/PRRSV_sanger
samtools depth /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/batch2_PRRSV.bam > batch2_PRRSV.depth
```

```{r}
PRRSV.depth <- read.table("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/PRRSV_sanger/batch2_PRRSV.depth")
colnames(PRRSV.depth) <- c("chr", "position", "depth")
mean(PRRSV.depth$depth)

# 选择区域 13941:15340
selected_region <- subset(PRRSV.depth, position >= 13941 & position <= 15340)
mean(selected_region$depth)

p1 <- ggplot(PRRSV.depth, aes(x=depth)) +
  geom_density(fill="blue", alpha=0.5) + 
  labs(subtitle ="Full Genomic Sequencing Depth Density",
       x="Sequencing Depth",
       y="Density") + theme_minimal()
p2 <- ggplot(selected_region, aes(x=depth)) +
  geom_density(fill="blue", alpha=0.5) + 
  labs(subtitle="Subset Genomic Sequencing Depth Density",
       x="Sequencing Depth",
       y="Density") + theme_minimal()
pden <- p1+p2
pden
ggsave("PRRSV_sanger/Seqencing_Depth_Density.pdf",plot = pden,width = 8,height = 4)

p3 <- ggplot(PRRSV.depth, aes(x=position, y=depth)) +
  geom_bar(stat="identity", width=1, fill="blue", alpha=0.5) +
  labs(subtitle ="Full Genomic Sequencing Depth per Position",
       x="Genomic Position",
       y="Sequencing Depth") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p4 <- ggplot(selected_region, aes(x=position, y=depth)) +
  geom_bar(stat="identity", width=1, fill="blue", alpha=0.5) +
  labs(subtitle="Subset Genomic Sequencing Depth per Position",
       x="Genomic Position",
       y="Sequencing Depth") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3+p4
pbar <- p3+p4
ggsave("PRRSV_sanger/Seqencing_Depth_per_Position.pdf",plot = pbar,width = 8,height = 4)
```



