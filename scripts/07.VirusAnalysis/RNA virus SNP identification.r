---
title: "RNA virus SNP identification"
author: "junwei"
date: "2023/6/30"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(vcfR)
library(data.table)
library(ggvenn)
```

#call snps
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

#snps valid
```{r}
DecodeR_snp <- read.vcfR("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/PRRSV_sanger/PRRSV_DecodeR_SNP.vcf")
Sanger_snp <- read.vcfR("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/PRRSV_sanger/PRRSV_Sanger_SNP.vcf")
DecodeR_snp1 <- read.vcfR("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/00analysis/PRRSV_sanger/PRRSV_DecodeR_SNP.vcf1")
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

