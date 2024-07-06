---
title: "Untitled"
author: "junwei"
date: "2024/1/18"
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
library(ggpubr)
library(ggrepel)
library(vcfR)
library(dplyr)
library(ggseqlogo)
library(ggplot2)
library(reshape2)
library(readxl)
library(tidyr)
library(forcats)

setwd("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/")

source("/mnt/raid61/Personal_data/songjunwei/AllScript/00Function/DIyplot.R")

accuracy_to_df <- function(file,epoch){
  content <- readLines(file)
  # 提取 Total, Median accuracy, Average accuracy, std
  total <- as.numeric(gsub("Mapping rate: ([0-9]+).*", "\\1", content[1]))
  median_accuracy <- as.numeric(gsub(".*Median accuracy: ([0-9.]+).*", "\\1", content[1]))
  average_accuracy <- as.numeric(gsub(".*Average accuracy: ([0-9.]+).*", "\\1", content[1]))
  std <- as.numeric(gsub(".*std: ([0-9.]+)", "\\1", content[1]))
  
  # 提取 Mismatch, Deletions, Insertions 的 Median 和 Average 值
  mismatch_median <- as.numeric(gsub(".*Mismatch: ([0-9.]+).*", "\\1", content[2]))
  deletions_median <- as.numeric(gsub(".*Deletions: ([0-9.]+).*", "\\1", content[2]))
  insertions_median <- as.numeric(gsub(".*Insertions: ([0-9.]+).*", "\\1", content[2]))
  length_median <- as.numeric(gsub(".*tlen: ([0-9.]+)", "\\1", content[2]))
  
  mismatch_average <- as.numeric(gsub(".*Mismatch: ([0-9.]+).*", "\\1", content[3]))
  deletions_average <- as.numeric(gsub(".*Deletions: ([0-9.]+).*", "\\1", content[3]))
  insertions_average <- as.numeric(gsub(".*Insertions: ([0-9.]+).*", "\\1", content[3]))
  length_average <- as.numeric(gsub(".*tlen: ([0-9.]+)", "\\1", content[3]))
  
   res <- data.frame(
        Epoch = epoch,
        Mapp_reads = total,
        Accuracy_Med = median_accuracy,
        Accuracy_Ave = average_accuracy,
        Std = std,
        Mismatch_Med = mismatch_median,
        Deletions_Med = deletions_median,
        Insertions_Med = insertions_median,
        Length_Med = length_median,
        Mismatch_Ave = mismatch_average,
        Deletions_Ave = deletions_average,
        Insertions_Ave = insertions_average,
        Length_Ave = length_average
      )
   return(res)
}
```

# fastq
```{r}
"/mnt/raid61/Personal_data/songjunwei/DRS_RTA/Rodan/outfile.fasta" 
# 331,500

library(ShortRead)
cbf <- readFasta("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata/230808_24sam_comb_filter.fa")
cbffq <- readFastq("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata/230808_24sam_comb.fq.gz")

#230808_24sam_comb.fq.gz #295747
#230808_24sam_comb_filter.fa #374204
```

#1 Rodan
```{linux}
## basecalling
rodan="/mnt/raid61/Personal_data/songjunwei/software/RODAN/"
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/
python $rodan/basecall.py -m $rodan/rna.torch SARS_fast5 > SARS_rodon.fa

```

```{linux}
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata
minimap2 -ax map-ont -t 30  /mnt/raid61/Personal_data/songjunwei/reference/SARS_CoV2/nCoV-2019.reference.fasta \
       230808_24sam_comb_filter.fq |  samtools sort -@ 30 | samtools view -b > 02align_filter/230808_24sam_comb_filter_sars.bam


cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/Rodan/
seqkit grep -f RTA_09_readid_filter.txt outfile.fasta > RTA09_Rodan_ft.fq
minimap2 -ax map-ont -t 30  /mnt/raid61/Personal_data/songjunwei/reference/SARS_CoV2/nCoV-2019.reference.fasta \
       RTA09_Rodan_ft.fq | samtools sort -@ 30 | samtools view -b > RTA09_Rodan_sars.bam
bcftools mpileup -f /mnt/raid61/Personal_data/songjunwei/reference/SARS_CoV2/nCoV-2019.reference.fasta \
   RTA09_Rodan_sars.bam  -Q 0 | bcftools call -mv -o RTA09_Rodan_sars.vcf
```


##QC
```{r }
Nano_res <- fread("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/Rodan/QC/bam_NanoComp-data.tsv.gz",header = T) 
Nano_res$dataset <- basename(Nano_res$dataset)
table(Nano_res$dataset)
#reads number
Nano_res$dataset <- gsub("RTA_09_sars.bam","Guppy",Nano_res$dataset)
Nano_res$dataset <- gsub("RTA09_Rodan_sars.bam","Rodan",Nano_res$dataset)
table(Nano_res$length_filter)
Reads_number <- Nano_res[,.N,by=dataset]
Reads_number$dataset <- factor(Reads_number$dataset,levels =Reads_number[order(N,decreasing = T)]$dataset )
ggplot(Reads_number,aes(x=dataset,y=N,fill=dataset))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")+
  geom_bar(stat = 'identity',position = 'stack')+ mytheme+xlab("Samples")+
   geom_text(aes(label = N),
              position = position_dodge(width = 0.9),
              vjust = -0.5,size = 4) + 
  ylab("Reads number")

#quality
ggplot(data=Nano_res,mapping = aes(x=dataset,y=mapQ,fill=dataset))+
       geom_violin(trim = FALSE,color="white") +
       geom_boxplot(outlier.colour = NA,position = position_dodge(0.9))+
      xlab("Samples")+ ylab("map Quality")+ mytheme +
      theme(
            panel.grid.minor = element_blank(),
            text = element_text(size = 18),
            axis.text.x= element_text(size=12,family="ArialMT",angle =  45,hjust = 1),
            legend.position = "none")

#length
mean(Nano_res$lengths)
min(Nano_res$lengths)
max(Nano_res$lengths)
ggplot(data=Nano_res,mapping = aes(x=lengths)) + 
  geom_density(aes(fill = dataset),alpha=0.3)+ scale_x_log10()+
  geom_vline(data = Nano_res, aes(xintercept = mean(lengths), color=dataset),linetype="dashed")+
  mytheme+ 
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank(),
        text = element_text(size = 12),
        axis.text.x= element_text(size=12,family="ArialMT"))+
    xlab("Length")+ ylab("Denisty")

ggplot(data=Nano_res,mapping = aes(x=dataset,y=lengths,fill=dataset))+
    #   geom_violin(trim = FALSE,color="white") +
       geom_boxplot(outlier.colour = NA,position = position_dodge(0.9))+
      xlab("Samples")+ ylab("Lengths")+ mytheme +
      theme(
            panel.grid.minor = element_blank(),
            text = element_text(size = 18),
            axis.text.x= element_text(size=12,family="ArialMT",angle =  45,hjust = 1),
            legend.position = "none")
sta <- cbind(Nano_res[,mean(lengths),by=dataset],Nano_res[,median(lengths),by=dataset])
colnames(sta) <- c("dataset","mean","dataset1","median")
sta$dataset1 <- NULL
```

##VCF
```{r}
vcf_files <- c("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/02split_SARS_align/RTA_09_sars.vcf",
"/mnt/raid61/Personal_data/songjunwei/DRS_RTA/Rodan/RTA09_Rodan_sars.vcf",
"/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/03NGS_res/20230421-247_sars.vcf")
file.exists(vcf_files)

snp_tab <- rbindlist(lapply(vcf_files, function(x){
  tmp1 <- read.vcfR(x)
  tmp2 <- data.table(data.frame(tmp1@fix))
  tmp2$sample <- gsub("_sars.vcf","",basename(x))
  tmp2$POS <- as.numeric(tmp2$POS)
tmp2$QUAL <- as.numeric(tmp2$QUAL)
  return(tmp2)
}))
table(snp_tab$sample)
snp_tab$ID <- paste0(snp_tab$CHROM,":",snp_tab$POS,":",snp_tab$REF,">",snp_tab$ALT)


snp_list <- split(snp_tab[QUAL>20],f = snp_tab[QUAL>20]$sample)
hist(snp_list$`20230421-247`$QUAL)
hist(snp_list$RTA_09$QUAL)
hist(snp_list$RTA09_Rodan$QUAL)

ggvenn(list(NGS=snp_list$`20230421-247`$ID,
            Guppy= snp_list$RTA_09$ID,
            Rodan= snp_list$RTA09_Rodan$ID))

snp_list$RTA_09$inNGS <- (snp_list$RTA_09$ID %in% snp_list$`20230421-247`$ID)
snp_list$RTA09_Rodan$inNGS <- (snp_list$RTA09_Rodan$ID %in% snp_list$`20230421-247`$ID)

data <- snp_list$RTA09_Rodan[,.N,by=inNGS]
data$precent <- round(data$N/sum(data$N)*100,2)
ggplot(data, aes(x="", y=precent, fill=inNGS)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  geom_text(aes(label = paste0(precent, "%")), position = position_stack(vjust=0.5)) +
  labs(x = NULL, y = NULL, fill = NULL)+ theme_bw() + theme(legend.text = element_text(size = 14),axis.text.y = element_blank(),axis.ticks=element_blank())
```

```{r}
intersect(intersect(snp_list$`20230421-247`$ID,
            snp_list$RTA_09$ID),
            snp_list$RTA09_Rodan$ID)

setdiff(intersect(snp_list$RTA09_Rodan$ID,
            snp_list$RTA_09$ID),snp_list$`20230421-247`$ID)

setdiff(intersect(snp_list$RTA09_Rodan$ID,
            snp_list$`20230421-247`$ID),snp_list$RTA_09$ID)
```

# 2 re-training
## sub sars f5
```{r}
sars_24sam_sub <- readRDS("01filter_rawdata/Combin_QC_align_filter.rds")
table(sars_24sam_sub$align)

sars_24sam_sub_sars <- sars_24sam_sub[align=="SARS"]


write.table(sars_24sam_sub_sars$readid,"08SARS_rodan/SARS_readid.txt",row.names = F,col.names = F,quote = F)
write.table(sars_24sam_sub_sars$readid[1:320],"08SARS_rodan/0Tran_model/Train_SARS_readid.txt",row.names = F,col.names = F,quote = F)
write.table(sars_24sam_sub_sars$readid[321:406],"08SARS_rodan/0Tran_model/Valid_SARS_readid.txt",row.names = F,col.names = F,quote = F)


# tsv <- fread("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/tayaki_res/SARS_readid.tsv")
# tsv_train <- tsv[UUID %in% sars_24sam_sub_sars$readid[1:320]]
# fwrite(tsv_train,"/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/0Tran_model/Train_SARS_readid.tsv",row.names = F,col.names = T,sep = "\t")
```

## tayaki
```{linux}
## fast5/fq process
conda activate base #2101

cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/0Tran_model/
sars_248_ref="/mnt/raid64/COVID19/Personal/lhkp/HX_NGS/20230426/04results/04fasta/20230421-248.fa" #RTA09 fa

for i in Train  Valid; 
do
fast5_subset -i /mnt/raid61/Personal_data/songjunwei/DRS_RTA/01Rawdata/batch5_SARS/230808_24sam \
            -s ${i}_SARS_fast5 \
            --batch_size  1000 \
            --read_id_list ${i}_SARS_readid.txt \
            --recursive  --ignore_symlinks \
            -t 10
seqkit grep -f ${i}_SARS_readid.txt /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/01filter_rawdata/230808_24sam_comb_filter.fq > ${i}_SARS_guppy.fq
minimap2 -ax map-ont -t 30  $sars_248_ref ${i}_SARS_guppy.fq |  samtools sort -@ 30 | samtools view -b > ${i}_SARS_guppy_248.bam
samtools index ${i}_SARS_guppy_248.bam
done

#### tayaki ####
conda activate base #2048
## install 
# /mnt/raid5/Personal/minion/miniconda3/bin/pip install -r requirements.txt
# python setup.py install

cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/0Tran_model/
tayaki="/mnt/raid5/Personal/minion/software/taiyaki/bin/"
sars_248_ref="/mnt/raid64/COVID19/Personal/lhkp/HX_NGS/20230426/04results/04fasta/20230421-248.fa" #RTA09 fa

for i in Train Valid ;
do
python $tayaki/generate_per_read_params.py ${i}_SARS_fast5 --output ${i}_SARS_readid.tsv
python $tayaki/get_refs_from_sam.py $sars_248_ref  ${i}_SARS_guppy_248.bam  --output ${i}_SARS_tayaki_248.fa
python $tayaki/prepare_mapped_reads.py ${i}_SARS_fast5 ${i}_SARS_readid.tsv /mnt/raid5/Personal/minion/Tmp/model_tayaki/${i}_tayaki_248.h5 \
        /mnt/raid5/Personal/minion/software/taiyaki/models/r941_rna_minion_upgrae.checkpoint ${i}_SARS_tayaki_248.fa --jobs 30 --overwrite 
cp /mnt/raid5/Personal/minion/Tmp/model_tayaki/${i}_tayaki_248.h5 ${i}_tayaki_248.h5
done

#### rodan ####
rodan="/mnt/raid61/Personal_data/songjunwei/software/RODAN/"
cd /mnt/raid5/Personal/minion/Tmp/model_tayaki/

## training
for i in Train Valid ;
do
mkdir rodan_${i}
mkdir rodan_${i}
python  $rodan/gendata.py -i ${i}_tayaki_248.h5 --outdir rodan_${i}
mv rodan_${i}/train.hdf5 rna-${i}.hdf5
done

mkdir runs
python  $rodan/model.py -c ./rna.config -a $rodan/rnaarch -n model248 --savedir ./runs --rna -w 25 -l 
# Validation loss: 1.4489204569866783
# Train losses: [2.225710514344667, 1.480465634873039, 1.458074877136632, 1.4474307173176815, 1.4449025863095333, 1.441629670168224, 1.441551054778852, 1.439997058165701, 1.4385978799117238, 1.4390837361938076, 1.4378141086352498, 1.4366237317260944, 1.436683454011616, 1.4363966286182404, 1.4366595007871326, 1.43619903922081, 1.4361280108753003, 1.4358080609848625, 1.4359178950912075, 1.435471955098604, 1.4352239401716935, 1.435779598198439, 1.4356051087379456, 1.4352151014302905, 1.4348549544811249, 1.4353058683244806, 1.4351538372667212, 1.4354981184005737, 1.4349316154655658, 1.4352964602018659]
# Valid losses: [1.5135465609399896, 1.4560558858670687, 1.4521880087099577, 1.4629587123268528, 1.4442975583829378, 1.4833475288591886, 1.4533569122615613, 1.445829529511301, 1.4501944278415881, 1.4501856126283343, 1.448550349787662, 1.4484174188814665, 1.4472031342355829, 1.4506466576927586, 1.447545609976116, 1.4492699158819098, 1.450933105067203, 1.4499160239571018, 1.4519202960164923, 1.4511395128149736, 1.4495273765764738, 1.4498613821832758, 1.44954346355639, 1.4478764408513118, 1.4495381556059186, 1.4515444542232312, 1.4496706034007825, 1.4486065663789447, 1.4480378063101518, 1.4489204569866783]
# Learning rate: [0.002, 0.002, 0.002, 0.001, 0.001, 0.0005, 0.0005, 0.00025, 0.00025, 0.000125, 0.000125, 6.25e-05, 6.25e-05, 3.125e-05, 3.125e-05, 1.5625e-05, 1.5625e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05]
cp ./runs/model248-epoch29.torch $rodan/model248-epoch29.torch

## basecalling
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/
f5="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/SARS_fast5"
python $rodan/basecall.py -m $rodan/model248-epoch29.torch $f5 --arch $rodan/rnaarch > SARS_rodon_t1Mod.fa
python $rodan/basecall.py -m $rodan/rna.torch $f5 --arch $rodan/rnaarch > SARS_rodon_rawMod.fa

python $rodan/basecall.py -m "/mnt/raid5/Personal/minion/Tmp/model_tayaki/torch/model248-ext.torch" $f5  --arch $rodan/rnaarch > SARS_rodon_t1Mod.fa

mkdir 02accuracy
cd 02accuracy
cp /mnt/raid64/COVID19/Personal/lhkp/HX_NGS/20230426/04results/04fasta/20230421-248.fa ref.fa
sars_248_ref="/mnt/raid64/COVID19/Personal/lhkp/HX_NGS/20230426/04results/04fasta/20230421-248.fa" #RTA09 fa
for i in SARS_guppy.fq SARS_rodon_rawMod.fa ;
do
minimap2 --secondary=no -ax map-ont -t 32 --cs ref.fa /${i} > 01accuracy/samfiles/${i}.sam
done

cd 01accuracy
python $rodan/accuracy.py samfiles/SARS_guppy.fq.sam ref.fa > accufiles/accuracy_guppy.txt
# Total: 406 
# Median accuracy: 0.8825791668779474 
# Average accuracy: 0.8796905447563035 
# std: 0.04241556862567157
# Median Mismatch: 0.019852758842105783 
# Deletions: 0.07780597856739989 
# Insertions: 0.01830385062092379
# Average Mismatch: 0.02255497659763202 
# Deletions: 0.07712021645371726 
# Insertions: 0.020634262192347225

python $rodan/accuracy.py samfiles/SARS_rodon_rawMod.fa.sam ref.fa > accufiles/accuracy_rodanraw.txt
# Total: 393 Median accuracy: 0.8682926829268293 Average accuracy: 0.8653626097729714 std: 0.0401615971300731
# Median  - Mismatch: 0.02540220152413209 Deletions: 0.07191011235955057 Insertions: 0.029508196721311476
# Average - Mismatch: 0.02748112539877543 Deletions: 0.07317707102388787 Insertions: 0.03397919380436529
```

#3 Vero-Infected
```{linux}
Vero-Infected="/mnt/raid65/PublicData/Zenodo/Cell_2020_SARS-Vero/Vero-Infected"
cd /mnt/raid65/PublicData/Zenodo/Cell_2020_SARS-Vero/Vero-Infected/20200222_0650_MN26136_FAN04901_ada1e2bf/
zstd -d VeroInf24h.all.fastq.zst

sars_2019_ref="/mnt/raid61/Personal_data/songjunwei/reference/SARS_CoV2/nCoV-2019.reference.fasta"
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/02Tran_model_CellData
minimap2 -ax map-ont -t 30  $sars_2019_ref /mnt/raid65/PublicData/Zenodo/Cell_2020_SARS-Vero/Vero-Infected/20200222_0650_MN26136_FAN04901_ada1e2bf/VeroInf24h.all.fastq |  samtools sort -@ 30 | samtools view -b >  ./01rawdata/SARS_CellD.bam
```

```{r}
sars_align <-
  readGAlignments(
    "/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/02Tran_model_CellData/01rawdata/SARS_CellD.bam",
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

write.table(names(sars_align)[1:454660],"/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/02Tran_model_CellData/01rawdata/Train_SARS_readid.txt",row.names = F,col.names = F,quote = F)
write.table(names(sars_align)[454661:568327],"/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/02Tran_model_CellData/01rawdata/Valid_SARS_readid.txt",row.names = F,col.names = F,quote = F)
```


```{linux}
## fast5/fq process
conda activate base #2101
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/02Tran_model_CellData
cd /mnt/raid5/Personal/minion/DRS_mul/model_train_CellData #2048
sars_2019_ref="/mnt/raid61/Personal_data/songjunwei/reference/SARS_CoV2/nCoV-2019.reference.fasta"
for i in Train  Valid; 
do
fast5_subset -i /mnt/raid65/PublicData/Zenodo/Cell_2020_SARS-Vero/Vero-Infected/20200222_0650_MN26136_FAN04901_ada1e2bf/fast5 \
            -s 01rawdata/${i}_SARS_fast5 \
            --batch_size  2000 \
            --read_id_list 01rawdata/${i}_SARS_readid.txt \
            --recursive  --ignore_symlinks \
            -t 30
seqkit grep -f 01rawdata/${i}_SARS_readid.txt /mnt/raid65/PublicData/Zenodo/Cell_2020_SARS-Vero/Vero-Infected/20200222_0650_MN26136_FAN04901_ada1e2bf/VeroInf24h.all.fastq > 01rawdata/${i}_SARS_guppy.fq
minimap2 -ax map-ont -t 30  $sars_2019_ref 01rawdata/${i}_SARS_guppy.fq |  samtools sort -@ 30 | samtools view -b > 01rawdata/${i}_SARS_guppy.bam
samtools index 01rawdata/${i}_SARS_guppy.bam
done

#### tayaki ####
conda activate base #2048
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/02Tran_model_CellData/
tayaki="/mnt/raid5/Personal/minion/software/taiyaki/bin/"
sars_2019_ref="/mnt/raid61/Personal_data/songjunwei/reference/SARS_CoV2/nCoV-2019.reference.fasta"

for i in Train Valid ;
do
python $tayaki/generate_per_read_params.py 01rawdata/${i}_SARS_fast5 --output 02tayaki/${i}_SARS_readid_raw.tsv
done
#Train从中提取 20000条，Valid从中提取2000条
cat Train_SARS_readid_raw.tsv|head -n 20001 >Train_SARS_readid.tsv
cat Valid_SARS_readid_raw.tsv|head -n 2001 > Valid_SARS_readid.tsv

for i in Train Valid ;
do
python $tayaki/get_refs_from_sam.py $sars_2019_ref  01rawdata/${i}_SARS_guppy.bam  --output 02tayaki/${i}_SARS_tayaki.fa --reverse
python $tayaki/prepare_mapped_reads.py 01rawdata/${i}_SARS_fast5 02tayaki/${i}_SARS_readid.tsv /mnt/raid5/Personal/minion/DRS_mul/model_train_CellData/02tayaki/${i}_tayaki.h5 \
        /mnt/raid5/Personal/minion/software/taiyaki/models/r941_rna_minion_upgrae.checkpoint     02tayaki/${i}_SARS_tayaki.fa --jobs 30 --overwrite 
done

rodan="/mnt/raid61/Personal_data/songjunwei/software/RODAN/"
cd /mnt/raid5/Personal/minion/DRS_mul/model_train_CellData/
for i in Train Valid ;
do
mkdir 02tayaki/rodan_${i}
python  $rodan/gendata.py -i 02tayaki/${i}_tayaki.h5 --outdir 02tayaki/rodan_${i}
mv 02tayaki/rodan_${i}/train.hdf5 02tayaki/rna-${i}.hdf5
rm -rf 02tayaki/rodan_${i}
done

cd /mnt/raid5/Personal/minion/DRS_mul/model_train_CellData/02tayaki
python  $rodan/model.py -c ./rna.config -a $rodan/rnaarch -n M_CellD --arch $rodan/rnaarch --savedir ./runs --rna -w 10 -l

cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/
f5="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/SARS_fast5"
rodan="/mnt/raid61/Personal_data/songjunwei/software/RODAN/"
python $rodan/basecall.py -m $rodan/M_CellD-epoch7.torch $f5 --arch $rodan/rnaarch -b 200 > SARS_rodon_t2CellD_ep7.fa
cd 01accuracy
minimap2 --secondary=no -ax map-ont -t 30 --cs ref.fa ../SARS_rodon_t2CellD_ep7.fa > SARS_rodon_t2CellD_ep7.sam
python $rodan/accuracy.py SARS_rodon_t2CellD_ep7.sam ref.fa
```


```{linux}
f5="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/SARS_fast5"
rodan="/mnt/raid61/Personal_data/songjunwei/software/RODAN/"
accuracy_dir="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/01accuracy/"

# 遍历所有的 M_CellD-epoch**.torch 文件
for model in /mnt/raid5/Personal/minion/DRS_mul/model_train_CellData/02tayaki/runs/M_CellD-epoch*.torch; do
    epoch=$(basename "$model" | sed 's/M_CellD-epoch\([0-9]*\).torch/\1/')
    output_fa="SARS_rodon_t2CellD_ep${epoch}.fa"
    output_sam="SARS_rodon_t2CellD_ep${epoch}.sam"
    accuracy_result_file="${accuracy_dir}/accuracy_epoch${epoch}.txt"

    # 运行 basecalling
    python "$rodan/basecall.py" -m "$model" "$f5" --arch "$rodan/rnaarch" -b 10 > /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/$output_fa

    # 运行 accuracy 测试并将结果保存到文本文件
    cd "$accuracy_dir"
    minimap2 --secondary=no -ax map-ont -t 32 --cs ref.fa "../$output_fa" > "$output_sam"
    python "$rodan/accuracy.py" "$output_sam" ref.fa > "$accuracy_result_file"
done

#QC
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/
NanoComp --fasta  SARS_guppy.fq   SARS_rodon_rawMod.fa  SARS_rodon_t2CellD_ep29.fa --name guppy rodan_raw rodan_new -o ./QC -t 30  -p fa_ --raw
```


```{r}
files <- list.files(path = "/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/01accuracy/accufiles/",pattern = "accuracy_.*\\.txt$",full.names = T)

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
write.csv(results_df,"08SARS_rodan/01accuracy/Epochs_res.csv",row.names = T)
```

```{r}
results_com <- results_df[c("8","accuracy_guppy.txt","accuracy_rodanraw.txt"),]
results_com <- data.table(results_com[,-1])
results_com$type <- c("Rodan_v1.0_Train","Guppy_v6.5.7","Rodan_v1.0_Raw")


pdata_com <- melt(results_com[,-c("Total","Std")])
pdata_com$type <- factor(pdata_com$type,levels = c("Rodan_v1.0_Train","Guppy_v6.5.7","Rodan_v1.0_Raw"))
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
ggsave("08SARS_rodan/01accuracy/Comparison_res.pdf",width = 9,height = 5)
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
ggsave("08SARS_rodan/01accuracy/Epochs_total.pdf",p1,width = 4,height = 4)
ggsave("08SARS_rodan/01accuracy/Epochs_accuracy_res.pdf",pall,width = 8,height = 6)
```


#4 m6A basecaller
```{linux}
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/
#3102
/mnt/raid62/Personal_data/zengwanqin/soft/ont-guppy/bin/guppy_basecaller \
-i /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/SARS_fast5/ \
--save_path 03m6aCaller/ --device cuda:0 \
-c /mnt/raid61/Personal_data/songjunwei/software/m6ABasecaller-main/basecalling_model/rna_r9.4.1_70bps_m6A_hac.cfg \
--min_qscore 7  --verbose_logs  \
--chunks_per_runner 256  --chunk_size 100 \
--num_callers 4 --gpu_runners_per_device 8

#2048
guppy_basecaller \
-i /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/SARS_fast5/ \
--save_path 05.1SUP/ --device cuda:0 \
-c rna_r9.4.1_70bps_sup.cfg~ \
--min_qscore 7  --verbose_logs  \
--chunks_per_runner 256  --chunk_size 100 \
--num_callers 4 --gpu_runners_per_device 8



minimap2 --secondary=no -ax map-ont -t 32 --cs ./01accuracy/ref.fa SARS_m6aCaller5.0.11.fq > ./01accuracy/SARS_m6aCaller5.0.11.sam
minimap2 --secondary=no -ax map-ont -t 32 --cs ./01accuracy/ref.fa SARS_guppy5.0.11.fq > ./01accuracy/SARS_guppy5.0.11.sam
minimap2 --secondary=no -ax map-ont -t 32 --cs ./01accuracy/ref.fa SARS_guppy657_SUP.fq > ./01accuracy/SARS_guppy657_SUP.sam

cd 01accuracy
rodan="/mnt/raid61/Personal_data/songjunwei/software/RODAN/"

python $rodan/accuracy.py SARS_m6aCaller5.0.11.sam ref.fa
# Total: 106 Median accuracy: 0.7359335199984445 Average accuracy: 0.7405692505048377 std: 0.027978806695137696
# Median  - Mismatch: 0.07334969961637755 Deletions: 0.15218061120588194 Insertions: 0.0334473997436367
# Average - Mismatch: 0.07268077543658238 Deletions: 0.1531301812678763 Insertions: 0.033619792790703484

python $rodan/accuracy.py SARS_guppy5.0.11.sam ref.fa
# Total: 346 Median accuracy: 0.8513248458988953 Average accuracy: 0.8457969296303429 std: 0.03991511879004343
# Median  - Mismatch: 0.03536529564288983 Deletions: 0.07404749309452956 Insertions: 0.03915810187452881
# Average - Mismatch: 0.03844432533262349 Deletions: 0.07386676179267904 Insertions: 0.041891983244354364

python $rodan/accuracy.py SARS_guppy657_SUP.sam ref.fa
# Total: 300 Median accuracy: 0.8293283155416131 Average accuracy: 0.82833685507421 std: 0.03819554401566258
# Median  - Mismatch: 0.04406830299687443 Deletions: 0.08436127852562945 Insertions: 0.03890811226451365
# Average - Mismatch: 0.045934074375060314 Deletions: 0.08497070122487656 Insertions: 0.04075836932585308

```

# 5 resquiggle
```{linux}
# Vero-infected 
input="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/02Tran_model_CellData/01rawdata/"
sars_2019_ref="/mnt/raid61/Personal_data/songjunwei/reference/SARS_CoV2/nCoV-2019.reference.fasta"

conda activate base #2101
for i in Train Valid ;
do
multi_to_single_fast5 --input_path $input/${i}_SARS_fast5  --save_path /mnt/data3/songjunwei/data/model_SARS/${i}_SARS_f5Sg --recursive -t 10 #线程不能太大 卡住跑不过去
done

cd /mnt/data3/songjunwei/data/model_SARS/
for i in  Valid ;
do
conda activate base 
tombo preprocess annotate_raw_with_fastqs --fast5-basedir ${i}_SARS_f5Sg --fastq-filenames $input/${i}_SARS_guppy.fq --processes 10 --overwrite
conda activate tom #2101
tombo resquiggle ${i}_SARS_f5Sg $sars_2019_ref \
 --q-score 0  --signal-matching-score 4  \
 --rna --overwrite --processes 20  --ignore-read-locks --num-most-common-errors 5
done

## clinical
sars_248_ref="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/SARS_clinical_248.fa"
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan
conda activate base #2101
multi_to_single_fast5 --input_path SARS_fast5  --save_path SARS_f5Sg --recursive -t 5 #线程不能太大 卡住跑不过去
tombo preprocess annotate_raw_with_fastqs --fast5-basedir SARS_f5Sg --fastq-filenames SARS_guppy.fq --processes 10 --overwrite
conda activate tom #2101
tombo resquiggle SARS_f5Sg $sars_248_ref \
 --q-score 0  --signal-matching-score 4  \
 --rna --overwrite --processes 20  --ignore-read-locks --num-most-common-errors 5

#zd download
cd /mnt/raid66/Personal_data/zhangdan/Basecalling/00.downloaddata/
multi_to_single_fast5 --input_path fast5_new  --save_path f5sg_new --recursive -t 5 
tombo preprocess annotate_raw_with_fastqs --fast5-basedir f5sg_new --fastq-filenames guppy/guppy_bc.fastq --processes 10 --overwrite
conda activate tom #2101
sars_2019_ref="/mnt/raid61/Personal_data/songjunwei/reference/SARS_CoV2/nCoV-2019.reference.fasta"
tombo resquiggle f5sg_new $sars_2019_ref \
 --q-score 0  --signal-matching-score 4  \
 --rna --overwrite --processes 20  --ignore-read-locks --num-most-common-errors 5

```

#6 RNAvirus model
## reads pre
```{linux}
# /mnt/raid61/Personal_data/songjunwei/DRS_RTA/08viurs_model/ #存储路径
cd /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/01rawdata/
guppy_basecaller \
  -i ./ -r \
  --save_path ./guppy --device cuda:0 \
  -c rna_r9.4.1_70bps_hac.cfg \
  --min_qscore 7  --verbose_logs  \
  --chunks_per_runner 256  --chunk_size 100 \
  --num_callers 4 --gpu_runners_per_device 8
cat ./guppy/pass/* > ./guppy/merge.fq
for i in PEDV PRRSV SVV; do
virus_ref="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/"
minimap2 -ax map-ont -t 30  $virus_ref/${i}.fasta .//guppy/merge.fq |  samtools sort -@ 30 | samtools view -b >  .//guppy/merge_${i}.bam
done
```

```{r}
readid_list <- 
lapply(list.files("/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/01rawdata/guppy/",pattern = "bam",full.names = T), function(bam){
  sars_align <-
    readGAlignments(bam,
                    param = ScanBamParam(
                      flag = scanBamFlag(
                        isSupplementaryAlignment = FALSE,
                        isPaired = FALSE,
                        isSecondaryAlignment = FALSE
                      ),
                      mapqFilter = 60
                    ),
                    use.names = T)
  names(sars_align)
})
names(readid_list) <- c("PEDV","PRRSV","SVV")

PRRSV_read <- fread("/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/Train_PRRSV_readid.txt",header = F)
Train <- intersect(PRRSV_read$V1,readid_list$PRRSV)
Valid_Test <- setdiff(readid_list$PRRSV,PRRSV_read$V1)
Valid <- Valid_Test[1:157]
Test <- Valid_Test[158:307]
write.table(Train,"/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/Train_PRRSV_readid.txt",row.names = F,col.names = F,quote = F)
write.table(Valid,"/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/Valid_PRRSV_readid.txt",row.names = F,col.names = F,quote = F)
write.table(Test,"/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/Test_PRRSV_readid.txt",row.names = F,col.names = F,quote = F)

SVV_read <- fread("/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/Train_SVV_readid.txt",header = F)
Train <- intersect(SVV_read$V1,readid_list$SVV)
Valid_Test <- setdiff(readid_list$SVV,SVV_read$V1)
Valid <- Valid_Test[1:1000]
Test <- Valid_Test[1001:2000]
write.table(Train,"/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/Train_SVV_readid.txt",row.names = F,col.names = F,quote = F)
write.table(Valid,"/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/Valid_SVV_readid.txt",row.names = F,col.names = F,quote = F)
write.table(Test,"/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/Test_SVV_readid.txt",row.names = F,col.names = F,quote = F)

sars_align <-
  readGAlignments(
    "/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/02Tran_model_CellData/01rawdata/SARS_CellD.bam",
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
SARS_read <- names(sars_align)
Train <- SARS_read[1:15000]
Valid <- SARS_read[15001:18000]
Test <- SARS_read[18001:21000]
write.table(Train,"/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/Train_SARS_readid.txt",row.names = F,col.names = F,quote = F)
write.table(Valid,"/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/Valid_SARS_readid.txt",row.names = F,col.names = F,quote = F)
write.table(Test,"/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/Test_SARS_readid.txt",row.names = F,col.names = F,quote = F)

write.table(readid_list$PEDV,"/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/Test_PEDV_readid.txt",row.names = F,col.names = F,quote = F)
```

```{linux}
cd /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w
cat Train_SARS_readid.txt Valid_SARS_readid.txt Test_SARS_readid.txt > /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/01rawdata/SARS_readid.txt
virus_ref="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/"
cat /mnt/raid61/Personal_data/songjunwei/reference/SARS_CoV2/nCoV-2019.reference.fasta $virus_ref/PRRSV.fasta  $virus_ref/SVV.fasta > virus_ref.fa

cat Train_*txt > Train_virus_readid.txt
cat Valid_*txt > Valid_virus_readid.txt


cd /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/01rawdata
fast5_subset -i /mnt/raid65/PublicData/Zenodo/Cell_2020_SARS-Vero/Vero-Infected/20200222_0650_MN26136_FAN04901_ada1e2bf/fast5 \
            -s SARS_fast5 \
            --batch_size  2000 \
            --read_id_list SARS_readid.txt \
            --recursive  --ignore_symlinks \
            -t 30
fast5_subset -i PRRSV_samples_fast5 \
            -s PRRSV_fast5 \
            --batch_size  4000 \
            --read_id_list PRRSV_readid.txt \
            --recursive  --ignore_symlinks \
            -t 30

for i in Train Valid Test; 
do
seqkit grep -f /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/${i}_SARS_readid.txt /mnt/raid65/PublicData/Zenodo/Cell_2020_SARS-Vero/Vero-Infected/20200222_0650_MN26136_FAN04901_ada1e2bf/VeroInf24h.all.fastq > /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/${i}_SARS_guppy.fq
done

cd /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/
mkdir QC
NanoComp --fastq  Train_SARS_guppy.fq  Valid_SARS_guppy.fq Test_SARS_guppy.fq -o ./QC -t 30  -p SARS_Cell_ --raw
```


## train
```{linux}
cd /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w
for i in Train Valid ;  
do
fast5_subset -i /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/01rawdata/ \
            -s ${i}_virus_fast5 \
            --batch_size  2000 \
            --read_id_list ${i}_virus_readid.txt \
            --recursive  --ignore_symlinks \
            -t 15
done

for i in Train Valid ;  
do
guppy_basecaller \
  -i ${i}_virus_fast5 -r \
  --save_path ${i}_virus_fast5/ --device cuda:0 \
  -c rna_r9.4.1_70bps_hac.cfg \
  --min_qscore 7  --verbose_logs  \
  --chunks_per_runner 256  --chunk_size 100 \
  --num_callers 1 --gpu_runners_per_device 1
cat ${i}_virus_fast5/pass/* > ${i}_virus_guppy.fq
minimap2 -ax map-ont -t 30  virus_ref.fa  ${i}_virus_guppy.fq |  samtools sort -@ 30 | samtools view -b > ${i}_virus_guppy.bam
done

# tayaki
conda activate base #2048
cd /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w
tayaki="/mnt/raid5/Personal/minion/software/taiyaki/bin/"
for i in Train Valid ;
do
python $tayaki/generate_per_read_params.py ${i}_virus_fast5 --output 01tayaki/${i}_virus_readid.tsv
done
for i in Train Valid ;
do
python $tayaki/get_refs_from_sam.py virus_ref.fa  ${i}_virus_guppy.bam  --output 01tayaki/${i}_virus_tayaki.fa --reverse
python $tayaki/prepare_mapped_reads.py ${i}_virus_fast5 01tayaki/${i}_virus_readid.tsv 01tayaki/${i}_virus_tayaki.h5 \
        /mnt/raid5/Personal/minion/software/taiyaki/models/r941_rna_minion_upgrae.checkpoint     01tayaki/${i}_virus_tayaki.fa --jobs 30 --overwrite 
done

# rodan
rodan="/mnt/raid61/Personal_data/songjunwei/software/RODAN/"
cd 01tayaki
for i in Train Valid ;
do
mkdir -p rodan_${i}
python  $rodan/gendata.py -i ${i}_virus_tayaki.h5 --outdir rodan_${i}
mv rodan_${i}/train.hdf5 rna-${i}.hdf5
rm -rf rodan_${i}
done
# cp $rodan/rna.config ./
python  $rodan/model.py -c ./rna.config -a $rodan/rnaarch -n Virus --arch $rodan/rnaarch --savedir ./runs --rna -w 10 -l


rodan="/mnt/raid5/Personal/minion/software/RODAN/"
cd /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/01tayaki
python  $rodan/model_mod.py -c ./rna.config -a $rodan/rnaarch -n Virus --arch $rodan/rnaarch --savedir ./runs_withRodan --rna -w 10 -m $rodan/rna.torch -l
```


## test
```{linux}
f5="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/SARS_fast5"
test_dir="/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/03Test/SARS/"
rodan="/mnt/raid61/Personal_data/songjunwei/software/RODAN/"
mkdir $test_dir/fa $test_dir/sam $test_dir/accuracy
for model in /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/01tayaki/runs/Virus-epoch*.torch; do
    epoch=$(basename "$model" | sed 's/M_CellD-epoch\([0-9]*\).torch/\1/')
    output_fa=$test_dir/fa/SARS_virus_ep${epoch}.fa
    output_sam=$test_dir/sam/SARS_virus_ep${epoch}.sam
    accuracy_result_file=$test_dir/accuracy/SARS_virus_ep${epoch}.txt

    # 运行 basecalling
    python "$rodan/basecall.py" -m "$model" "$f5" --arch "$rodan/rnaarch" -b 100 > $output_fa
    minimap2 --secondary=no -ax map-ont -t 32 --cs /mnt/raid61/Personal_data/songjunwei/reference/SARS_CoV2/nCoV-2019.reference.fasta $output_fa > $output_sam
    python "$rodan/accuracy.py" "$output_sam" /mnt/raid61/Personal_data/songjunwei/reference/SARS_CoV2/nCoV-2019.reference.fasta  > "$accuracy_result_file"
done

# RodanNew
cd /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/
for i in  SARS ; do #PEDV  PRRSV SVV
fast5_subset -i 01rawdata/${i}_fast5 \
            -s  03Test/Test_${i}_fast5 \
            --batch_size  2000 \
            --read_id_list 02Model_use1.5w/Test_${i}_readid.txt \
            --recursive  --ignore_symlinks \
            -t 15
test_dir="03Test/$i"
mkdir -p $test_dir/fa $test_dir/sam $test_dir/accuracy
rodan="/mnt/raid61/Personal_data/songjunwei/software/RODAN/"
virus_ref="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/"
for model in /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/01tayaki/epoch2n/Virus-epoch26.torch; do
    epoch=$(basename "$model" | sed 's/\([0-9]*\).torch/\1/')
    output_fa=$test_dir/fa/${epoch}.fa
    output_sam=$test_dir/sam/${epoch}.sam
    accuracy_result_file=$test_dir/accuracy/${epoch}.txt

    # 运行 basecalling
    #python "$rodan/basecall.py" -m "$model" 03Test/Test_${i}_fast5 --arch "$rodan/rnaarch" -b 80 > $output_fa
    minimap2 --secondary=no -ax map-ont -t 32 --cs $virus_ref/${i}.fasta $output_fa > $output_sam
    python "$rodan/accuracy.py" $output_sam $virus_ref/${i}.fasta  > $accuracy_result_file
done
done


# guppy and rodan
cd /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/03Test
virus_ref="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/"
rodan="/mnt/raid61/Personal_data/songjunwei/software/RODAN/"
for i in SARS ; do #PEDV PRRSV SVV
  # guppy_basecaller \
  #   -i Test_${i}_fast5 \
  #   --save_path ${i}/guppy --device cuda:0 \
  #   -c rna_r9.4.1_70bps_hac.cfg \
  #   --min_qscore 7  --verbose_logs  \
  #   --chunks_per_runner 256  --chunk_size 100 \
  #   --num_callers 8 --gpu_runners_per_device 16
  # cat ${i}/guppy/pass/* > ${i}/${i}_guppy.fq
  # python $rodan/basecall.py -m $rodan/rna.torch Test_${i}_fast5 --arch $rodan/rnaarch -b 80 > ${i}/${i}_roadnRaw.fa

  minimap2 --secondary=no -ax map-ont -t 32 --cs $virus_ref/${i}.fasta ${i}/${i}_guppy.fq > ${i}/${i}_guppy.sam
  minimap2 --secondary=no -ax map-ont -t 32 --cs $virus_ref/${i}.fasta ${i}/${i}_roadnRaw.fa > ${i}/${i}_roadnRaw.sam
  python "$rodan/accuracy.py" ${i}/${i}_guppy.sam $virus_ref/${i}.fasta  > ${i}/${i}_guppy_accuracy.txt
  python "$rodan/accuracy.py" ${i}/${i}_roadnRaw.sam $virus_ref/${i}.fasta  > ${i}/${i}_rodanRaw_accuracy.txt
done
```


```{r}
files <- list.files(path = "/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/03Test/",pattern = "Virus-epoch.*\\.txt",full.names = T,recursive = T)
files <- grep("old",files,invert = T,value = T)
results_df <- data.frame()
# 遍历每个文件，读取内容，并提取所需信息
for (file in files) {
  epoch <- gsub("SARS_virus_epVirus-epoch([0-9]+).*\\.txt", "\\1", basename(file))  # 提取 epoch 编号
  res <- accuracy_to_df(file,epoch=epoch)
  results_df <- rbind(results_df, res)
}

results_df$speice <- lapply(c("PEDV","PRRSV","SARS","SVV"), function(x) rep(x,29)) %>% unlist()
results_df$idex <- paste0(results_df$speice,results_df$Epoch)
rownames(results_df) <-results_df$idex 
results_df$Epoch <- as.numeric(as.character(results_df$Epoch))
results_df <- data.table(results_df)[order(speice,Epoch)]

write.csv(results_df,"../08viurs_model/Virus_Epochs_res.csv",row.names = F)


results_df <- na.omit(results_df)
# 1. Total 随 Epoch 变化的图
p1 <- ggarrange(plotlist = lapply(split(results_df,f = results_df$speice), function(pdata){
  ggplot(data = pdata, aes(x = Epoch, y = Total,color=speice)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Total Over Epochs", x = "Epoch", y = "Total")
}))


# 2. Median 和 Average Accuracy 随 Epoch 变化的图
accuracy_df <- melt(results_df[, c("speice","Epoch", "MedianAccuracy", "AverageAccuracy")], id.vars = c("Epoch","speice"))
p2 <- ggplot(data = accuracy_df, aes(x = Epoch, y = value, color = speice,linetype = variable)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Accuracy", x = "Epoch", y = "Accuracy") +
  scale_color_brewer(palette = "Set1")

# 3. Median 和 Average Mismatch 随 Epoch 变化的图
mismatch_df <- melt(results_df[, c("speice","Epoch", "MismatchMedian", "MismatchAverage")], id.vars = c("Epoch","speice"))
p3 <- ggplot(data = mismatch_df, aes(x = Epoch, y = value, color = speice,linetype = variable)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Mismatch", x = "Epoch", y = "Mismatch") +
  scale_color_brewer(palette = "Set1")

# 4. Median 和 Average Deletions 随 Epoch 变化的图
deletions_df <- melt(results_df[, c("speice","Epoch", "DeletionsMedian", "DeletionsAverage")], id.vars = c("Epoch","speice"))
p4 <- ggplot(data = deletions_df, aes(x = Epoch, y = value, color = speice,linetype = variable)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Deletions", x = "Epoch", y = "Deletions") +
  scale_color_brewer(palette = "Set1")

# 5. Median 和 Average Insertions 随 Epoch 变化的图
insertions_df <- melt(results_df[, c("speice","Epoch", "InsertionsMedian", "InsertionsAverage")], id.vars = c("Epoch","speice"))
p5 <- ggplot(data = insertions_df, aes(x = Epoch, y = value, color = speice,linetype = variable)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Insertions", x = "Epoch", y = "Insertions") +
  scale_color_brewer(palette = "Set1")

pall <- ggarrange(p2,p3,p4,p5)

ggsave("../08viurs_model/Virus_Epochs_total.pdf",p1,width = 6,height = 5)
ggsave("../08viurs_model/Virus_Epochs_accuracy_res.pdf",pall,width = 8,height = 6)
```

## compare
```{r}
#guppy和Rodan原始
files <- list.files("/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/03Test/",pattern = "accuracy.txt",recursive = T,full.names = T)
files <- grep("old",files,invert = T,value = T)[7:14]
files <- c(files,"/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/01accuracy/accufiles/accuracy_guppy.txt",
            "/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/01accuracy/accufiles/accuracy_rodanraw.txt",
           "/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/03Test/old/SARS/accuracy/SARS_virus_epVirus-epoch26.torch.txt")
df1 <- lapply(files, function(file){
  df <- accuracy_to_df(file,epoch = strsplit(basename(file),"_")[[1]][2])
  df$species <- strsplit(basename(file),"_")[[1]][1]
  return(df)
}) %>% rbindlist()
df1[9:11,"species"] <- "SARS_lab406"
df1[5:6,"species"] <- "SARS_CellData"
df1[9:11,"Epoch"] <- c("guppy","rodanRaw","Virus-epoch26.txt")


#混和模型选用26
select <-list.files("/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/03Test/",pattern = "Virus-epoch26.txt",recursive = T,full.names = T)
df2 <- lapply(select, function(file){
  df <- accuracy_to_df(file,epoch = gsub("SARS_virus_epVirus-epoch([0-9]+).*\\.txt", "\\1", basename(file)))
  return(df)
}) %>% rbindlist()
df2$species <- c("PEDV","PRRSV","SARS_CellData","SVV")
```

#7 Human Pb
```{linux}
cd /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/03Test
virus_ref="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/"
rodan="/mnt/raid61/Personal_data/songjunwei/software/RODAN/"

python $rodan/basecall.py -m /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/01tayaki/epoch2n/Virus-epoch26.torch /mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/04f5sg --arch $rodan/rnaarch -b 80 > /mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/03accuracy/fafiles/rodanep26.fa
cd /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/03Test/Pberghei
ln -s /mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/03accuracy/fafiles/rodanep26.fa Pberghei_rodanep26.fa
ln -s /mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/03accuracy/fafiles/Pb_rodan_raw.fa Pberghei_rodanRaw.fa
ln -s /mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/03accuracy/fafiles/Pb_guppy.fastq Pberghei_guppy.fq


i="HomSap"
fast5="/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/01rawdata/HomSap_fast5"
multi_to_single_fast5 --input_path /mnt/raid61/Personal_data/songjunwei/DRS_RTA/01Rawdata/Tumor4/raw_f5  --save_path /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/01rawdata/HomSap_fast5 --recursive -t 8  --ignore_symlinks

guppy_basecaller \
  -i $fast5 -r \
  --save_path ${i}/guppy --device cuda:0 \
  -c rna_r9.4.1_70bps_hac.cfg \
  --min_qscore 7  --verbose_logs  \
  --chunks_per_runner 256  --chunk_size 100 \
  --num_callers 16 --gpu_runners_per_device 16
cat ${i}/guppy/pass/* > ${i}/${i}_guppy.fq

rodan="/mnt/raid5/Personal/minion/software/RODAN"
python $rodan/basecall.py -m $rodan/rna.torch $fast5 --arch $rodan/rnaarch -b 50 > ${i}/${i}_rodanRaw.fa
python $rodan/basecall.py -m /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/01tayaki/epoch2n/Virus-epoch26.torch $fast5 --arch $rodan/rnaarch -b 50 > ${i}/${i}_rodanep26.fa

cd /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/03Test
virus_ref="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/"
for i in Pberghei HomSap;do
for sof in guppy.fq rodanRaw.fa rodanep26.fa;do
  minimap2 --secondary=no -ax map-ont -t 32 --cs $virus_ref/${i}.fasta ${i}/${i}_${sof} > ${i}/${i}_${sof}.sam
  python "$rodan/accuracy.py" ${i}/${i}_${sof}.sam $virus_ref/${i}.fasta  > ${i}/${i}_${sof}_accuracy.txt
done
done

# human bam
cd /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/03Test/HomSap
virus_ref="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/"
minimap2 -ax map-ont -t 30  $virus_ref/HomSap.fasta \
       HomSap_guppy.fq |  samtools sort -@ 30 | samtools view -b > HomSap_guppy.bam
```

```{r}
Hum_Pb <-c(list.files("/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/03Test/HomSap/",pattern = "_accuracy.txt",recursive = T,full.names = T),list.files("/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/03Test/Pberghei//",pattern = "_accuracy.txt",recursive = T,full.names = T))
df3 <- lapply(Hum_Pb, function(file){
  df <- accuracy_to_df(file,epoch =  str_split_fixed(basename(file),"\\.",2)[,1])
  return(df)
}) %>% rbindlist()

df3$species <- str_split_fixed(df3$Epoch,"_",2)[,1]
df3$Epoch <- str_split_fixed(df3$Epoch,"_",2)[,2]

df_res <- rbind(df1,df2,df3)
df_res$Epoch <- gsub("Virus-epoch26.txt","rodanep26",df_res$Epoch)
df_res$Epoch <- factor(df_res$Epoch,levels = c("guppy","rodanRaw","rodanep26"))
df_res <- df_res[order(df_res$species,df_res$Epoch),]
df_res <- data.frame(df_res)[,c("species",colnames(df_res)[1:13])]
df_res[,4:14] <- round(df_res[,4:14],3)
write.csv(df_res,"../08viurs_model/Virus_compared_res.csv",row.names = F)
```

#8 lla model
```{linux}
su lla
conda activate basecall

cd /mnt/raid5/lla/RNA/src/
model='/mnt/raid5/lla/RNA/src/model/pb+virus+rodan/weights_10_98.00.tar'
output="/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/04lla_Test/fafiles"
batchsize=64
time python basecall_rodan.py  --checkpoint-file ${model} --output-file $output/Pberghei_lla10.fa /mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/04f5sg --max-reads 3000 --batchsize ${batchsize} --reverse  --beam-size 5
time python basecall_rodan.py  --checkpoint-file ${model} --output-file   $output/HomSap_lla10.fa /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/01rawdata/Tumor4_f5sg --max-reads 3000 --batchsize ${batchsize} --reverse  --beam-size 5
time python basecall_rodan.py  --checkpoint-file ${model} --output-file $output/LabSARS_lla10.fa /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/SARS_fast5 --max-reads 3000 --batchsize ${batchsize} --reverse  --beam-size 5
for i in  PEDV PRRSV SVV SARS;do #
time python basecall_rodan.py  --checkpoint-file ${model} --output-file $output/${i}_lla10.fa /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/03Test/Test_${i}_fast5 --max-reads 3000 --batchsize ${batchsize} --reverse  --beam-size 5
done
```

```{linux}
cd /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/04lla_Test/
mkdir samfiles
mkdir accfiles
ref="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/"
i="LabSARS"
minimap2 --secondary=no -ax map-ont -t 32 --cs $ref/SARS.fasta  fafiles/${i}_lla10.fa > samfiles/${i}_lla10.sam
python /mnt/raid5/Personal/minion/software/RODAN/accuracy.py samfiles/${i}_lla10.sam $ref/SARS.fasta  > accfiles/${i}_lla10_accuracy.txt

for i in Pberghei HomSap PEDV PRRSV SVV SARS ;do
  minimap2 --secondary=no -ax map-ont -t 32 --cs $ref/${i}.fasta  fafiles/${i}_lla10.fa > samfiles/${i}_lla10.sam
  python /mnt/raid5/Personal/minion/software/RODAN/accuracy.py samfiles/${i}_lla10.sam $ref/${i}.fasta  > accfiles/${i}_lla10_accuracy.txt
done
```

```{r}
#lla 模型
select <-list.files("/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/04lla_Test/accfiles/",pattern = "_accuracy.txt",recursive = T,full.names = T)
df4 <- lapply(select, function(file){
  df <- accuracy_to_df(file,epoch = gsub("_lla10_accuracy.txt", "\\1", basename(file)))
  return(df)
}) %>% rbindlist()
df4$species <- df4$Epoch
df4$Epoch <- "lla_model"

df4 <- data.frame(df4)[,c("species",colnames(df4)[1:13])]
df4[,4:14] <- round(df4[,4:14],3)
df4$species[2] <- "SARS_lab406"

df_res <- read.csv("../08viurs_model/Virus_compared_res.csv")
df_res1 <- rbind(df_res,df4)
df_res1$Epoch <- factor(df_res1$Epoch,levels = c("guppy","rodanRaw","rodanep26","lla_model"))
df_res1 <- df_res1[order(df_res1$species,df_res1$Epoch),]
write.csv(df_res1,"../08viurs_model/Virus_compared_res1.csv",row.names = F)
```

```{r}
lla_test <- list.files("/mnt/raid5/lla/RNA/src/","txt",full.names = T)[c(1,5,6,7,8,11)]
df5 <- lapply(lla_test, function(file){
  df <- accuracy_to_df(file,epoch = gsub(".txt", "\\1", basename(file)))
  return(df)
}) %>% rbindlist()
write.csv(df5,"../08viurs_model/Rodan_Species_compared_reslla.csv")
```


#9 genernal model with pre-trained(rodan)
```{linux}
cd /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/01rawdata
cp -r /mnt/raid61/Personal_data/songjunwei/DRS_RTA/07Pb_model/01rawdata/*fast5 ./Pberghei_fast5/
cp -r /mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/SARS_f5Sg ./LabSARS_fast5
```

## reads pre
```{r}
#input
virs <- list.files("/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/01rawdata/guppy/",pattern = "bam",full.names = T)
bamfiles <-
  c(virs,
    "/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/03Test/HomSap/HomSap_guppy.bam",
    "/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/02Tran_model_CellData/01rawdata/SARS_CellD.bam",
    "/mnt/raid61/Personal_data/songjunwei/DRS_RTA/02align/batch3_9Pb/PbergheiANKA.bam"
  )

readid_list <- 
lapply(bamfiles, function(bam){
  align <- readGAlignments(bam,param = ScanBamParam(
                      flag = scanBamFlag(
                        isSupplementaryAlignment = FALSE,
                        isPaired = FALSE,
                        isSecondaryAlignment = FALSE),mapqFilter = 60),use.names = T)
  names(align)
})
names(readid_list) <- c("PEDV","PRRSV","SVV","HomSap","SARS","Pberghei")
SARS_lab <- fread("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/08SARS_rodan/SARS_readid.txt",header = F)
readid_list$LabSARS <- SARS_lab$V1
lapply(readid_list, length)

#output
data_list <- data.frame(
  names = c("SARS", "HomSap", "Pberghei", "SVV", "PRRSV", "LabSARS"),
  train = c(2000, 5000, 5000, 2000, 1272, 300),
  valid = c(500, 500, 500, 500, 150, 56),
  test = c(2000, 2000, 2000, 2000, 150, 50)
)
output_dir="/mnt/raid5/Personal/minion/DRS_mul/09gneral_model/01read_process/"
for (i in 1:nrow(data_list)) {
  name <- data_list$names[i]
  read_ids <- unique(readid_list[[name]])
  train_ids <- read_ids[1:data_list$train[i]]
  valid_ids <- read_ids[(data_list$train[i] + 1):(data_list$train[i] + data_list$valid[i])]
  test_ids <- read_ids[(data_list$train[i] + data_list$valid[i] + 1):(data_list$train[i] + data_list$valid[i] + data_list$test[i])]
  
  # 将读取ID写入到文本文件中
  write.table(train_ids, paste0(output_dir,"/Train_", name, "_readid.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(valid_ids, paste0(output_dir,"Valid_", name, "_readid.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(test_ids, paste0(output_dir,"Test_", name, "_readid.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}
write.table(readid_list$PEDV, paste0(output_dir,"/Test_PEDV_readid.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

```

train
```{linux}
virus_ref="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/"
cat /mnt/raid61/Personal_data/songjunwei/reference/SARS_CoV2/nCoV-2019.reference.fasta $virus_ref/PRRSV.fasta  $virus_ref/SVV.fasta $virus_ref/HomSap.fasta $virus_ref/Pberghei.fasta > /mnt/raid5/Personal/minion/DRS_mul/09gneral_model/01read_process/general_ref.fa

cd /mnt/raid5/Personal/minion/DRS_mul/09gneral_model/01read_process/
cat Train_*txt > Train_general_readid.txt
cat Valid_*txt > Valid_general_readid.txt
for i in Train Valid ;  
do
fast5_subset -i /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/01rawdata/ \
            -s ${i}_general_fast5 \
            --batch_size  2000 \
            --read_id_list ${i}_general_readid.txt \
            --recursive  --ignore_symlinks \
            -t 15
done

for i in  LabSARS PEDV SARS SVV PRRSV HomSap Pberghei  ; do 
fast5_subset -i /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/01rawdata/${i}_fast5 \
            -s  /mnt/raid5/Personal/minion/DRS_mul/09gneral_model/02test_data/Test_${i}_fast5 \
            --batch_size  4000 \
            --read_id_list /mnt/raid5/Personal/minion/DRS_mul/09gneral_model/01read_process/Test_${i}_readid.txt \
            --recursive  --ignore_symlinks \
            -t 30
done

for i in Train Valid ;  
do
guppy_basecaller \
  -i ${i}_general_fast5 -r \
  --save_path ${i}_general_fast5/ --device cuda:0 \
  -c rna_r9.4.1_70bps_hac.cfg \
  --min_qscore 7  --verbose_logs  \
  --chunks_per_runner 256  --chunk_size 100 \
  --num_callers 1 --gpu_runners_per_device 1
cat ${i}_general_fast5/pass/* > ${i}_general_guppy.fq
minimap2 -ax map-ont -t 30  general_ref.fa  ${i}_general_guppy.fq |  samtools sort -@ 30 | samtools view -b > ${i}_general_guppy.bam
done

# tayaki
tayaki="/mnt/raid5/Personal/minion/software/taiyaki/bin/"
for i in Train Valid ;
do
python $tayaki/generate_per_read_params.py ${i}_general_fast5 --output 01tayaki/${i}_general_readid.tsv
python $tayaki/get_refs_from_sam.py general_ref.fa  ${i}_general_guppy.bam  --output 01tayaki/${i}_general_tayaki.fa --reverse
python $tayaki/prepare_mapped_reads.py ${i}_general_fast5 01tayaki/${i}_general_readid.tsv 01tayaki/${i}_general_tayaki.h5 \
        /mnt/raid5/Personal/minion/software/taiyaki/models/r941_rna_minion_upgrae.checkpoint     01tayaki/${i}_general_tayaki.fa --jobs 30 --overwrite 
done

# rodan
rodan="/mnt/raid5/Personal/minion/software/RODAN/"
cd /mnt/raid5/Personal/minion/DRS_mul/09gneral_model/01read_process/01tayaki
for i in Train Valid ;
do
mkdir -p rodan_${i}
python  $rodan/gendata.py -i ${i}_general_tayaki.h5 --outdir rodan_${i}
mv rodan_${i}/train.hdf5 rna-${i}.hdf5
rm -rf rodan_${i}
done

cp "/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/01tayaki/rna.config" ./
rodan="/mnt/raid5/Personal/minion/software/RODAN/"
python  $rodan/model_mod.py -c ./rna.config -a $rodan/rnaarch -n General --arch $rodan/rnaarch --savedir ./runs_withRodan --rna -w 10 -m $rodan/rna.torch -l

```

##test
```{linux}
cd /mnt/raid5/Personal/minion/DRS_mul/09gneral_model/03try
rodan="/mnt/raid5/Personal/minion/software/RODAN/"
rodan="/mnt/raid61/Personal_data/songjunwei/software/RODAN/"
for i in  SARS PEDV PRRSV SVV HomSap Pberghei; do  #LabSARS
  python $rodan/basecall.py -m /mnt/raid5/Personal/minion/DRS_mul/08viurs_model/02Model_use1.5w/01tayaki/runs_withRodan/Virus-epoch5.torch \
      /mnt/raid5/Personal/minion/DRS_mul/09gneral_model/02test_data/Test_${i}_fast5 --arch $rodan/rnaarch -b 100 > fafiles/${i}_virusWRon.fa  # batch size 100:30% memory
  minimap2 --secondary=no -ax map-ont -t 32 --cs $virus_ref/${i}.fasta fafiles/${i}_virusWRon.fa > samfiles/${i}_virusWRon.sam
  python $rodan/accuracy.py samfiles/${i}_virusWRon.sam $virus_ref/${i}.fasta  > accufiles/${i}_virusWRon_accuracy.txt
done



#10 densecall ⭐️
## Basecalling
```{linux}
rodan="/mnt/raid5/Personal/minion/software/RODAN/"
output="/mnt/raid5/Personal/minion/DRS_mul/09gneral_model/04try_lla/fafiles"
## data1
for i in SARS PEDV  SVV HomSap Pberghei; do  
test_data="/mnt/raid5/Personal/minion/DRS_mul/09gneral_model/02test_data/"
cdup densecall
densecall basecaller --amp --output-file $output/${i}.fq rna_r9.4.1_hac@v1.0  \
   $test_data/Test_${i}_fast5 -b 300
guppy_basecaller \
    -i $test_data/Test_${i}_fast5 \
    --save_path $output/${i}_guppy --device cuda:0 \
    -c rna_r9.4.1_70bps_hac.cfg \
    --min_qscore 7  --verbose_logs  \
    --chunks_per_runner 256  --chunk_size 100 \
    --num_callers 8 --gpu_runners_per_device 16
  cat $output/${i}_guppy/pass/* > $output/${i}_guppy.fq
cdup base
python $rodan/basecall.py -m $rodan/rna.torch \
  $test_data/Test_${i}_fast5 \
  --arch $rodan/rnaarch -b 300 > $output/${i}_roadn.fa
done


## data2
for i in LabSARS PRRSV; do  
test_data="/mnt/raid5/Personal/minion/DRS_mul/08viurs_model/01rawdata/"
cdup densecall
densecall basecaller --amp --output-file $output/${i}.fq rna_r9.4.1_hac@v1.0  \
   $test_data/${i}_fast5 -b 300
guppy_basecaller \
    -i $test_data/${i}_fast5 -r \
    --save_path $output/${i}_guppy --device cuda:0 \
    -c rna_r9.4.1_70bps_hac.cfg \
    --min_qscore 7  --verbose_logs  \
    --chunks_per_runner 256  --chunk_size 100 \
    --num_callers 32 --gpu_runners_per_device 32
cat $output/${i}_guppy/pass/* > $output/${i}_guppy.fq
cdup base
python $rodan/basecall.py -m $rodan/rna.torch \
  $test_data/${i}_fast5 \
  --arch $rodan/rnaarch -b 300 > $output/${i}_roadn.fa
done

## data3
for i in arabidopsis human mouse poplar yeast;do
  mkdir /mnt/raid5/Personal/minion/DRS_mul/09gneral_model/04try_lla/rodan_data/${i}/
  cp /mnt/raid5/lla/RNA/data/${i}/dataset/0/* /mnt/raid5/Personal/minion/DRS_mul/09gneral_model/04try_lla/rodan_data/${i}/
done
cp /mnt/raid5/lla/RNA/data/arabidopsis/0/* /mnt/raid5/Personal/minion/DRS_mul/09gneral_model/04try_lla/rodan_data/arabidopsis/

rodan="/mnt/raid5/Personal/minion/software/RODAN/"
output="/mnt/raid5/Personal/minion/DRS_mul/09gneral_model/04try_lla/fafiles"
for i in arabidopsis human mouse poplar yeast; do  
test_data="/mnt/raid5/Personal/minion/DRS_mul/09gneral_model/04try_lla/rodan_data/"
cdup densecall
densecall basecaller --amp --output-file $output/${i}.fq rna_r9.4.1_hac@v1.0  \
   $test_data/${i} -b 300
guppy_basecaller \
    -i $test_data/${i} -r \
    --save_path $output/${i}_guppy --device cuda:0 \
    -c rna_r9.4.1_70bps_hac.cfg \
    --min_qscore 7  --verbose_logs  \
    --chunks_per_runner 256  --chunk_size 100 \
    --num_callers 32 --gpu_runners_per_device 32
cat $output/${i}_guppy/pass/* > $output/${i}_guppy.fq
cdup base
python $rodan/basecall.py -m $rodan/rna.torch $test_data/${i} \
  --arch $rodan/rnaarch -b 300 > $output/${i}_roadn.fa
done
```

## all:align&accu
```{linux}
cd /mnt/raid5/Personal/minion/DRS_mul/09gneral_model/04try_lla
virus_ref="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/"
rodan="/mnt/raid5/Personal/minion/software/RODAN/"
for i in LabSARS PRRSV PEDV SVV HomSap Pberghei ; do  #SARS
  ref_name=`echo $i | sed 's/LabSARS/SARS/'`
  for type in ".fq" "_guppy.fq" "_roadn.fa"; do
    minimap2 --secondary=no -ax map-ont -t 32 --cs $virus_ref/${ref_name}.fasta fafiles/${i}${type} > sam_full_files/${i}${type}.sam
    python "$rodan/accuracy.py" sam_full_files/${i}${type}.sam $virus_ref/${ref_name}.fasta  > accufiles/${i}${type}_accuracy.txt
  done
done

cdup base
cd /mnt/raid5/Personal/minion/DRS_mul/09gneral_model/04try_lla
virus_ref="/mnt/raid5/lla/RNA/data/ref/"
rodan="/mnt/raid5/Personal/minion/software/RODAN/"
for i in arabidopsis human mouse yeast poplar; do 
  for type in ".fq" "_guppy.fq" "_roadn.fa"; do
    minimap2 --secondary=no -ax map-ont -t 32 --cs $virus_ref/${i}_reference.fasta fafiles/${i}${type} > sam_full_files/${i}${type}.sam
    python "$rodan/accuracy.py" sam_full_files/${i}${type}.sam $virus_ref/${i}_reference.fasta > accufiles/${i}${type}_accuracy.txt
  done
done

## stat
for i in LabSARS PRRSV PEDV SVV HomSap Pberghei arabidopsis human mouse yeast poplar; do 
  for type in ".fq" "_guppy.fq" "_roadn.fa"; do
    # 定义SAM文件路径
    sam_file="sam_full_files/${i}${type}.sam"
    # 使用samtools flagstat，并提取所需的信息
    stats=$(samtools flagstat "$sam_file")
    total_reads=$(echo "$stats" | grep -oP '^\d+' | head -n 1)
    mapped_reads=$(echo "$stats" | grep -oP '^\d+ \+ \d+ mapped' | grep -oP '^\d+')
    
    # 计算映射百分比
    mapped_percentage=$(echo "scale=2; $mapped_reads*100/$total_reads" | bc)
    
    # 将文件名、总reads数和映射百分比写入到输出文件，格式为三列
    echo -e "${i}${type}\t$total_reads\t$mapped_percentage%" >> sam_full_files/flagstat_summary.txt
  done
done
```


## 拆分500条 accuracy
```{linux}
#### 每个样本每500条
## lab
cd /mnt/raid5/Personal/minion/DRS_mul/09gneral_model/04try_lla
virus_ref="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/ref/" 
rodan="/mnt/raid5/Personal/minion/software/RODAN/"

for i in  LabSARS PRRSV  PEDV SVV HomSap Pberghei; do 
 ref_name=$(echo $i | sed 's/LabSARS/SARS/')
  for type in ".fq" "_guppy.fq" "_roadn.fa"; do
    sam_file="sam_full_files/${i}${type}.sam"
    output_prefix="accu10files/${i}${type}"
    split_dir="sam_split_files/${i}${type}_split"  # 指定分割后文件的存储目录
    mkdir -p $split_dir  # 确保每个组合的分割目录都存在
    header_file="${sam_file}_header"
    grep "^@" $sam_file > $header_file
    total_lines=$(grep -v "^@" $sam_file | wc -l)
    ((lines_per_file = 500)) # 每个文件500条reads
    # 使用 split 命令按照500条reads分割文件，并将分割后的文件保存到指定目录
    grep -v "^@" $sam_file | split -l $lines_per_file -d -a 3 - ${split_dir}/${i}${type}_part_
    # 获取实际分割的文件数量，用于后续循环
    num_parts=$(ls ${split_dir}/${i}${type}_part_* | wc -l)
    # 对每个分割后的文件添加头部并运行accuracy.py
    for part in $(seq -f "%03g" 0 $((num_parts - 1))); do
      part_file="${split_dir}/${i}${type}_part_$part"
      cat $header_file $part_file > "${part_file}_with_header"
      python "$rodan/accuracy.py" "${part_file}_with_header" $virus_ref/${ref_name}.fasta > ${output_prefix}_part_${part}_accuracy.txt
    done
  done
done

## rodan data
cd /mnt/raid5/Personal/minion/DRS_mul/09gneral_model/04try_lla
virus_ref="/mnt/raid5/lla/RNA/data/ref/"
rodan="/mnt/raid5/Personal/minion/software/RODAN/"

for i in arabidopsis human mouse poplar yeast ; do 
  for type in ".fq" "_guppy.fq" "_roadn.fa"; do
    sam_file="sam_full_files/${i}${type}.sam"
    output_prefix="accu10files/${i}${type}"
    split_dir="sam_split_files/${i}${type}_split"  # 指定分割后文件的存储目录
    mkdir -p $split_dir  # 确保每个组合的分割目录都存在
    header_file="${sam_file}_header"
    grep "^@" $sam_file > $header_file
    total_lines=$(grep -v "^@" $sam_file | wc -l)
    ((lines_per_file = 500)) # 每个文件500条reads

    # 使用 split 命令按照500条reads分割文件，并将分割后的文件保存到指定目录
    grep -v "^@" $sam_file | split -l $lines_per_file -d -a 3 - ${split_dir}/${i}${type}_part_
    # 获取实际分割的文件数量，用于后续循环
    num_parts=$(ls ${split_dir}/${i}${type}_part_* | wc -l)

    # 对每个分割后的文件添加头部并运行accuracy.py
    for part in $(seq -f "%03g" 0 $((num_parts - 1))); do
      part_file="${split_dir}/${i}${type}_part_$part"
      cat $header_file $part_file > "${part_file}_with_header"
      python "$rodan/accuracy.py" "${part_file}_with_header" $virus_ref/${i}_reference.fasta > ${output_prefix}_part_${part}_accuracy.txt
    done
  done
done
```


## 特定物种训练


