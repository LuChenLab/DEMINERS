---
title: "Untitled"
author: "junwei"
date: "2024-03-26"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/")
library(vcfR)
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
library(UpSetR)
source("/mnt/raid61/Personal_data/songjunwei/AllScript/00Function/DIyplot.R")
setwd("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/06SARS_24sam/")

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


#4 mouse spec model
```{linux}
##01 fast5 make
conda activate base #2048
cd /mnt/raid5/Personal/minion/DRS_mul/09gneral_model/05mouse_model/01rawdata
for i in Train  Valid Test; do #
    mkdir -p /mnt/raid5/Personal/minion/DRS_mul/09gneral_model/05mouse_model/01rawdata/${i}_mouse_fast5/
done
    
for i in {1..5}; do
    #chmod 755 $i/*
    cp -r /mnt/raid5/lla/RNA/data/mouse/dataset/$i/ Train_mouse_fast5/
done
    cp -r /mnt/raid5/lla/RNA/data/mouse/dataset/6 Valid_mouse_fast5/
    cp -r /mnt/raid5/lla/RNA/data/mouse/dataset/0 Test_mouse_fast5/
    cp -r /mnt/raid5/lla/RNA/data/mouse/dataset/7 Test_mouse_fast5/
      cp -r /mnt/raid5/lla/RNA/data/mouse/dataset/8 Test_mouse_fast5/
      
##02 call
for i in Train  Valid ; do #
guppy_basecaller \
  -i ${i}_mouse_fast5 -r \
  --save_path ${i}_mouse_fast5/ --device cuda:0 \
  -c rna_r9.4.1_70bps_hac.cfg \
  --min_qscore 7  --verbose_logs  \
  --chunks_per_runner 256  --chunk_size 100 \
  --num_callers 32 --gpu_runners_per_device 32
cat ${i}_mouse_fast5/pass/* > ${i}_mouse_guppy.fq
cdup base
minimap2 -ax map-ont -t 30  /mnt/raid5/lla/RNA/data/ref/mouse_reference.fasta ${i}_mouse_guppy.fq |  samtools sort -@ 30 | samtools view -b > ${i}_mouse_guppy.bam
done

##03 tayaki
mkdir 0tayaki
tayaki="/mnt/raid5/Personal/minion/software/taiyaki/bin/"
for i in Train Valid ; do
python $tayaki/generate_per_read_params.py ${i}_mouse_fast5 --output 0tayaki/${i}_mouse_readid.tsv
python $tayaki/get_refs_from_sam.py /mnt/raid5/lla/RNA/data/ref/mouse_reference.fasta  ${i}_mouse_guppy.bam  --output 0tayaki/${i}_mouse_tayaki.fa --reverse
python $tayaki/prepare_mapped_reads.py ${i}_mouse_fast5 0tayaki/${i}_mouse_readid.tsv 0tayaki/${i}_mouse_tayaki.h5 \
        /mnt/raid5/Personal/minion/software/taiyaki/models/r941_rna_minion_upgrae.checkpoint     0tayaki/${i}_mouse_tayaki.fa --jobs 30 --overwrite 
done

##04 rodan
rodan="/mnt/raid5/Personal/minion/software/RODAN/"
cd 0tayaki
for i in Train Valid ;
do
mkdir -p rodan_${i}
python  $rodan/gendata.py -i ${i}_mouse_tayaki.h5 --outdir rodan_${i}
done
mv rodan_Train h5data_use
mv rodan_Valid/train.hdf5 h5data_use/valid.hdf5
rm -rf rodan_Valid

##05 species-specific model
### densecall 
cdup densecall
densecall train --config /mnt/raid5/Personal/minion/software/densecall/densecall/models/configs/rna_r9.4.1@v1.toml --data_dir h5data_use --workdir densecall_model
### rodan
#rodan="/mnt/raid5/Personal/minion/software/RODAN/"
#python  $rodan/model.py -c ./rodan_model/rna.config -a $rodan/rnaarch -n humRodan --arch $rodan/rnaarch --savedir ./rodan_model --rna -w 30 -l

##06 basecall
rodan="/mnt/raid5/Personal/minion/software/RODAN/"
#output="/mnt/raid61/Personal_data/songjunwei/DRS_RTA/10basecall_new/01fastq/"
output="/mnt/raid5/Personal/minion/DRS_mul/09gneral_model/04try_lla/fafiles"
for i in mouse; do
test_data="/mnt/raid5/Personal/minion/DRS_mul/09gneral_model/04try_lla/rodan_data/"
guppy_basecaller \
    -i $test_data/${i} -r \
    --save_path $output/${i}_guppy --device cuda:0 \
    -c rna_r9.4.1_70bps_hac.cfg \
    --min_qscore 7  --verbose_logs  \
    --chunks_per_runner 256  --chunk_size 100 \
    --num_callers 32 --gpu_runners_per_device 32
cat $output/${i}_guppy/pass/* > $output/${i}_guppy.fq
cdup base
python $rodan/basecall.py -m $rodan/rna.torch $test_data/${i} --arch $rodan/rnaarch -b 300 > $output/${i}_rodan.fq
cdup densecall
densecall basecaller --amp --output-file $output/${i}_densecall.fq rna_r9.4.1_hac@v1.0  $test_data/${i} -b 300
densecall basecaller --amp --output-file $output/${i}_densecall-spec1.fq rna_mouse_r9.4.1@v1.0  $test_data/${i} -b 300
densecall basecaller --amp --output-file $output/${i}_densecall-spec2.fq rna_mouse_r9.4.1@v2.0  $test_data/${i} -b 300
done

##07 accuracy
cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/10basecall_new/04accuracy/0class_genome
ref="/mnt/raid5/lla/RNA/data/ref/mouse_reference.fasta" ## 标准参考
for i in guppy rodan densecall densecall-spec1 densecall-spec2; do
  minimap2 --secondary=no -ax map-ont -t 32 --cs $ref /mnt/raid5/Personal/minion/DRS_mul/09gneral_model/04try_lla/fafiles/mouse_${i}.fq > mouse_${i}.sam
  python "/mnt/raid5/Personal/minion/software/RODAN/accuracy.py" mouse_${i}.sam $ref > mouse_${i}_accuracy.txt
done

cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/10basecall_new/04accuracy/1mouse_bam
for i in guppy rodan densecall densecall-spec1 densecall-spec2; do
  minimap2 --secondary=no -ax map-ont -t 32 --cs /mnt/raid61/Personal_data/songjunwei/reference/mouse/GRCm38/release101/Mus_musculus.GRCm38.dna.primary_assembly.fa \
  mouse_${i}.fq | samtools sort -o mouse_${i}_sorted.bam
  samtools index mouse_${i}_sorted.bam
done

#06 basecall 4000*3=1.2w

rodan="/mnt/raid5/Personal/minion/software/RODAN/"
output="/mnt/raid5/Personal/minion/DRS_mul/09gneral_model/04try_lla/mouse_spe_fa"
test_data="/mnt/raid5/Personal/minion/DRS_mul/09gneral_model/05mouse_model/01rawdata/Test_mouse_fast5"
i="mouse"
guppy_basecaller \
    -i $test_data/ -r \
    --save_path $output/${i}_guppy --device cuda:0 \
    -c rna_r9.4.1_70bps_hac.cfg \
    --min_qscore 7  --verbose_logs  \
    --chunks_per_runner 256  --chunk_size 100 \
    --num_callers 32 --gpu_runners_per_device 32
cat $output/${i}_guppy/pass/* > $output/${i}_guppy.fq
cdup base
python $rodan/basecall.py -m $rodan/rna.torch $test_data/ --arch $rodan/rnaarch -b 300 > $output/${i}_rodan.fq
cdup densecall
densecall basecaller --amp --output-file $output/${i}_densecall.fq rna_r9.4.1_hac@v1.0  $test_data/ -b 300
densecall basecaller --amp --output-file $output/${i}_densecall-spec.fq rna_mouse_r9.4.1@v1.0  $test_data/ -b 300

cdup base
ref="/mnt/raid5/lla/RNA/data/ref/mouse_reference.fasta" ## 标准参考
for i in guppy rodan densecall densecall-spec; do
  minimap2 --secondary=no -ax map-ont -t 32 --cs $ref $output/mouse_${i}.fq > $output/mouse_${i}.sam
  python "/mnt/raid5/Personal/minion/software/RODAN/accuracy.py" $output/mouse_${i}.sam $ref > $output/mouse_${i}_accuracy.txt
done
```


准确度分析
```{r fig.width=10,fig.height=6}
df_accu <- lapply(list.files("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/10basecall_new/04accuracy/0class_genome/",pattern = "_accuracy.txt",recursive = T,full.names = T), function(file){
  df <- accuracy_to_df(file,epoch = gsub("_accuracy.txt", "\\1", basename(file)))
  return(df)
}) %>% rbindlist()
df_accu <- df_accu[grep("mouse",df_accu$Epoch),]
df_accu$Mapping.Rate <- df_accu$Mapp_reads/4000
df_accu$Data.name <- "mouse"
df_accu$Basecaller <- c("DEMINERS","DEMINERS-specific1","DEMINERS-specific2","Guppy","RODAN")
df_accu_res <- df_accu[,c("Data.name","Basecaller","Mapping.Rate","Accuracy_Med","Length_Med","Std","Mismatch_Med","Deletions_Med","Insertions_Med")]
df_accu_res[,c(3:4,6:9)] <- round(df_accu_res[,c(3:4,6:9)]*100,2)
openxlsx::write.xlsx(df_accu_res,"/mnt/raid61/Personal_data/songjunwei/DRS_RTA/09general_model/densecall_mouse_specific_0410.xlsx")

df_accu_res <- df_accu_res[-3,]
df_accu_res[2,2] <- "DEMINERS-specific"
df_accu_res$Basecaller <- factor(df_accu_res$Basecaller,levels = c("DEMINERS-specific","DEMINERS","RODAN","Guppy"))
df_accu_pdata <- df_accu_res %>% gather(key = "Parameter", value = "Value", -c(`Data.name`, Basecaller)) %>% data.table()

basecaller_colors <- c("DEMINERS-specific" = "#984EA3", "DEMINERS" = "#E41A1C", "RODAN" = "#4DAF4A", "Guppy" = "#377EB8")

metrics <- c("Mapping.Rate", "Accuracy_Med", "Length_Med", "Mismatch_Med", "Deletions_Med", "Insertions_Med")

plot_list <- lapply(metrics, function(metric) {
  gg <- ggplot(df_accu_res, aes(x = Basecaller, y = .data[[metric]], fill = Basecaller)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(label = sprintf("%.2f", .data[[metric]])), vjust = -0.3, position = position_dodge(0.9)) +
    scale_fill_manual(values = basecaller_colors) +
    theme_minimal() + mytheme + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = paste0(gsub("_Med","",metric)," (%)"), x = NULL)
  return(gg)
})
ggarrange(plotlist = plot_list,common.legend = T)

ggsave("../09general_model/desecall_compare_bar.pdf",width = 10,height = 6)
```



# 11 densecall 统计
## 所有样本 accu table
```{r fig.width=6,fig.height=5}
select <-list.files("/mnt/raid5/Personal/minion/DRS_mul/09gneral_model/04try_lla/accufiles",
                    pattern = "_accuracy.txt",recursive = T,full.names = T)
df <- lapply(select, function(file){
  df <- accuracy_to_df(file,epoch = gsub("_accuracy.txt", "\\1", basename(file)))
  return(df)
}) %>% rbindlist()
tmp <- str_split_fixed(gsub(".fa","",gsub(".fq","",df$Epoch)),"_",2)
df$Epoch <- tmp[,2]
df$species <- tmp[,1]
df$Epoch[df$Epoch==""] <- "densecall"

df$Data_from <- ifelse(df$species %in% c("arabidopsis","human","mouse","yeast","poplar"),"Rodan","lab")
df$mapping_rate <- df$Mapp_reads/4000
write.csv(df,"../09general_model/rodan_5species.csv",row.names = F)


virus <-list.files("/mnt/raid5/Personal/minion/DRS_mul/09gneral_model/05virus_densecall/02res/",
                    pattern = "_accuracy.txt",recursive = T,full.names = T)
df_virus <- lapply(virus, function(file){
  df <- accuracy_to_df(file,epoch = gsub("_accuracy.txt", "\\1", basename(file)))
  return(df)
}) %>% rbindlist()
df_virus$species <- df_virus$Epoch
df_virus$Epoch <- "densecall_virus"
df_m <- rbind(df,df_virus)

test_num <- data.frame(
  species = c("SARS", "HomSap", "Pberghei", "SVV", "PRRSV", "LabSARS","PEDV"),
  test = c(2000, 2000, 2000, 1983, 1563, 406,11)
)
new <- merge(df,test_num,by="species")
new$mapping_rate <- new$Mapp_reads/new$test
write.csv(new,"../09general_model/LabSARS_guppy_rodan_res.csv",row.names = F)
```


```{r fig.width=10,fig.height=3}
# save 后 input
densecall <- data.table::data.table(openxlsx::read.xlsx("../09general_model/densecall_compared_resV2_0325.xlsx"))

densecall$Basecaller <- gsub("densecall","DEMINERS",densecall$Basecaller)
densecall$Basecaller <- gsub("guppy","Guppy",densecall$Basecaller)
densecall$Basecaller <- gsub("rodan","RODAN",densecall$Basecaller)
densecall[,c(3,4,6:9)] <- round(densecall[,c(3,4,6:9)]*100,2)
densecall$Basecaller <- factor(densecall$Basecaller,levels = c("DEMINERS_virus","DEMINERS","RODAN","Guppy"))

densecall <- densecall[order(densecall$Data.from,densecall$Data.name,densecall$Basecaller)]
openxlsx::write.xlsx(densecall,"../09general_model/densecall_compared_resV2_0325_new.xlsx")

```

## 拆分10个 
```{r fig.width=6,fig.height=5}
split10 <- lapply(list.files("/mnt/raid61/Personal_data/songjunwei/DRS_RTA/09general_model/accu10files/",pattern = "_accuracy.txt",recursive = T,full.names = T), function(file){
  df <- accuracy_to_df(file,epoch = gsub("_accuracy.txt", "\\1", basename(file)))
  return(df)
}) %>% rbindlist()
split10 <- na.omit(split10)
split10$Epoch <- gsub(".fa","",split10$Epoch)
split10$Epoch <- gsub(".fq","",split10$Epoch)

tmp <- str_split_fixed(split10$Epoch,"_",3)
split10$species <- tmp[,1]
split10$Epoch <- tmp[,2]
split10$part <- tmp[,3]

split10$part <- gsub("part_","",split10$part)
split10$Epoch[grep("guppy",split10$Epoch)] <- "Guppy"
split10$Epoch[grep("roadn",split10$Epoch)] <- "RODAN"
split10$Epoch[split10$Epoch=="part"] <- "DEMINERS"
write.csv(split10,"../09general_model/general_compare_split10.csv",row.names = F)


```

## plot
### boxplot 主图
```{r}
densecall <- data.table::data.table(openxlsx::read.xlsx("../09general_model/densecall_compared_resV2_0325_new.xlsx"))
data_long <- densecall %>% gather(key = "Parameter", value = "Value", -c(`Data.name`, Basecaller, `Data.from`)) %>% data.table()
prarmeter <- c("Accuracy","Mismatch","Deletions","Insertions") #c("Accuracy","Mapping.Rate","Length","Mismatch","Deletions","Insertions","Std")
data_long <- data_long[Parameter %in% prarmeter]
data_long <- data_long[Basecaller!="DEMINERS_virus"]
data_long$Parameter <- factor(data_long$Parameter,levels = prarmeter)
data_long$Basecaller <- factor(data_long$Basecaller,levels = c("DEMINERS","RODAN","Guppy"))

medians <- data_long[,.(median=median(Value)),by=c("Parameter", "Basecaller")]

# Map the basecallers to the extracted colors
boxtheme <- theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
basecaller_colors <- c("DEMINERS_virus" = "#984EA3", "DEMINERS" = "#E41A1C", "RODAN" = "#4DAF4A", "Guppy" = "#377EB8")


ggplot(data_long, aes(x = Basecaller, y = Value, fill = Basecaller)) +
  geom_boxplot() +
  geom_text(data = medians, aes(label = sprintf("%.2f", median), y = median + 0.0), 
            position = position_dodge(width = 0.75), vjust = 0, size = 3) +
  facet_wrap(~ Parameter, scales = "free", ncol = 4) +
  theme(axis.text.x= element_text(size=12,colour = "black",angle = 45,hjust = 1))+
  theme_minimal()+ boxtheme + stat_compare_means(method = 'wilcox.test', 
                             comparisons = list(c("Guppy", "RODAN"), c("DEMINERS", "Guppy")))+
  labs(title = "Basecaller Performance Comparison", x = "Basecaller", y = "Value")
ggsave("../09general_model/desecall_compare_box_4pic_v2.pdf",width = 10,height = 3)
# ggsave("../09general_model/desecall_compare_box_SARS-.pdf",width = 9,height = 8)


## 物种特异性bar图
# Generate the plot with the specified colors and add value labels on top of the bars
ggplot(densecall[Data.name %in% c("SARS","PRRSV","SVV")], aes(x = Basecaller, y = Accuracy, fill = Basecaller)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = sprintf("%.2f", Accuracy)), vjust = -0.3, position = position_dodge(0.9)) +
  scale_fill_manual(values = basecaller_colors) +
  facet_wrap(~ Data.name, scales = "free_x") +
  theme_minimal() + mytheme + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Accuracy (%)", x = "", title = "Accuracy Comparison by Basecaller and Data.name")
ggsave("../09general_model/desecall_compare_bar.pdf",width = 10,height = 3)
```


```{r fig.width=10,fig.height=10}
ggplot(densecall, aes(x = Basecaller, y = Accuracy, fill = Basecaller)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Data.name, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Accuracy", x = "", title = "Accuracy Comparison by Basecaller and Data.name")
ggsave("../09general_model/desecall_compare_bar.pdf",width = 10,height = 10)
```


```{r fig.width=10,fig.height=3}
split10 <- data.table(read.csv("../09general_model/general_compare_split10.csv"))
#split10_lab <- split10[split10$species %in% c("HomSap","LabSARS","Pberghei","PRRSV","SVV","PEDV")]
#split10_rodan <- split10[!split10$species %in% c("HomSap","LabSARS","Pberghei","PRRSV","SVV")]

data_long <- split10 %>% gather(key = "Parameter", value = "Value", -c(species,Epoch,part)) %>% data.table()
table(data_long$Epoch)
data_long <- data_long[Parameter %in% c("Accuracy_Med","Mismatch_Med","Deletions_Med","Insertions_Med")]
data_long[Parameter!="Length"]$Value <- data_long[Parameter!="Length"]$Value *100
data_long$Parameter <- factor(data_long$Parameter,levels = c("Accuracy_Med","Mismatch_Med","Deletions_Med","Insertions_Med"))
data_long$Epoch <- factor(data_long$Epoch,levels = c("DEMINERS","RODAN","Guppy"))


medians <- data_long[,.(median=median(Value)),by=c("Parameter", "Epoch")]
ggplot(data_long, aes(x = Epoch, y = Value, fill = Epoch)) +
  geom_boxplot() +
  geom_text(data = medians, aes(label = sprintf("%.2f", median), y = median + 0.0), 
            position = position_dodge(width = 0.75), vjust = 0, size = 3) +
  facet_wrap(~Parameter, scales = "free", ncol = 4) +
  theme(axis.text.x= element_text(size=16,colour = "black",angle = 45,hjust = 1),guide_legend(title.position =  "none"))+ 
 scale_fill_manual(values = basecaller_colors) +
  theme_minimal() + boxtheme + stat_compare_means(method = 'wilcox.test', 
                             comparisons = list(c("DEMINERS", "RODAN"), c("DEMINERS", "Guppy")))+
  labs(title = "", x = "", y = "Performance%")
ggsave("../09general_model//desecall_compare_box_split10.pdf",width = 10,height = 3)
```


```{r fig.width=10,fig.height=3.5}
rodan_spe <- c("arabidopsis","human","poplar","yeast","mouse")
Accu_med <- data_long[Parameter %in% "Accuracy_Med"][ species %in% rodan_spe] #c("arabidopsis","human","poplar","yeast","HomSap","PRRSV","SVV","mouse","Pberghei")
Accu_med$species <- factor(Accu_med$species,rodan_spe)

medians <-Accu_med[,.(median=median(Value)),by=c("species", "Epoch")]
ggplot(Accu_med, aes(x = Epoch, y = Value, fill = Epoch)) +
  geom_boxplot() +
  geom_text(data = medians, aes(label = sprintf("%.2f", median), y = median ), 
            position = position_dodge(width = 0.75), vjust = 0, size = 3) +
  facet_wrap(~ species, scales = "free", ncol = 5) +
  scale_fill_manual(values = basecaller_colors) + boxtheme +  theme(axis.text.x= element_text(size=8,colour = "black",angle = 45,hjust = 1))+ 
  stat_compare_means(method = 'wilcox.test',  label = "p.format",label.x = 1.5, comparisons = list(c("DEMINERS", "RODAN"), c("DEMINERS", "Guppy")))+
  labs(title = "", x = "Basecaller", y = "Accuracy percentage%")
ggsave("../09general_model//desecall_compare_box_split10_speices.pdf",width = 10,height = 3.5)
```





