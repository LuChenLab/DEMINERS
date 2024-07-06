# Alignmennt:
minimap2 -ax map-ont -t 10 /mnt/raid61/Personal_data/tangchao/Document/gencode/human/GRCh37/GRCh37.primary_assembly.genome.fa /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/Tumor4/DMG_1.fq | samtools sort -@ 4 | samtools view -b > /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/DMG_1.sorted.bam

minimap2 -ax map-ont -t 10 /mnt/raid61/Personal_data/tangchao/Document/gencode/human/GRCh37/GRCh37.primary_assembly.genome.fa /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/Tumor4/DMG_2.fq | samtools sort -@ 4 | samtools view -b > /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/DMG_2.sorted.bam

minimap2 -ax map-ont -t 10 /mnt/raid61/Personal_data/tangchao/Document/gencode/human/GRCh37/GRCh37.primary_assembly.genome.fa /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/Tumor4/DMG_3.fq | samtools sort -@ 4 | samtools view -b > /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/DMG_3.sorted.bam

minimap2 -ax map-ont -t 10 /mnt/raid61/Personal_data/tangchao/Document/gencode/human/GRCh37/GRCh37.primary_assembly.genome.fa /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/Tumor4/GBM_6.fq | samtools sort -@ 4 | samtools view -b > /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/GBM_6.sorted.bam

samtools index /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/DMG_1.sorted.bam
samtools index /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/DMG_2.sorted.bam
samtools index /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/DMG_3.sorted.bam
samtools index /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/GBM_6.sorted.bam


# Mutation calling:

samtools mpileup -g -f /mnt/raid61/Personal_data/tangchao/Document/gencode/human/GRCh37/GRCh37.primary_assembly.genome.fa /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/DMG_1.sorted.bam | bcftools call -mv -Ob -o /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/DMG_1.vcf

samtools mpileup -g -f /mnt/raid61/Personal_data/tangchao/Document/gencode/human/GRCh37/GRCh37.primary_assembly.genome.fa /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/DMG_2.sorted.bam | bcftools call -mv -Ob -o /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/DMG_2.vcf

samtools mpileup -g -f /mnt/raid61/Personal_data/tangchao/Document/gencode/human/GRCh37/GRCh37.primary_assembly.genome.fa /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/DMG_3.sorted.bam | bcftools call -mv -Ob -o /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/DMG_3.vcf

samtools mpileup -g -f /mnt/raid61/Personal_data/tangchao/Document/gencode/human/GRCh37/GRCh37.primary_assembly.genome.fa /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/GBM_6.sorted.bam | bcftools call -mv -Ob -o /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/GBM_6.vcf


bcftools view DMG_1.vcf > DMG_1_plain.vcf
bcftools view DMG_2.vcf > DMG_2_plain.vcf
bcftools view DMG_3.vcf > DMG_3_plain.vcf
bcftools view GBM_6.vcf > GBM_6_plain.vcf

cd /mnt/raid61/Personal_data/songjunwei/DRS_RTA/03DecodeR_res/align/
samtools mpileup -g -f /mnt/raid61/Personal_data/tangchao/Document/gencode/human/GRCh37/GRCh37.primary_assembly.genome.fa DMG_1.sorted.bam DMG_2.sorted.bam DMG_3.sorted.bam GBM_6.sorted.bam | bcftools call -mv -Ob -o Tumor.vcf

