minimap2 -V
# 2.17-r941

samtools --version
# samtools 1.11

cd /mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/BaseCall/20201127/guppy
cat fast5_pass/*fastq > fast5_pass.fastq

minimap2 -ax map-ont -t 10 /mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/Reference/RNAsequence.fa fast5_pass.fastq | samtools sort -@ 4 | samtools view -b > /mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/alignment/20201127/minimap2/20201127.sorted.bam

cd /mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/alignment/20201127/minimap2
samtools index 20201127.sorted.bam
