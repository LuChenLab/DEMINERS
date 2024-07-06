minimap2 -V
# 2.17-r941

samtools --version
# samtools 1.11

cd /mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/BaseCall

for i in 20200605 20200620 20200902;
do 
  minimap2 -ax map-ont -t 10 /mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/Reference/RNAsequence_$i.fa ./$i/guppy/fast5_pass/*fastq | samtools sort -@ 4 | samtools view -b > /mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/alignment/$i/minimap2/$i.sorted.bam
done

for i in 20200605 20200620 20200902;
do 
  samtools index /mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/alignment/$i/minimap2/$i.sorted.bam
done


# 20210703

cd /mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex
minimap2 -ax map-ont -t 10 ./data/Reference/RNAsequence_20210703.fa ./data/BaseCall/20210703/guppy/fast5_pass/pass/*fastq | samtools sort -@ 4 | samtools view -b > ./data/alignment/20210703/minimap2/20210703.sorted.bam
