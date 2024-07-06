# Severe:1099
guppy_basecaller -v
# : Guppy Basecalling Software, (C) Oxford Nanopore Technologies, Limited. Version 4.2.3+8aca2af8

guppy_basecaller -i /mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/Rawdata/20201127/20201126/no_sample/20201126_1357_MN26652_FAO61922_e7fbd202/fast5_pass \
                 -s /mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/BaseCall/20201127/guppy/fast5_pass \
                 --u_substitution off \
                 --device "cuda:0" \
                 -c rna_r9.4.1_70bps_hac.cfg

# Severe:1099
guppy_basecaller -v
# : Guppy Basecalling Software, (C) Oxford Nanopore Technologies, Limited. Version 4.2.3+8aca2af8

guppy_basecaller -i /mnt/raid62/BetaCoV/RawData/20200605/20200605_0941_MN29097_FAM92974_9ad8351b/fast5_pass \
                 -s /mnt/raid64/Nanopore_RNA_directly_Sequencing/analysis/00.GuppyBaseCall/20200605/guppy/fast5_pass \
                 --u_substitution off \
                 --device "cuda:0" \
                 -c rna_r9.4.1_70bps_hac.cfg
mv -R /mnt/raid64/Nanopore_RNA_directly_Sequencing/analysis/00.GuppyBaseCall/20200605/guppy/fast5_pass \
   /mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/BaseCall/20200605/guppy/fast5_pass

# Severe:1099
guppy_basecaller -v
# : Guppy Basecalling Software, (C) Oxford Nanopore Technologies, Limited. Version 4.2.3+8aca2af8

guppy_basecaller -i /mnt/raid62/BetaCoV/RawData/20200620/20200620_0916_MN27328_FAM92901_ea10d670/fast5_pass \
                 -s /mnt/raid64/Nanopore_RNA_directly_Sequencing/analysis/00.GuppyBaseCall/20200620/guppy/fast5_pass \
                 --u_substitution off \
                 --device "cuda:0" \
                 -c rna_r9.4.1_70bps_hac.cfg
mv /mnt/raid64/Nanopore_RNA_directly_Sequencing/analysis/00.GuppyBaseCall/20200620/guppy/fast5_pass \
   /mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/BaseCall/20200620/guppy/fast5_pass


# Severe:1099
guppy_basecaller -v
# : Guppy Basecalling Software, (C) Oxford Nanopore Technologies, Limited. Version 4.2.3+8aca2af8

guppy_basecaller -i /mnt/raid64/Nanopore_RNA_directly_Sequencing/data/RawData/20200902RNA/20200902_1234_MN27328_FAL79918_bc368336/fast5_pass \
                 -s /mnt/raid64/Nanopore_RNA_directly_Sequencing/analysis/00.GuppyBaseCall/20200902/guppy/fast5_pass \
                 --u_substitution off \
                 --device "cuda:0" \
                 -c rna_r9.4.1_70bps_hac.cfg
mv /mnt/raid64/Nanopore_RNA_directly_Sequencing/analysis/00.GuppyBaseCall/20200902/guppy/fast5_pass \
   /mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/BaseCall/20200902/guppy/fast5_pass

# Severe:2048
guppy_basecaller -v
# : Guppy Basecalling Software, (C) Oxford Nanopore Technologies, Limited. Version 5.0.11+2b6dbff

guppy_basecaller -i /mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/Rawdata/20210703/no_sample/20210703_1132_MN26652_FAP34145_f2607660/fast5_pass/ \
                 -s /mnt/raid5/Personal/tangchao/project/Nanopore/BarcodeDecomplex/data/BaseCall/20210703/guppy/fast5_pass \
                 --u_substitution off \
                 --device "cuda:0" \
                 -c rna_r9.4.1_70bps_hac.cfg

