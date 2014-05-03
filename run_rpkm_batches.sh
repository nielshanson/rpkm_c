#!/bin/bash

while read line
do

folder=$line\_combined_unique
assembly=$line\_combined_unique.fasta
sam_se=$line\_combined_unique_aln.se.sam
sam_pe=$line\_combined_unique_aln.pe.sam
gff=$line\_combined_unique.unannot.gff
pgdb=e$line\_combined_unique.basepathways.txt
output=$line.basepathways_rpkm.txt
stdout=$line\_rpkm_stdout.txt
stderr=$line\_rpkm_stderr.txt

./rpkm -c /Users/mayab3/Desktop/Hallam_Projects/NESAP/jgi_illumina/MetaPathways_2_0/output/Mar27_2014/$folder/preprocessed/$assembly --r1 /Volumes/3TB/Hallam_Projects2/NESAP/jgi_illumina/Metagenomes/BWA_aligned_SAM_files_Apr2014/$sam_se --r2 /Volumes/3TB/Hallam_Projects2/NESAP/jgi_illumina/Metagenomes/BWA_aligned_SAM_files_Apr2014/$sam_pe -O /Users/mayab3/Desktop/Hallam_Projects/NESAP/jgi_illumina/MetaPathways_2_0/output/Mar27_2014/$folder/orf_prediction/$gff -p /Users/mayab3/Desktop/Hallam_Projects/NESAP/jgi_illumina/MetaPathways_2_0/output/Mar27_2014/pgdbs/$pgdb --multireads --format sam-2 -o /Users/mayab3/Desktop/Hallam_Projects/NESAP/jgi_illumina/MetaPathways_2_0/output/Mar27_2014/pgdbs/RPKM_tables/$output $1>/Users/mayab3/Desktop/Hallam_Projects/NESAP/jgi_illumina/MetaPathways_2_0/output/Mar27_2014/pgdbs/RPKM_tables/$stdout $2>/Users/mayab3/Desktop/Hallam_Projects/NESAP/jgi_illumina/MetaPathways_2_0/output/Mar27_2014/pgdbs/RPKM_tables/$stderr

done<completed_bwa.list
