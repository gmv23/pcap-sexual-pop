#!/usr/bin/bash

#This script runs Tassel pipeline to call SNPs on P. cap. sexual population

#ARGUMENTS
#1) tassel 5 standalone directory path
#2) genome directory path
#3) keyfile path
#4) path to directory with fastq files

tassel_path=$1
genome_dir_path=$2
keyfile_path=$3
fastq_dir_path=$4

#Copy key file to working directory

#First run GBSSeqToTagDBPlugin
$tassel_path/run_pipeline.pl -Xms50G -Xmx150G -fork1 -GBSSeqToTagDBPlugin \
-db capsici_gbsv2.db -e ApeKI -i $fastq_dir_path -k $keyfile_path \
-endPlugin -runfork1

#Run TagExportToFastqPlugin
$tassel_path/run_pipeline.pl -Xms50G -Xmx150G -fork1 -TagExportToFastqPlugin \
-db capsici_gbsv2.db -o capsici_tags_to_align.fa.gz -endPlugin -runfork1

#rename scaffolds on reference genome
cat $genome_dir_path/Phyca11_unmasked_genomic_scaffolds.fasta | sed 's/PHYCAscaffold_//g' > Phyca11_renamed.fasta

#rename scaffold in  Frank's mitochondrial assembly to 1000
cat $genome_dir_path/Mt_genome_from_Frank_martin_from_email.fasta | sed "s/>.*/>1000/" > mt_renamed.fasta

#combine files and remove all line breaks in sequences so line lengths match up
cat Phyca11_renamed.fasta mt_renamed.fasta | awk '!/^>/{printf "%s",$0;next}; 1' | \
sed '2,$s/>/\n>/' | sed '$s/$/\n/' > Pcap_complete_genome.fasta

#Index genome
bwa index -a bwtsw Pcap_complete_genome.fasta

#Align tags to genome
bwa aln -t 4 Pcap_complete_genome.fasta capsici_tags_to_align.fa.gz > capsici_tags_aligned.sai

#Convert to sam file
bwa samse Pcap_complete_genome.fasta capsici_tags_aligned.sai capsici_tags_to_align.fa.gz > capsici_tags_aligned.sam

#print out info on number of tags aligned per chromosome
samtools view -S -b capsici_tags_aligned.sam > capsici_tags_aligned.bam
samtools sort capsici_tags_aligned.bam sorted_capsici_tags_aligned
samtools index sorted_capsici_tags_aligned.bam
samtools idxstats sorted_capsici_tags_aligned.bam > number_tags_aligned_per_scaffold.txt

#Update database with tag alignment information (SAMTOGBSdbPlugin)
$tassel_path/run_pipeline.pl -fork1 -SAMToGBSdbPlugin -i \
capsici_tags_aligned.sam -db capsici_gbsv2.db -minMAPQ 30 -endPlugin -runfork1

#Run Discovery SNP Caller Plugin
$tassel_path/run_pipeline.pl -Xms20G -Xmx150G -fork1 -DiscoverySNPCallerPluginV2 -db \
capsici_gbsv2.db -endPlugin -runfork1

#Get taxa file that eliminates all water blanks to get quality info only with actual biological samples
grep -i -E -v "blank|\sH20\s" $keyfile_path | cut -f19 | \
tail -n +2  > non_blank_taxa.txt

#Run Quality Profiler on all biological samples
$tassel_path/run_pipeline.pl -fork1 -SNPQualityProfilerPlugin -db capsici_gbsv2.db \
-taxa non_blank_taxa.txt -tname "non_blank" -statFile capsici_snps_quality.txt -endPlugin -runfork1

#Run Production pipeline
$tassel_path/run_pipeline.pl -fork1 -ProductionSNPCallerPluginV2 -Xms50G -Xmx150G \
-db capsici_gbsv2.db -e ApeKI -i $fastq_dir_path -k $keyfile_path \
-o capsici_geno.vcf -endPlugin -runfork1


