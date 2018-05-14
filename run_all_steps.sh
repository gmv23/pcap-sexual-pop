#!/usr/bin/bash

# This script will make a folder with the date and
#1) run tassel pipeline
#2) run filtering steps
#3) get allele depth ratio file and estimate ploidies, removing non-diploid individuals
#4) get pairwise IBS matrix and perform clone correction
#5) Final filtering

#Or you can provide the date as an argument and it will supply what is missing from a folder started on that date

#Stop on errors
set -e

SCRIPTS='/home/greg/sexual_pop/scripts'

if [ -z $1 ]; then
	DATE=$(date +%F)
else
	DATE=$1
	echo "Date received"
fi

if [ ! -d pipeline_$DATE ]; then
	mkdir pipeline_$DATE
fi

cd pipeline_$DATE
#Copy current version of this script to directory
cp "$SCRIPTS"/run_all_steps.sh .

if [ ! -d tassel_run ]; then
	mkdir tassel_run
	cd tassel_run
	cp "$SCRIPTS"/run_discovery_pipeline.sh .
	cp "$SCRIPTS"/combine_key_files.py .
	cp "$SCRIPTS"/merge_sample_names.py .
	#Combine key files
	python combine_key_files.py /home/greg/sexual_pop/sequencekeyfiles/ combined_key_file.txt
	#Hard code to fix INCORRECTLY labeled samples (12PF15A and 16A))}1
	awk '$1 == "C4B1JACXX" && $4 == "12PF_16A" {gsub(/12PF_16A/,"12PF_15A")}1' combined_key_file.txt > combined_key_file_fix1.txt
        awk '$1 == "C6P86ANXX" && $4 == "12PF_15A" {gsub(/12PF_15A/,"12PF_16A")}1' combined_key_file_fix1.txt > combined_key_file_fix2.txt
	#Rename duplicate samples so TASSEL will combined their reads
	python merge_sample_names.py combined_key_file_fix2.txt merged_key_file.txt
	bash run_discovery_pipeline.sh /home/greg/software/tassel-5-standalone \
	/home/greg/sexual_pop/genome merged_key_file.txt /home/greg/sexual_pop/raw_reads > run_discovery_pipeline.log
	cd ..
fi

if [ ! -d filtering ]; then
	mkdir filtering
	cd filtering
	cp "$SCRIPTS"/filter_vcf.sh .
	cp ../tassel_run/capsici_geno.vcf .
	bash filter_vcf.sh > filter_vcf.log
	cd ..
fi

if [ ! -d clone_correction ]; then
        mkdir clone_correction
        cd clone_correction
        cp "$SCRIPTS"/get_pairwise_IBS.R .
        cp ../filtering/capPF_filtered* .
        Rscript get_pairwise_IBS.R > get_pairwise_IBS.log
        cp "$SCRIPTS"/Clone_correction_and_trends.R .
        mkdir plots
        mkdir filtered_data
        Rscript Clone_correction_and_trends.R > Clone_correction.log
        cd ..
fi

if [ ! -d ploidy_variation ]; then
	mkdir ploidy_variation
	cd ploidy_variation
	cp "$SCRIPTS"/get_genotype_depth.py .
	cp ~/sexual_pop/lg_info/* .
	cp "$SCRIPTS"/Ploidy_by_lg.R .
	cp ../clone_correction/filtered_data/Clonal_groups_including_0664_reintros.csv . 
	cp ../filtering/capPF_filtered.012.indv .
	mkdir plots
	mkdir plots/all_isolates
	mkdir plots/all_isolates_by_lg
	python get_genotype_depth.py ../filtering/genomic_pass_initial_filter.recode.vcf depth_ratios.txt 15 > depth_ratios.log
	Rscript Ploidy_by_lg.R > Ploidy_by_lg.log

	cd ..
fi

if [ ! -d final_filtering ]; then
	mkdir final_filtering
	cd final_filtering
	cp "$SCRIPTS"/get_genotype_depth.py .
	cp ../filtering/genomic_pass_initial_filter.recode.vcf .

	#Just keep clone corrected progeny and all parental reps for now and remove non diplod isolates
	cat <( grep -E "664|6180" ../filtering/capPF_filtered.012.indv ) ../clone_correction/filtered_data/capPF_cc.012.indv > isolates_to_keep.txt 
	vcftools --vcf genomic_pass_initial_filter.recode.vcf --keep isolates_to_keep.txt --recode --out capPF_cc_parental_reps
	vcftools --vcf capPF_cc_parental_reps.recode.vcf --remove ../ploidy_variation/non_diploid_isolates.txt --recode --out capPF_cc_parental_reps_dip
	vcftools --vcf capPF_cc_parental_reps_dip.recode.vcf --012 --out capPF_cc_parental_reps_dip

	#Get depth ratios on new vcf file
	python get_genotype_depth.py capPF_cc_parental_reps_dip.recode.vcf cc_diploid_depth_ratios.txt 4

	#Make file of scaffold lengths
	paste <(grep ">" ../tassel_run/Pcap_complete_genome.fasta) <(cat ../tassel_run/Pcap_complete_genome.fasta | awk '/^[T|A|G|C|N]/ {print length($0)}') \
	> scaffold_lengths.txt
	sed -i 's/>//g' scaffold_lengths.txt

	#Run script
	cp "$SCRIPTS"/Final_filtering.R .
	Rscript Final_filtering.R > Final_filtering.log

	cd ..
fi

cd ..
