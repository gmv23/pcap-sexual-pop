#!/usr/bin/bash

#Get depth by individuals so we can filter out failed samples
vcftools --vcf capsici_geno.vcf --depth --out full_depth_by_individual

#Which samples have <1 mean read depth from my plate?
grep "CB7GTANXX" full_depth_by_individual.idepth | grep -E "\s0\."

#It appears that 7 samples failed (from plate CB7GTANXX):
# The H20 blank
# 15PF78A
# 15PF77B
# 15PF13C
# 15FP01A
# 15PF33A
# 15PF48A

#Just keep blight farm isolates and parents 
cat full_depth_by_individual.idepth | awk -F: '$1 ~ /PF|664|6180/ {print}' | cut -f1 > pf_and_parents.txt
vcftools --vcf capsici_geno.vcf --keep pf_and_parents.txt --recode --out step1_pf_and_parents

#Look at missingness by individual
vcftools --vcf step1_pf_and_parents.recode.vcf --missing-indv --out pf_and_parents

#Print average %missingness
awk '{ total += $5 } END { print total/NR }' pf_and_parents.imiss
#How many samples have greater than 40% missingness?
awk '$5 > 0.4 { count++ } END { print count }' pf_and_parents.imiss
#What samples are they?
awk '$5 > 0.4 {print}' pf_and_parents.imiss | sort -k 5

#Get rid of samples with more than %40 missingness (is 40 a good number)?
vcftools --vcf step1_pf_and_parents.recode.vcf --remove <(tail pf_and_parents.imiss -n +2 | awk '$5 > 0.4 {print $1}' \
pf_and_parents.imiss) --recode --out step2_low_missing_inds

#Get rid of indels
vcftools --vcf step2_low_missing_inds.recode.vcf --remove-indels --recode --out step3_no_indels

#Get rid of non biallelic SNPs 
vcftools --vcf step3_no_indels.recode.vcf --min-alleles 2 --max-alleles 2 --recode --out step4_biallelic

#Exclude genotypes with fewer than 5 reads
vcftools --vcf step4_biallelic.recode.vcf --minDP 5 --recode --out step5_high_depth

#Look at missingness on a per site basis
vcftools --vcf step5_high_depth.recode.vcf --missing-site --out site_missingness

#Exclude SNPs with >20% missing data
vcftools --vcf step5_high_depth.recode.vcf --max-missing 0.8 --recode --out step6_low_missing

#Look at mean site depth
vcftools --vcf step6_low_missing.recode.vcf --site-mean-depth --out site_mean_depth

#Get rid of sites with mean read depth <10 or >50
vcftools --vcf step6_low_missing.recode.vcf --min-meanDP 10 --max-meanDP 50 --recode --out step7_good_read_depth

#Look at MAF distribution
vcftools --vcf step7_good_read_depth.recode.vcf --freq --out maf

#Filter based on MAF > 0.05
vcftools --vcf step7_good_read_depth.recode.vcf --maf 0.05 --recode --out step8_good_maf

#Get rid of mitochondrial chromosome (for now)
vcftools --vcf step8_good_maf.recode.vcf --not-chr 1000 --recode --out genomic_pass_initial_filter 

#Get 012 files
vcftools --vcf genomic_pass_initial_filter.recode.vcf --012 --out capPF_filtered 
