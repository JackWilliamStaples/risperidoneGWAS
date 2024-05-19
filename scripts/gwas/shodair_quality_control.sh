# make sure the path points to the bioinformatics software directory
PATH=$PATH:~/bioinformatics_software/

# create a working directory for the quality control if it does not exist
if [[ -d ~/Psychiatric_GWAS/Psychiatric_GWAS/quality_control_working_directory ]]
then
  echo "Quality control working directory exists"
else
  echo "Created working directory for quality control"
  mkdir ~/Psychiatric_GWAS/Psychiatric_GWAS/quality_control_working_directory/
fi

# The majority of the commands and R scripts for this QC analysis were adapted from the manuscript: 
  # Marees, et al. (2018). International journal of methods in psychiatric research, 
  # 27(2), e1608. https://doi.org/10.1002/mpr.1608 (see Maurees github page: https://github.com/MareesAT/GWA_tutorial/)

# working directory for this is
cd ~/Psychiatric_GWAS/Psychiatric_GWAS/quality_control_working_directory/

# create a text file for removing participant SCH-615 because they were not taking risperidone
if [[ -f ./no_risperidone.txt ]]
then
  echo "no_risperidone.txt file already exists"
else 
  echo "no_risperidone.txt file does not exist"
  echo "8 SCH-615" > no_risperidone.txt
  echo "no_risperidone.txt file created"
fi

# create a text file for updating the gender for participant SCH-937 and SCH-254 
# (correct gender for both participants is male = 1)
if [[ -f ./update_sex.txt ]]
then
  echo "update_sex.txt file already exists"
else
  echo "update_sex.txt file does not exist"
  echo "93 SCH-937 1" > update_sex.txt
  echo "105 SCH-254 1" >> update_sex.txt
  echo "update_sex.txt file created"
fi

# load initial shodair genotype .ped file into plink 
  # remove participant ID SCH-615 (no_risperidone.txt) because they were not taking risperidone
  # update the gender of SCH-937 and SCH-254 to male
  # retain autosomal and X chromosome variants only
plink --file ../input_plink_genotype_data/GSA2016_125_025 --keep-allele-order --remove no_risperidone.txt --update-sex update_sex.txt --chr 1-22 X --make-bed --out variant_QC_1
# 732163 variants and 167 people pass filters and QC.
# This will generate a warning message: "het. haploid genotypes present (see variant_QC_1.hh); many commands treat these as missing."

# check to see if the warning message is due to variants being present in pseudo-autosomal (XY) regions
# the warning message should disappear if so
plink --bfile variant_QC_1 --keep-allele-order --split-x b37 no-fail --make-bed --out variant_QC_2

# Since the warning message persists, check for any remaining sex discrepancies that could be a cause for the warning message
plink --bfile variant_QC_2 --keep-allele-order --check-sex
# Identify individuals with sex discrepancy.
grep "PROBLEM" plink.sexcheck | awk '{print$1,$2}'> sex_discrepancy.txt
# there are NO individuals with a sex discrepancy here

# Since there are no sex discrepancies, set the problematic heterozygous haploid genotypes to missing with the --set-hh-missing flag 
# to get rid of the warning message in the commands following this one
plink --bfile variant_QC_2 --set-hh-missing --keep-allele-order --make-bed --out variant_QC_3

# Reduce the genotyping data down to just single nucleotide polymorphisms (snps) only (INDELs are not compatible with imputation software!),
# and remove snps with non-traditional allele codes (just-acgt) because these variants are not compatible with imputation software either.
# For more info on input data formats for the TOPMED imputation server, see: https://topmedimpute.readthedocs.io/en/latest/pipeline.html
plink --bfile variant_QC_3 --keep-allele-order --snps-only just-acgt --make-bed --out variant_QC_4
# 722669 variants and 167 people pass filters and QC.

# convert the plink file to a VCF file for updating variant IDs with BCFtools
plink --bfile variant_QC_4 --keep-allele-order --recode vcf --out variant_QC_5

# convert all the existing variant IDs to unique variant IDs with the format chromosome number, position, reference allele, alternate allele 
# (variant IDs must be unique for PLINK to perform calculations properly)
bcftools annotate variant_QC_5.vcf -x ID -I +'%CHROM:%POS:%REF:%ALT' -o variant_QC_6.vcf

# output only the first occurance of duplicate variants (same chrom,pos,ref,alt) if a variant is present multiple times
bcftools norm --rm-dup exact variant_QC_6.vcf -o variant_QC_7.vcf
# total/split/realigned/skipped:	722669/0/0/0

# since variants may be duplicated more than once, run the same command again
bcftools norm --rm-dup exact variant_QC_7.vcf -o variant_QC_8.vcf
# total/split/realigned/skipped:	716234/0/0/0

# run the same command again, until the total number of variants remains constant (no duplicates remain)
bcftools norm --rm-dup exact variant_QC_8.vcf -o variant_QC_9.vcf
# total/split/realigned/skipped:	716234/0/0/0

# create a new file of participant IDs and genders from the variant_QC_4.fam file
rm update_sex.txt
awk '{print$1" "$2" "$5}' variant_QC_4.fam > update_sex.txt

# convert the vcf back to a plink file and update the participant gender
plink --vcf variant_QC_9.vcf --keep-allele-order --update-sex update_sex.txt --make-bed --out variant_QC_10
# 716234 variants and 167 people pass filters and QC.

# check to see that all variant IDs are unique (this step should print nothing to the console)
cut -f2 variant_QC_10.bim | sort | uniq -c | awk '$1>1'

# Investigate missingness per individual and per SNP
plink --bfile variant_QC_10 --keep-allele-order --missing    
# output: plink.imiss and plink.lmiss

# Delete SNPs and individuals with high levels of missingness 
# Delete SNPs with missingness >0.2.
plink --bfile variant_QC_10 --keep-allele-order --geno 0.2 --make-bed --out variant_QC_11
# 294 variants removed due to missing genotype data (--geno).
# 715940 variants and 167 people pass filters and QC.

# Delete individuals with missingness >0.2.
plink --bfile variant_QC_11 --keep-allele-order --mind 0.2 --make-bed --out variant_QC_12
# 0 people removed due to missing genotype data

# Delete SNPs with missingness >0.02.
plink --bfile variant_QC_12 --keep-allele-order --geno 0.02 --make-bed --out variant_QC_13
# 10340 variants removed due to missing genotype data (--geno).
# 705600 variants and 167 people pass filters and QC.

# Delete individuals with missingness >0.02.
plink --bfile variant_QC_13 --keep-allele-order --mind 0.02 --make-bed --out variant_QC_14
# 0 people removed due to missing genotype data (--mind)

#### According to Psychiatric Genomics Consortium (PGC) paper supplementary methods, 
    # reference: Ripke, et al., Nature, 511, 421-427, 2014 (https://www.nature.com/articles/nature13595)
    # remove X chromosome markers from all of the analysis based on snp missingness and deviation from HWE in female individuals only 
    # excerpt from paper:
        # "Chromosome X imputation was conducted for subjects passing quality control 
        # for the autosomal analysis with the additional exclusions of chrX SNPs with 
        # missingness >= 0.05 or HWE p < 10^-6 in females"
        
# isolate a plink file for female participants only for chromosome 23 (X chromosome)
plink --bfile variant_QC_14 --keep-allele-order --chr 23 --filter-females --make-bed --out variant_QC_14_X_chrom_females
# 132 people removed due to gender filter (--filter-females).
# 27035 variants and 35 people pass filters and QC.

#### delete snps with missingness >= 0.05 and with HWE p < 10^-6 in the female cohort only

# delete snps with missingness >=0.05
plink --bfile variant_QC_14_X_chrom_females --keep-allele-order --geno 0.05 --make-bed --out variant_QC_14_X_chrom_females_geno0.05
# 36 variants removed due to missing genotype data (--geno).
# 26999 variants and 35 people pass filters and QC.

# and with HWE p < 10^-6 in the female cohort only
plink --bfile variant_QC_14_X_chrom_females_geno0.05 --keep-allele-order --hwe include-nonctrl 1e-6 --make-bed --out variant_QC_14_X_chrom_females_geno0.05_hwe1e-6
# --hwe: 0 variants removed due to Hardy-Weinberg exact test.
# 26999 variants and 35 people pass filters and QC.

# identify the variants that were removed from the females only by identifying the variant IDs that 
# are in the variant_QC_14_X_chrom_females but not in the variant_QC_14_X_chrom_females_geno0.05_hwe1e-6 file
# the tidyverse package must be installed for this script to run
# this script will return a file called VariantIDsToRemove.txt to be used in the next command
Rscript --no-save ../scripts/gwas/identify_female_X_chromosome_snps_to_remove.R

# remove these snps from all of the data, not just the data for females
plink --bfile variant_QC_14 --keep-allele-order --exclude VariantIDsToRemove.txt --make-bed --out variant_QC_15
# 705564 variants and 167 people pass filters and QC.

# Remove SNPs with a MAF < 0.05
plink --bfile variant_QC_15 --keep-allele-order --maf 0.05 --make-bed --out variant_QC_16
# 398956 variants removed due to minor allele threshold(s)
# 306608 variants and 167 people pass filters and QC.

# Delete SNPs which are not in Hardy-Weinberg equilibrium [ HWE p-value < 10^-6 ]
plink --bfile variant_QC_16 --keep-allele-order --hwe include-nonctrl 1e-6 --make-bed --out variant_QC_17
# --hwe: 21 variants removed due to Hardy-Weinberg exact test.
# 306587 variants and 167 people pass filters and QC.

# Calculate heterozygosity rate of each individual using linkage disequilibrium pruned snps. 
# Individuals with heterozygosity rate +/- 3 standard deviations from the population mean heterozygosity rate need to be removed from the analysis

# obtain a set of linkage disequilibrium pruned snps with window size of 200kb, step size of 100kb, and linkage disequilibrium r-squared cutoff of 0.2
# while excluding regions of the genome with extremely high linkage disequilibrium (inversion.txt file)
plink --bfile variant_QC_17 --keep-allele-order --exclude range ../region_to_exclude_from_LD_pruning/inversion.txt --indep-pairwise 200 100 0.2 --out indepSNP

# extract the linkage disequilibrium pruned variants from the microarray genotyping data
plink --bfile variant_QC_17 --keep-allele-order --extract indepSNP.prune.in --het --out R_check
# results are written to R_check.het file

# generate a list of individuals deviating by more than 3sd from heterozygosity rate mean 
# with the R-script for the Maurees, et al. manuscript
Rscript --no-save ../scripts/gwas/heterozygosity_outliers_list.R
# this will create a file called fail-het-qc.txt

# Adapt this file to make it compatible for PLINK, by removing all quotation marks from the file and selecting only the first two columns.
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt

# Remove heterozygosity rate outliers.
plink --bfile variant_QC_17 --keep-allele-order --remove het_fail_ind.txt --make-bed --out variant_QC_18
# 167 people (132 males, 35 females) loaded from .fam.
# --remove: 164 people remaining.

# checking for cryptic relatedness, removing people above a pi-hat of 0.2 (i.e., related by the second degree or closer)
# check for relatedness of people of 2nd degree or closer with linkage disequilibrium pruned variants
plink --bfile variant_QC_18 --keep-allele-order --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2
# results written to pihat_min0.2.genome

# get individual missingness reports for all individuals in the study,
# one participant in each pair of participants related by 2nd degree or more in the pihat_min0.2.genome need to be removed
# based on greater proportion of missing genotypes
plink --bfile variant_QC_18 --keep-allele-order --missing
# report written to plink.imiss

# identify individuals that need to be removed from each related pair
Rscript --no-save ../scripts/gwas/identify_related_individuals_to_remove.R

# remove the individuals that have the lower call rates from the related pairs
plink --bfile variant_QC_18 --keep-allele-order --remove 0.2_low_call_rate_pihat.txt --make-bed --out variant_QC_19
# 164 people (130 males, 34 females) loaded from .fam.
# --remove: 161 people remaining.

# save the final quality controlled plink file to the file system as a vcf file
mkdir ../quality_control_output/
plink --bfile variant_QC_19 --keep-allele-order --recode vcf --out variant_QC_output_unphased
# 161 people (129 males, 32 females) loaded from .fam.
# 306587 variants and 161 people pass filters and QC.
mv variant_QC_output_unphased.vcf ../quality_control_output/

#If needed, download 1000 genomes variant calling data with the command below
if [[ -f ../Genome_1000_Data/ALL.2of4intersection.20100804.genotypes.vcf.gz ]]
then
  echo "1000 genomes variant calling data already exists"
else
  echo "1000 genomes variant calling data does not exist"
  mkdir ../Genome_1000_Data/
  wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz
  mv ALL.2of4intersection.20100804.genotypes.vcf.gz ../Genome_1000_Data/
  echo "1000 genomes variant calling data was downloaded"
fi

# create plink files with the 1000 genomes vcf file
# for autosomal and X chromosome variants only
plink --vcf ../Genome_1000_Data/ALL.2of4intersection.20100804.genotypes.vcf.gz --keep-allele-order --chr 1-22 X --make-bed --out 1K_genomes_1

# Reduce the genotyping data down to just single nucleotide polymorphisms (snps) only,
# and remove snps with non-traditional allele codes (just-acgt) 
plink --bfile 1K_genomes_1 --keep-allele-order --snps-only just-acgt --make-bed --out 1K_genomes_2
rm 1K_genomes_1*

# convert the plink file to a VCF file for updating variant IDs with BCFtools
plink --bfile 1K_genomes_2 --keep-allele-order --recode vcf --out 1K_genomes_3
rm 1K_genomes_2*

# convert all the existing variant IDs to unique variant IDs with the format 
# chromosome number, position, reference allele, alternate allele 
# (variant IDs must be unique for PLINK to perform calculations properly)
bcftools annotate 1K_genomes_3.vcf -x ID -I +'%CHROM:%POS:%REF:%ALT' -o 1K_genomes_4.vcf
rm 1K_genomes_3*

# output only the first occurance of duplicate variants (same chrom,pos,ref,alt) if a variant is present multiple times
bcftools norm --rm-dup exact 1K_genomes_4.vcf -o 1K_genomes_5.vcf
rm 1K_genomes_4*

# convert the VCF file back to a plink file
plink --vcf 1K_genomes_5.vcf --keep-allele-order --make-bed --out 1K_genomes_6
rm 1K_genomes_5*

# check to see that all variant IDs are unique (this step should print nothing to the console)
cut -f2 1K_genomes_6.bim | sort | uniq -c | awk '$1>1'

## QC on 1000 Genomes data.
# Remove variants based on missing genotype data.
plink --bfile 1K_genomes_6 --keep-allele-order --geno 0.2 --make-bed --out 1K_genomes_7

# Remove individuals based on missing genotype data.
plink --bfile 1K_genomes_7 --keep-allele-order --mind 0.2 --make-bed --out 1K_genomes_8

# Remove variants based on missing genotype data.
plink --bfile 1K_genomes_8 --keep-allele-order --geno 0.02 --make-bed --out 1K_genomes_9

# Remove individuals based on missing genotype data.
plink --bfile 1K_genomes_9 --keep-allele-order --mind 0.02 --make-bed --out 1K_genomes_10

# Remove variants based on MAF less than 0.05
plink --bfile 1K_genomes_10 --keep-allele-order --maf 0.05 --make-bed --out 1K_genomes_11

# Extract the variants present in the shodair dataset from the 1000 genomes dataset.
awk '{print$2}' variant_QC_19.bim > shodair_SNPs.txt
plink --bfile 1K_genomes_11 --keep-allele-order --extract shodair_SNPs.txt --make-bed --out 1K_genomes_12

# Extract the variants present in 1000 Genomes dataset from the shodair dataset.
awk '{print$2}' 1K_genomes_12.bim > 1K_genomes_12_SNPs.txt
plink --bfile variant_QC_19 --keep-allele-order --extract 1K_genomes_12_SNPs.txt --recode --make-bed --out shodair_MDS
# The datasets now contain the exact same variants.

## The datasets must have the same build. Change the 1000 Genomes data build.
awk '{print$2,$4}' shodair_MDS.map > build_shodair.txt
# build_shodair.txt contains one SNP-id and physical position per line.

plink --bfile 1K_genomes_12 --keep-allele-order --update-map build_shodair.txt --make-bed --out 1K_genomes_13
# the 1000 genomes data and shodair data now have the same build

## Merge the shodair and 1000 Genomes data sets
  # Prior to merging 1000 Genomes data with the shodair data we want to make sure that the files are mergeable, for this we conduct 3 steps:
    # 1) Make sure the reference genome is similar in the shodair and the 1000 Genomes Project datasets.
    # 2) Resolve strand issues.
    # 3) Remove the SNPs which after the previous two steps still differ between datasets.

# 1) set reference genome 
awk '{print$2,$5}' 1K_genomes_13.bim > 1kg_ref-list.txt
plink --bfile shodair_MDS --keep-allele-order --reference-allele 1kg_ref-list.txt --make-bed --out shodair-adj

# 2) Resolve strand issues.
# Check for potential strand issues.
awk '{print$2,$5,$6}' 1K_genomes_13.bim > 1K_genomes_13_tmp
awk '{print$2,$5,$6}' shodair-adj.bim > shodair-adj_tmp
sort 1K_genomes_13_tmp shodair-adj_tmp |uniq -u > all_differences.txt
# there are no strand issues if the all_differences.txt is empty

## Flip SNPs for resolving strand issues. (there are none)
# Print SNP-identifier and remove duplicates.
awk '{print$1}' all_differences.txt | sort -u > flip_list.txt
# Flip the non-corresponding snps (there are none) 
plink --bfile shodair-adj --keep-allele-order --flip flip_list.txt --reference-allele 1kg_ref-list.txt --make-bed --out corrected_shodair
#--flip: 0 SNPs flipped.

# Check for SNPs which are still problematic after they have been flipped. (there are none)
awk '{print$2,$5,$6}' corrected_shodair.bim > corrected_shodair_tmp
sort 1K_genomes_13_tmp corrected_shodair_tmp |uniq -u  > uncorresponding_SNPs.txt
# there are no uncorresponding snps

# 3) Remove problematic SNPs from shodair data and 1000 Genomes.
awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exlusion.txt
# there are zero snps here

# Remove the problematic snps (there are none in this case)
plink --bfile corrected_shodair --keep-allele-order --exclude SNPs_for_exlusion.txt --make-bed --out shodair_MDS2
plink --bfile 1K_genomes_13 --keep-allele-order --exclude SNPs_for_exlusion.txt --make-bed --out 1K_genomes_14

# Merge shodair with 1000 Genomes Data.
plink --bfile shodair_MDS2 --keep-allele-order --bmerge 1K_genomes_14.bed 1K_genomes_14.bim 1K_genomes_14.fam --make-bed --out MDS_merge2

## Perform multidimensional scaling (MDS) on the shodair data anchored by 1000 Genomes data.
# Using a set of linkage-disequilibrium pruned SNPs
plink --bfile MDS_merge2 --keep-allele-order --extract indepSNP.prune.in --genome --out MDS_merge2

# perform the MDS with the data anchored by 1000 genomes data
plink --bfile MDS_merge2 --keep-allele-order --read-genome MDS_merge2.genome --cluster --mds-plot 10 --out MDS_merge2

# Download the file with population information for the 1000 genomes dataset.
# The file 20100804.ALL.panel contains continental population codes of the individuals of 1000 genomes.
if [[ -f ../Genome_1000_Data/20100804.ALL.panel ]]
then
  echo "1000 genomes population code data already exists"
else
  echo "1000 genomes population code data does not exist"
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/20100804.ALL.panel
  mv 20100804.ALL.panel ../Genome_1000_Data/
  echo "1000 genomes population code data was downloaded"
fi

# Convert 1000 genomes subpopulation codes into superpopulation codes (i.e., AFR, AMR, ASN, and EUR).
awk '{print$1,$1,$2}' ../Genome_1000_Data/20100804.ALL.panel > populationCode_1kG.txt
sed 's/JPT/ASN/g' populationCode_1kG.txt > populationCode_1kG2.txt
sed 's/ASW/AFR/g' populationCode_1kG2.txt > populationCode_1kG3.txt
sed 's/CEU/EUR/g' populationCode_1kG3.txt > populationCode_1kG4.txt
sed 's/CHB/ASN/g' populationCode_1kG4.txt > populationCode_1kG5.txt
sed 's/CHD/ASN/g' populationCode_1kG5.txt > populationCode_1kG6.txt
sed 's/YRI/AFR/g' populationCode_1kG6.txt > populationCode_1kG7.txt
sed 's/LWK/AFR/g' populationCode_1kG7.txt > populationCode_1kG8.txt
sed 's/TSI/EUR/g' populationCode_1kG8.txt > populationCode_1kG9.txt
sed 's/MXL/AMR/g' populationCode_1kG9.txt > populationCode_1kG10.txt
sed 's/GBR/EUR/g' populationCode_1kG10.txt > populationCode_1kG11.txt
sed 's/FIN/EUR/g' populationCode_1kG11.txt > populationCode_1kG12.txt
sed 's/CHS/ASN/g' populationCode_1kG12.txt > populationCode_1kG13.txt
sed 's/PUR/AMR/g' populationCode_1kG13.txt > populationCode_1kG14.txt

# Create a populationCodefile of the shodair data
awk '{print$1,$2,"OWN"}' shodair_MDS.fam > populationCodefile_own.txt

# Concatenate populationCodefiles and add a header to the file
echo "FID IID populationCode" > fileHeader.txt
cat fileHeader.txt populationCode_1kG14.txt populationCodefile_own.txt > populationCodefile.txt

# Generate population stratification plot.
Rscript --no-save ../scripts/gwas/MDS_merged.R
# this will create MDS.pdf file with a multi-dimensional scaling plot

## Create covariates based on multidimensional scaling components using the shodair data only (no longer anchored by the 1000 genomes data).
  # Perform an MDS on the shodair data with all participants included on a linkage disequilibrium pruned set of snps.
  # The values of the 10 MDS components are used as covariates in the regression association analysis to adjust 
  # for any potential confounding due to population stratification
plink --bfile shodair_MDS --keep-allele-order --extract indepSNP.prune.in --genome --out shodair_MDS
plink --bfile shodair_MDS --keep-allele-order --read-genome shodair_MDS.genome --cluster --mds-plot 10 --out shodair_MDS_results

# convert the shodair_MDS_results.mds file to a plink covariate file (covar_mds.txt) for performing regression analysis
Rscript --no-save ../scripts/gwas/prepare_MDS_covariate_file.R

### delete intermediate files to save memory
rm *bed
rm *bim
rm *fam
rm *ped
rm *vcf



###############################
