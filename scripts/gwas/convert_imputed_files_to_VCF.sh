# this is a script for merging dosage VCF files output from the TOPMED imputation server,
# performing post-imputation quality control, reformatting VCF files, and converting to PLINK files
  # (some steps for this can be found at this link: https://github.com/huw-morris-lab/imputation )

# set the path to the bioinformatics software directory
PATH=$PATH:~/bioinformatics_software/

# set the working directory
cd ~/Psychiatric_GWAS/Psychiatric_GWAS/TOPMEDimputationResults/

# concatenate the imputed VCF files together with BCFtools
bcftools concat chr_1/chr1.dose.vcf.gz chr_2/chr2.dose.vcf.gz chr_3/chr3.dose.vcf.gz chr_4/chr4.dose.vcf.gz chr_5/chr5.dose.vcf.gz chr_6/chr6.dose.vcf.gz chr_7/chr7.dose.vcf.gz chr_8/chr8.dose.vcf.gz chr_9/chr9.dose.vcf.gz chr_10/chr10.dose.vcf.gz chr_11/chr11.dose.vcf.gz chr_12/chr12.dose.vcf.gz chr_13/chr13.dose.vcf.gz chr_14/chr14.dose.vcf.gz chr_15/chr15.dose.vcf.gz chr_16/chr16.dose.vcf.gz chr_17/chr17.dose.vcf.gz chr_18/chr18.dose.vcf.gz chr_19/chr19.dose.vcf.gz chr_20/chr20.dose.vcf.gz chr_21/chr21.dose.vcf.gz chr_22/chr22.dose.vcf.gz chr_x/chrX.dose.vcf.gz -Oz -o temp1.vcf.gz
# filter to variants with an INFO score ≥ 0.3, this is a common threshold for GWAS
bcftools view -i 'R2>=0.30' temp1.vcf.gz -Oz -o temp2.vcf.gz
# split the rows containing multi-allelic sites if there are any present in the file
bcftools norm -m-any temp2.vcf.gz -Oz -o temp3.vcf.gz

# create a file for updating the chromosome labels in the vcf file (the "chr" prefix needs to be removed)
for chromosome in {1..22} X
do
  echo "chr${chromosome} ${chromosome}" >> new_chromosome_codes.txt
done

# update the chromosome codes in the vcf file
bcftools annotate --rename-chrs new_chromosome_codes.txt temp3.vcf.gz -Oz -o temp4.vcf.gz
# convert all the existing variant IDs to unique variant IDs with the format chromosome number, position, reference allele, alternate allele. the variant IDs MUST be unique before loading into PLINK
bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' temp4.vcf.gz -Oz -o temp5.vcf.gz
# filter to variants with a minor allele frequency ≥ 0.05 and with a non-missing genotype in at least 20 individuals (successfully imputed in ≥ 20 individuals)
vcftools --gzvcf temp5.vcf.gz --maf 0.05 --max-missing-count 20 --recode --recode-INFO-all --out imputed_variants_allCHROMs_INFO_0.3_MAF_0.05_maxMissing20.dose
# delete temporary files
rm temp*


