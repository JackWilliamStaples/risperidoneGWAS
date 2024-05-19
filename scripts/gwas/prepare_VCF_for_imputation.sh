# this is a script for preparing VCF files that are output from quality control for phasing and imputation 
# with the TOPMED imputation server (https://topmedimpute.readthedocs.io/en/latest/prepare-your-data.html)

# set the path to the bioinformatics software directory
PATH=$PATH:~/bioinformatics_software/

# create a working directory for this if it does not exist
if [[ -d ~/Psychiatric_GWAS/Psychiatric_GWAS/VCF_files_for_imputation ]]
then
  echo "directory for vcf files for imputation exists"
else
  echo "created directory for vcf files for imputation"
  mkdir ~/Psychiatric_GWAS/Psychiatric_GWAS/VCF_files_for_imputation/
fi

# set the working directory
cd ~/Psychiatric_GWAS/Psychiatric_GWAS/VCF_files_for_imputation/

# sort the variants based on genomic position
bcftools sort ../quality_control_output/variant_QC_output_unphased.vcf -o variant_QC_output_unphased_sorted.vcf

# for each chromosome in the VCF file
for currentChromosome in {1..23}
do
  # split the vcf file by chromosome
  plink --vcf variant_QC_output_unphased_sorted.vcf --keep-allele-order --chr ${currentChromosome} --recode vcf --out temp
  # compress and index the resulting file
  bgzip -c temp.vcf > "variant_QC_output_unphased_sorted_chr${currentChromosome}.vcf.gz"
  tabix -p vcf "variant_QC_output_unphased_sorted_chr${currentChromosome}.vcf.gz"
  rm temp*
done

# submit each zipped VCF file for each chromosome to the TOPMED imputation server (https://imputation.biodatacatalyst.nhlbi.nih.gov/#!pages/home):
    # For the imputation, the parameters and flags selected were:
      # reference panel: TOPMED r2
      # array build: GRCh37/hg19
      # rsq filter: 0.3
      # phasing: Eagle v2.4 (phased output)
      # population: Skip
      # Mode: Quality Control & Imputation
      
# the TOPMED imputation server will initially perform additional quality control steps before imputation
# if there are too many variants with strand flip issues, the imputation will fail and a snps-excluded.txt 
# can be downloaded to see which variants had strand flip issues

# obtain a .bim file of snps in the variant_QC_output_unphased_sorted.vcf file
plink --vcf variant_QC_output_unphased_sorted.vcf --keep-allele-order --make-bed --out variant_QC_output_unphased_sorted

# obtain the variants that had strand flip issues from the snps-excluded.txt file from the initial imputation attempt
# the snps-excluded.txt file is stored in the ../QCreport_TOPMEDimputationServer/ directory
Rscript --no-save ../scripts/gwas/identify_snps_with_strand_flip_issues.R
# this will create a file called snps_to_flip.txt

# use plink to flip the variants that had strand issues in the initial TOPMED imputation run
plink --vcf variant_QC_output_unphased_sorted.vcf --keep-allele-order --flip snps_to_flip.txt --recode vcf --out variant_QC_output_unphased_sorted_flipped

# make sure the vcf file is sorted based on genomic position again
bcftools sort variant_QC_output_unphased_sorted_flipped.vcf -o variant_QC_output_unphased_sorted_flipped_sorted.vcf

# delete the old vcf files that were split by chromosome
rm variant_QC_output_unphased_sorted_chr*

# split the vcf file again by chromosome and compress and index with bgzip and tabix, respectively
for currentChromosome in {1..23}
do
  # split the vcf file by chromosome
  plink --vcf variant_QC_output_unphased_sorted_flipped_sorted.vcf --keep-allele-order --chr ${currentChromosome} --recode vcf --out temp
  # compress and index the resulting file
  bgzip -c temp.vcf > "variant_QC_output_unphased_sorted_chr${currentChromosome}.vcf.gz"
  tabix -p vcf "variant_QC_output_unphased_sorted_chr${currentChromosome}.vcf.gz"
  rm temp*
done

# Submit the VCF files to the TOPMED imputation server with the parameters specified in the comments above.
# Submit the files as separate jobs with only 4 chromosomes at a time or else the job will fail for too many snps needing
# to be excluded from the analysis.

# You will receive emails from the TOPMED imputation server with the account you registered with when each job is finished.
# When the jobs are finished, store the results in the newly created ../TOPMEDimputationResults/ directory
if [[ -d ~/Psychiatric_GWAS/Psychiatric_GWAS/TOPMEDimputationResults ]]
then
  echo "directory for imputation results exists"
else
  echo "created directory imputation results"
  mkdir ~/Psychiatric_GWAS/Psychiatric_GWAS/TOPMEDimputationResults/
fi

### delete intermediate files to save memory
rm *bed
rm *bim
rm *fam
rm *vcf*

