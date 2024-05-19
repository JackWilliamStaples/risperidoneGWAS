# this is a script for performing genome-wide association analysis of microarray genotyping data and clinical risperidone treatment outcomes

# make sure the path points to the bioinformatics software directory
PATH=$PATH:~/bioinformatics_software/

# create a working directory for the gwas analysis if it does not exist
if [[ -d ~/Psychiatric_GWAS/Psychiatric_GWAS/gwas_analysis ]]
then
  echo "GWAS analysis working directory exists"
else
  echo "Created working directory for GWAS analysis"
  mkdir ~/Psychiatric_GWAS/Psychiatric_GWAS/gwas_analysis/
fi

# working directory for this is
cd ~/Psychiatric_GWAS/Psychiatric_GWAS/gwas_analysis/

# run the script for preparing continuous and binary PLINK phenotype files with clinical outcomes data
Rscript --no-save ../scripts/gwas/prepare_plink_phenotype_files.R

# update the participant IDs to IIDs in the VCF file so they are compatible with the next step
plink --vcf ../TOPMEDimputationResults/imputed_variants_allCHROMs_INFO_0.3_MAF_0.05_maxMissing20.dose.recode.vcf.gz --keep-allele-order --recode vcf-iid --out temp
# test individual alleles for association with continuous phenotypes with an additive allele dose linear regression model
plink --vcf temp.vcf --keep-allele-order --double-id --pheno PLINKphenotypeFile_continuous.txt --allow-no-sex --all-pheno --linear --ci 0.95 --out continuous_phenotypes
# test individual alleles for association with binary phenotypes with an allele dose logistic regression model
plink --vcf temp.vcf --keep-allele-order --double-id --pheno PLINKphenotypeFile_binary.txt --1 --allow-no-sex --all-pheno --logistic --ci 0.95 --seed 777 --out binary_phenotypes

# Perform SNP clumping with PLINK for continuous phenotypes (for more information on clumping, see: https://www.cog-genomics.org/plink/1.9/postproc#clump)
# to reduce the SNPs down to the most associated SNPs with phenotype (i.e., removing any SNPs that are strongly correlated with the causative SNP based on linkage disequilibrium).
# Do this with the association results for each phenotype
for currentPhenotype in {ChangeInGAF,DaysOnRisperidone,DaysToReadmittance,GAF_admit,GAF_discharge,LengthOfStay,MaxDoseRisperidone_mg}
do
  plink --vcf temp.vcf --keep-allele-order --double-id --allow-no-sex --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump "continuous_phenotypes.${currentPhenotype}.assoc.linear" --clump-snp-field SNP --clump-field P --out "clump_file_additive_assoc.${currentPhenotype}"
done

# Perform SNP clumping for binary phenotypes as well
for currentPhenotype in {ReadmittedOnRisperidone,Risperidone_discontinued,PatientReadmitted}
do
  plink --vcf temp.vcf --keep-allele-order --double-id --allow-no-sex --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump "binary_phenotypes.${currentPhenotype}.assoc.logistic" --clump-snp-field SNP --clump-field P --out "clump_file_logistic_assoc.${currentPhenotype}"
done

# filter the clump files to variants that were significant after Bonferroni adjustment for the number of
# independent, causative variants
Rscript --no-save ../scripts/gwas/identify_variants_with_GWAS_significance.R

# create manhattan and qq-plots to visualize the association plots
Rscript --no-save ../scripts/gwas/manhattan_and_qqplots.R

# create a PLINK file of the variants that were significant after Bonferroni adjustment for the number of independent, causative variants
# create genotype matrices with these variants with genotypes as 0, 1, or 2 copies of the variant allele
plink --vcf temp.vcf --keep-allele-order --double-id --allow-no-sex --extract SignificantVariantsFromGWAS_forFilteringPLINKfiles.txt --recode A --out SignificantVariantsFromGWAS_geno_matrix

# compute summary statistics of clinical and demographic data
# and perform regression and survival analysis of study outcomes versus genetic data and clinical/demographic data
Rscript --no-save ../scripts/gwas/PhenotypeRegressionAnalysis.R

# delete temporary files
rm *nosex
rm *log
rm *linear
rm *clumped
rm *logistic
rm temp*
