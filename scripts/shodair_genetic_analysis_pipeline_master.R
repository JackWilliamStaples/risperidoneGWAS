######## this is a wrapper script for running the entire shodair psychiatric GWAS analysis pipeline

# set the working directory
setwd(dir = "~/Psychiatric_GWAS/Psychiatric_GWAS/")

# print current time at start to the screen
system(command = 'date +"%r"')

### grant permission to execute all scripts in the pipeline
system(command = "chmod +x ./scripts/gwas/shodair_quality_control.sh")
system(command = "chmod +x ./scripts/gwas/prepare_VCF_for_imputation.sh")
system(command = "chmod +x ./scripts/gwas/convert_imputed_files_to_VCF.sh")
system(command = "chmod +x ./scripts/variant_annotation/ANNOVAR_annotation.sh")
system(command = "chmod +x ./scripts/variant_annotation/call_star_alleles.sh")
system(command = "chmod +x ./scripts/gwas/association_analysis.sh")
### execute each script

# perform quality control of microarray genotype data with PLINK
system(command = "./scripts/gwas/shodair_quality_control.sh")

# prepare vcf files for phasing and imputation with the TOPMED imputation server
system(command = "./scripts/gwas/prepare_VCF_for_imputation.sh")

# perform quality control and filtering of imputed variant call data
system(command = "./scripts/gwas/convert_imputed_files_to_VCF.sh")

# annotate variant calling data with ANNOVAR, the variant effect predictor, pharmVar, and pharmGKB data
system(command = "./scripts/variant_annotation/ANNOVAR_annotation.sh")

# call star alleles on the variant calling data with starGazer and PGxPOP for all star alleles
system(command = "./scripts/variant_annotation/call_star_alleles.sh")

# perform genome-wide association study analysis
system(command = "./scripts/gwas/association_analysis.sh")

# print current time at start to the screen
system(command = 'date +"%r"')


