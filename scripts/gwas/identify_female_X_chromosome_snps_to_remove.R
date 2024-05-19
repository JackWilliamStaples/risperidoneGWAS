# this is a script for identifying variant IDs that need to be removed from the analysis 
# based on quality control of the X chromosome in only the female participants

# load necessary packages
source(file = "../scripts/load_necessary_packages.R")

# read in the variant_QC_14_X_chrom_females.bim and variant_QC_14_X_chrom_females_geno0.05_hwe1e-6.bim files from the file system 
variant_QC_14_X_chrom_females <-
  read.delim(
    file = "variant_QC_14_X_chrom_females.bim",
    stringsAsFactors = FALSE,
    header = FALSE
    )
  
variant_QC_14_X_chrom_females_geno0.05_hwe1e_6 <-
  read.delim(
    file = "variant_QC_14_X_chrom_females_geno0.05_hwe1e-6.bim",
    stringsAsFactors = FALSE,
    header = FALSE
    )

# select the variant ID columns from both the variant_QC_14_X_chrom_females.bim and variant_QC_14_X_chrom_females_geno0.05_hwe1e-6.bim files
variant_QC_14_X_chrom_females_variantIDs <-
  variant_QC_14_X_chrom_females %>%
  # select the variant ID column
  select(
    .data = .,
   "variantID" = V2
    ) %>%
  pull(
    .data = .,
    variantID
    )

variant_QC_14_X_chrom_females_geno0.05_hwe1e_6_variantIDs <-
  variant_QC_14_X_chrom_females_geno0.05_hwe1e_6 %>%
  # select the variant ID column
  select(
    .data = .,
    "variantID" = V2
  ) %>%
  pull(
    .data = .,
    variantID
    )

# find the difference in the two variant ID sets
VariantIDsToRemove <-
  setdiff(
    x = variant_QC_14_X_chrom_females_variantIDs,
    y = variant_QC_14_X_chrom_females_geno0.05_hwe1e_6_variantIDs
    ) %>%
  # convert to a tibble 
  tibble()
  
# save the single column tibble of variant IDs to remove to the file system
write.table(
  x = VariantIDsToRemove,
  file = "VariantIDsToRemove.txt",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
  )