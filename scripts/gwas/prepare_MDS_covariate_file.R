# this is a script for converting a plink multi-dimensional scaling (MDS) results file
# to a plink covariates file

# load necessary packages
source(file = "../scripts/load_necessary_packages.R")

# load the MDS file
MDSfile <-
  read.table(
    file = "./shodair_MDS_results.mds",
    header = TRUE
    )

# save a covariates file to the file system with participant IDs and 10 MDS components
MDSfile %>%
  mutate(
    .data = .,
    FID = IID
    ) %>%
  select(
    .data = .,
    -SOL
    ) %>%
  write_delim(
    x = .,
    file = "./covar_mds.txt"
    )
