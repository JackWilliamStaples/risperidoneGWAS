# load necessary packages
source(file = "../scripts/load_necessary_packages.R")

# this is a script for identifying participants that need to be 
# removed from 2nd degree relative pairs based on greater proportion of missing 
# genotype calls

# load the file with individual genotype missingness
MissingnessData <-
  read.table(
    file = "./plink.imiss",
    header = TRUE
    )

# load the file with related pairs of individuals
RelatedPairs <-
  read.table(
    file = "./pihat_min0.2.genome",
    header = TRUE
    )

RelatedPairs <-
  RelatedPairs %>%
    # add a column with proportion of missing genotypes for IID1 participant IDs
    mutate(
      .data = .,
      F_MISS1 = MissingnessData %>%
                select(
                  .data = .,
                  IID,
                  F_MISS
                  ) %>%
                left_join(
                  x = RelatedPairs,
                  y = .,
                  by = c("IID1" = "IID")
                  ) %>%
                pull(
                  .data = .,
                  F_MISS
                  )
      ) %>%
    # add a column with proportion of missing genotypes for IID2 participant IDs
    mutate(
      .data = .,
      F_MISS2 = MissingnessData %>%
                select(
                  .data = .,
                  IID,
                  F_MISS
                  ) %>%
                left_join(
                  x = RelatedPairs,
                  y = .,
                  by = c("IID2" = "IID")
                  ) %>%
                pull(
                  .data = .,
                  F_MISS
                  )
      ) %>%
  # identify the individual with the greater proportion of missing genotypes
  mutate(
    .data = .,
    ParticipantToRemove = if_else(
                                  condition = F_MISS1 > F_MISS2,
                                  true = paste(FID1,IID1),
                                  false = paste(FID2,IID2)
                                    )
    )


# save a file of the individuals to remove to the file system
RelatedPairs %>%
  select(
    .data = .,
    ParticipantToRemove
    ) %>%
  write.table(
    x = .,
    file = "./0.2_low_call_rate_pihat.txt",
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE
      )

