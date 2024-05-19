# this is a script for identifying variants with strand-flip issues from the TOPMED imputation server
# quality control report (snps-excluded.txt)

# load necessary packages
source(file = "../scripts/load_necessary_packages.R")

# load the snps-excluded.txt file from the initial imputation attempt
SNPsExcluded <-
  read.delim(
    file = "../QCreport_TOPMEDimputationServer/snps-excluded.txt",
    header = TRUE
    ) 

SNPsExcluded <-
  # identify snps with strand flip issues and save a list of them to the file system
  SNPsExcluded %>%
    tibble() %>%
    filter(
      .data = .,
      (FilterType == "Strand flip") | (FilterType == "Strand flip and Allele switch")
      ) %>%
    mutate(
      .data = .,
      # create a column with the chromosome number only only
      chrom = str_extract(
                          # isolate the "chr" and chromosome number string
                          string = X.Position,
                          pattern = "chr[:alnum:]*"
                          ) %>%
              # remove the "chr" from the string
              gsub(
                pattern = "chr",
                replacement = "",
                x = .
                  ) %>%
              # replace the "X" in X chromosome with 23
              gsub(
                pattern = "X",
                replacement = 23,
                x = .
                ),
      # create a genomic position column
      position = str_split(
                           string = X.Position,
                           pattern = ":"
                           ) %>%
                  lapply(
                    X = .,
                    FUN = function(currentArray)
                    {
                      ValueToReturn <-
                        currentArray[2]
                      
                      return(ValueToReturn)
                    }
                      ) %>%
                    unlist(x = .),
      # remove the "chr" from the X.Position column of variant IDs
      X.Position = gsub(
                        pattern = "chr",
                        replacement = "",
                        x = X.Position
                        ) %>%
                   # replace X chromosome with 23
                   gsub(
                     pattern = "X",
                     replacement = 23,
                     x = .
                     )
      ) %>%
    mutate(
      .data = .,
      # create a column of the chromosome and position concatenated
      chrom.position = paste0(chrom,":",position)
      ) 

# load the .bim file of variants that will be used for imputation
VariantsForImputation <-
  read.table(
    file = "./variant_QC_output_unphased_sorted.bim",
    header = FALSE
    )

# update the column names of the VariantsForImputation
names(x = VariantsForImputation) <- c("chrom","snp","position_cM","position","A1","A2")

VariantsForImputation <-
  VariantsForImputation %>%
    # concatenate the chromosome and genomic position columns together
    mutate(
      .data = .,
      chrom.position = paste0(chrom,":",position)
    )

# filter the sequencing data to variants that need to be flipped from the SNPsExcluded file
VariantsForImputation %>%
  filter(
    .data = .,
    chrom.position %in% SNPsExcluded$chrom.position
    ) %>%
  select(
    .data = .,
    snp
    ) %>%
  write.table(
    x = .,
    file = "./snps_to_flip.txt",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
    )
