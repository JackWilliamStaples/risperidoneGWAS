# this is a script for filtering clump files to variants that were statistically
# significant after Bonferroni adjustment for the number of independent, causative
# variants and combining with variant annotation data

# load necessary R packages 
source(file = "../scripts/load_necessary_packages.R")

# identify clump files to load for all phenotypes by identifying files 
# that end with "clumped"
FilesToLoad <-
    list.files(path = "./") %>%
      grepl(
        pattern = "*clumped$",
        x = .
      ) %>%
      list.files(path = "./")[.]

# load each of the clump files and filter to variants with genome-wide
# significance
SignificantVariantsFromGWAS <-
  FilesToLoad %>%
  mclapply(
    X = .,
    FUN = function(currentFile)
    {
      # identify the current phenotype by parsing the clump file name
      currentPhenotype <-
        currentFile %>%
        gsub(
          pattern = "(clump_file_additive_assoc.|clump_file_logistic_assoc.)",
          replacement = "",
          x = .
        ) %>%
        gsub(
          pattern = ".clumped",
          replacement = "",
          x = .
        )

      GWASresultsFileToLoad <-
        # load the GWAS results file for the current phenotype to obtain
        # regression summary statistics
        list.files(path = "./") %>%
          grepl(
            pattern = "(*linear$|*logistic$)",
            x = .
            ) %>%
          list.files(path = "./")[.] %>%
          grepl(
            pattern = currentPhenotype,
            x = .
            ) %>%
          list.files(path = "./")[grepl(pattern = "(*linear$|*logistic$)",x = list.files(path = "./"))][.]

       # load the clump file for the current phenotype
       LoadedClumpFile <-
          read.table(
            file = glue("./{currentFile}"),
            header = TRUE,
            quote = ""
            )

       FilteredClumpFile <-
           LoadedClumpFile %>%
           mutate(
             .data = .,
             P = P %>% as.numeric(x = .)
             ) %>%
           # filter the loaded file to P-values < 3x10^-7, estimated from the number of independent snps (number of rows in clump file)
           filter(
             .data = .,
             P < as.numeric(x = 0.05/nrow(x = LoadedClumpFile))
             ) %>%
           # select columns with the variant ID and P-value
           select(
             .data = .,
             SNP,
             P
             ) %>%
           # create a column to label the current phenotype
           mutate(
             .data = .,
             Phenotype = currentPhenotype
             ) %>%
           # create a column with the P-value significance threshold
           mutate(
             .data = .,
             alpha = as.numeric(x = 0.05/nrow(x = LoadedClumpFile))
             )

       # load the GWAS results file for the current phenotype
       LoadedGWASresultsFile <-
         read.table(
           file = paste0("./",GWASresultsFileToLoad),
           header = TRUE,
           sep = ""
           ) %>%
         mutate(
           .data = .,
           P = P %>% as.numeric(x = .)
           ) %>%
         # filter the loaded file to P-values < 3x10^-7, estimated from the number of independent snps (number of rows in clump file)
         filter(
           .data = .,
           P < as.numeric(x = 0.05/nrow(x = LoadedClumpFile))
         )

       # if there are any genome-wide significant variants present,
       # join the clump file data with the GWAS summary statistics data
       if(
         nrow(x = LoadedGWASresultsFile)>0 &
         nrow(x = FilteredClumpFile)>0)
       {
         FileToReturn <-
           FilteredClumpFile %>%
               left_join(
                 x = .,
                 y = LoadedGWASresultsFile %>%
                     select(
                       .data = .,
                       SNP,
                       A1,
                       TEST,
                       NMISS,
                       BETA,
                       SE,
                       L95,
                       U95,
                       STAT
                     ),
                 by = "SNP"
                 )
       } else {
         FileToReturn <- NULL
       }


       return(FileToReturn)

    },mc.cores = detectCores()
      ) %>%
  do.call(
    what = "rbind",
    args = .
    ) %>%
  # join the significant variants with variant annotation data
  left_join(
    x = .,
    y = read.delim(
                   file = "../annotation_output/allVariants_annovar_pharmGKB_and_VEP_data.tsv",
                   header = TRUE
                   ),
    by = c("SNP"="VariantID_fromInputFile")
    )

# save the significant variants with annotation data added to the file system
SignificantVariantsFromGWAS %>%
  write_delim(
    x = .,
    file = "./SignificantVariantsFromGWAS.txt",
    delim = "\t"
    )

# save a file of variant IDs only for the significant variants for filtering PLINK files
SignificantVariantsFromGWAS %>%
  select(
    .data = .,
    SNP
  ) %>%
  write_delim(
    x = .,
    file = "./SignificantVariantsFromGWAS_forFilteringPLINKfiles.txt",
    col_names = FALSE
    )


# load each of the clump files and filter to variants with the top 10 most significant P-values 
# for each individual phenotype that was tested
Top10MostSignificantVariantsForEachPhenotypeFromGWAS <-
  FilesToLoad %>%
  mclapply(
    X = .,
    FUN = function(currentFile)
    {
      print(x = currentFile)
      # identify the current phenotype by parsing the clump file name
      currentPhenotype <-
        currentFile %>%
        gsub(
          pattern = "(clump_file_additive_assoc.|clump_file_logistic_assoc.)",
          replacement = "",
          x = .
        ) %>%
        gsub(
          pattern = ".clumped",
          replacement = "",
          x = .
        )

      GWASresultsFileToLoad <-
        # load the GWAS results file for the current phenotype to obtain
        # regression summary statistics
        list.files(path = "./") %>%
          grepl(
            pattern = "(*linear$|*logistic$)",
            x = .
            ) %>%
          list.files(path = "./")[.] %>%
          grepl(
            pattern = currentPhenotype,
            x = .
            ) %>%
          list.files(path = "./")[grepl(pattern = "(*linear$|*logistic$)",x = list.files(path = "./"))][.]

       # load the clump file for the current phenotype
       LoadedClumpFile <-
          read.table(
            file = glue("./{currentFile}"),
            header = TRUE,
            quote = ""
            ) 

       FilteredClumpFile <-
           LoadedClumpFile %>%
           mutate(
             .data = .,
             P = P %>% as.numeric(x = .)
             ) %>%
           # sort the results by P-value in ascending order
           arrange(
             .data = .,
             P
             ) %>%
           # select columns with the variant ID and P-value
           select(
             .data = .,
             SNP,
             P
             ) %>%
           # create a column to label the current phenotype
           mutate(
             .data = .,
             Phenotype = currentPhenotype
             ) %>%
           # create a column with the P-value significance threshold
           mutate(
             .data = .,
             alpha = as.numeric(x = 0.05/nrow(x = LoadedClumpFile))
             ) %>%
           # filter to the top 10 most significant P-values
           head(x = .,10)

       # load the GWAS results file for the current phenotype
       LoadedGWASresultsFile <-
         read.table(
           file = paste0("./",GWASresultsFileToLoad),
           header = TRUE,
           sep = ""
           ) %>%
         mutate(
           .data = .,
           P = P %>% as.numeric(x = .)
           ) 
       
       print(x = GWASresultsFileToLoad)
       
       # if the files are linear association files,
       # select specific linear columns from the LoadedGWASresultsFile
       if(grepl(pattern = "*linear$",x = GWASresultsFileToLoad)==TRUE)
       {
           # join the clump file data with the GWAS summary statistics data
           # this will return a file of the 10 most significant SNPs after clumping with
           # GWAS summary statistics added
           FileToReturn <-
             FilteredClumpFile %>%
             left_join(
               x = .,
               y = LoadedGWASresultsFile %>%
                 select(
                   .data = .,
                   SNP,
                   A1,
                   TEST,
                   NMISS,
                   "BETA_or_OR" = BETA,
                   SE,
                   L95,
                   U95,
                   STAT
                 ),
               by = "SNP"
             )
       }
       
       # if the files are logistic association files,
       # select specific logistic columns from the LoadedGWASresultsFile
       if(grepl(pattern = "*logistic$",x = GWASresultsFileToLoad)==TRUE)
       {
           # join the clump file data with the GWAS summary statistics data
           # this will return a file of the 10 most significant SNPs after clumping with
           # GWAS summary statistics added
           FileToReturn <-
             FilteredClumpFile %>%
             left_join(
               x = .,
               y = LoadedGWASresultsFile %>%
                 select(
                   .data = .,
                   SNP,
                   A1,
                   TEST,
                   NMISS,
                   "BETA_or_OR" = OR,
                   SE,
                   L95,
                   U95,
                   STAT
                 ),
               by = "SNP"
             )
       }

       return(FileToReturn)

    },mc.cores = detectCores()
      ) %>%
  do.call(
    what = "rbind",
    args = .
    ) %>%
  # join the top 10 most significant variants with variant annotation data
  left_join(
    x = .,
    y = read.delim(
                   file = "../annotation_output/allVariants_annovar_pharmGKB_and_VEP_data.tsv",
                   header = TRUE
                   ),
    by = c("SNP"="VariantID_fromInputFile")
    )

# save the top 10 most significant variants for each phenotype to the file system
Top10MostSignificantVariantsForEachPhenotypeFromGWAS %>%
write_delim(
  x = .,
  file = "./Top10MostSignificantVariantsForEachPhenotypeFromGWAS.txt",
  delim = "\t"
  )

