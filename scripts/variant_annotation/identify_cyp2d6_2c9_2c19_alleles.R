# This is a script for filtering ANNOVAR and VEP variant annotation data to 
# variants in the cyp2d6, cyp2c9, and cyp2c19 loci

# load necessary packages
source("../scripts/load_necessary_packages.R")
source("../scripts/pipeline_functions.R")


  # load the variant annotation data
  read.delim(
    file = "../annotation_output/allVariants_annovar_pharmGKB_and_VEP_data.tsv",
    header = TRUE,
    sep = "\t"
      ) %>%
    # filter the annotation data to variants in cyp2d6, 2c9, and 2c19 gene regions
  filter(
    .data = .,
    grepl(pattern = "(CYP2D6|CYP2C9|CYP2C19)",x = Gene.refGeneWithVer)==TRUE |
    grepl(pattern = "(CYP2D6|CYP2C9|CYP2C19)",x = Gene_VEP)==TRUE
    ) %>%
  # select relevant columns
  select(
    .data = .,
    rsNumber_VEP,
    Chr,
    Start,
    End,
    Ref,
    Alt,
    Func.refGeneWithVer,
    Gene.refGeneWithVer,
    Gene_VEP,
    starAlleles_pharmVar,
    pharmGKB_EvidenceLevel
    ) %>%
    # save the results to the file system
    write_delim(
      x = .,
      file = "../annotation_output/cyp2d6_2c9_2c19_alleles.txt",
      delim = "\t"
        )
  
  