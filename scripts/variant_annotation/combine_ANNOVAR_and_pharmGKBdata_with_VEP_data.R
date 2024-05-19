# this is a script for combining ANNOVAR and pharmGKB annotated data with useful data from the VEP annotation output
# this script adds pharmGKB levels of evidence to the ANNOVAR annotated data based on rs-number

# load necessary packages
source("../scripts/load_necessary_packages.R")
source("../scripts/pipeline_functions.R")

# load the ANNOVAR annotated data with ADME optimized prediction framework scores computed
ANNOVARfileNames <-
  list(
    "allVariants" = "../annotation_output/annovar_and_ADMEoptimizedPredictions.tsv"
  )

ANNOVARdata <-
  ANNOVARfileNames %>%
    lapply(
      X = .,
      FUN = function(currentFile)
      {
        LoadedFile <-
          read.delim(
            file = currentFile,
            header = TRUE
            )

        return(LoadedFile)
      }
        )

VEPannotationData <-
  # load the VEP annotation data for each chromosome
  c(1:23) %>%
  lapply(
    X = .,
    FUN = function(currentChromosome)
    {
      FileName <-
         glue("../VEP_annotation_output/chr{currentChromosome}.txt")

      return(FileName)
    }
      ) %>%
  lapply(
    X = .,
    FUN = function(currentFile)
    {
      
          # load the VEP file
          LoadedFile <-
            read.delim(
              file = currentFile
              ) %>%
              # select useful columns
              select(
                .data = .,
                "VariantID_VEP" = X.Uploaded_variation,
                "Location_VEP" = Location,
                GIVEN_REF,
                USED_REF,
                "Allele_VEP" = Allele,
                "Consequence_VEP" = Consequence,
                "IMPACT_rating_VEP" = IMPACT,
                "Gene_VEP" = SYMBOL,
                "Gene_code_VEP" = Gene,
                Feature_type:BIOTYPE,
                "exon_number_VEP" = EXON,
                "intron_number_VEP" = INTRON,
                HGVSc:CDS_position,
                "Protein_position_VEP" = Protein_position,
                "Amino_acids_VEP" = Amino_acids,
                Codons,
                "existing_variant_VEP" = Existing_variation,
                STRAND,
                SOURCE,
                "VAR_SYNONYMS_VEP" = VAR_SYNONYMS,
                "Pubmed_ID" = PUBMED,
                "PHENOTYPES_VEP" = PHENOTYPES,
                MOTIF_NAME,
                DOMAINS,
                "TRANSCRIPTION_FACTORS_VEP" = TRANSCRIPTION_FACTORS,
                "CADD_RAW_VEP" = CADD_RAW,
                "CADD_PHRED_VEP" = CADD_PHRED,
                "LoFtool_VEP" = LoFtool
                ) %>%
              # obtain one row for each variant ID
              distinct(
                .data = .,
                VariantID_VEP,
                .keep_all = TRUE
                ) %>%
              # create a duplicate copy of the existing variant ID column
              mutate(
                .data = .,
                rsNumber_VEP = existing_variant_VEP
                )

          # if there is at least one variant in the annotation data for the current chromosome,
          # remove alternate variant IDs and replace missing rs-numbers
          if(nrow(x = LoadedFile)>0)
          {
            LoadedFile <-
              LoadedFile %>%
              # but only with the rs-numbers in the existing variant ID column, if an rs-number exists
              RemoveAlternateVariantIDsAndReplaceMissingRSnumbers(
                InputDataFrame = .,
                rsIDcolumnName = "rsNumber_VEP"
                )
          }

          return(LoadedFile)

    }#,mc.cores = detectCores()-1
  ) %>%
  # bind each file for each chromosome together by row
  do.call(
    what = "rbind",
    args = .
    ) 

# load files downloaded from pharmVar with rs-numbers and star alleles
StarAlleleDataFromPharmVar <-
  # list the files in the pharmVar star allele data directory
  list.files(path = "../pharmVar_starAlleleData/") %>%
  # loop through each file in the pharmVar star allele data directory
  lapply(
    X = .,
    FUN = function(currentFolder)
    {
      FileToReturn <-
        glue("../pharmVar_starAlleleData/{currentFolder}/GRCh37/") %>%
        list.files(path = .) %>%
        paste0("../pharmVar_starAlleleData/",currentFolder,"/GRCh37/",.) %>%
        read.delim(
          file = .,
          header = TRUE,
          skip = 1
          )
        
      return(FileToReturn)
    }
      ) %>%
  # bind each file together by row
  do.call(
    what = "rbind",
    args = .
    ) 
  

# load the table of clinical pharmGKB annotations with clinical levels of evidence for a variant influencing drug metabolism or drug response
pharmGKBevidenceData <-
  read.delim(
    file = "./pharmGKB_annotation_file/clinical_annotations.tsv",
    header = TRUE
      ) %>%
  # separate multiple entries in the Variant.Haplotypes column into multiple rows based on a comma
  separate_rows(
    data = .,
    Variant.Haplotypes,
    sep = ", "
    ) %>%
  # select the Variant.Haplotypes and Level.of.Evidence columns only
  select(
    .data = .,
    Variant.Haplotypes,
    Level.of.Evidence
    ) %>%
    # group by the Variant.Haplotypes and split into a list of dataframes based on Variant.Haplotypes
    group_by(
      .data = .,
      Variant.Haplotypes
    ) %>%
    group_split(.tbl = .) %>%
    # collapse the Variant.Haplotypes and levels of evidence into a single row
    # after sorting the level(s) of evidence in numeric order
    lapply(
      X = .,
      FUN = function(currentDataFrame)
      {
        TableToReturn <-
          tibble(
            "Variant.Haplotypes" = unique(x = currentDataFrame$Variant.Haplotypes),
            "pharmGKB_EvidenceLevel" = currentDataFrame$Level.of.Evidence %>%
              unique(x = .) %>%
              sort(x = .) %>%
              paste(.,collapse = " | ")
          )
        
        return(TableToReturn)
      }
    ) %>%
    # bind the list of dataframes back together by row
    do.call(
      what = "rbind",
      args = .
    ) %>%
    # replace any star alleles in the Variant.Haplotypes with an rs-number if the star allele is present in the StarAlleleDataFromPharmVar table
    mutate(
      .data = .,
      Variant.Haplotypes = Variant.Haplotypes %>%
                           lapply(
                             X = .,
                             FUN = function(currentVariant)
                             {
                                 # if the current variant is not an rs-number, meaning it is a star allele,
                                 # replace the star allele with an rs-number, if the star allele is present in the StarAlleleTable,
                                 # if the current variant is an rs-number, leave it as is
                                 if(grepl(pattern = "rs",x = currentVariant)==FALSE)
                                 {
                                     # if the current variant is in the star allele table, replace the star allele with the rs-number,
                                     # otherwise, leave the star allele as is
                                     if(any(StarAlleleDataFromPharmVar$Haplotype.Name == currentVariant)==TRUE)
                                     {
                                         VariantToReturn <-
                                           StarAlleleDataFromPharmVar %>%
                                           filter(
                                             .data = .,
                                             Haplotype.Name == currentVariant
                                             ) %>%
                                           head(x = .,1) %>%
                                           pull(
                                             .data = .,
                                             rsID
                                             )
  
                                          # if the rs-number is missing, leave the star allele as is
                                          if(grepl(pattern = "rs",x = VariantToReturn)==FALSE)
                                          {
                                              VariantToReturn <- currentVariant
                                          }
  
                                     } else {
                                       VariantToReturn <- currentVariant
                                     }
  
                                 } else {
                                   VariantToReturn <- currentVariant
                                 }
  
                               return(VariantToReturn)
                             }
                             ) %>%
                            unlist(x = .)
      ) 


# join the ANNOVAR, pharmGKB, and pharmVar star allele files with the VEP annotations
# save the resulting file to the file system
mapply(
  FUN = function(currentANNOVARfile,currentVEPfile,currentFileName)
    {
        # create a directory for storing data if it doesn't exist
        if(!dir.exists("../annotation_output/"))
        {
          dir.create(path = "../annotation_output/")
        }
    
      FileToSave <-
        currentANNOVARfile %>%
          # join the annovar and vep files, keeping all of the variant IDs in the annovar file
          left_join(
            x = .,
            y = currentVEPfile,
            by = c("VariantID_fromInputFile" = "VariantID_VEP")
            ) %>%
          # join the resulting file with the pharmGKB evidence table based on the rs-numbers present in the VEP data
          left_join(
            x = .,
            y = pharmGKBevidenceData,
            by = c("rsNumber_VEP" = "Variant.Haplotypes")
            ) %>%
          # join the resulting file with the pharmVar star allele data based on rs-number to get all possible star alleles that correspond to the rs-number
          left_join(
            x = .,
            y = StarAlleleDataFromPharmVar %>%
                # filter out missing rs-IDs
                filter(
                  .data = .,
                  grepl(pattern = "rs",x = rsID)
                  ) %>%
                # select useful columns
                select(
                  .data = .,
                  Haplotype.Name,
                  rsID
                  ) %>%
                # group by the rsID
                group_by(
                  .data = .,
                  rsID
                  ) %>%
                # split into a list of data frames by rsID
                group_split(.tbl = .) %>%
                # loop through each dataframe
                lapply(
                  X = .,
                  FUN = function(currentDataFrame)
                  {

                    TableToReturn <-
                      tibble(
                        "rsID" = currentDataFrame$rsID %>% unique(x = .),
                        "starAlleles_pharmVar" = currentDataFrame$Haplotype.Name %>% unique(x = .) %>% paste(.,collapse = " | ")
                      )

                    return(TableToReturn)
                  }
                    ) %>%
                do.call(
                  what = "rbind",
                  args = .
                  ),
            by = c("rsNumber_VEP" = "rsID")
            )


      FileToSave %>%
          write_tsv(
            x = .,
            file = paste0("../annotation_output/",currentFileName,"_","annovar_pharmGKB_and_VEP_data.tsv")
              )
      
      # save a table of rs-numbers and variant IDs to be used as an easy lookup table
      FileToSave %>%
        select(
          .data = .,
          VariantID_fromInputFile,
          rsNumber_VEP
          ) %>%
        write_tsv(
          x = .,
          file = "../annotation_output/variantID_rsNumber_map.tsv"
            )
      
      return(NULL)

    },
  ANNOVARdata,
  list(VEPannotationData),
  names(x = ANNOVARdata),
  SIMPLIFY = FALSE
  )

