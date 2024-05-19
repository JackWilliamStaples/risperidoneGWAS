# this is a script for preparing manhattan and qq-plots of linear and logistic GWAS association results with the qqman package

# load necessary R packages 
source(file = "../scripts/load_necessary_packages.R")

# load the table of rs-numbers and variant IDs for annotating manhattan plots
VariantIDrsNumberMap <-
    read.delim(
      file = "../annotation_output/variantID_rsNumber_map.tsv",
      header = TRUE
      )

# load the table of variants that were significant in the GWAS
# after Bonferroni adjustment for the number of causative, independent snps
SignificantVariantsFromGWAS <-
  read.delim(
    file = "./SignificantVariantsFromGWAS.txt",
    header = TRUE
    )

# identify linear and logistic gwas results files
# by identifying files that end with "linear" or "logistic"
FilesToLoad <-
    list.files(path = "./") %>%
      grepl(
        pattern = "(*linear$|*logistic$)",
        x = .
      ) %>%
      list.files(path = "./")[.]

# create a directory for storing manhattan plots if it doesn't exist
if(!dir.exists(paths = "./manhattanPlots/"))
{
    dir.create(path = "./manhattanPlots/")
}

# create a directory for storing qq-plots if it doesn't exist
if(!dir.exists(paths = "./qqPlots/"))
{
    dir.create(path = "./qqPlots/")
}

  FilesToLoad %>%
  # load each GWAS results file
  mclapply(
    X = .,
    FUN = function(currentFile)
    {
        LoadedFile <-
          read.table(
            file = paste0("./",currentFile),
            header = TRUE,
            quote = ""
            )

        return(LoadedFile)
    },mc.cores = detectCores()
      ) %>%
  # create a manhattan plot for each phenotype
  mcmapply(
    FUN = function(currentGWASresult,currentPhenotype)
    {
            # identify the current phenotype by parsing the phenotype association file name
            currentPhenotype <-
              currentPhenotype %>%
              gsub(
                x = .,
                pattern = "(binary_phenotypes.|continuous_phenotypes.)",
                replacement = ""
              ) %>%
              gsub(
                x = .,
                pattern = "(.assoc.logistic|.assoc.linear)",
                replacement = ""
              )
            
            # identify the alpha value for the significance cutoff for the current phenotype
            AlphaValue <-
              SignificantVariantsFromGWAS %>%
                # filter to the current phenotype
                filter(
                  .data = .,
                  Phenotype == currentPhenotype
                  ) %>%
                pull(
                  .data = .,
                  alpha
                  ) %>%
                # extract the alpha value
                unique(x = .) %>%
                as.numeric(x = .)
      
            # create a tiff manhattan plot for the current phenotype
            tiff(
              filename = glue("./manhattanPlots/{currentPhenotype}.tiff"),
              compression = "none",
              units = "in",
              width = 7,
              height = 5,
              res = 300
              )
  
  
            currentGWASresult <-
                currentGWASresult %>%
                # select required columns for making a manhattan plot
                select(
                  .data = .,
                  CHR,
                  BP,
                  SNP,
                  P
                  ) %>%
              # make sure the P-value column is a numeric
              mutate(
                .data = .,
                P = P %>% as.numeric(x = .)
                ) %>%
              # remove variants with missing P-values
              filter(
                .data = .,
                !is.na(x = P)
                ) %>%
              # join the currentGWASresult table with the variant ID rs-number map
              left_join(
                x = .,
                y = VariantIDrsNumberMap,
                by = c("SNP" = "VariantID_fromInputFile")
                  ) %>%
              # replace the variant ID with the rs-number if the rs-number is not missing
              mutate(
                .data = .,
                SNP = if_else(
                              condition = (rsNumber_VEP=="*rsNA") | (is.na(x = rsNumber_VEP)==TRUE),
                              true = SNP,
                              false = rsNumber_VEP
                                )
                ) %>%
              # remove the rsNumber_VEP column now that it is no longer needed
              select(
                .data = .,
                -rsNumber_VEP
                )
            
                currentGWASresult %>%
              # create the manhattan plot
                manhattan(
                  x = .,
                  chr = "CHR",
                  bp = "BP",
                  snp = "SNP",
                  p = "P",
                  # alternate colors for each chromosome
                  col = c("mediumvioletred","darkcyan"),
                  # set the suggestive significance line as the alpha value
                  # if there are any significant variants for the current phenotype
                  suggestiveline = if(length(x = AlphaValue)>0)
                                   {
                                       as.numeric(x = -log10(x = AlphaValue))
                                   } else {
                                       FALSE
                                   },
                  # annotate data points with variants that are significant
                  # if there are any significant variants for the current phenotype
                  annotatePval = 
                                  # if(length(x = AlphaValue)>0)
                                  # {
                                  #    as.numeric(x = AlphaValue)
                                  # } else {
                                    FALSE,
                                  #},
                  annotateTop = FALSE,
                  # label the title of the graph with the current phenotype
                  main = currentPhenotype,
                  # adjust y-axis limits
                  ylim = c(0,9)
                    )
  
            dev.off()
  
  
            # create a tiff qq-plot
            tiff(
              filename = glue("./qqPlots/{currentPhenotype}.tiff"),
              compression = "none",
              units = "in",
              width = 7,
              height = 5,
              res = 300
              )
  
            # create the qq-plot with the P-values
            qq(
              pvector = currentGWASresult$P,
              # label the title of the graph with the current phenotype
              main = currentPhenotype
              )
  
            dev.off()
        
    },
    .,
    FilesToLoad,
    SIMPLIFY = FALSE,mc.cores = detectCores()
      )



