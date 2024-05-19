# this is a script for storing functions that are used in the psychiatric GWAS pipeline

##################################################################################################################################
  ### define a function for adding levels of significance to P-value columns with astericks
    # this function is designed to be used with tidyverse mutate on an input column vector
  AddAstericksToPvalues <-
    function(
      columnVector = NULL # a column vector from a dataframe
      )
    {
        VectorToReturn <-
          columnVector %>%
          # convert to a numeric so the if-else statements below can be evaluated
          as.numeric(x = .) %>%
          lapply(
            X = .,
            FUN = function(currentPvalue)
            {
                # if the current P-value is not NA, add astericks for level of significance,
                # otherwise, return an NA
                if(!is.na(x = currentPvalue))
                {
                     # if the P-value is less than 0.001, add 3 astericks
                      if (currentPvalue < 0.001) {
                        ValueToReturn <-
                          paste0(currentPvalue,"***")
                      } # if the P-value is less than 0.01, add 2 astericks
                      else if (currentPvalue < 0.01) {
                        ValueToReturn <-
                          paste0(currentPvalue,"**")
                      } # if the P-value is less than 0.05, add 1 astericks
                      else if (currentPvalue < 0.05)
                    {
                      ValueToReturn <-
                        paste0(currentPvalue,"*")
                    } else {
                      ValueToReturn <- as.character(x = currentPvalue)
                    }

                } else {
                  ValueToReturn <- as.character(x = NA)
                }

                 return(ValueToReturn)
          }
          ) %>%
          unlist(x = .)

        return(VectorToReturn)
    }

#########################################################################################################################################################################

FilterToExonicVariantsFromVEPannotation <-
  function(
          InputDataFrame = NULL, # a dataframe that is the result of annotating a VCF file with the variant ensemble effect predictor (VEP)
          NameOfVEPconsequenceColumn = NULL # the name of the VEP consequence column as a string
          ) 
  {
    
    ### this is a function for filtering a dataframe that is the result of annotating a VCF file with the variant ensemble effect predictor (VEP)
    ### to exonic variants only based on the VEP consequences, for more information on VEP consequences, see: https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html 
   
    DataFrameToReturn <-
    InputDataFrame %>%
      # filter the data to exonic variants only
      filter(
        .data = .,
          (eval(expr = parse(text = NameOfVEPconsequenceColumn)) == "frameshift_variant") |
          (eval(expr = parse(text = NameOfVEPconsequenceColumn)) == "frameshift_variant,splice_region_variant") |
          (eval(expr = parse(text = NameOfVEPconsequenceColumn)) == "inframe_deletion") |
          (eval(expr = parse(text = NameOfVEPconsequenceColumn)) == "missense_variant") |
          (eval(expr = parse(text = NameOfVEPconsequenceColumn)) == "missense_variant,splice_region_variant") |
          (eval(expr = parse(text = NameOfVEPconsequenceColumn)) == "splice_region_variant,synonymous_variant") |
          (eval(expr = parse(text = NameOfVEPconsequenceColumn)) == "start_lost") |
          (eval(expr = parse(text = NameOfVEPconsequenceColumn)) == "stop_gained") |
          (eval(expr = parse(text = NameOfVEPconsequenceColumn)) == "synonymous_variant")
          )
    
    return(DataFrameToReturn)
  }

#########################################################################################################################################################################

RemoveAlternateVariantIDsAndReplaceMissingRSnumbers <-
  function(
    InputDataFrame = NULL, # a dataframe that is the result of annotating a VCF file with the variant ensemble effect predictor (VEP))
    rsIDcolumnName = NULL # a string with the column name of the rs-ids to update
          )
  {
    
          ###### this is a function for removing alternate variant identifiers from the VEP annotation column of existing variant identifiers, leaving only rs-ids
          ###### this function also replaces any missing rs-identifiers with an "*rsNA"
          DataFrameToReturn <-
          InputDataFrame %>%
          # remove any cosmic variant IDs (COSV) or CM variant IDs if there is an rs-number present
          mutate(
            .data = .,
            !!rsIDcolumnName := eval(expr = parse(text = rsIDcolumnName)) %>%
                                lapply(
                                  X = .,
                                  FUN = function(currentString)
                                  {
                                    # if the current string contains "rs", remove "COSV" or "CM" or "CS" variant IDs and any commas
                                    if(grepl(pattern = "rs",x = currentString)==TRUE)
                                    {
                                      StringToReturn <-
                                        currentString %>%
                                        str_remove_all(string = .,pattern = "COSV[:digit:]+") %>%
                                        str_remove_all(string = .,pattern = "CM[:digit:]+") %>%
                                        str_remove_all(string = .,pattern = "CS[:digit:]+") %>%
                                        str_remove_all(string = .,pattern = "CR[:digit:]+") %>%
                                        str_remove_all(string = .,pattern = "[:punct:]")
                                    } else {
                                      StringToReturn <- currentString
                                    }

                                    return(StringToReturn)
                                  }
                                ) %>%
                                unlist(x = .)
          ) %>%
          # if there is no rs-id in the existing_variant_VEP, change it to a "*rsNA" for no rs-number
          mutate(
            .data = .,
            !!rsIDcolumnName := if_else(
                                        condition = eval(expr = parse(text = rsIDcolumnName)) == "-",
                                        true = "*rsNA",
                                        false = eval(expr = parse(text = rsIDcolumnName))
                                        )
          )
          
          return(DataFrameToReturn)
  }


#########################################################################################################################################################################

### create a function for creating a tidy summary of the linear robust standard errors regression results
# that are output from the R lm_robust() function from the estimatrr package

TidyLinearRobustStandardErrorsRegressionResults <-
  function(
    DataForRegression = NULL, # the dataframe used for the linear regression with observations with missing values removed
    RegressionModel = NULL # the lm_robust() regression model object
  )
  {
      # summarize the regression model
      ModelSummary <-
        RegressionModel %>%
        summary(object = .)

      ResultsToReturn <-
        RegressionModel %>%
        # tidy the regression model with the etimatr tidy() function
        estimatr::tidy(x = .) %>%
        mutate(
          .data = .,
          # add the sample size
          N = DataForRegression %>% na.omit(object = .) %>% nrow(x = .),
          # add the unadjusted and adjusted R-squared
          R2 = ModelSummary$r.squared %>% round(x = .,digits = 3),
          R2_adj = ModelSummary$adj.r.squared %>% round(x = .,digits = 3),
          # add a column with the estimate and standard error pasted together,
          estimate_std.error = paste(
                                    signif(x = estimate,digits = 3),
                                    "±",
                                    signif(x = std.error,digits = 3)
                                    ),
          # round the p.value and add astericks for level of significance
          p.value = p.value %>% signif(x = .,digits = 3) %>% AddAstericksToPvalues(columnVector = .),
          # add the F-statistic
          F_stat = ModelSummary$fstatistic[["value"]] %>% signif(x = .,digits = 3),
          # add the F-statistic p-value by using the pf() function and add astericks for level of significance
          F_stat_p_value = pf(
                              q = ModelSummary$fstatistic["value"],
                              df1 = ModelSummary$fstatistic["numdf"],
                              df2 = ModelSummary$fstatistic["dendf"],
                              lower.tail = FALSE
                              ) %>%
                            signif(
                              x = .,
                              digits = 3
                              ) %>%
                            AddAstericksToPvalues(columnVector = .)
        )


      return(ResultsToReturn)
  }


#########################################################################################################################################################################

    ### create a function for creating a tidy summary of the logistic regression results that are output from the R glm(family = "binomial") function

    TidyLogisticRegressionResults <-
      function(
        DataForRegression = NULL, # the dataframe used for the logistic regression with observations with missing values removed
        RegressionModel = NULL # the glm() regression model object
      )
      {
        
          # summarize the regression
          RegressionSummary <-
            RegressionModel %>%
            summary(object = .)

          ResultsToReturn <-
              # convert the regression coefficients to a dataframe
              RegressionSummary$coefficients %>%
              data.frame(
                .,
                row.names = NULL
                ) %>%
              # convert the rownames to an actual column
              mutate(
                .data = .,
                term = rownames(x = RegressionSummary$coefficients)
                ) %>%
              # rearrange and rename columns
              select(
                .data = .,
                term,
                "p.value" = Pr...z..,
                everything()
                ) %>%
              mutate(
                .data = .,
                # add a column with the sample size
                N = DataForRegression %>% na.omit(object = .) %>% nrow(x = ),
                # add a column with the null deviance
                Null.deviance = RegressionSummary$null.deviance,
                # add a column with the residual deviance
                Deviance = RegressionSummary$deviance,
                # add a column with the Akaike Information Criterion (AIC)
                AIC = RegressionSummary$aic,
                # add a column with the estimate and standard error pasted together
                estimate_std.error = paste0(signif(x = Estimate,digits = 3)," ± ",signif(x = Std..Error,digits = 3)),
                # round the p-value and add astericks with levels of significance
                p.value = p.value %>% signif(x = .,digits = 3) %>% AddAstericksToPvalues(columnVector = .),
                # Compute the P-value for significance of the overall model by comparing the model with predictors to the null regression model
                # with a distributed Chi-squared test with degrees of freedom = degrees of freedom of model with predictors - degrees of freedom on NULL model.
                # See: https://stats.idre.ucla.edu/mplus/dae/logit-regression/ for more details on this.
                ChiSq_p_value = with(
                                    RegressionModel,
                                    pchisq(
                                      null.deviance - deviance,
                                      df = df.null - df.residual,
                                      lower.tail = FALSE
                                      )
                                    ) %>%
                                  signif(x = .,digits = 3) %>%
                                  AddAstericksToPvalues(columnVector = .)
               ) %>%
              # calculate the odds-ratio and 95% confidence interval by exponentiating the model coefficients and the confidence interval
                  # bind the odds-ratio and 95% confidence intervals with the regression summary by column
              cbind(
                  .,
                  cbind(
                    "OR" = coef(object = RegressionModel),
                    confint(object = RegressionModel)
                    ) %>%
                  exp(x = .) %>%
                  data.frame(
                    .,
                    row.names = NULL
                    )
                ) %>%
              # rename the 2.5 % and 97.5 % confidence interval columns
              dplyr::rename(
                .data = .,
                CI_2.5 = X2.5..,
                CI_97.5 = X97.5..
                ) %>%
              # add a column with the odds-ratio and confidence intervals pasted together
              mutate(
                .data = .,
                OR_CI = paste0(signif(x = OR,digits = 3)," (",signif(x = CI_2.5,digits = 3),"-",signif(x = CI_97.5,digits = 3),")")
                ) %>%
              # filter out the intercept term
             filter(
               .data = .,
               term != "(Intercept)"
               )

          return(ResultsToReturn)
      }

#########################################################################################################################################################################

UpdateUGT1AfamilyGeneNames <-
  function(
    InputDataFrame = NULL, # a dataframe that is the result of annotating a VCF file with the variant ensemble effect predictor (VEP)
    VEPgeneColumnName = NULL # a string with the name of the VEP gene symbol column
    )
  {

    # replace any occurrence of "UGT1A1,3-10" in the gene column with "UGT1A" since UGT1A1 was not truly sequenced
    # UGT1A4 was actually sequenced but the UGT1A genes (UGT1A1,3-10) have overlapping transcripts so the whole region is annotated as UGT1A1
    DataFrameToReturn <-
      InputDataFrame %>%
      mutate(
            .data = .,
            !!VEPgeneColumnName := if_else(
                                          condition = grepl(
                                                            pattern = "UGT1A",
                                                            x = eval(expr = parse(text = VEPgeneColumnName))
                                                            )==TRUE,
                                          true = "UGT1A",
                                          false = eval(expr = parse(text = VEPgeneColumnName))
                                         )
            )

     return(DataFrameToReturn)
  }

#########################################################################################################################################################################
