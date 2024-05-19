
# load necesssary packages
source("../scripts/load_necessary_packages.R")

# load the de-identified phenotype data obtained from electronic health records (EHRs)
PhenotypeData <-
  read_excel(
    path = "../PhenotypeData/EHR_phenotypes.xlsx",
    sheet = 1,
    na = c("NA")
    )

PLINKphenotypeFile <-
  # select the columns that will be used as the major endpoints of the study
  PhenotypeData %>%
    select(
      .data = .,
      RandomID,
      GAF_admit,
      GAF_discharge,
      Risperidone_discontinued,
      PatientReadmitted,
      DaysToReadmittance,
      ReadmittedOnRisperidone,
      DaysOnRisperidone,
      LengthOfStay,
      "MaxDoseRisperidone_mg" = MaxDoseRisperidone_mgPerDay
      # "Rx_acute_dose_mg_risperidone" = Rx_acute_dose_mg_risperidone_35636,
      # "Rx_residential_dose_mg_risperidone" = Rx_residential_dose_mg_risperidone_35636
    ) %>%
  # loop through each individual column
  lapply(
    X = .,
    FUN = function(currentColumn)
    {
      ColumnToReturn <-
        # loop through each element of each column
        currentColumn %>%
        lapply(
          X = .,
          FUN = function(currentElement)
          {
            
            # if the currentElement is a "y" for yes, convert to a character 1,
            # if the currentElement is a "n" for no, convert to a character 0,
            # otherwise, leave as the currentElement
            if(!is.na(x = currentElement))
            {
              if(currentElement == "y")
              {
                ValueToReturn <- as.character(x = 1)
              } else if (currentElement == "n"){
                ValueToReturn <- as.character(x = 0)
              } else {
                ValueToReturn <- currentElement
              }
            } else {
              ValueToReturn <- currentElement
            }
            
              return(ValueToReturn)
          }
            ) %>%
          unlist(x = .) 
      
      return(ColumnToReturn)
    }
      ) %>%
  as.data.frame(x = .) %>%
  # # identify the max dose of risperidone in acute and residential stays
  # mutate(
  #   .data = .,
  #   max_acute_risperidone_dose_mg = Rx_acute_dose_mg_risperidone %>%
  #                                   lapply(
  #                                     X = .,
  #                                     FUN = function(currentDose)
  #                                     {
  #                                       # if the current dose is missing,
  #                                       # convert to a zero character value for 0 mg of the dose
  #                                       if(is.na(x = currentDose))
  #                                       {
  #                                         currentDose <-
  #                                           as.character(x = 0)
  #                                       }
  # 
  #                                       ValueToReturn <-
  #                                         currentDose %>%
  #                                           # split the sequence of doses based on the semi-colon
  #                                           str_split(string = .,pattern = ";") %>%
  #                                           # unlist the result
  #                                           unlist(x = .) %>%
  #                                           # convert to a numeric
  #                                           as.numeric(x = .) %>%
  #                                           # calculate the max dose
  #                                           max(.,na.rm = TRUE)
  # 
  #                                       return(ValueToReturn)
  #                                     }
  #                                       ) %>%
  #                                    unlist(x = .),
  #   max_residential_risperidone_dose_mg = Rx_residential_dose_mg_risperidone %>%
  #                                         lapply(
  #                                           X = .,
  #                                           FUN = function(currentDose)
  #                                           {
  #                                             # if the current dose is missing,
  #                                             # convert to a zero character value for 0 mg of the dose
  #                                             if(is.na(x = currentDose))
  #                                             {
  #                                               currentDose <-
  #                                                 as.character(x = 0)
  #                                             }
  #                                   
  #                                             ValueToReturn <-
  #                                               currentDose %>%
  #                                                 # split the sequence of doses based on the semi-colon
  #                                                 str_split(string = .,pattern = ";") %>%
  #                                                 # unlist the result
  #                                                 unlist(x = .) %>%
  #                                                 # convert to a numeric
  #                                                 as.numeric(x = .) %>%
  #                                                 # calculate the max dose
  #                                                 max(.,na.rm = TRUE)
  #                                   
  #                                             return(ValueToReturn)
  #                                           }
  #                                             ) %>%
  #                                          unlist(x = .)
  #   ) %>%
  mutate(
    .data = .,
    # # create a max dose of risperidone column by finding the pair-wise max (pmax) of the max doses over acute and residential stays
    # MaxDoseRisperidone_mg = pmax(
    #                             max_acute_risperidone_dose_mg,
    #                             max_residential_risperidone_dose_mg,
    #                             na.rm = TRUE
    #                             ),
    # create a change in GAF column as a percentage over baseline by 
    # subtracting discharge and admit gaf and dividing by the admit GAF
    ChangeInGAF = ((as.numeric(x = GAF_discharge) - as.numeric(x = GAF_admit))/as.numeric(x = GAF_admit))*100
    ) %>%
  # if the ChangeInGAF is a negative number, sensor the number as a zero since PLINK cannot handle 
  # negative phenotype values
  mutate(
    .data = .,
    ChangeInGAF = ChangeInGAF %>%
                  as.numeric(x = .) %>%
                  lapply(
                    X = .,
                    FUN = function(currentValue)
                    {
                      
                        if(is.na(x = currentValue)==FALSE)
                        {
                            if(currentValue < 0)
                            {
                              ValueToReturn <-
                                as.numeric(x = 0)
                            } else {
                              ValueToReturn <- as.numeric(x = currentValue)
                            }
                        } else {
                          ValueToReturn <- as.numeric(x = currentValue)
                        }
                      
                      return(ValueToReturn)
                    }
                      ) %>%
                   unlist(x = .)
    ) %>%
  # # deselect columns that are no longer needed
  # select(
  #   .data = .,
  #   -Rx_acute_dose_mg_risperidone,
  #   -Rx_residential_dose_mg_risperidone,
  #   -max_acute_risperidone_dose_mg,
  #   -max_residential_risperidone_dose_mg
  #   )  %>%
  # convert each column to a numeric and replace missing values with a -9,
  # which is the missing value indicator for PLINK
  lapply(
    X = .,
    FUN = function(currentColumn)
    {
          ArrayToReturn <-
            currentColumn %>%
            # convert the current column to a numeric
            as.numeric(x = .) %>%
            # convert each missing value in the column to a -9
            sapply(
              X = .,
              FUN = function(currentValue)
              {
                if(is.na(x = currentValue))
                {
                  ValueToReturn <-
                    as.numeric(x = -9)
                } else {
                  ValueToReturn <- currentValue
                }

                return(ValueToReturn)
              },
              simplify = TRUE
                )

            return(ArrayToReturn)
    }
      ) %>%
  as.data.frame(x = .) %>%
  # create and FID and IID column with the random participant IDs as required for PLINK phenotype files
  mutate(
    .data = .,
    FID = paste0("SCH-",RandomID),
    IID = paste0("SCH-",RandomID)
    ) %>%
  select(
    .data = .,
    -RandomID
    ) %>%
  select(
    .data = .,
    FID,
    IID,
    everything()
    )

# save a plink phenotype file of continuous outcomes to the file system
PLINKphenotypeFile %>%
  select(
    .data = .,
    FID,
    IID,
    GAF_admit,
    GAF_discharge,
    ChangeInGAF,
    DaysToReadmittance,
    MaxDoseRisperidone_mg,
    LengthOfStay,
    DaysOnRisperidone
    ) %>%
  write_tsv(
    x = .,
    file = "./PLINKphenotypeFile_continuous.txt"
    )

# save a plink phenotype file of binary outcomes to the file system
PLINKphenotypeFile %>%
  select(
    .data = .,
    FID,
    IID,
    Risperidone_discontinued,
    PatientReadmitted,
    ReadmittedOnRisperidone
    ) %>%
  write_tsv(
    x = .,
    file = "./PLINKphenotypeFile_binary.txt"
    )

