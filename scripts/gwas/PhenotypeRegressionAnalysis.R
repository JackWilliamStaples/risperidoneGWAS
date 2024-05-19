# This is a script for:
   # (1) computing basic summary statistics of risperidone treatment outcomes and demographic data
   # (2) creating histograms and box and whisker plots of continuous study outcomes to visualize outlying data points
   # (3) performing primary and multiple linear and logistic regressions of risperidone treatment outcomes
   #     versus demographic covariates and variant genotypes that were significant in GWAS analysis
   # (4) creating boxplots of continuous phenotypes stratified by variant genotype for variants that were significant in GWAS analysis
   # (5) performing primary Kaplan-Meier survival analysis of risperidone discontinuation
   # (6) performing primary and multiple Cox proportional hazard regression analysis of risperidone discontinuation

# load necessary R packages
source(file = "../scripts/load_necessary_packages.R")
# load pipeline functions
source(file = "../scripts/pipeline_functions.R")

# load the significant variants from the GWAS analysis
SignificantVariantsFromGWAS <-
  read.delim(
    file = "./SignificantVariantsFromGWAS.txt",
    header = TRUE
    ) %>%
  # create a column with the variant ID in the GTEx format
  mutate(
    .data = .,
    VariantID_GTEx_format = paste0("chr",CHROM_fromInputFile,"_",POS_fromInputFile,"_",Ref,"_",Alt,"_b37")
    ) %>%
  # select useful columns
  select(
    .data = .,
    SNP:STAT,
    Ref,
    Alt,
    Gene.refGeneWithVer,
    VariantID_GTEx_format,
    GTEx_V6p_gene,
    GTEx_V6p_tissue,
    dbscSNV_ADA_SCORE,
    dbscSNV_RF_SCORE,
    HGVSc,
    HGVSp,
    starts_with(match = "X1000g2015aug")
    ) %>%
  # create a column with the beta and standard error pasted together
  mutate(
    .data = .,
    BETA_SE = paste(BETA,"±",SE)
    )

# load the variant ID rs-number map for labeling variant IDs
VariantIDrsNumberMap <-
  read.delim(
    file ="../annotation_output/variantID_rsNumber_map.tsv",
    header = TRUE
    ) %>%
  # filter to variant IDs that were significant in the GWAS analysis
  filter(
    .data = .,
    VariantID_fromInputFile %in% SignificantVariantsFromGWAS$SNP
      )

# load the genotype matrix of the significant variants from the GWAS
GenotypeMatrix <-
  read.table(
    file = "./SignificantVariantsFromGWAS_geno_matrix.raw",
    header = TRUE,
    quote = ""
    ) %>%
  # remove unneccessary columns
  select(
    .data = .,
    -PAT,
    -MAT,
    -SEX,
    -PHENOTYPE,
    -IID,
    "ParticipantID" = FID
    )

# load the multidimensional scaling (MDS) components (population structure estimate) from
# the GWAS analysis
MDScomponents <-
  read.table(
    file = "../quality_control_working_directory/covar_mds.txt",
    header = TRUE,
    quote = ""
    ) %>%
   select(
     .data = .,
     -IID,
     "ParticipantID" = FID
     )

# load the continuous PLINK phenotype file that was used for association analysis
PhenotypeData_continuous <-
  read.table(
    file = "./PLINKphenotypeFile_continuous.txt",
    header = TRUE,
    quote = ""
    ) %>%
  select(
    .data = .,
    -IID,
    "ParticipantID" = FID
    ) %>%
  # change values that are equal to -9 (missing value code for PLINK) to a numeric NA
  mapply(
    FUN = function(currentColumn,currentColumnName)
    {
                # if the current column is not the ParticipantID column,
                # change values equal to -9 to missing
                if(currentColumnName!="ParticipantID")
                {
                  ValuesToReturn <-
                    currentColumn %>%
                    sapply(
                      X = .,
                      FUN = function(currentValue)
                      {
                          if(currentValue==as.numeric(x = -9))
                          {
                              ValueToReturn <-
                                as.numeric(x = NA)
                          } else {
                            ValueToReturn <-
                              as.numeric(x = currentValue)
                          }

                        return(ValueToReturn)

                      },
                      simplify = FALSE
                        ) %>%
                     unlist(x = .)

                } else {

                  ValuesToReturn <-currentColumn

                }

              return(ValuesToReturn)
    },
    .,
    names(x = .),
    SIMPLIFY = FALSE
      ) %>%
  data.frame()

# load the binary PLINK phenotype file that was used for association analysis
PhenotypeData_binary <-
  read.table(
    file = "./PLINKphenotypeFile_binary.txt",
    header = TRUE,
    quote = ""
    ) %>%
  select(
    .data = .,
    -IID,
    "ParticipantID" = FID
    ) %>%
  # change values that are equal to -9 (missing value code for PLINK) to a numeric NA
  mapply(
    FUN = function(currentColumn,currentColumnName)
    {
                # if the current column is not the ParticipantID column,
                # change values equal to -9 to missing
                if(currentColumnName!="ParticipantID")
                {
                  ValuesToReturn <-
                    currentColumn %>%
                    sapply(
                      X = .,
                      FUN = function(currentValue)
                      {
                          if(currentValue==as.numeric(x = -9))
                          {
                              ValueToReturn <-
                                as.numeric(x = NA)
                          } else {
                            ValueToReturn <-
                              as.numeric(x = currentValue)
                          }

                        return(ValueToReturn)

                      },
                      simplify = FALSE
                        ) %>%
                     unlist(x = .)

                } else {

                  ValuesToReturn <-currentColumn

                }

              return(ValuesToReturn)
    },
    .,
    names(x = .),
    SIMPLIFY = FALSE
      ) %>%
  data.frame()

# join the clinical data with the genotype matrix of significant variants
GenotypeMatrixAndClinicalData <-
  GenotypeMatrix %>%
  # join the genotype matrix with the continuous phenotype data by participant ID
    left_join(
      x = .,
      y = PhenotypeData_continuous,
      by = "ParticipantID"
      ) %>%
  # join the genotype matrix with the binary phenotype data by participant ID
    left_join(
      x = .,
      y = PhenotypeData_binary,
      by = "ParticipantID"
      ) %>%
   # join the result with the rest of the clinical EHR data
   left_join(
     x = .,
     y = read_excel(
                    path = "../PhenotypeData/EHR_phenotypes.xlsx",
                    sheet = 1,
                    na = c("NA")
                    ) %>%
          mutate(
            .data = .,
            "ParticipantID" = paste0("SCH-",RandomID)
            ) %>%
          select(
            .data = .,
            -RandomID,
            -GAF_admit,
            -GAF_discharge,
            -DaysToReadmittance,
            -LengthOfStay,
            -DaysOnRisperidone,
            -Risperidone_discontinued,
            -PatientReadmitted,
            -ReadmittedOnRisperidone
          ),
       by = "ParticipantID"
       ) %>%
  # join the results with the multidimensional scaling components
  # (estimates of population structure) from the GWAS analysis
  left_join(
    x = .,
    y = MDScomponents,
    by = "ParticipantID"
      ) %>%
   # select useful columns for analysis
   select(
     .data = .,
  ##### study participant IDs
        ParticipantID,
  ##### signficant variant genotypes from the GWAS analysis
        starts_with(match = "X"),
  ##### main outcomes
       ChangeInGAF,
       DaysToReadmittance,
       GAF_admit,
       GAF_discharge,
       LengthOfStay,
       MaxDoseRisperidone_mg,
       ReadmittedOnRisperidone,
       Risperidone_discontinued,
       PatientReadmitted,
  ##### max risperidone dose dosing information
      MaxDoseRisperidone_mgPerDose,
      MaxDoseRisperidone_frequency,
      MaxDoseRisperidone_mgPerDay,
  ##### days in acute and residential treatment
       acute_days,
       residential_days,
  ##### unit transfers
       UnitTransfers,
  ##### discharge unit
       Discharge_unit,
  ##### residential visit
       Residential_visit,
  ##### readmittance
     Readmitted_unit,
     ReadmittedRisperidoneDose_mg,
###### antipsychotic drugs
  ##### aripiprazole
       AripiprazoleAfterDischarge,
       DaysOnAripiprazole,
  ##### chlorpromazine
       ChlorpromazineAfterDischarge,
       DaysOnChlorpromazine,
  ##### haloperidol
       HaloperidolAfterDischarge,
       DaysOnHaloperidol,
  ##### olanzapine
       OlanzapineAfterDischarge,
       DaysOnOlanzapine,
  ##### paliperidone
       PaliperidoneAfterDischarge,
       DaysOnPaliperidone,
  ##### quetiapine
       QuetiapineAfterDischarge,
       DaysOnQuetiapine,
  ##### risperidone
       RisperidoneBeforeAdmit,
       DaysOnRisperidone,
  ##### ziprasidone
       ZiprasidoneAfterDischarge,
       DaysOnZiprasidone,
  ##### reason for risperidone discontinuation
       starts_with(match = "ReasonForRispDiscontinuation"),
  ##### reason for readmittance
       starts_with(match = "ReasonForReadmittance"),
    ##### baseline characteristics
       Age_years,
       Weight_kilograms,
       Height_centimeters,
       Gender_MorF,
       Insurance_medicaid,
       Insurance_private,
       GAF_admit,
       Admit_unit,
    ##### psychiatric history (inpatient, outpatient, and family)
       starts_with(match = "PsychiatricHistory"),
    ##### neglect and abuse
       NeglectAndAbuseHistory,
    ##### substance abuse
       DSMdischargeAx1_AlcoholAbuse,
       DSMdischargeAx1_MarijuanaDependence,
       DSMdischargeAx1_NarcoticsDependence_HO,
       DSMdischargeAx1_OpiateUseDisorder,
       DSMdischargeAx1_PolysubstanceAbuse,
       SubstanceAbuseHistory,
    ##### psychiatric diagnosis
       DSMdischargeAx1_ADHD_combinedType,
       DSMdischargeAx1_AnxietyDisorder_NOS,
       DSM5discharge_GeneralizedAnxietyDisorder,
       DSMdischargeAx1_AutismSpectrumDisorder,
       DSMdischargeAx1_BipolarDisorder,
       DSM5discharge_BipolarDisorder,
       DSMdischargeAx1_IntermittentExplosiveDisorder,
       DSMdischargeAx1_MajorDepressiveDisorder,
       DSMdischargeAx1_MoodDisorder_NOS,
       DSMdischargeAx1_ObsessiveCompulsiveDisorder,
       DSMdischargeAx1_OppositionalDefiantDisorder,
       DSMdischargeAx1_PostraumaticStressDisorder,
       DSM5discharge_PosttraumaticStressDisorder,
       DSMdischargeAx1_PsychoticSymptoms,
       DSM5discharge_Psychosis,
    ##### antipsychotic medications
         # aripiprazole
       "Rx_acute_aripiprazole" = Rx_acute_aripiprazole_89013,
       "Rx_residential_aripiprazole" = Rx_residential_aripiprazole_89013,
         # chlorpromazine
       "Rx_acute_chlorpromazine" = Rx_acute_chlorpromazine_2403,
       "Rx_acute_chlorpromazineHCl" = Rx_acute_chlorpromazinehydrochloride_104728,
       "Rx_residential_chlorpromazine" = Rx_residential_chlorpromazine_2403,
       "Rx_residential_chlorpromazineHCl" = Rx_residential_chlorpromazinehydrochloride_104728,
         # haloperidol
       "Rx_acute_haloperidol" = Rx_acute_haloperidol_5093,
       "Rx_residential_haloperidol" = Rx_residential_haloperidol_5093,
         # olanzapine
       "Rx_acute_olanzapine" = Rx_acute_olanzapine_61381,
       "Rx_residential_olanzapine" = Rx_residential_olanzapine_61381,
         # paliperidone
       "Rx_acute_paliperidone" = Rx_acute_paliperidone_679314,
       "Rx_residential_paliperidone" = Rx_residential_paliperidone_679314,
         # quetiapine
       "Rx_acute_quetiapine" = Rx_acute_quetiapine_51272,
       "Rx_acute_quetiapineFumarate" = Rx_acute_quetiapinefumarate_221153,
       "Rx_residential_quetiapine" = Rx_residential_quetiapine_51272,
       "Rx_residential_quetiapineFumarate" = Rx_residential_quetiapinefumarate_221153,
         # risperidone
       "Rx_acute_risperidone" = Rx_acute_risperidone_35636,
       "Rx_residential_risperidone" = Rx_residential_risperidone_35636,
         # ziprasidone
       "Rx_acute_ziprasidone" = Rx_acute_ziprasidone_115698,
       "Rx_residential_ziprasidone" = Rx_residential_ziprasidone_115698,
    ##### CYP2D6 substrates and inhibitors
        # acetaminophen
       "Rx_acute_acetaminophen" = Rx_acute_acetaminophen_161,
       "Rx_acute_acetaminophenButalbitalCaffeine" = Rx_acute_acetaminophenbutalbitalcaffeine_214130,
       "Rx_residential_acetaminophen" = Rx_residential_acetaminophen_161,
       "Rx_residential_acetaminophenAspirinCaffeine" = Rx_residential_acetaminophenaspirincaffeine_466584,
       "Rx_residential_acetaminophenHydrocodone" = Rx_residential_acetaminophenhydrocodone_214182,
        # amphetamine
       "Rx_acute_amphetamine" = Rx_acute_amphetamineaspartateamphetaminesulfatedextroamphetaminesaccharatedextroamphetaminesulfate_822929,
       "Rx_residential_amphetamine" = Rx_residential_amphetamineaspartateamphetaminesulfatedextroamphetaminesaccharatedextroamphetaminesulfate_822929,
        # lisdexamfetamine (pro-drug for amphetamine)
       "Rx_acute_lisdexamfetamine" = Rx_acute_lisdexamfetamine_700810,
       "Rx_acute_lisdexamfetamineDimesylate" = Rx_acute_lisdexamfetaminedimesylate_673579,
       "Rx_residential_lisdexamfetamineDimesylate" = Rx_residential_lisdexamfetaminedimesylate_673579,
        # amitriptyline
       "Rx_acute_amitriptyline" = Rx_acute_amitriptyline_704,
        # aripiprazole (see antipsychotic medications above)
        # atomoxetine
       "Rx_acute_atomoxetine" = Rx_acute_atomoxetine_38400,
       "Rx_acute_atomoxetineHCl" = Rx_acute_atomoxetinehydrochloride_353103,
       "Rx_residential_atomoxetine" = Rx_residential_atomoxetine_38400,
       "Rx_residential_atomoxetineHCl" = Rx_residential_atomoxetinehydrochloride_353103,
        # benzocaine
       "Rx_residential_benzocaine" = Rx_residential_benzocaine_1399,
        # bupropion
       "Rx_acute_bupropion" = Rx_acute_bupropion_42347,
       "Rx_acute_bupropionHCl" = Rx_acute_bupropionhydrochloride_203204,
       "Rx_residential_bupropion" = Rx_residential_bupropion_42347,
       "Rx_residential_bupropionHCl" = Rx_residential_bupropionhydrochloride_203204,
        # buspirone
       "Rx_acute_buspironeHCl" = Rx_acute_buspironehydrochloride_203116,
       "Rx_residential_buspironeHCl" = Rx_residential_buspironehydrochloride_203116,
        # chlorpheniramine
       "Rx_acute_chlorpheniraminePhenylephrine" = Rx_acute_chlorpheniraminephenylephrine_214393,
       "Rx_residential_chlorpheniraminePhenylephrine" = Rx_residential_chlorpheniraminephenylephrine_214393,
       "Rx_residential_chlorpheniraminePseudoephedrine" = Rx_residential_chlorpheniraminepseudoephedrine_214395,
        # chlorpromazine (see antipsychotic medications above)
        # cimetidine
       "Rx_residential_cimetidine" = Rx_residential_cimetidine_2541,
        # citalopram
       "Rx_acute_citalopram" = Rx_acute_citalopram_2556,
       "Rx_acute_citalopramHydrobromide" = Rx_acute_citalopramhydrobromide_221078,
       "Rx_residential_citalopramHydrobromide" = Rx_residential_citalopramhydrobromide_221078,
        # clonidine
       "Rx_acute_clonidine" = Rx_acute_clonidine_2599,
       "Rx_acute_clonidineHCl" = Rx_acute_clonidinehydrochloride_142432,
       "Rx_residential_clonidine" = Rx_residential_clonidine_2599,
       "Rx_residential_clonidineHCl" = Rx_residential_clonidinehydrochloride_142432,
        # clotrimazole
       "Rx_residential_betamethazoneClotrimazole" = Rx_residential_betamethasoneclotrimazole_106928,
        # codeine
       "Rx_acute_codeineGuaifenesin" = Rx_acute_codeineguaifenesin_214442,
        # dextromethorphan
       "Rx_acute_dextromethorphan" = Rx_acute_dextromethorphan_3289,
       "Rx_acute_dextromethorphanGuaifenesin" = Rx_acute_dextromethorphanguaifenesin_214488,
       "Rx_residential_dextromethorphan" = Rx_residential_dextromethorphan_3289,
       "Rx_residential_dextromethorphanGuaifenesin" = Rx_residential_dextromethorphanguaifenesin_214488,
        # diphenhydramine
       "Rx_acute_diphenhydramine" = Rx_acute_diphenhydramine_3498,
       "Rx_acute_diphenhydramineHCl" = Rx_acute_diphenhydraminehydrochloride_1362,
       "Rx_residential_diphenhydramine" = Rx_residential_diphenhydramine_3498,
       "Rx_residential_diphenhydramineHCl" = Rx_residential_diphenhydraminehydrochloride_1362,
        # duloxetine
       "Rx_acute_duloxetine" = Rx_acute_duloxetine_72625,
       "Rx_acute_duloxetineHCl" = Rx_acute_duloxetinehydrochloride_476250,
       "Rx_residential_duloxetine" = Rx_residential_duloxetine_72625,
        # escitalopram
       "Rx_acute_escitalopram" = Rx_acute_escitalopram_321988,
        # fluoxetine
       "Rx_acute_fluoxetine" = Rx_acute_fluoxetine_4493,
       "Rx_acute_fluoxetineHCl" = Rx_acute_fluoxetinehydrochloride_227224,
       "Rx_residential_fluoxetine" = Rx_residential_fluoxetine_4493,
       "Rx_residential_fluoxetineHCl" = Rx_residential_fluoxetinehydrochloride_227224,
        # fluvoxamine
       "Rx_acute_fluvoxamine" = Rx_acute_fluvoxamine_42355,
        # haloperidol (see antipsychotic medications above)
        # hydroxyzine
       "Rx_acute_hydroxyzine" = Rx_acute_hydroxyzine_5553,
       "Rx_acute_hydroxyzinePamoate" = Rx_acute_hydroxyzinepamoate_203182,
       "Rx_residential_hydroxyzine" = Rx_residential_hydroxyzine_5553,
       "Rx_residential_hydroxyzineHCl" = Rx_residential_hydroxyzinehydrochloride_154987,
       "Rx_residential_hydroxyzinePamoate" = Rx_residential_hydroxyzinepamoate_203182,
        # imipramine
       "Rx_acute_imipramine" = Rx_acute_imipramine_5691,
        # ketoconazole
       "Rx_acute_ketoconazole" = Rx_acute_ketoconazole_6135,
       "Rx_residential_ketoconazole" = Rx_residential_ketoconazole_6135,
        # lidocaine
       "Rx_residential_lidocainePrilocaine" = Rx_residential_lidocaineprilocaine_166283,
        # loperamide
       "Rx_acute_loperamide" = Rx_acute_loperamide_6468,
       "Rx_residential_loperamide" = Rx_residential_loperamide_6468,
        # loratadine
       "Rx_acute_loratadine" = Rx_acute_loratadine_28889,
       "Rx_residential_loratadine" = Rx_residential_loratadine_28889,
       "Rx_residential_loratadinePseudoephedrine" = Rx_residential_loratadinepseudoephedrine_214682,
        # miconazole
       "Rx_residential_miconazole" = Rx_residential_miconazole_6932,
       "Rx_residential_miconazoleNitrate" = Rx_residential_miconazolenitrate_42790,
        # mirtazapine
       "Rx_acute_mirtazapine" = Rx_acute_mirtazapine_15996,
        # nicotine
       "Rx_acute_nicotine" = Rx_acute_nicotine_7407,
       "Rx_residential_nicotine" = Rx_residential_nicotine_7407,
        # olanzapine (see antipsychotic medications above)
        # omeprazole
       "Rx_acute_omeprazole" = Rx_acute_omeprazole_7646,
       "Rx_residential_omeprazole" = Rx_residential_omeprazole_7646,
        # ondansetron
       "Rx_acute_ondansetronHCl" = Rx_acute_ondansetronhydrochloride_203148,
       "Rx_residential_ondansetron" = Rx_residential_ondansetron_26225,
        # oxybutynin
       "Rx_acute_oxybutyninChloride" = Rx_acute_oxybutyninchloride_54251,
       "Rx_residential_oxybutynin" = Rx_residential_oxybutynin_32675,
        # paliperidone (see antipsychotic medications above)
        # propanolol
       "Rx_acute_propranolol" = Rx_acute_propranolol_8787,
        # promethazine
       "Rx_residential_promethazineHCl" = Rx_residential_promethazinehydrochloride_142445,
        # quetiapine (see antipsychotic medications above)
        # ranitidine
       "Rx_acute_ranitidine" = Rx_acute_ranitidine_9143,
        # risperidone (see antipsychotic medications above)
        # sertraline
       "Rx_acute_sertraline" = Rx_acute_sertraline_36437,
       "Rx_acute_sertralineHCl" = Rx_acute_sertralinehydrochloride_155137,
       "Rx_residential_sertraline" = Rx_residential_sertraline_36437,
        # St. Johns Wort
       "Rx_residential_stJohnsWort" = Rx_residential_stjohnswortextract_258326,
        # terbinafine
       "Rx_acute_terbinafineHCl" = Rx_acute_terbinafinehydrochloride_235838,
       "Rx_residential_terbinafineHCl" = Rx_residential_terbinafinehydrochloride_235838,
        # tramadol
       "Rx_residential_tramadol" = Rx_residential_tramadol_10689,
        # trazodone
       "Rx_acute_trazodone" = Rx_acute_trazodone_10737,
       "Rx_acute_trazodoneHCl" = Rx_acute_trazodonehydrochloride_82112,
       "Rx_residential_trazodone" = Rx_residential_trazodone_10737,
       "Rx_residential_trazodoneHCl" = Rx_residential_trazodonehydrochloride_82112,
   ##### dopaminergic modulators
       # amphetamine and lisdexamfetamine (amphetamine pro-drug) (see cyp2d6 modulators above)
       # amantadine
       "Rx_acute_amantadine" = Rx_acute_amantadine_620,
       "Rx_residential_amantadine" = Rx_residential_amantadine_620,
       # aripiprazole (see cyp2d6 modulators above)
       # atomoxetine (see cyp2d6 modulators above)
       # benztropine
       "Rx_acute_benztropine" = Rx_acute_benztropine_1424,
       "Rx_acute_benztropineMesylate" = Rx_acute_benztropinemesylate_18927,
       "Rx_residential_benztropine" = Rx_residential_benztropine_1424,
       "Rx_residential_benztropineMesylate" = Rx_residential_benztropinemesylate_18927,
       # bupropion (see cyp2d6 modulators above)
       # buspirone (see cyp2d6 modulators above)
       # carbamazepine
       "Rx_acute_carbamazepine" = Rx_acute_carbamazepine_2002,
       # carbidopa levodopa
       "Rx_acute_carbidopaLevodopa" = Rx_acute_carbidopalevodopa_103990,
       # chlorpheniramine phenylephrine (see cyp2d6 modulators above)
       # chlorpheniramine pseudophedrine (see cyp2d6 modulators above)
       # chlorpromazine (see cyp2d6 modulators above)
       # duloxetine (see cyp2d6 modulators above)
       # escitalopram (see cyp2d6 modulators above)
       # haloperidol (see cyp2d6 modulators above)
       # imipramine (see cyp2d6 modulators above)
       # lamotrigine
       "Rx_acute_lamotrigine" = Rx_acute_lamotrigine_28439,
       "Rx_residential_lamotrigine" = Rx_residential_lamotrigine_28439,
       # methylphenidate/dexmethylphenidate
       "Rx_acute_dexmethylphenidateHCl" = Rx_acute_dexmethylphenidatehydrochloride_353105,
       "Rx_residential_dexmethylphenidateHCl" = Rx_residential_dexmethylphenidatehydrochloride_353105,
       "Rx_residential_dexmethylphenidate" = Rx_residential_dexmethylphenidate_352372,
       "Rx_acute_methylphenidate" = Rx_acute_methylphenidate_6901,
       "Rx_residential_methylphenidate" = Rx_residential_methylphenidate_6901,
       "Rx_acute_methylphenidateHCl" = Rx_acute_methylphenidatehydrochloride_203188,
       "Rx_residential_methylphenidateHCl" = Rx_residential_methylphenidatehydrochloride_203188,
       # nicotine (see cyp2d6 modulators above)
       # olanzapine (see cyp2d6 modulators above)
       # paliperidone (see cyp2d6 modulators above)
       # promethazine (see cyp2d6 modulators above)
       # pseudoephedrine (see chlorpheniraminePseudophedrine and loratidinePseudoephedrine in cyp2d6 modulators above)
       # quetiapine (see cyp2d6 modulators above)
       # risperidone (see cyp2d6 modulators above)
       # sertraline (see cyp2d6 modulators above)
       # ziprasidone (see antipsychotic medications above)
   ##### serotonergic modulators
       # amphetamine and lisdexamfetamine (pro-drug for amphetamine) (see cyp2d6 modulators above)
       # amantadine (see dopaminergic modulators above)
       # amitriptyline (see cyp2d6 modulators above)
       # aripiprazole (see cyp2d6 modulators above)
       # atomoxetine (see cyp2d6 modulators above)
       # benztropine (see dopaminergic modulators above)
       # bupropion (see dopaminergic modulators above)
       # buspirone (see dopaminergic modulators above)
       # carbidopia levodopa (see dopaminergic modulators above)
       # chlorpheniramine phenylephrine (see dopaminergic modulators above)
       # chlorpheniramine pseudoephedrine (see dopaminergic modulators above)
       # chlorpromazine (see dopaminergic modulators above)
       # citalopram (see cyp2d6 modulators above)
       # cyproheptadine
       "Rx_acute_cyproheptadine" = Rx_acute_cyproheptadine_3013,
       # dextromethorphan (see cyp2d6 modulators above)
       # duloxetine (see dopaminergic modulators above)
       # escitalopram (see dopaminergic modulators above)
       # fluoxetine (see cyp2d6 modulators above)
       # fluvoxamine (see cyp2d6 modulators above)
       # guaifenesin Dextromethorphan (dextromethorphan -- see cyp2d6 modulators above)
       # haloperidol (see dopaminergic modulators above)
       # imipramine (see dopaminergic modulators above)
       # lamotrigine (see dopaminergic modulators above)
       # methylphenidate/dexmethylphenidate (see dopaminergic modulators above)
       # mirtazapine (see cyp2d6 modulators above)
       # olanzapine (see cyp2d6 modulators above)
       # ondansetron (see cyp2d6 modulators above)
       # oxymetazoline
       "Rx_residential_oxymetazoline" = Rx_residential_oxymetazoline_7812,
       # paliperidone (see dopaminergic modulators above)
       # propranolol (see cyp2d6 modulators above)
       # pseudoephedrine (see dopaminergic modulators above)
       # quetiapine (see dopaminergic modulators above)
       # rizatriptan
       "Rx_residential_rizatriptan" = Rx_residential_rizatriptan_88014,
       # rocuronium bromide
       "Rx_residential_rocuroniumBromide" = Rx_residential_rocuroniumbromide_32521,
       # sertraline (see dopaminergic modulators above)
       # St. John's Wort (see cyp2d6 modulators above)
       # sumatriptan
       "Rx_residential_sumatriptan" = Rx_residential_sumatriptan_37418,
       # tramadol (see cyp2d6 modulators above)
       # trazodone (see cyp2d6 modulators above)
       # ziprasidone (see dopaminergic modulators above)
   ##### multidimensional scaling (MDS) components (estimates of population structure from genomic data)
       C1:C10
   ) %>%
   # convert the DSM diagnostic code columns to one-hot encoded vectors
   mutate_at(
     .tbl = .,
     .vars = vars(starts_with(match = "DSM")),
     .funs = function(currentColumn)
             {

               ArrayToReturn <-
                 currentColumn %>%
                 sapply(
                   X = .,
                   FUN = function(currentValue){

                       # if the current value is missing, set as a numeric 0
                       if(is.na(x = currentValue))
                       {
                         ValueToReturn <- as.numeric(x = 0)
                       }

                       # if the current value is not a "0", set as a numeric 1
                       if(currentValue!="0")
                       {
                         ValueToReturn <- as.numeric(x = 1)
                       # otherwise, set as a 0
                       } else {
                         ValueToReturn <- as.numeric(x = 0)
                       }

                       return(ValueToReturn)
                   },
                   simplify = TRUE
                   ) %>%
                  as.numeric(x = .)

                 return(ArrayToReturn)
             }
       ) %>%
   # convert NAs in the medication columns to zero
   mutate_at(
     .tbl = .,
     .vars = vars(starts_with(match = "Rx_")),
     .funs = function(currentColumn)
             {
               ArrayToReturn <-
                 currentColumn %>%
                   sapply(
                     X = .,
                     FUN = function(currentValue)
                     {
                         # if the value is missing, convert to zero
                         # otherwise, leave as is
                         if(is.na(x = currentValue))
                         {
                           ValueToReturn <- as.numeric(x = 0)
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
   # combine DSM4 and DSM5 columns for the same diagnosis
   mutate(
     .data = .,
     AnxietyDisorder = ((DSMdischargeAx1_AnxietyDisorder_NOS+DSM5discharge_GeneralizedAnxietyDisorder)>0) %>%
                       as.numeric(x = .),
     BipolarDisorder = ((DSMdischargeAx1_BipolarDisorder+DSM5discharge_BipolarDisorder)>0) %>%
                       as.numeric(x = .),
     PTSD = ((DSMdischargeAx1_PostraumaticStressDisorder+DSM5discharge_PosttraumaticStressDisorder)>0) %>%
            as.numeric(x = .),
     PsychoticSymptoms = ((DSMdischargeAx1_PsychoticSymptoms+DSM5discharge_Psychosis)>0) %>%
                         as.numeric(x = .)
     ) %>%
   # remove the DSM4 and DSM5 columns that were combined
   select(
     .data = .,
     -DSMdischargeAx1_AnxietyDisorder_NOS,
     -DSM5discharge_GeneralizedAnxietyDisorder,
     -DSMdischargeAx1_BipolarDisorder,
     -DSM5discharge_BipolarDisorder,
     -DSMdischargeAx1_PostraumaticStressDisorder,
     -DSM5discharge_PosttraumaticStressDisorder,
     -DSMdischargeAx1_PsychoticSymptoms,
     -DSM5discharge_Psychosis
     ) %>%
   # convert acute and residential antipsychotic medication columns into one one-hot encoded columns
   mutate(
     .data = .,
     ##### antipsychotic medications
          # aripiprazole
          aripiprazole = ((Rx_acute_aripiprazole+Rx_residential_aripiprazole)>0) %>% as.numeric(x = .),
          # chlorpromazine
          chlorpromazine = ((Rx_acute_chlorpromazine+Rx_acute_chlorpromazineHCl+Rx_residential_chlorpromazine+Rx_residential_chlorpromazineHCl)>0) %>%
                           as.numeric(x = .),
          # haloperidol
          haloperidol = ((Rx_acute_haloperidol+Rx_residential_haloperidol)>0) %>% as.numeric(x = .),
          # olanzapine
          olanzapine = ((Rx_acute_olanzapine+Rx_residential_olanzapine)>0) %>% as.numeric(x = .),
          # paliperidone
          paliperidone = ((Rx_acute_paliperidone+Rx_residential_paliperidone)>0) %>% as.numeric(x = .),
          # quetiapine
          quetiapine = ((Rx_acute_quetiapine+Rx_acute_quetiapineFumarate+Rx_residential_quetiapine+Rx_residential_quetiapineFumarate)>0) %>%
                       as.numeric(x = .),
          # risperidone
          risperidone = ((Rx_acute_risperidone+Rx_residential_risperidone)>0) %>% as.numeric(x = .),
          # ziprasidone
          ziprasidone = ((Rx_acute_ziprasidone+Rx_residential_ziprasidone)>0) %>% as.numeric(x = .)
     ) %>%
   # convert acute and residential CYP2D6 substrate and inhibitor columns into one-hot encoded column vectors
   mutate(
     .data = .,
     ##### CYP2D6 substrates and inhibitors
         # acetaminophen
         acetaminophen = ((Rx_acute_acetaminophen+Rx_acute_acetaminophenButalbitalCaffeine+Rx_residential_acetaminophen+Rx_residential_acetaminophenAspirinCaffeine+Rx_residential_acetaminophenHydrocodone)>0) %>%
                         as.numeric(x = .),
         # amphetamine and lisdexamfetamine (amphetamine pro-drug)
         amphetamine = ((Rx_acute_amphetamine+Rx_residential_amphetamine+Rx_acute_lisdexamfetamine+Rx_acute_lisdexamfetamineDimesylate+Rx_residential_lisdexamfetamineDimesylate)>0) %>%
                       as.numeric(x = .),
         # amitriptyline
         amitriptyline = (Rx_acute_amitriptyline>0) %>% as.numeric(x = .),
         # aripiprazole (see antipsychotic medications above)
         # atomoxetine
         atomoxetine = ((Rx_acute_atomoxetine+Rx_acute_atomoxetineHCl+Rx_residential_atomoxetine+Rx_residential_atomoxetineHCl)>0) %>%
                       as.numeric(x = .),
         # benzocaine
         benzocaine = (Rx_residential_benzocaine>0) %>% as.numeric(x = .),
         # bupropion
         bupropion = ((Rx_acute_bupropion+Rx_acute_bupropionHCl+Rx_residential_bupropion+Rx_residential_bupropionHCl)>0) %>%
                     as.numeric(x = .),
         # buspirone
         buspirone = ((Rx_acute_buspironeHCl+Rx_residential_buspironeHCl)>0) %>% as.numeric(x = .),
         # chlopheniramine
         chlopheniramine = ((Rx_acute_chlorpheniraminePhenylephrine+Rx_residential_chlorpheniraminePhenylephrine+Rx_residential_chlorpheniraminePseudoephedrine)>0) %>%
                          as.numeric(x = .),
         # chlorpromazine (see antipsychotic medications above)
         # cimetidine
         cimetidine = (Rx_residential_cimetidine>0) %>% as.numeric(x = .),
         # citalopram
         citalopram = ((Rx_acute_citalopram+Rx_acute_citalopramHydrobromide+Rx_residential_citalopramHydrobromide)>0) %>%
                      as.numeric(x = .),
         # clonidine
         clonidine = ((Rx_acute_clonidine+Rx_acute_clonidineHCl+Rx_residential_clonidine+Rx_residential_clonidineHCl)>0) %>%
                     as.numeric(x = .),
         # clotrimazole
         betamethazoneClotrimazole = (Rx_residential_betamethazoneClotrimazole>0) %>% as.numeric(x = .),
         # codeine
         codeine = (Rx_acute_codeineGuaifenesin>0) %>% as.numeric(x = .),
         # dextromethorphan
         dextromethorphan = ((Rx_acute_dextromethorphan+Rx_acute_dextromethorphanGuaifenesin+Rx_residential_dextromethorphan+Rx_residential_dextromethorphanGuaifenesin)>0) %>%
                            as.numeric(x = .),
         # diphenhydramine
         diphenhydramine = ((Rx_acute_diphenhydramine+Rx_acute_diphenhydramineHCl+Rx_residential_diphenhydramine+Rx_residential_diphenhydramineHCl)>0) %>%
                           as.numeric(x = .),
         # duloxetine
         duloxetine = ((Rx_acute_duloxetine+Rx_acute_duloxetineHCl+Rx_residential_duloxetine)>0) %>% as.numeric(x = .),
         # escitalopram
         escitalopram = (Rx_acute_escitalopram>0) %>% as.numeric(x = .),
         # fluoxetine
         fluoxetine = ((Rx_acute_fluoxetine+Rx_acute_fluoxetineHCl+Rx_residential_fluoxetine+Rx_residential_fluoxetineHCl)>0) %>%
                      as.numeric(x = .),
         # fluvoxamine
         fluvoxamine = (Rx_acute_fluvoxamine>0) %>% as.numeric(x = .),
         # haloperidol (see antipsychotic medications above)
         # hydroxyzine
         hydroxyzine = ((Rx_acute_hydroxyzine+Rx_acute_hydroxyzinePamoate+Rx_residential_hydroxyzine+Rx_residential_hydroxyzineHCl+Rx_residential_hydroxyzinePamoate)>0) %>%
                       as.numeric(x = .),
         # imipramine
         imipramine = (Rx_acute_imipramine>0) %>% as.numeric(x = .),
         # ketoconazole
         ketoconazole = ((Rx_acute_ketoconazole+Rx_residential_ketoconazole)>0) %>% as.numeric(x = .),
         # lidocaine
         lidocaine = (Rx_residential_lidocainePrilocaine>0) %>% as.numeric(x = .),
         # loperamide
         loperamide = ((Rx_acute_loperamide+Rx_residential_loperamide)>0) %>% as.numeric(x = .),
         # loratadine
         loratadine = ((Rx_acute_loratadine+Rx_residential_loratadine+Rx_residential_loratadinePseudoephedrine)>0) %>%
                      as.numeric(x = .),
         # miconazole
         miconazole = ((Rx_residential_miconazole+Rx_residential_miconazoleNitrate)>0) %>% as.numeric(x = .),
         # mirtazapine
         mirtazapine = (Rx_acute_mirtazapine>0) %>% as.numeric(x = .),
         # nicotine
         nicotine = ((Rx_acute_nicotine+Rx_residential_nicotine)>0) %>% as.numeric(x = .),
         # olanzapine (see antipsychotic medications above)
         # omeprazole
         omeprazole = ((Rx_acute_omeprazole+Rx_residential_omeprazole)>0) %>% as.numeric(x = .),
         # ondansetron
         ondansetron = ((Rx_acute_ondansetronHCl+Rx_residential_ondansetron)>0) %>% as.numeric(x = .),
         # oxybutynin
         oxybutynin = ((Rx_acute_oxybutyninChloride+Rx_residential_oxybutynin)>0) %>% as.numeric(x = .),
         # paliperidone (see antipsychotic medications above)
         # propanolol
         propanolol = (Rx_acute_propranolol>0) %>% as.numeric(x = .),
         # promethazine
         promethazine = (Rx_residential_promethazineHCl>0) %>% as.numeric(x = .),
         # quetiapine (see antipsychotic medications above)
         # ranitidine
         ranitidine = (Rx_acute_ranitidine>0) %>% as.numeric(x = .),
         # risperidone (see antipsychotic medications above)
         # sertraline
         sertraline = ((Rx_acute_sertraline+Rx_acute_sertralineHCl+Rx_residential_sertraline)>0) %>% as.numeric(x = .),
         # St. Johns Wort
         stJohnsWort = (Rx_residential_stJohnsWort>0) %>% as.numeric(x = .),
         # terbinafine
         terbinafine = ((Rx_acute_terbinafineHCl+Rx_residential_terbinafineHCl)>0) %>% as.numeric(x = .),
         # tramadol
         tramadol = (Rx_residential_tramadol>0) %>% as.numeric(x = .),
         # trazodone
         trazodone = ((Rx_acute_trazodone+Rx_acute_trazodoneHCl+Rx_residential_trazodone+Rx_residential_trazodoneHCl)>0) %>%
                     as.numeric(x = .)
     ) %>%
  # convert acute and residential dopaminergic modulators into one-hot encoded column vectors
   mutate(
     .data = .,
     ##### dopaminergic modulators
         # amphetamine and lisdexamfetamine (pro-drug for amphetamine) (see cyp2d6 modulators above)
         # amantadine
         amantadine = ((Rx_acute_amantadine+Rx_residential_amantadine)>0) %>% as.numeric(x = .),
         # aripiprazole (see cyp2d6 modulators above)
         # atomoxetine (see cyp2d6 modulators above)
         # benztropine
         benztropine = ((Rx_acute_benztropine+Rx_acute_benztropineMesylate+Rx_residential_benztropine+Rx_residential_benztropineMesylate)>0) %>%
                       as.numeric(x = .),
         # bupropion (see cyp2d6 modulators above)
         # buspirone (see cyp2d6 modulators above)
         # carbamazepine
         carbamazepine = (Rx_acute_carbamazepine>0) %>% as.numeric(x = .),
         # carbidopia levodopa
         carbidopiaLevodopa = (Rx_acute_carbidopaLevodopa>0) %>% as.numeric(x = .),
         # chlopheniramine phenylephrine (see cyp2d6 modulators above)
         # chlorpheniramine psudophedrine (see cyp2d6 modulators above)
         # chlorpromazine (see cyp2d6 modulators above)
         # duloxetine (see cyp2d6 modulators above)
         # escitalopram (see cyp2d6 modulators above)
         # haloperidol (see cyp2d6 modulators above)
         # imipramine (see cyp2d6 modulators above)
         # lamotrigine
         lamotrigine = ((Rx_acute_lamotrigine+Rx_residential_lamotrigine)>0) %>% as.numeric(x = .),
         # methylphenidate/dexmethylphenidate
         methylphenidate = ((Rx_acute_dexmethylphenidateHCl+Rx_residential_dexmethylphenidateHCl+Rx_residential_dexmethylphenidate+Rx_acute_methylphenidate+Rx_residential_methylphenidate+Rx_acute_methylphenidateHCl+Rx_residential_methylphenidateHCl)>0) %>%
                           as.numeric(x = .),
         # nicotine (see cyp2d6 modulators above)
         # olanzapine (see cyp2d6 modulators above)
         # paliperidone (see cyp2d6 modulators above)
         # promethazine (see cyp2d6 modulators above)
         # pseudoephedrine
         pseudoephedrine = ((Rx_residential_chlorpheniraminePseudoephedrine+Rx_residential_loratadinePseudoephedrine)>0) %>%
                           as.numeric(x = .)
         # quetiapine (see cyp2d6 modulators above)
         # risperidone (see cyp2d6 modulators above)
         # sertraline (see cyp2d6 modulators above)
         # ziprasidone (see antipsychotic medications above)
     ) %>%
    # convert acute and residential serotonergic modulators into one-hot encoded column vectors
   mutate(
     .data = .,
     ##### serotonergic modulators
         # amphetamine/lisdexamfetamine (see cyp2d6 modulators above)
         # amantadine (see dopaminergic modulators above)
         # amitriptyline (see cyp2d6 modulators above)
         # aripiprazole (see cyp2d6 modulators above)
         # atomoxetine (see cyp2d6 modulators above)
         # benztropine (see dopaminergic modulators above)
         # bupropion (see dopaminergic modulators above)
         # buspirone (see dopaminergic modulators above)
         # carbidopia levodopa (see dopaminergic modulators above)
         # chlopheniramine phenylephrine (see dopaminergic modulators above)
         # chlorpheniramine pseudophedrine (see dopaminergic modulators above)
         # chlorpromazine (see dopaminergic modulators above)
         # citalopram (see cyp2d6 modulators above)
         # cyproheptadine
         cyproheptadine = (Rx_acute_cyproheptadine>0) %>% as.numeric(x = .),
         # dextromethorphan (see cyp2d6 modulators above)
         # duloxetine (see dopaminergic modulators above)
         # escitalopram (see dopaminergic modulators above)
         # fluoxetine (see cyp2d6 modulators above)
         # fluvoxamine (see cyp2d6 modulators above)
         # guaifenesin Dextromethorphan (dextromethorphan -- see cyp2d6 modulators above)
         # haloperidol (see dopaminergic modulators above)
         # imipramine (see dopaminergic modulators above)
         # lamotrigine (see dopaminergic modulators above)
         # methylphenidate/dexmethylphenidate (see dopaminergic modulators above)
         # mirtazapine (see cyp2d6 modulators above)
         # olanzapine (see cyp2d6 modulators above)
         # ondansetron (see cyp2d6 modulators above)
         # oxymetazoline
         oxymetazoline = (Rx_residential_oxymetazoline>0) %>% as.numeric(x = .),
         # paliperidone (see dopaminergic modulators above)
         # propanolol (see cyp2d6 modulators above)
         # pseudoephedrine (see dopaminergic modulators above)
         # quetiapine (see dopaminergic modulators above)
         # risperidone (see dopaminergic modulators above)
         # rizatriptan
         rizatriptan = (Rx_residential_rizatriptan>0) %>% as.numeric(x = .),
         # rocuronium bromide
         rocuroniumBromide = (Rx_residential_rocuroniumBromide>0) %>% as.numeric(x = .),
         # sertraline (see dopaminergic modulators above)
         # stJohnsWort (see cyp2d6 modulators above)
         # sumatriptan
         sumatriptan = (Rx_residential_sumatriptan>0) %>% as.numeric(x = .)
         # tramadol (see cyp2d6 modulators above)
         # trazodone (see cyp2d6 modulators above)
         # ziprasidone (see dopaminergic modulators above)
     ) %>%
   # remove the medication columns that were combined and one-hot encoded
   select(
     .data = .,
     -Rx_acute_aripiprazole,
     -Rx_residential_aripiprazole,
     -Rx_acute_chlorpromazine,
     -Rx_acute_chlorpromazineHCl,
     -Rx_residential_chlorpromazine,
     -Rx_residential_chlorpromazineHCl,
     -Rx_acute_haloperidol,
     -Rx_residential_haloperidol,
     -Rx_acute_olanzapine,
     -Rx_residential_olanzapine,
     -Rx_acute_paliperidone,
     -Rx_residential_paliperidone,
     -Rx_acute_quetiapine,
     -Rx_acute_quetiapineFumarate,
     -Rx_residential_quetiapine,
     -Rx_residential_quetiapineFumarate,
     -Rx_acute_risperidone,
     -Rx_residential_risperidone,
     -Rx_acute_ziprasidone,
     -Rx_residential_ziprasidone,
     -Rx_acute_acetaminophen,
     -Rx_acute_acetaminophenButalbitalCaffeine,
     -Rx_residential_acetaminophen,
     -Rx_residential_acetaminophenAspirinCaffeine,
     -Rx_residential_acetaminophenHydrocodone,
     -Rx_acute_amphetamine,
     -Rx_residential_amphetamine,
     -Rx_acute_lisdexamfetamine,
     -Rx_acute_lisdexamfetamineDimesylate,
     -Rx_residential_lisdexamfetamineDimesylate,
     -Rx_acute_amitriptyline,
     -Rx_acute_atomoxetine,
     -Rx_acute_atomoxetineHCl,
     -Rx_residential_atomoxetine,
     -Rx_residential_atomoxetineHCl,
     -Rx_residential_benzocaine,
     -Rx_acute_bupropion,
     -Rx_acute_bupropionHCl,
     -Rx_residential_bupropion,
     -Rx_residential_bupropionHCl,
     -Rx_acute_buspironeHCl,
     -Rx_residential_buspironeHCl,
     -Rx_acute_chlorpheniraminePhenylephrine,
     -Rx_residential_chlorpheniraminePhenylephrine,
     -Rx_residential_chlorpheniraminePseudoephedrine,
     -Rx_residential_cimetidine,
     -Rx_acute_citalopram,
     -Rx_acute_citalopramHydrobromide,
     -Rx_residential_citalopramHydrobromide,
     -Rx_acute_clonidine,
     -Rx_acute_clonidineHCl,
     -Rx_residential_clonidine,
     -Rx_residential_clonidineHCl,
     -Rx_residential_betamethazoneClotrimazole,
     -Rx_acute_codeineGuaifenesin,
     -Rx_acute_dextromethorphan,
     -Rx_acute_dextromethorphanGuaifenesin,
     -Rx_residential_dextromethorphan,
     -Rx_residential_dextromethorphanGuaifenesin,
     -Rx_acute_diphenhydramine,
     -Rx_acute_diphenhydramineHCl,
     -Rx_residential_diphenhydramine,
     -Rx_residential_diphenhydramineHCl,
     -Rx_acute_duloxetine,
     -Rx_acute_duloxetineHCl,
     -Rx_residential_duloxetine,
     -Rx_acute_escitalopram,
     -Rx_acute_fluoxetine,
     -Rx_acute_fluoxetineHCl,
     -Rx_residential_fluoxetine,
     -Rx_residential_fluoxetineHCl,
     -Rx_acute_fluvoxamine,
     -Rx_acute_hydroxyzine,
     -Rx_acute_hydroxyzinePamoate,
     -Rx_residential_hydroxyzine,
     -Rx_residential_hydroxyzineHCl,
     -Rx_residential_hydroxyzinePamoate,
     -Rx_acute_imipramine,
     -Rx_acute_ketoconazole,
     -Rx_residential_ketoconazole,
     -Rx_residential_lidocainePrilocaine,
     -Rx_acute_loperamide,
     -Rx_residential_loperamide,
     -Rx_acute_loratadine,
     -Rx_residential_loratadine,
     -Rx_residential_loratadinePseudoephedrine,
     -Rx_residential_miconazole,
     -Rx_residential_miconazoleNitrate,
     -Rx_acute_mirtazapine,
     -Rx_acute_nicotine,
     -Rx_residential_nicotine,
     -Rx_acute_omeprazole,
     -Rx_residential_omeprazole,
     -Rx_acute_ondansetronHCl,
     -Rx_residential_ondansetron,
     -Rx_acute_oxybutyninChloride,
     -Rx_residential_oxybutynin,
     -Rx_acute_propranolol,
     -Rx_residential_promethazineHCl,
     -Rx_acute_ranitidine,
     -Rx_acute_sertraline,
     -Rx_acute_sertralineHCl,
     -Rx_residential_sertraline,
     -Rx_residential_stJohnsWort,
     -Rx_acute_terbinafineHCl,
     -Rx_residential_terbinafineHCl,
     -Rx_residential_tramadol,
     -Rx_acute_trazodone,
     -Rx_acute_trazodoneHCl,
     -Rx_residential_trazodone,
     -Rx_residential_trazodoneHCl,
     -Rx_acute_amantadine,
     -Rx_residential_amantadine,
     -Rx_acute_benztropine,
     -Rx_acute_benztropineMesylate,
     -Rx_residential_benztropine,
     -Rx_residential_benztropineMesylate,
     -Rx_acute_carbamazepine,
     -Rx_acute_carbidopaLevodopa,
     -Rx_acute_lamotrigine,
     -Rx_residential_lamotrigine,
     -Rx_acute_dexmethylphenidateHCl,
     -Rx_residential_dexmethylphenidateHCl,
     -Rx_residential_dexmethylphenidate,
     -Rx_acute_methylphenidate,
     -Rx_residential_methylphenidate,
     -Rx_acute_methylphenidateHCl,
     -Rx_residential_methylphenidateHCl,
     -Rx_acute_cyproheptadine,
     -Rx_residential_oxymetazoline,
     -Rx_residential_rizatriptan,
     -Rx_residential_rocuroniumBromide,
     -Rx_residential_sumatriptan
     ) %>%
   # create one-hot encoded cyp2d6, dopaminergic, and serotonergic modulator columns
   # without acetaminophen or diphenhydramine included since they were used only as needed or without risperidone included since everyone received risperidone
   mutate(
     .data = .,
     cyp2d6_modulator = ((amphetamine+amitriptyline+aripiprazole+atomoxetine+benzocaine+bupropion+buspirone+chlopheniramine+chlorpromazine+cimetidine+citalopram+clonidine+betamethazoneClotrimazole+codeine+dextromethorphan+duloxetine+escitalopram+fluoxetine+fluvoxamine+haloperidol+hydroxyzine+imipramine+ketoconazole+lidocaine+loperamide+loratadine+miconazole+mirtazapine+nicotine+olanzapine+omeprazole+ondansetron+oxybutynin+paliperidone+propanolol+promethazine+quetiapine+ranitidine+sertraline+stJohnsWort+terbinafine+tramadol+trazodone)>0) %>%
                        as.numeric(x = .),
     dopaminergic_modulator = ((amphetamine+amantadine+aripiprazole+atomoxetine+benztropine+bupropion+buspirone+carbamazepine+carbidopiaLevodopa+chlopheniramine+chlorpromazine+duloxetine+escitalopram+haloperidol+imipramine+lamotrigine+methylphenidate+nicotine+olanzapine+paliperidone+promethazine+pseudoephedrine+quetiapine+sertraline+ziprasidone)>0) %>%
                              as.numeric(x = .),
     serotonergic_modulator = ((amphetamine+amantadine+amitriptyline+aripiprazole+atomoxetine+benztropine+bupropion+buspirone+carbidopiaLevodopa+chlopheniramine+chlorpromazine+citalopram+cyproheptadine+dextromethorphan+duloxetine+escitalopram+fluoxetine+fluvoxamine+haloperidol+imipramine+lamotrigine+methylphenidate+mirtazapine+olanzapine+ondansetron+oxymetazoline+paliperidone+propanolol+pseudoephedrine+quetiapine+rizatriptan+rocuroniumBromide+sertraline+stJohnsWort+sumatriptan+tramadol+trazodone+ziprasidone)>0) %>%
                              as.numeric(x = .)
     ) %>%
   mutate(
     .data = .,
     # convert weight in kilograms to pounds
     Weight_pounds = Weight_kilograms*2.20462,
     # convert height in centimeters to inches
     Height_inches = Height_centimeters/2.54
     ) %>%
   # calculate body mass index = ( weight (lb) / [height (in)]^2 ) x 703
   mutate(
     .data = .,
     BMI_kgPerSqM = (Weight_pounds/(Height_inches*Height_inches))*703
     ) %>%
   select(
     .data = .,
     -Weight_kilograms,
     -Height_centimeters
     ) %>%
   # one-hot encode the Insurance_medicaid and Insurance_private columns then
   # combine the Insurance_medicaid and Insurance_private columns into a single insurance column
   # If the patient has both private and medicaid insurance, set the insurance to private
  mutate(
    .data = .,
    Insurance_medicaid = if_else(
                                 condition = Insurance_medicaid=="0",
                                 true = as.numeric(x = 0),
                                 false = as.numeric(x = 1)
                                   ),
    Insurance_private = if_else(
                                condition = Insurance_private=="0",
                                true = as.numeric(x = 0),
                                false = as.numeric(x = 1)
                              )
    ) %>%
  mutate(
    .data = .,
    Insurance_medicaidAndPrivate = if_else(
                                           condition = Insurance_medicaid==as.numeric(x = 1) & Insurance_private==as.numeric(x = 1),
                                           true = as.numeric(x = 1),
                                           false = as.numeric(x = 0)
                                             )
    ) %>%
  mutate(
    .data = .,
    Insurance = mapply(
                      FUN = function(medicaidInsuranceColumn,privateInsuranceColumn)
                      {
                            # if the patient has both medicaid and private insurance,
                            # set the insurance as private,
                            # otherwise, set as medicaid or private insurance
                            if((medicaidInsuranceColumn+privateInsuranceColumn)==as.numeric(x = 2))
                            {
                                  ValueToReturn <- "private"
                            } else {

                                if(medicaidInsuranceColumn==as.numeric(x = 1) & privateInsuranceColumn==as.numeric(x = 0))
                                {
                                  ValueToReturn <- "medicaid"
                                } else if(medicaidInsuranceColumn==as.numeric(x = 0) & privateInsuranceColumn==as.numeric(x = 1))
                                {
                                  ValueToReturn <- "private"
                                } else {
                                  ValueToReturn <- as.character(x = NA)
                                }

                            }

                          return(ValueToReturn)
                      },
                      Insurance_medicaid,
                      Insurance_private,
                      SIMPLIFY = FALSE
                        ) %>%
                      unlist(x = .)
    ) %>%
   # Create a PsychiatricHistory column by combining inpatient and outpatient
   # psychiatric history columns. If the patient has either type of psychiatric history, code the value as "y" for yes.
   mutate(
     .data = .,
     PsychiatricHistory = if_else(
                                  condition = (PsychiatricHistory_inpatient=="y") |
                                              (PsychiatricHistory_outpatient=="y"),
                                  true = "y",
                                  false = "n"
                                    )
     ) %>%
    # one-hot encode the reason for risperidone discotinuation columns,
    # with missing values left as missing
    mutate_at(
      .tbl = .,
      .vars = vars(starts_with(match = "ReasonForRispDiscontinuation_")),
      .funs = function(currentColumn)
              {
                 ArrayToReturn <-
                   currentColumn %>%
                   sapply(
                     X = .,
                     FUN = function(currentElement)
                     {
                          # if the current value is missing, leave it as missing
                         if(is.na(x = currentElement))
                         {
                           ValueToReturn <-
                             as.numeric(x = NA)
                         }
                         # if the current element is not a character zero,
                         # meaning there was a reason for discontinuation,
                         # set as a numeric 1
                         else if(currentElement!="0")
                         {
                           ValueToReturn <-
                             as.numeric(x = 1)
                         # otherwise, set as a numeric 0
                         } else {
                           ValueToReturn <-
                             as.numeric(x = 0)
                         }

                         return(ValueToReturn)
                     },
                     simplify = TRUE
                       )

                 return(ArrayToReturn)
              }
        ) %>%
    # one-hot encode the reason for readmittance columns,
    # with missing values left as missing
    mutate_at(
      .tbl = .,
      .vars = vars(starts_with(match = "ReasonForReadmittance_")),
      .funs = function(currentColumn)
              {
                 ArrayToReturn <-
                   currentColumn %>%
                   sapply(
                     X = .,
                     FUN = function(currentElement)
                     {
                          # if the current value is missing, leave it as missing
                         if(is.na(x = currentElement))
                         {
                           ValueToReturn <-
                             as.numeric(x = NA)
                         }
                         # if the current element is not a character zero,
                         # meaning there was a reason for readmittance,
                         # set as a numeric 1
                         else if(currentElement!="0")
                         {
                           ValueToReturn <-
                             as.numeric(x = 1)
                         # otherwise, set as a numeric 0
                         } else {
                           ValueToReturn <-
                             as.numeric(x = 0)
                         }

                         return(ValueToReturn)
                     },
                     simplify = TRUE
                       )

                 return(ArrayToReturn)
              }
        ) %>%
    # Combine reason for risperidone discontinuation columns into similar categories
    # (lack of efficacy, side effects, or other reason)
    mutate(
      .data = .,
      # create a general side effect category for risperidone discontinuation
      ReasonForRispDiscontinuation_GeneralSideEffectCategory =
        as.numeric(x =
                      (
                      as.numeric(x = ReasonForRispDiscontinuation_aggression==as.numeric(x = 1))+
                      as.numeric(x = ReasonForRispDiscontinuation_AttentionalIssues==as.numeric(x = 1))+
                      as.numeric(x = ReasonForRispDiscontinuation_drowsiness==as.numeric(x = 1))+
                      as.numeric(x = ReasonForRispDiscontinuation_ExtrapyramidalSymptoms==as.numeric(x = 1))+
                      as.numeric(x = ReasonForRispDiscontinuation_Insomnia==as.numeric(x = 1))+
                      as.numeric(x = ReasonForRispDiscontinuation_irritability==as.numeric(x = 1))+
                      as.numeric(x = ReasonForRispDiscontinuation_LactationRelated==as.numeric(x = 1))+
                      as.numeric(x = ReasonForRispDiscontinuation_Nausea==as.numeric(x = 1))+
                      as.numeric(x = ReasonForRispDiscontinuation_SideEffects==as.numeric(x = 1))+
                      as.numeric(x = ReasonForRispDiscontinuation_SuicidalIdeation==as.numeric(x = 1))+
                      as.numeric(x = ReasonForRispDiscontinuation_WeightGain==as.numeric(x = 1))
                      )>0
                   ),
      # create a general lack of efficacy category for risperidone discontinuation
      ReasonForRispDiscontinuation_GeneralLackOfEfficacy =
        as.numeric(x =
                     (
                     as.numeric(x = ReasonForRispDiscontinuation_Efficacy==as.numeric(x = 1))+
                     as.numeric(x = ReasonForRispDiscontinuation_NotHelpful==as.numeric(x = 1))
                     )>0
                  ),
      # create a general other reason category for risperidone discontinuation
      ReasonForRispDiscontinuation_GeneralOtherReason =
        as.numeric(x =
                      (
                        as.numeric(x = ReasonForRispDiscontinuation_DiagnosisChange==as.numeric(x = 1))+
                        as.numeric(x = ReasonForRispDiscontinuation_MedicationWash==as.numeric(x = 1))+
                        as.numeric(x = ReasonForRispDiscontinuation_NotAvailable==as.numeric(x = 1))+
                        as.numeric(x = ReasonForRispDiscontinuation_PatientRefusal==as.numeric(x = 1))+
                        as.numeric(x = ReasonForRispDiscontinuation_PrescriberDecision==as.numeric(x = 1))+
                        as.numeric(x = ReasonForRispDiscontinuation_PrescriberFearOfSideEffects==as.numeric(x = 1))
                      )>0
                  )
      ) %>%
     # create a single risperidone discontinuation column with one of the general three discontinuation categories
     # (side effects, lack of efficacy, other)
    mutate(
      .data = .,
      ReasonForRisperidoneDiscontinuation = if_else(
                                                    condition = !is.na(x = ReasonForRispDiscontinuation_GeneralSideEffectCategory) & ReasonForRispDiscontinuation_GeneralSideEffectCategory==as.numeric(x = 1),
                                                    true = "side_effects",
                                                    false = if_else(
                                                                    condition = !is.na(x = ReasonForRispDiscontinuation_GeneralLackOfEfficacy) & ReasonForRispDiscontinuation_GeneralLackOfEfficacy==as.numeric(x = 1),
                                                                    true = "lack_of_efficacy",
                                                                    false = if_else(
                                                                                    condition = !is.na(x = ReasonForRispDiscontinuation_GeneralOtherReason) & ReasonForRispDiscontinuation_GeneralOtherReason==as.numeric(x = 1),
                                                                                    true = "other_reason",
                                                                                    false = as.character(x = NA)
                                                                                      )
                                                                      )
                                                      )
      ) %>%
  # make sure that all variables that are included in the analysis are the appropriate data type for performing regressions
  mutate(
    .data = .,
    ChangeInGAF = ChangeInGAF %>% as.numeric(x = .),
    DaysOnRisperidone = DaysOnRisperidone %>% as.numeric(x = .),
    DaysToReadmittance = DaysToReadmittance %>% as.numeric(x = .),
    GAF_admit = GAF_admit %>% as.numeric(x = .),
    GAF_discharge = GAF_discharge %>% as.numeric(x = .),
    LengthOfStay = LengthOfStay %>% as.numeric(x = .),
    MaxDoseRisperidone_mg = MaxDoseRisperidone_mg %>% as.numeric(x = .),
    ReadmittedOnRisperidone = ReadmittedOnRisperidone %>% factor(x = .),
    PatientReadmitted = PatientReadmitted %>% factor(x = .),
    MaxDoseRisperidone_mgPerDose = MaxDoseRisperidone_mgPerDose %>% as.numeric(x = .),
    MaxDoseRisperidone_frequency = MaxDoseRisperidone_frequency %>% factor(x = .),
    MaxDoseRisperidone_mgPerDay = MaxDoseRisperidone_mgPerDay %>% as.numeric(x = .),
    Age_years = Age_years %>% as.numeric(x = .),
    Gender_MorF = Gender_MorF %>% factor(x = .),
    BMI_kgPerSqM = BMI_kgPerSqM %>% as.numeric(x = .),
    Residential_visit = Residential_visit %>% factor(x = .),
    Insurance = Insurance %>% factor(x = .),
    PsychiatricHistory_family = PsychiatricHistory_family %>% factor(x = .),
    NeglectAndAbuseHistory = NeglectAndAbuseHistory %>% factor(x = .),
    SubstanceAbuseHistory = SubstanceAbuseHistory %>% factor(x = .),
    DSMdischargeAx1_ADHD_combinedType = DSMdischargeAx1_ADHD_combinedType %>% factor(x = .),
    AnxietyDisorder = AnxietyDisorder %>% factor(x = .),
    DSMdischargeAx1_AutismSpectrumDisorder = DSMdischargeAx1_AutismSpectrumDisorder %>% factor(x = .),
    BipolarDisorder = BipolarDisorder %>% factor(x = .),
    DSMdischargeAx1_IntermittentExplosiveDisorder = DSMdischargeAx1_IntermittentExplosiveDisorder %>% factor(x = .),
    DSMdischargeAx1_MajorDepressiveDisorder = DSMdischargeAx1_MajorDepressiveDisorder %>% factor(x = .),
    DSMdischargeAx1_MoodDisorder_NOS = DSMdischargeAx1_MoodDisorder_NOS %>% factor(x = .),
    DSMdischargeAx1_OppositionalDefiantDisorder = DSMdischargeAx1_OppositionalDefiantDisorder %>% factor(x = .),
    PTSD = PTSD %>% factor(x = .),
    PsychoticSymptoms = PsychoticSymptoms %>% factor(x = .),
    cyp2d6_modulator = cyp2d6_modulator %>% factor(x = .),
    dopaminergic_modulator = dopaminergic_modulator %>% factor(x = .),
    serotonergic_modulator = serotonergic_modulator %>% factor(x = .),
    Risperidone_discontinued = Risperidone_discontinued %>% factor(x = .)
  ) %>%
  mutate_at(
    .tbl = .,
    .vars = vars(starts_with(match = "ReasonForRispDiscontinuation_")),
    .funs = function(currentColumn)
    {
      ArrayToReturn <-
        currentColumn %>%
        factor(x = .)

      return(ArrayToReturn)
    }
      )

############ create a summary table of baseline clinical and demographic characteristics

BaselineClinicalAndDemographicCharacteristics <-
  GenotypeMatrixAndClinicalData %>%
  summarise(
    .data = .,
    # baseline characteristics
    #"Baseline characteristics" = paste0(""),
    "Age (y)" = paste(signif(x = mean(x = Age_years,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = Age_years,na.rm = TRUE),digits = 3)),
    "Age range (y)" = paste0("(",min(Age_years,na.rm = TRUE)," - ",max(Age_years,na.rm = TRUE),")"),
    #"Weight (lb)" = paste(signif(x = mean(x = Weight_pounds,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = Weight_pounds,na.rm = TRUE),digits = 3)),
    #"Height (in)" = paste(signif(x = mean(x = Height_inches,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = Height_inches,na.rm = TRUE),digits = 3)),
    "BMI (kg/sq. m)" = paste(signif(x = mean(x = BMI_kgPerSqM,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = BMI_kgPerSqM,na.rm = TRUE),digits = 3)),
    "BMI range (kg/ sq. m)" = paste0("(",signif(x = min(BMI_kgPerSqM,na.rm = TRUE),digits = 3)," - ",signif(x = max(BMI_kgPerSqM,na.rm = TRUE),digits = 3),")"),
    "Gender (N)" = paste0(""),
    "Gender (F) (N)" = paste0(sum(as.numeric(x = Gender_MorF=="F"),na.rm = TRUE)," (",round(x = (sum(as.numeric(x = Gender_MorF=="F"),na.rm = TRUE)/nrow(x = .))*100,digits = 1),"%)"),
    "Gender (M) (N)" = paste0(sum(as.numeric(x = Gender_MorF=="M"),na.rm = TRUE)," (",round(x = (sum(as.numeric(x = Gender_MorF=="M"),na.rm = TRUE)/nrow(x = .))*100,digits = 1),"%)"),
    "Insurance (N)" = paste0(""),
    "Insurance (Medicaid) (N)" = paste0(sum(as.numeric(x = Insurance_medicaid!="0" & Insurance_private=="0"),na.rm = TRUE)," (",round(x = (sum(as.numeric(x = Insurance_medicaid!="0" & Insurance_private=="0"),na.rm = TRUE)/nrow(x = .))*100,digits = 1),"%)"),
    "Insurance (private) (N)" = paste0(sum(as.numeric(x = Insurance_private!="0" & Insurance_medicaid=="0"),na.rm = TRUE)," (",round(x = (sum(as.numeric(x = Insurance_private!="0" & Insurance_medicaid=="0"),na.rm = TRUE)/nrow(x = .))*100,digits = 1),"%)"),
    "Insurance (Medicaid and private) (N)" = paste0(sum(as.numeric(x = Insurance_medicaidAndPrivate==1),na.rm = TRUE)," (",round(x = (sum(as.numeric(x = Insurance_medicaidAndPrivate==1),na.rm = TRUE)/nrow(x = .))*100,digits = 1),"%)"),
    "Insurance (none) (N)" = paste0(sum(as.numeric(x = Insurance_medicaid=="0" & Insurance_private=="0"),na.rm = TRUE)," (",round(x = (sum(as.numeric(x = Insurance_medicaid=="0" & Insurance_private=="0"),na.rm = TRUE)/nrow(x = .))*100,digits = 1),"%)"),
    # "Admit unit (acute) (N)" = sum(as.numeric(x = Admit_unit=="a"),na.rm = TRUE),
    # "Admit unit (residential) (N)" = sum(as.numeric(x = Admit_unit=="r"),na.rm = TRUE),
    "Acute visit only (N)" = paste0(sum(as.numeric(x = Residential_visit=="n"),na.rm = TRUE)," (",round(x = (sum(as.numeric(x = Residential_visit=="n"),na.rm = TRUE)/nrow(x = .))*100,digits = 1),"%)"),
    "Acute length of stay (days)" = paste(signif(x = mean(x = acute_days,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = acute_days,na.rm = TRUE),digits = 3)),
    "Acute days range" = paste0("(",min(acute_days,na.rm = TRUE)," - ",max(acute_days,na.rm = TRUE),")"),
    "Residential visit only (N)" = paste0(sum(as.numeric(x = Admit_unit=="r" & Discharge_unit=="r" & UnitTransfers==as.numeric(x = 0)),na.rm = TRUE)," (",round(x = (sum(as.numeric(x = Admit_unit=="r" & Discharge_unit=="r" & UnitTransfers==as.numeric(x = 0)),na.rm = TRUE)/nrow(x = .))*100,digits = 1),"%)"),
    "Residential length of stay (days)" = paste(signif(x = mean(x = residential_days,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = residential_days,na.rm = TRUE),digits = 3)),
    "Residential days range" = paste0("(",min(residential_days,na.rm = TRUE)," - ",max(residential_days,na.rm = TRUE),")"),
    "Acute AND residential stay (both) (N)" = paste0(sum(UnitTransfers>0,na.rm = TRUE)," (",round(x = (sum(UnitTransfers>0,na.rm = TRUE)/nrow(x = .))*100,digits = 1),"%)"),
    "GAF score (admission)" = paste(signif(x = mean(x = GAF_admit,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = GAF_admit,na.rm = TRUE),digits = 3)),
    "GAF score (admission) range" = paste0("(",min(GAF_admit,na.rm = TRUE)," - ",max(GAF_admit,na.rm = TRUE),")"),
    # Risperidone dosing information
    "Risperidone dosing" = paste0(""),
    "Dose range (mg/dose)" = paste0("(",min(MaxDoseRisperidone_mgPerDose,na.rm = TRUE)," - ",max(MaxDoseRisperidone_mgPerDose,na.rm = TRUE),")"),
    "Dose (mg/dose)" = paste(signif(x = mean(x = MaxDoseRisperidone_mgPerDose,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = MaxDoseRisperidone_mgPerDose,na.rm = TRUE),digits = 3)),
    "Dose range (mg/day)" = paste0("(",min(MaxDoseRisperidone_mgPerDay,na.rm = TRUE)," - ",max(MaxDoseRisperidone_mgPerDay,na.rm = TRUE),")"),
    "Dose (mg/day)" = paste(signif(x = mean(x = MaxDoseRisperidone_mgPerDay,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = MaxDoseRisperidone_mgPerDay,na.rm = TRUE),digits = 3)),
    "Risperidone dose frequency" = paste0(""), 
    "qd" = paste0(sum(as.numeric(x = MaxDoseRisperidone_frequency=="qd"),na.rm = TRUE)," (",round(x = (sum(as.numeric(x = MaxDoseRisperidone_frequency=="qd"),na.rm = TRUE)/nrow(x = .))*100,digits = 1),"%)"),
    "bid" = paste0(sum(as.numeric(x = MaxDoseRisperidone_frequency=="bid"),na.rm = TRUE)," (",round(x = (sum(as.numeric(x = MaxDoseRisperidone_frequency=="bid"),na.rm = TRUE)/nrow(x = .))*100,digits = 1),"%)"),
    "tid" = paste0(sum(as.numeric(x = MaxDoseRisperidone_frequency=="tid"),na.rm = TRUE)," (",round(x = (sum(as.numeric(x = MaxDoseRisperidone_frequency=="tid"),na.rm = TRUE)/nrow(x = .))*100,digits = 1),"%)"),
    "qid" = paste0(sum(as.numeric(x = MaxDoseRisperidone_frequency=="qid"),na.rm = TRUE)," (",round(x = (sum(as.numeric(x = MaxDoseRisperidone_frequency=="qid"),na.rm = TRUE)/nrow(x = .))*100,digits = 1),"%)"),
    # psychiatric history
    #"Psychiatric history (N)" = paste0(""),
    #"Psych. hist. (inpatient or outpatient)" = sum(as.numeric(x = PsychiatricHistory=="y"),na.rm = TRUE),
    #"Psych. hist. (inpatient)" = sum(as.numeric(x = PsychiatricHistory_inpatient=="y"),na.rm = TRUE),
    #"Psych. hist. (outpatient)" = sum(as.numeric(x = PsychiatricHistory_outpatient=="y"),na.rm = TRUE),
    "Psych. hist. (family) (N)" = paste0(sum(as.numeric(x = PsychiatricHistory_family=="y"),na.rm = TRUE)," (",round(x = (sum(as.numeric(x = PsychiatricHistory_family=="y"),na.rm = TRUE)/nrow(x = .))*100,digits = 1),"%)"),
    # neglect and abuse
    #"Neglect and abuse history (N)" = paste0(""),
    "Hist. (neglect and abuse) (N)" = paste0(sum(as.numeric(x = NeglectAndAbuseHistory=="y"),na.rm = TRUE)," (",round(x = (sum(as.numeric(x = NeglectAndAbuseHistory=="y"),na.rm = TRUE)/nrow(x = .))*100,digits = 1),"%)"),
    # substance abuse
    # "Substance abuse (N)" = paste0(""),
    # "Alcohol abuse" = sum(DSMdischargeAx1_AlcoholAbuse, na.rm = TRUE),
    # "Marijuana dependence" = sum(DSMdischargeAx1_MarijuanaDependence, na.rm = TRUE),
    # "Narcotics dependence" = sum(DSMdischargeAx1_NarcoticsDependence_HO, na.rm = TRUE),
    # "Opiate use disorder" = sum(DSMdischargeAx1_OpiateUseDisorder, na.rm = TRUE),
    # "Polysubstance abuse" = sum(DSMdischargeAx1_PolysubstanceAbuse, na.rm = TRUE),
    "Substance abuse history (N)" = paste0(sum(as.numeric(x = SubstanceAbuseHistory=="y"),na.rm = TRUE)," (",round(x = (sum(as.numeric(x = SubstanceAbuseHistory=="y"),na.rm = TRUE)/nrow(x = .))*100,digits = 1),"%)"),
    # antipsychotic medications
    # "Antipsychotic medications (N)" = paste0(""),
    # ##### antipsychotic medications
    #  # aripiprazole
    #  aripiprazole = aripiprazole %>% sum(.,na.rm = TRUE),
    #  # chlorpromazine
    #  chlorpromazine = chlorpromazine %>% sum(.,na.rm = TRUE),
    #  # haloperidol
    #  haloperidol = haloperidol %>% sum(.,na.rm = TRUE),
    #  # olanzapine
    #  olanzapine = olanzapine %>% sum(.,na.rm = TRUE),
    #  # paliperidone
    #  paliperidone = paliperidone %>% sum(.,na.rm = TRUE),
    #  # quetiapine
    #  quetiapine = quetiapine %>% sum(.,na.rm = TRUE),
    #  # risperidone
    #  risperidone = risperidone %>% sum(.,na.rm = TRUE),
    #  # ziprasidone
    #  ziprasidone = ziprasidone %>% sum(.,na.rm = TRUE),
    # cyp2d6 substrates and inhibitors
    #"CYP2D6 substrates and inhibitors (N)" = paste0(""),
    "CYP2D6 modulators (N) total w/o acetaminophen, diphenhydramine, or risperidone" = cyp2d6_modulator %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(., na.rm = TRUE),
    "CYP2D6 modulators (%)" = round(x = (sum(as.numeric(x = as.character(x = cyp2d6_modulator)),na.rm = TRUE)/161)*100,digits = 1),
    #   # acetaminophen
    #   acetaminophen = acetaminophen %>% sum(.,na.rm = TRUE),
    #   # amphetamine and lisdexamfetamine (pro-drug for amphetamine)
    #   amphetamine = amphetamine %>% sum(.,na.rm = TRUE),
    #   # amitriptyline
    #   amitriptyline = amitriptyline %>% sum(.,na.rm = TRUE),
    #   # aripiprazole (see antipsychotic medications above)
    #   aripiprazole_cyp2d6 = aripiprazole %>% sum(.,na.rm = TRUE),
    #   # atomoxetine
    #   atomoxetine = atomoxetine %>% sum(.,na.rm = TRUE),
    #   # benzocaine
    #   benzocaine = benzocaine %>% sum(.,na.rm = TRUE),
    #   # bupropion
    #   bupropion = bupropion %>% sum(.,na.rm = TRUE),
    #   # buspirone
    #   buspirone = buspirone %>% sum(.,na.rm = TRUE),
    #   # chlopheniramine
    #   chlopheniramine = chlopheniramine %>% sum(.,na.rm = TRUE),
    #   # chlorpromazine (see antipsychotic medications above)
    #   chlorpromazine_cyp2d6 = chlorpromazine %>% sum(.,na.rm = TRUE),
    #   # cimetidine
    #   cimetidine = cimetidine %>% sum(.,na.rm = TRUE),
    #   # citalopram
    #   citalopram = citalopram %>% sum(.,na.rm = TRUE),
    #   # clonidine
    #   clonidine = clonidine %>% sum(.,na.rm = TRUE),
    #   # clotrimazole
    #   betamethazoneClotrimazole = betamethazoneClotrimazole %>% sum(.,na.rm = TRUE),
    #   # codeine
    #   codeine = codeine %>% sum(.,na.rm = TRUE),
    #   # dextromethorphan
    #   dextromethorphan = dextromethorphan %>% sum(.,na.rm = TRUE),
    #   # diphenhydramine
    #   diphenhydramine = diphenhydramine %>% sum(.,na.rm = TRUE),
    #   # duloxetine
    #   duloxetine = duloxetine %>% sum(.,na.rm = TRUE),
    #   # escitalopram
    #   escitalopram = escitalopram %>% sum(.,na.rm = TRUE),
    #   # fluoxetine
    #   fluoxetine = fluoxetine %>% sum(.,na.rm = TRUE),
    #   # fluvoxamine
    #   fluvoxamine = fluvoxamine %>% sum(.,na.rm = TRUE),
    #   # haloperidol (see antipsychotic medications above)
    #   haloperidol_cyp2d6 = haloperidol %>% sum(.,na.rm = TRUE),
    #   # hydroxyzine
    #   hydroxyzine = hydroxyzine %>% sum(.,na.rm = TRUE),
    #   # imipramine
    #   imipramine = imipramine %>% sum(.,na.rm = TRUE),
    #   # ketoconazole
    #   ketoconazole = ketoconazole %>% sum(.,na.rm = TRUE),
    #   # lidocaine
    #   lidocaine = lidocaine %>% sum(.,na.rm = TRUE),
    #   # loperamide
    #   loperamide = loperamide %>% sum(.,na.rm = TRUE),
    #   # loratadine
    #   loratadine = loratadine %>% sum(.,na.rm = TRUE),
    #   # miconazole
    #   miconazole = miconazole %>% sum(.,na.rm = TRUE),
    #   # mirtazapine
    #   mirtazapine = mirtazapine %>% sum(.,na.rm = TRUE),
    #   # nicotine
    #   nicotine = nicotine %>% sum(.,na.rm = TRUE),
    #   # olanzapine (see antipsychotic medications above)
    #   olanzapine_cyp2d6 = olanzapine %>% sum(.,na.rm = TRUE),
    #   # omeprazole
    #   omeprazole = omeprazole %>% sum(.,na.rm = TRUE),
    #   # ondansetron
    #   ondansetron = ondansetron %>% sum(.,na.rm = TRUE),
    #   # oxybutynin
    #   oxybutynin = oxybutynin %>% sum(.,na.rm = TRUE),
    #   # paliperidone (see antipsychotic medications above)
    #   paliperidone_cyp2d6 = paliperidone %>% sum(.,na.rm = TRUE),
    #   # propanolol
    #   propanolol = propanolol %>% sum(.,na.rm = TRUE),
    #   # promethazine
    #   promethazine = promethazine %>% sum(.,na.rm = TRUE),
    #   # quetiapine (see antipsychotic medications above)
    #   quetiapine_cyp2d6 = quetiapine %>% sum(.,na.rm = TRUE),
    #   # ranitidine
    #   ranitidine = ranitidine %>% sum(.,na.rm = TRUE),
    #   # risperidone (see antipsychotic medications above)
    #   risperidone_cyp2d6 = risperidone %>% sum(.,na.rm = TRUE),
    #   # sertraline
    #   sertraline = sertraline %>% sum(.,na.rm = TRUE),
    #   # St. Johns Wort
    #   stJohnsWort = stJohnsWort %>% sum(.,na.rm = TRUE),
    #   # terbinafine
    #   terbinafine = terbinafine %>% sum(.,na.rm = TRUE),
    #   # tramadol
    #   tramadol = tramadol %>% sum(.,na.rm = TRUE),
    #   # trazodone
    #   trazodone = trazodone %>% sum(.,na.rm = TRUE),
    # # dopaminergic modulators
    # "Dopaminergic modulators (N)" = paste0(""),
    # "total w/o risperidone" = dopaminergic_modulator %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(., na.rm = TRUE),
    #   # amphetamine (see cyp2d6 modulators above)
    #   amphetamine_dopaminergic = amphetamine %>% sum(.,na.rm = TRUE),
    #   # amantadine
    #   amantadine = amantadine %>% sum(.,na.rm = TRUE),
    #   # aripiprazole (see cyp2d6 modulators above)
    #   aripiprazole_dopaminergic = aripiprazole %>% sum(.,na.rm = TRUE),
    #   # atomoxetine (see cyp2d6 modulators above)
    #   atomoxetine_dopaminergic = atomoxetine %>% sum(.,na.rm = TRUE),
    #   # benztropine
    #   benztropine = benztropine %>% sum(.,na.rm = TRUE),
    #   # bupropion (see cyp2d6 modulators above)
    #   bupropion_dopaminergic = bupropion %>% sum(.,na.rm = TRUE),
    #   # buspirone (see cyp2d6 modulators above)
    #   buspirone_dopaminergic = buspirone %>% sum(.,na.rm = TRUE),
    #   # carbamazepine
    #   carbamazepine = carbamazepine %>% sum(.,na.rm = TRUE),
    #   # carbidopia levodopa
    #   carbidopiaLevodopa = carbidopiaLevodopa %>% sum(.,na.rm = TRUE),
    #   # chlopheniramine (see cyp2d6 modulators above)
    #   chlopheniramine_dopaminergic = chlopheniramine %>% sum(.,na.rm = TRUE),
    #   # chlorpromazine (see cyp2d6 modulators above)
    #   chlorpromazine_dopaminergic = chlorpromazine %>% sum(.,na.rm = TRUE),
    #   # duloxetine (see cyp2d6 modulators above)
    #   duloxetine_dopaminergic = duloxetine %>% sum(.,na.rm = TRUE),
    #   # escitalopram (see cyp2d6 modulators above)
    #   escitalopram_dopaminergic = escitalopram %>% sum(.,na.rm = TRUE),
    #   # haloperidol (see cyp2d6 modulators above)
    #   haloperidol_dopaminergic = haloperidol %>% sum(.,na.rm = TRUE),
    #   # imipramine (see cyp2d6 modulators above)
    #   imipramine_dopaminergic = imipramine %>% sum(.,na.rm = TRUE),
    #   # lamotrigine
    #   lamotrigine = lamotrigine %>% sum(.,na.rm = TRUE),
    #   # methylphenidate/dexmethylphenidate
    #   methylphenidate = methylphenidate %>% sum(.,na.rm = TRUE),
    #   # nicotine (see cyp2d6 modulators above)
    #   nicotine_dopaminergic = nicotine %>% sum(.,na.rm = TRUE),
    #   # olanzapine (see cyp2d6 modulators above)
    #   olanzapine_dopaminergic = olanzapine %>% sum(.,na.rm = TRUE),
    #   # paliperidone (see cyp2d6 modulators above)
    #   paliperidone_dopaminergic = paliperidone %>% sum(.,na.rm = TRUE),
    #   # promethazine (see cyp2d6 modulators above)
    #   promethazine_dopaminergic = promethazine %>% sum(.,na.rm = TRUE),
    #   # pseudoephedrine
    #   pseudoephedrine = pseudoephedrine %>% sum(.,na.rm = TRUE),
    #   # quetiapine (see cyp2d6 modulators above)
    #   quetiapine_dopaminergic = quetiapine %>% sum(.,na.rm = TRUE),
    #   # risperidone (see cyp2d6 modulators above)
    #   risperidone_dopaminergic = risperidone %>% sum(.,na.rm = TRUE),
    #   # sertraline (see cyp2d6 modulators above)
    #   sertraline_dopaminergic = sertraline %>% sum(.,na.rm = TRUE),
    #   # ziprasidone (see antipsychotic medications above)
    #   ziprasidone_dopaminergic = ziprasidone %>% sum(.,na.rm = TRUE),
    # # serotonergic modulators
    # "Serotonergic modulators (N)" = paste0(""),
    #   "total without risperidone" = serotonergic_modulator %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(., na.rm = TRUE),
    #   # amphetamine and lisdexamfetamine (amphetamine pro-drug) (see cyp2d6 modulators above)
    #   amphetamine_serotonergic = amphetamine %>% sum(.,na.rm = TRUE),
    #   # amantadine (see dopaminergic modulators above)
    #   amantadine_serotonergic = amantadine %>% sum(.,na.rm = TRUE),
    #   # amitriptyline (see cyp2d6 modulators above)
    #   amitriptyline_serotonergic = amitriptyline %>% sum(.,na.rm = TRUE),
    #   # aripiprazole (see cyp2d6 modulators above)
    #   aripiprazole_serotonergic = aripiprazole %>% sum(.,na.rm = TRUE),
    #   # atomoxetine (see cyp2d6 modulators above)
    #   atomoxetine_serotonergic = atomoxetine %>% sum(.,na.rm = TRUE),
    #   # benztropine (see dopaminergic modulators above)
    #   benztropine_serotonergic = benztropine %>% sum(.,na.rm = TRUE),
    #   # bupropion (see dopaminergic modulators above)
    #   bupropion_serotonergic = bupropion %>% sum(.,na.rm = TRUE),
    #   # buspirone (see dopaminergic modulators above)
    #   buspirone_serotonergic = buspirone %>% sum(.,na.rm = TRUE),
    #   # carbidopia levodopa (see dopaminergic modulators above)
    #   carbidopiaLevodopa_serotonergic = carbidopiaLevodopa %>% sum(.,na.rm = TRUE),
    #   # chlopheniramine (see dopaminergic modulators above)
    #   chlopheniramine_serotonergic = chlopheniramine %>% sum(.,na.rm = TRUE),
    #   # chlorpromazine (see dopaminergic modulators above)
    #   chlorpromazine_serotonergic = chlorpromazine %>% sum(.,na.rm = TRUE),
    #   # citalopram (see cyp2d6 modulators above)
    #   citalopram_serotonergic = citalopram %>% sum(.,na.rm = TRUE),
    #   # cyproheptadine
    #   cyproheptadine = cyproheptadine %>% sum(.,na.rm = TRUE),
    #   # dextromethorphan (see cyp2d6 modulators above)
    #   dextromethorphan_serotonergic = dextromethorphan %>% sum(.,na.rm = TRUE),
    #   # duloxetine (see dopaminergic modulators above)
    #   duloxetine_serotonergic = duloxetine %>% sum(.,na.rm = TRUE),
    #   # escitalopram (see dopaminergic modulators above)
    #   escitalopram_serotonergic = escitalopram %>% sum(.,na.rm = TRUE),
    #   # fluoxetine (see cyp2d6 modulators above)
    #   fluoxetine_serotonergic = fluoxetine %>% sum(.,na.rm = TRUE),
    #   # fluvoxamine (see cyp2d6 modulators above)
    #   fluvoxamine_serotonergic = fluvoxamine %>% sum(.,na.rm = TRUE),
    #   # haloperidol (see dopaminergic modulators above)
    #   haloperidol_serotonergic = haloperidol %>% sum(.,na.rm = TRUE),
    #   # imipramine (see dopaminergic modulators above)
    #   imipramine_serotonergic = imipramine %>% sum(.,na.rm = TRUE),
    #   # lamotrigine (see dopaminergic modulators above)
    #   lamotrigine_serotonergic = lamotrigine %>% sum(.,na.rm = TRUE),
    #   # methylphenidate/dexmethylphenidate (see dopaminergic modulators above)
    #   methylphenidate_serotonergic = methylphenidate %>% sum(.,na.rm = TRUE),
    #   # mirtazapine (see cyp2d6 modulators above)
    #   mirtazapine_serotonergic = mirtazapine %>% sum(.,na.rm = TRUE),
    #   # olanzapine (see cyp2d6 modulators above)
    #   olanzapine_serotonergic = olanzapine %>% sum(.,na.rm = TRUE),
    #   # ondansetron (see cyp2d6 modulators above)
    #   ondansetron_serotonergic = ondansetron %>% sum(.,na.rm = TRUE),
    #   # oxymetazoline
    #   oxymetazoline = oxymetazoline %>% sum(.,na.rm = TRUE),
    #   # paliperidone (see dopaminergic modulators above)
    #   paliperidone_serotonergic = paliperidone %>% sum(.,na.rm = TRUE),
    #   # propanolol (see cyp2d6 modulators above)
    #   propanolol_serotonergic = propanolol %>% sum(.,na.rm = TRUE),
    #   # pseudoephedrine (see dopaminergic modulators above)
    #   pseudoephedrine_serotonergic = pseudoephedrine %>% sum(.,na.rm = TRUE),
    #   # quetiapine (see dopaminergic modulators above)
    #   quetiapine_serotonergic = quetiapine %>% sum(.,na.rm = TRUE),
    #   # risperidone (see dopaminergic modulators above)
    #   risperidone_serotonergic = risperidone %>% sum(.,na.rm = TRUE),
    #   # rizatriptan
    #   rizatriptan = rizatriptan %>% sum(.,na.rm = TRUE),
    #   # rocuronium bromide
    #   rocuroniumBromide = rocuroniumBromide %>% sum(.,na.rm = TRUE),
    #   # sertraline (see dopaminergic modulators above)
    #   sertraline_serotonergic = sertraline %>% sum(.,na.rm = TRUE),
    #   # stJohnsWort (see cyp2d6 modulators above)
    #   stJohnsWort_serotonergic = stJohnsWort %>% sum(.,na.rm = TRUE),
    #   # sumatriptan
    #   sumatriptan = sumatriptan %>% sum(.,na.rm = TRUE),
    #   # tramadol (see cyp2d6 modulators above)
    #   tramadol_serotonergic = tramadol %>% sum(.,na.rm = TRUE),
    #   # trazodone (see cyp2d6 modulators above)
    #   trazodone_serotonergic = trazodone %>% sum(.,na.rm = TRUE),
    #   # ziprasidone (see dopaminergic modulators above)
    #   ziprasidone_serotonergic = ziprasidone %>% sum(.,na.rm = TRUE)
    # psychiatric diagnosis
    "Psychiatric diagnosis (N)" = paste0("")
    ) %>%
  # transpose the rows and columns
  t(x = .) %>%
  # convert to a dataframe
  data.frame(.) %>%
  # convert the rownames to a column
  rownames_to_column(
    .data = .,
    var = "Characteristic"
    ) %>%
  select(
    .data = .,
    Characteristic,
    "Value" = `.`
    ) %>%
  mutate(
    .data = .,
    Characteristic = Characteristic %>%
                     gsub(pattern = "_cyp2d6",replacement = "",x = .) %>%
                     gsub(pattern = "_dopaminergic",replacement = "",x = .) %>%
                     gsub(pattern = "_serotonergic",replacement = "",x = .)
    )

BaselineClinicalAndDemographicCharacteristics <- 
  GenotypeMatrixAndClinicalData %>%
    summarise(
      .data = .,
      "Attention deficit hyperactivity disorder" = DSMdischargeAx1_ADHD_combinedType %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(., na.rm = TRUE),
      "Anxiety disorder" = AnxietyDisorder %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(., na.rm = TRUE),
      "Autism spectrum disorder" = DSMdischargeAx1_AutismSpectrumDisorder %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(., na.rm = TRUE),
      "Bipolar disorder" = BipolarDisorder %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(., na.rm = TRUE),
      "Intermittent explosive disorder" = DSMdischargeAx1_IntermittentExplosiveDisorder %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(., na.rm = TRUE),
      "Major depressive disorder" = DSMdischargeAx1_MajorDepressiveDisorder %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(., na.rm = TRUE),
      "Mood disorder" = DSMdischargeAx1_MoodDisorder_NOS %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(., na.rm = TRUE),
      "Obsessive compulsive disorder" = DSMdischargeAx1_ObsessiveCompulsiveDisorder %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(., na.rm = TRUE),
      "Oppositional defiant disorder" = DSMdischargeAx1_OppositionalDefiantDisorder %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(., na.rm = TRUE),
      "Posttraumatic stress disorder" = PTSD %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(., na.rm = TRUE),
      "Psychotic disorder" = PsychoticSymptoms %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(., na.rm = TRUE),
     ) %>%
  # transpose the rows and columns
  t(x = .) %>%
  # convert to a dataframe
  data.frame(.) %>%
  # convert the rownames to a column
  rownames_to_column(
    .data = .,
    var = "Characteristic"
  ) %>%
  select(
    .data = .,
    Characteristic,
    "Value" = `.`
  ) %>%
  # arrange the psychiatric diagnosis count values in descending order
  arrange(
    .data = .,
    desc(x = Value)
    ) %>%
  # calculate the number of participants with each psychiatric illness as a percentage of the total participants (N=161)
  mutate(
    .data = .,
    Percentage = signif(x = (Value/161)*100,digits = 2)
  ) %>%
  # paste the value count and the percentage together
  mutate(
    .data = .,
    Value = paste0(Value," (",Percentage,"%)")
    ) %>%
  # deselect the percentage column
  select(
    .data = .,
    -Percentage
    ) %>%
  rbind(
    BaselineClinicalAndDemographicCharacteristics,
    .
  )

############ create a summary table of clinical outcomes

OutcomesSummaryTable <-
  GenotypeMatrixAndClinicalData %>%
  summarise(
    .data = .,
    ### main outcomes
      #"Main outcomes" = paste0(""),
      "Change in GAF (%)" = paste(signif(x = mean(x = ChangeInGAF,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = ChangeInGAF,na.rm = TRUE),digits = 3)),
      "Change in GAF range (%)" = paste0("(",min(ChangeInGAF,na.rm = TRUE)," - ",max(ChangeInGAF,na.rm = TRUE),")"),
      "Days to readmit." = paste(signif(x = mean(x = DaysToReadmittance,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = DaysToReadmittance,na.rm = TRUE),digits = 3)),
      "Days to readmit. range" = paste0("(",min(DaysToReadmittance,na.rm = TRUE)," - ",max(DaysToReadmittance,na.rm = TRUE),")"),
      "GAF (admit)" = paste(signif(x = mean(x = GAF_admit,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = GAF_admit,na.rm = TRUE),digits = 3)),
      "GAF (admit) range" = paste0("(",min(GAF_admit,na.rm = TRUE)," - ",max(GAF_admit,na.rm = TRUE),")"),
      "GAF (discharge)" = paste(signif(x = mean(x = GAF_discharge,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = GAF_discharge,na.rm = TRUE),digits = 3)),
      "GAF (discharge) range" = paste0("(",min(GAF_discharge,na.rm = TRUE)," - ",max(GAF_discharge,na.rm = TRUE),")"),
      "Length of stay (days)" = paste(signif(x = mean(x = LengthOfStay,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = LengthOfStay,na.rm = TRUE),digits = 3)),
      "Length of stay range (days)" = paste0("(",min(LengthOfStay,na.rm = TRUE)," - ",max(LengthOfStay,na.rm = TRUE),")"),
      "Max risp. dose (mg)" = paste(signif(x = mean(x = MaxDoseRisperidone_mg,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = MaxDoseRisperidone_mg,na.rm = TRUE),digits = 3)),
      "Max risp. dose range (%)" = paste0("(",min(MaxDoseRisperidone_mg,na.rm = TRUE)," - ",max(MaxDoseRisperidone_mg,na.rm = TRUE),")"),
      "Readmitted on risperidone (N)" = ReadmittedOnRisperidone %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
      "Readmitted on risperidone (%)" = round(x = ((ReadmittedOnRisperidone %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE))/nrow(x = GenotypeMatrixAndClinicalData))*100,digits = 1),
      "Patient readmit. (N)" = PatientReadmitted %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
      "Patient readmit. (%)" = round(x = ((PatientReadmitted %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE))/nrow(x = GenotypeMatrixAndClinicalData))*100,digits = 1),
    # ### days in acute and residential treatment
    #"Acute and residential length of stay" = paste0(""),
    #"Acute days" = paste(signif(x = mean(x = acute_days,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = acute_days,na.rm = TRUE),digits = 3)),
    #"Acute days range" = paste0("(",min(acute_days,na.rm = TRUE)," - ",max(acute_days,na.rm = TRUE),")"),
    #"Residential days" = paste(signif(x = mean(x = residential_days,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = residential_days,na.rm = TRUE),digits = 3)),
    #"Residential days range" = paste0("(",min(residential_days,na.rm = TRUE)," - ",max(residential_days,na.rm = TRUE),")"),
    # ### unit transfers
    #   "Unit transfers" = paste(signif(x = mean(x = UnitTransfers,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = UnitTransfers,na.rm = TRUE),digits = 3)),
    # ### discharge unit
    #   "Discharge unit (acute)" = sum(as.numeric(x = Discharge_unit=="a"),na.rm = TRUE),
    #   "Discharge unit (residential)" = sum(as.numeric(x = Discharge_unit=="r"),na.rm = TRUE),
    # ### number of residential and non-residential hospital visits
    #   "Non-residential visits" = sum(as.numeric(x = Residential_visit=="n"),na.rm = TRUE),
    #   "Residential visits" = sum(as.numeric(x = Residential_visit=="y"),na.rm = TRUE),
    # ### readmittance
    #   "Readmittance" = paste0(""),
    #   "Readmit. unit (acute)" = sum(as.numeric(x = Readmitted_unit=="a"),na.rm = TRUE),
    #   "Readmit. unit (residential)" = sum(as.numeric(x = Readmitted_unit=="r"),na.rm = TRUE),
    #   "Readmit. risperidone dose (mg)" = paste(signif(x = mean(x = ReadmittedRisperidoneDose_mg,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = ReadmittedRisperidoneDose_mg,na.rm = TRUE),digits = 3)),
    # "Antipsychotic medications" = paste0(""),
    # ### aripiprazole
    #   "Days on aripiprazole" = paste(signif(x = mean(x = DaysOnAripiprazole,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = DaysOnAripiprazole,na.rm = TRUE),digits = 3)),
    #   "Aripiprazole discontinued (N)" = sum(as.numeric(x = AripiprazoleAfterDischarge=="n"),na.rm = TRUE),
    # ### chlorpromazine
    #   "Days on chlorpromazine" = paste(signif(x = mean(x = DaysOnChlorpromazine,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = DaysOnChlorpromazine,na.rm = TRUE),digits = 3)),
    #   "Chlorpromazine discontinued (N)" = sum(as.numeric(x = ChlorpromazineAfterDischarge=="n"),na.rm = TRUE),
    # ### haloperidol
    #   "Days on haloperidol" = paste(signif(x = mean(x = DaysOnHaloperidol,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = DaysOnHaloperidol,na.rm = TRUE),digits = 3)),
    #   "Haloperidol discontinued (N)" = sum(as.numeric(x = HaloperidolAfterDischarge=="n"),na.rm = TRUE),
    # ### olanzapine
    #   "Days on olanzapine" = paste(signif(x = mean(x = DaysOnOlanzapine,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = DaysOnOlanzapine,na.rm = TRUE),digits = 3)),
    #   "Olanzapine discontinued (N)" = sum(as.numeric(x = OlanzapineAfterDischarge=="n"),na.rm = TRUE),
    # ### paliperidone
    #   "Days on paliperidone" = paste(signif(x = mean(x = DaysOnPaliperidone,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = DaysOnPaliperidone,na.rm = TRUE),digits = 3)),
    #   "Paliperidone discontinued (N)" = sum(as.numeric(x = PaliperidoneAfterDischarge=="n"),na.rm = TRUE),
    # ### quetiapine
    #   "Days on quetiapine" = paste(signif(x = mean(x = DaysOnQuetiapine,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = DaysOnQuetiapine,na.rm = TRUE),digits = 3)),
    #   "Quetiapine discontinued (N)" = sum(as.numeric(x = QuetiapineAfterDischarge=="n"),na.rm = TRUE),
    ### risperidone
      "Days treated with risperidone" = paste(signif(x = mean(x = DaysOnRisperidone,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = DaysOnRisperidone,na.rm = TRUE),digits = 3)),
      "Days treated with risperidone range" = paste0("(",min(DaysOnRisperidone,na.rm = TRUE)," - ",max(DaysOnRisperidone,na.rm = TRUE),")"),
      # "Risperidone before admit (y)" = sum(as.numeric(x = RisperidoneBeforeAdmit=="y"),na.rm = TRUE),
      # "Risperidone before admit (n)" = sum(as.numeric(x = RisperidoneBeforeAdmit=="n"),na.rm = TRUE),
      # "Risperidone before admit (NA)" = sum(as.numeric(x = is.na(x =RisperidoneBeforeAdmit)),na.rm = TRUE),
      "Risperidone discontinued (N)" = Risperidone_discontinued %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
      "Risperidone discontinued (%)" = round(x = ((Risperidone_discontinued %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE))/nrow(x = GenotypeMatrixAndClinicalData))*100,digits = 1),
    # ### ziprasidone
    #   "Days on ziprasidone" = paste(signif(x = mean(x = DaysOnZiprasidone,na.rm = TRUE),digits = 3),"±",signif(x = sd(x = DaysOnZiprasidone,na.rm = TRUE),digits = 3)),
    #   "Ziprasidone discontinued (N)" = sum(as.numeric(x = ZiprasidoneAfterDischarge=="n"),na.rm = TRUE),
    # ### reason for risperidone discontinuation
    # "Reason for risperidone discontinuation" = paste0(""),
    #   "General side effect category (N)" = ReasonForRispDiscontinuation_GeneralSideEffectCategory %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
    #       # "Reason for risp. discon. (aggression)" = ReasonForRispDiscontinuation_aggression %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
    #       # "Reason for risp. discon. (attentional issues)" = ReasonForRispDiscontinuation_AttentionalIssues %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
    #       # "Reason for risp. discon. (drowsiness)" = ReasonForRispDiscontinuation_drowsiness %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
    #       # "Reason for risp. discon. (extrapyramidal symptoms)" = ReasonForRispDiscontinuation_ExtrapyramidalSymptoms %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
    #       # "Reason for risp. discon. (insomnia)" = ReasonForRispDiscontinuation_Insomnia %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
    #       # "Reason for risp. discon. (irritability)" = ReasonForRispDiscontinuation_irritability %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
    #       # "Reason for risp. discon. (lactation side effect)" = ReasonForRispDiscontinuation_LactationRelated %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
    #       # "Reason for risp. discon. (nausea)" = ReasonForRispDiscontinuation_Nausea %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
    #       # "Reason for risp. discon. (side effects)" = ReasonForRispDiscontinuation_SideEffects %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
    #       # "Reason for risp. discon. (suicidal ideation)" = ReasonForRispDiscontinuation_SuicidalIdeation %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
    #       # "Reason for risp. discon. (weight gain)" = ReasonForRispDiscontinuation_WeightGain %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
    #   "General lack of efficacy category (N)" = ReasonForRispDiscontinuation_GeneralLackOfEfficacy %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
    #       # "Reason for risp. discon. (efficacy)" = ReasonForRispDiscontinuation_Efficacy %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
    #       # "Reason for risp. discon. (not helpful)" = ReasonForRispDiscontinuation_NotHelpful %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
    #   "General other reason category (N)" = ReasonForRispDiscontinuation_GeneralOtherReason %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
          # "Reason for risp. discon. (diagnosis change)" = ReasonForRispDiscontinuation_DiagnosisChange %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
          # "Reason for risp. discon. (medication washout)" = ReasonForRispDiscontinuation_MedicationWash %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
          # "Reason for risp. discon. (N/A)" = ReasonForRispDiscontinuation_NotAvailable %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
          # "Reason for risp. discon. (patient refusal)" = ReasonForRispDiscontinuation_PatientRefusal %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
          # "Reason for risp. discon. (prescriber decision)" = ReasonForRispDiscontinuation_PrescriberDecision %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
          # "Reason for risp. discon. (prescriber fear of side effects)" = ReasonForRispDiscontinuation_PrescriberFearOfSideEffects %>% as.character(x = .) %>% as.numeric(x = .) %>% sum(.,na.rm = TRUE),
    # ### reason for readmittance
    # "Reason for readmittance" = paste0(""),
    #   "Reason for readmit. (aggression)" = sum(as.numeric(x = ReasonForReadmittance_Aggression!="0"),na.rm = TRUE),
    #   "Reason for readmit. (behavioral issues)" = sum(as.numeric(x = ReasonForReadmittance_BehavioralIssues!="0"),na.rm = TRUE),
    #   "Reason for readmit. (evaluation and stabilization)" = sum(as.numeric(x = ReasonForReadmittance_EvaluationAndStabilization!="0"),na.rm = TRUE),
    #   "Reason for readmit. (homicidal ideation)" = sum(as.numeric(x = ReasonForReadmittance_HomicidalIdeation!="0"),na.rm = TRUE),
    #   "Reason for readmit. (increasing mood dysregulation)" = sum(as.numeric(x = ReasonForReadmittance_IncreasingMoodDysregulation!="0"),na.rm = TRUE),
    #   "Reason for readmit. (intoxication)" = sum(as.numeric(x = ReasonForReadmittance_Intoxication!="0"),na.rm = TRUE),
    #   "Reason for readmit. (psychosis)" = sum(as.numeric(x = ReasonForReadmittance_psychosis!="0"),na.rm = TRUE),
    #   "Reason for readmit. (school behavior)" = sum(as.numeric(x = ReasonForReadmittance_SchoolBehavior!="0"),na.rm = TRUE),
    #   "Reason for readmit. (self harm)" = sum(as.numeric(x = ReasonForReadmittance_SelfHarm!="0"),na.rm = TRUE),
    #   "Reason for readmit. (sexualized behavior)" = sum(as.numeric(x = ReasonForReadmittance_SexualizedBehavior!="0"),na.rm = TRUE),
    #   "Reason for readmit. (suicidal threats)" = sum(as.numeric(x = ReasonForReadmittance_SuicidalThreats!="0"),na.rm = TRUE),
    #   "Reason for readmit. (unknown)" = sum(as.numeric(x = ReasonForReadmittance_Unknown!="0"),na.rm = TRUE)
    ) %>%
    # transpose the rows and columns
    t(x = .) %>%
    # convert to a dataframe
    data.frame(.) %>%
    # convert the rownames to a column
    rownames_to_column(
      .data = .,
      var = "Characteristic"
      ) %>%
    select(
      .data = .,
      Characteristic,
      "Value" = `.`
      ) 

# save summary statistics tables to the file system
if(!dir.exists(paths = "./DemographicSummaryStatistics/"))
{
      dir.create(path = "./DemographicSummaryStatistics/")
}

list(
  "BaselineClinicalAndDemographicCharacteristics" = BaselineClinicalAndDemographicCharacteristics,
  "OutcomesSummaryTable" = OutcomesSummaryTable
  ) %>%
  mapply(
    FUN = function(currentTable,currentTableName)
    {
        currentTable %>%
        write_tsv(
          x = .,
          file = glue("./DemographicSummaryStatistics/{currentTableName}.tsv")
          )
    },
    .,
    names(x = .)
      )

############## create box and whisker plots and histograms of continuous outcomes

  c(
    "ChangeInGAF",
    "DaysOnRisperidone",
    "DaysToReadmittance",
    "GAF_admit",
    "GAF_discharge",
    "LengthOfStay",
    "MaxDoseRisperidone_mg"
    ) %>%
  lapply(
    X = .,
    FUN = function(currentOutcome)
    {
      DataToPlot <-
        GenotypeMatrixAndClinicalData %>%
        # select the current outcome
        select(
          .data = .,
          "Outcome" = all_of(x = currentOutcome)
          ) %>%
        # remove missing values
        na.omit(object = .)

      # create the box and whisker plot
      BoxAndWhiskerPlot <-
        DataToPlot %>%
          ggplot(
            data = .,
            mapping = aes(
                          x = "",
                          y = Outcome
                            )
              ) +
          geom_boxplot(outlier.shape = NA) +
          geom_point(
            data = DataToPlot,
            mapping = aes(x = "",y = Outcome),
            shape = 21,
            colour = "black",
            fill = "darkorchid1",
            size = 1.5,
            stroke = 2/3,
            alpha = 1,
            position = position_jitter()
            ) +
          theme_classic() +
          scale_y_continuous(
            labels = seq(
                         0,
                         max(DataToPlot$Outcome,na.rm = TRUE)+5,
                         max(DataToPlot$Outcome,na.rm = TRUE)/10
                         ) %>%
                      round(x = .,digits = 4),
            breaks = seq(
                         0,
                         max(DataToPlot$Outcome,na.rm = TRUE)+5,
                         max(DataToPlot$Outcome,na.rm = TRUE)/10
                         ) %>%
                     round(x = .,digits = 4)
              ) +
          xlab(label = currentOutcome) +
          ylab(label = "Value") +
          theme(
            axis.title = element_text(family = "Arial",face = "bold",size = 18),
            axis.text.x = element_text(family = "Arial",face = "bold",size = 18),
            axis.text.y = element_text(family = "Arial",face = "bold",size = 18)
          )

      # create a directory for saving the box and whisker plots if it does not exist
      if(!dir.exists(paths = "./BoxAndWhiskerPlotsContinuousOutcomes/"))
      {
        dir.create(path = "./BoxAndWhiskerPlotsContinuousOutcomes/")
      }

      tiff(
        filename = glue("./BoxAndWhiskerPlotsContinuousOutcomes/{currentOutcome}_BoxAndWhisker.tiff"),
        width = 5,
        height = 7,
        units = "in",
        compression = "none",
        res = 300
        )

      print(x = BoxAndWhiskerPlot)

      dev.off()

      # perform a Shapiro-Wilke's test for normality with the current outcome
      TestForNormalityResult <-
        DataToPlot$Outcome %>%
        shapiro.test(x = .)

      # calculate the binwidth required for 30 histogram bins
      HistogramBinWidth <-
        DataToPlot %>%
        # calculate value max, min, and the binwidth for 30 bins
        summarise(
          .data = .,
          "MaxValue" = max(Outcome,na.rm = TRUE),
          "MinValue" = min(Outcome,na.rm = TRUE),
          "BinWidth" = (max(Outcome,na.rm = TRUE) - min(Outcome,na.rm = TRUE))/30
        ) %>%
        pull(.data = .,BinWidth)

      HistogramPlot <-
          DataToPlot %>%
          # create a histogram
          ggplot(
            data = .,
            mapping = aes(x = Outcome)
          ) +
          geom_histogram(
            binwidth = HistogramBinWidth,
            fill="#69b3a2",
            color="black",
            alpha=0.9
          ) +
          annotate(
            geom = "text",
            x = max(DataToPlot$Outcome)*0.75,
            y = 40,
            label = paste0("Non-normality test (p = ",signif(x = TestForNormalityResult$p.value,digits = 3),")"),
            size = 5,
            family = "Arial"
            ) +
          theme_classic() +
          xlab(label = currentOutcome) +
          ylab(label = "Participant count") +
          ggtitle(label = currentOutcome) +
          theme(
            axis.title = element_text(family = "Arial",face = "bold",size = 18),
            axis.text.x = element_text(family = "Arial",face = "bold",size = 18),
            axis.text.y = element_text(family = "Arial",face = "bold",size = 18),
            plot.title = element_text(family = "Arial",face = "bold",size = 11,hjust = 0.5)
          )

      # save the final plot to the file system as a tiff
      if(!dir.exists(paths = "./HistogramsContinuousOutcomes/"))
      {
        dir.create(path = "./HistogramsContinuousOutcomes/")
      }

      tiff(
        filename = glue("./HistogramsContinuousOutcomes/{currentOutcome}_histogram.tiff"),
        width = 7,
        height = 5,
        units = "in",
        compression = "none",
        res = 300
        )

      print(x = HistogramPlot)



      dev.off()

      return(NULL)
    }
      )

##############

# define an array of the continuous study outcomes
ContinuousStudyOutcomes <-
  c(
    "ChangeInGAF",
    "DaysOnRisperidone",
    "DaysToReadmittance",
    "GAF_admit",
    "GAF_discharge",
    "LengthOfStay",
    "MaxDoseRisperidone_mg"
  )

# define an array of the binary study outcomes
BinaryStudyOutcomes <-
  c(
    "ReadmittedOnRisperidone",
    "PatientReadmitted",
    "Risperidone_discontinued"
  )

# define an array of covariates for each outcome
CovariatesForStudyOutcomes <-
  c(
    "Age_years",
    "Gender_MorF",
    "BMI_kgPerSqM",
    "GAF_admit",
    "Residential_visit",
    "Insurance",
    "PsychiatricHistory_family",
    "NeglectAndAbuseHistory",
    "SubstanceAbuseHistory",
    "cyp2d6_modulator"#,
    #"dopaminergic_modulator",
    #"serotonergic_modulator"
  )

#############

#### create a directory for saving miscellaneous plots if it does not exist
if(!dir.exists(paths = "./miscellaneousPlots/"))
{
  dir.create(path = "./miscellaneousPlots/")
}


  # create a scatter plot of days to readmittance versus participant age
  # and max risperidone dose versus participant age
  c("DaysToReadmittance","MaxDoseRisperidone_mg") %>%
  lapply(
    X = .,
    FUN = function(currentOutcome)
    {
        PlotToReturn <-
          GenotypeMatrixAndClinicalData %>%
          select(
            .data = .,
            Age_years,
            "Outcome" = all_of(x = currentOutcome)
            ) %>%
          na.omit(object = .) %>%
          ggplot(
            data = .,
            mapping = aes(x = Age_years,y = Outcome)
              ) +
        geom_point() +
        theme_classic() +
        xlab(label = "Age (years)") +
        ylab(label = currentOutcome) +
        theme(
          axis.title = element_text(family = "Arial",face = "bold",size = 11),
          axis.text.x = element_text(family = "Arial",face = "bold",size = 11),
          axis.text.y = element_text(family = "Arial",face = "bold",size = 11)
        )
        
        tiff(
          filename = glue("./miscellaneousPlots/{currentOutcome}_versus_age_scatterPlot.tiff"),
          units = "in",
          compression = "none",
          res = 300,
          width = 7,
          height = 5
        )
        
        print(x = PlotToReturn)
        dev.off()
    }
      )

#############
  
############## perform primary linear regression analysis of continuous outcomes versus individual clinical and demographic covariates
  
  # create a directory for saving primary regressions
  if(!dir.exists(paths = "./primaryRegressions/"))
  {
      dir.create(path = "./primaryRegressions/")
  }
  
  primaryLinearRegressionsContinuousOutcomes <-
    ContinuousStudyOutcomes %>%
    # loop through each continuous outcome
    lapply(
      X = .,
      FUN = function(currentOutcome)
      {
        
          RegressionResults_RSE <-
                # loop through each covariate
                CovariatesForStudyOutcomes %>%
                  lapply(
                    X = .,
                    FUN = function(currentPredictor)
                    {
                          # select the data for the regression and remove missing values
                          DataForRegression <-
                            GenotypeMatrixAndClinicalData %>%
                            select(
                              .data = .,
                              "Outcome" = all_of(x = currentOutcome),
                              "Predictor" = all_of(x = currentPredictor)
                              ) %>%
                            na.omit(object = .)
                      
                          # if the current outcome is not the current predictor,
                          # and the current predictor has more than one unique category or value, perform the regression
                          if(
                            (currentOutcome!=currentPredictor) &
                            (length(x = unique(x = DataForRegression$Predictor))>1)
                            )
                          {
                            
                              # perform a linear regression with robust standard errors of the outcome versus the predictor
                              RegressionModel <-
                                DataForRegression %>%
                                lm_robust(
                                  formula = Outcome ~ Predictor,
                                  data = .
                                  )

                              # perform an F-test for variance homogeneity to see if there is any heteroscedasticity with
                              # the lm() function; the lm_robust() function cannot be used in the ols_test_f() function
                              VarianceHomogeneityTestResult <-
                                DataForRegression %>%
                                lm(
                                  formula = Outcome ~ Predictor,
                                  data = .
                                ) %>%
                                ols_test_f(
                                  model = .
                                  )

                              RegressionResults_RSE <-
                                # create a tidy summary of the linear robust standard errors regression
                                TidyLinearRobustStandardErrorsRegressionResults(
                                  DataForRegression = DataForRegression,
                                  RegressionModel = RegressionModel
                                  ) %>%
                                 # add a label for the current outcome and predictor
                                 # and a label for the regression performed
                                mutate(
                                  .data = .,
                                  Outcome = currentOutcome,
                                  Predictor = currentPredictor,
                                  Regression = paste0(currentOutcome," ~ ",currentPredictor)
                                  ) %>%
                                # add the P-value from the test for Heteroscedastiticy with astericks for levels of significance
                                mutate(
                                  .data = .,
                                  Heteroscedastiticy_p.value = VarianceHomogeneityTestResult$p %>%
                                                               signif(x = .,digits = 3) %>%
                                                               AddAstericksToPvalues(columnVector = .)
                                  ) %>%
                                # filter out the intercept term
                                filter(
                                  .data = .,
                                  term != "(Intercept)"
                                  )
                              
                              return(RegressionResults_RSE)
                           }


                    }
                      ) %>%
                  do.call(
                    what = "rbind",
                    args = .
                    )

      return(RegressionResults_RSE)

      }
        ) %>%
      do.call(
        what = "rbind",
        args = .
        ) %>%
    # select useful columns
    # and update column names
    select(
      .data = .,
      "Outcome" = Outcome,
      Predictor,
      "Term" = term,
      "Beta ± S.E." = estimate_std.error,
      "P-value" = p.value,
      "F-stat. P-value" = F_stat_p_value,
      R2,
      "R2-adj." = R2_adj,
      N,
      Regression
      )
  
  primaryLinearRegressionsContinuousOutcomes %>%
    write_delim(
      x = .,
      file = "./primaryRegressions/primaryLinearRegressionsContinuousOutcomes.txt",
      delim = "\t"
        )
  
##############

############## perform multiple regressions with robust standard errors of each continuous phenotype versus 
             # 10 multidimensional scaling components that are used to adjust for genomic population stratification
             # to estimate the variability explained by MDS components alone for each continuous phenotype
  MultipleLinearRegresssionsContinuousOutcomesMDScomponents <-
    ContinuousStudyOutcomes %>%
    # loop through each continuous outcome
    lapply(
      X = .,
      FUN = function(currentOutcome)
      {
            DataForRegression <-
              GenotypeMatrixAndClinicalData %>%
              # select the current outcome and MDS components
              select(
                .data = .,
                "Outcome" = all_of(x = currentOutcome),
                C1:C10
              ) %>%
              # remove observations with missing values
              na.omit(object = .)
            
            # define the regression formula: current outcome ~ 10 MDS components
            RegressionFormula <-
              DataForRegression %>% 
              select(
                .data = .,
                -Outcome,
                C1:C10
                ) %>%
              names(x = .) %>%
              paste(collapse = " + ") %>%
              paste(
                "Outcome",
                "~",
                .
              )
              
            # perform a linear regression with robust standard errors of the outcome versus the 10 MDS components
            RegressionModel <-
              DataForRegression %>%
              lm_robust(
                formula = as.formula(RegressionFormula),
                data = .
                )
            
            # perform an F-test for variance homogeneity to see if there is any heteroscedasticity with
            # the lm() function; the lm_robust() function cannot be used in the ols_test_f() function
            VarianceHomogeneityTestResult <-
              DataForRegression %>%
              lm(
                formula = as.formula(RegressionFormula),
                data = .
              ) %>%
              ols_test_f(
                model = .
                )

            RegressionResults_RSE <-
              # create a tidy summary of the linear robust standard errors regression
              TidyLinearRobustStandardErrorsRegressionResults(
                DataForRegression = DataForRegression,
                RegressionModel = RegressionModel
                ) %>%
               # add a label for the current outcome
               # and a label for the regression performed
              mutate(
                .data = .,
                Outcome = currentOutcome,
                Regression = RegressionFormula
                ) %>%
              # add the P-value from the test for Heteroscedastiticy with astericks for levels of significance
              mutate(
                .data = .,
                Heteroscedastiticy_p.value = VarianceHomogeneityTestResult$p %>%
                                             signif(x = .,digits = 3) %>%
                                             AddAstericksToPvalues(columnVector = .)
                ) %>%
              # filter out the intercept term
              filter(
                .data = .,
                term != "(Intercept)"
                )

            return(RegressionResults_RSE)

      }
        ) %>%
      do.call(
        what = "rbind",
        args = .
        ) %>%
    # select useful columns
    # and update column names
    select(
      .data = .,
      "Outcome" = Outcome,
      "F-stat. P-value" = F_stat_p_value,
      R2,
      "R2-adj." = R2_adj,
      N,
      Regression
      ) %>%
    unique(x = .)
  
  # create a directory for saving MDS components regressions
  if(!dir.exists(paths = "./VariabilityExplainedByMDScomponentsOnly/"))
  {
      dir.create(path = "./VariabilityExplainedByMDScomponentsOnly/")
  }
  
  # save the MDS components regressions to the file system
  MultipleLinearRegresssionsContinuousOutcomesMDScomponents %>%
    write_delim(
      x = .,
      file = "./VariabilityExplainedByMDScomponentsOnly/VariabilityExplainedByMDScomponentsOnly.txt",
      delim = "\t"
        )
  
############## 
  
############## perform primary logistic regression analysis of binary outcomes versus individual clinical and demographic covariates
  
  # Temporarily suppress non-problematic warning messages that are produced from the primary logistic regressions
  # for cases with very low sample sizes
  options(warn = -1)   
  
  primaryLogisticRegressionsBinaryOutcomes <-
      BinaryStudyOutcomes %>%
        # loop through each binary outcome
       lapply(
         X = .,
         FUN = function(currentOutcome)
         {
             ResultsToReturn <-
                 # loop through each covariate
                   CovariatesForStudyOutcomes %>%
                   lapply(
                     X = .,
                     FUN = function(currentPredictor)
                     {
                       
                         # select data for the regression and remove missing values
                         DataForRegression <-
                           GenotypeMatrixAndClinicalData %>%
                           select(
                             .data = .,
                             "Outcome" = all_of(x = currentOutcome),
                             "Predictor" = all_of(x = currentPredictor)
                             ) %>%
                           na.omit(object = .)
                       
                           # if the current outcome is not the current predictor,
                           # and the current predictor has more than one unique category or values, 
                           # and the current outcome has more than one unique category or values,
                           # perform the regression
                           if(
                             (currentOutcome!=currentPredictor) &
                             (length(x = unique(x = DataForRegression$Predictor))>1) &
                             (length(x = unique(x = DataForRegression$Outcome))>1)
                             )
                           {
                             
                                 # set a random number seed in case one is required
                                 # for reproducibility of the logistic regression
                                 set.seed(seed = 777)
                                 
                                 # define the logistic regression model
                                 RegressionModel <-
                                   DataForRegression %>%
                                   glm(
                                     formula = Outcome ~ Predictor,
                                     family = "binomial",
                                     data = .
                                   )
  
                                 # tidy the regression summary and return the results
                                 ResultsToReturn <-
                                   TidyLogisticRegressionResults(
                                     DataForRegression = DataForRegression,
                                     RegressionModel = RegressionModel
                                     ) %>%
                                    mutate(
                                     .data = .,
                                     # add columns with the outcome and predictor
                                     Outcome = currentOutcome,
                                     Predictor = currentPredictor,
                                     # add a column for the regression
                                     Regression = paste0(currentOutcome," ~ ",currentPredictor)
                                     )

                                 return(ResultsToReturn)
                           }
                      
                     }) %>%
                    do.call(
                      what = "rbind",
                      args = .
                      )
             
             return(ResultsToReturn)
         }
       ) %>%
      do.call(
        what = "rbind",
        args = .
        )
  
  # turn warning messages back on 
  options(warn = 0)
  
  primaryLogisticRegressionsBinaryOutcomes %>%
    write_delim(
      x = .,
      file = "./primaryRegressions/primaryLogisticRegressionsBinaryOutcomes.txt",
      delim = "\t"
        )
  
############### 
  
# ############# Perform linear regressions with robust standard errors
# of continuous outcomes versus individual significant variant genotypes from the GWAS for the variants
# that had a significant association with a phenotype after correction for multiple comparisons.
# The only phenotypes that had a significant variant association in the GWAS were:
    # DaysOnRisperidone, DaysToReadmittance, LengthOfStay, and MaxDoseRisperidone_mg.      

GenotypeRegressionResults_RSE <-
    # loop through each of the significant variants from the GWAS
    SignificantVariantsFromGWAS %>%
      group_by(
        .data = .,
        SNP
        ) %>%
      group_split(.tbl = .) %>%
      lapply(
        X = .,
        FUN = function(currentVariantDataFrame)
        {
              # select the current variant and phenotype
              currentVariant <- currentVariantDataFrame$SNP %>%
                                  gsub(
                                    pattern = ":",
                                    replacement = ".",
                                    x = .
                                    )

              currentPhenotype <- currentVariantDataFrame$Phenotype
              
              # select the current alpha value for significance
              currentAlphaValue <- currentVariantDataFrame$alpha
              
              DataForGenotypeBoxplot <-
                # select the column from the genotype matrix with clinical data
                # for the current variant, current phenotype, and the participant IDs
                GenotypeMatrixAndClinicalData %>%
                  select(
                    .data = .,
                    ParticipantID,
                    "VariantGenotype" = contains(match = currentVariant),
                    "PhenotypeValue" = all_of(x = currentPhenotype)
                    ) %>%
                     # filter to values that are greater than zero so they can be plotted on a log-scale boxplot
                    filter(
                      .data = .,
                      PhenotypeValue > 0
                      ) 

              # create a table of observation counts for each genotype (0, 1, or 2)
              GenotypeCounts <-
                data.frame(
                  "Genotype_0" = sum(as.numeric(x = DataForGenotypeBoxplot$VariantGenotype==0),na.rm = TRUE),
                  "Genotype_1" = sum(as.numeric(x = DataForGenotypeBoxplot$VariantGenotype==1),na.rm = TRUE),
                  "Genotype_2" = sum(as.numeric(x = DataForGenotypeBoxplot$VariantGenotype==2),na.rm = TRUE),
                  row.names = NULL
                )

                # if the homozygous reference allele observation count is < 3,
                # update the data for regression for a recessive regression model
                if(GenotypeCounts$Genotype_0<as.numeric(x = 3))
                {
                  DataForGenotypeBoxplot <-
                      DataForGenotypeBoxplot %>%
                      mutate(
                        .data = .,
                        # if the genotype is 0 or 1 copies of the variant allele, set as 1 copy
                        # if the genotype is 2 copies of the variant allele, set as 2 copies
                        VariantGenotype = if_else(
                                                  condition = (VariantGenotype==as.numeric(x = 0)) | (VariantGenotype==as.numeric(x = 1)),
                                                  true = as.numeric(x = 1),
                                                  false = as.numeric(x = 2)
                                                    )
                        )
                }

                # if the homozygous alternate allele observation count is < 3,
                # update the data for regression for a dominant allele regression model
                if(GenotypeCounts$Genotype_2<as.numeric(x = 3))
                {
                  DataForGenotypeBoxplot <-
                        DataForGenotypeBoxplot %>%
                        mutate(
                          .data = .,
                          # if the genotype is 1 or 2 copies of the variant allele, set as 1 copy
                          # if the genotype is 0 copies of the variant allele, set as 0 copies
                          VariantGenotype = if_else(
                                                    condition = (VariantGenotype==as.numeric(x = 1)) | (VariantGenotype==as.numeric(x = 2)),
                                                    true = as.numeric(x = 1),
                                                    false = as.numeric(x = 0)
                                                      )
                          )
                }

                # select columns with the phenotype value and the variant genotype
                # for performing a regression with robust standard errors of the current phenotype versus genotype
                DataForRegression <-
                  DataForGenotypeBoxplot %>%
                  select(
                    .data = .,
                    "Y" = PhenotypeValue,
                    "X" = VariantGenotype
                    ) %>%
                  na.omit(object = .)

                # create a robust standard errors regression model of the phenotype value versus genotype
                RegressionModel <-
                  DataForRegression %>%
                  lm_robust(
                    formula = Y ~ X,
                    data = .
                    )

                # perform an F-test for variance homogeneity to see if there is any heteroscedasticity with
                # the lm() function; the lm_robust() function cannot be used in the ols_test_f() function
                VarianceHomogeneityTestResult <-
                  DataForRegression %>%
                  lm(
                    formula = Y ~ X,
                    data = .
                  ) %>%
                  ols_test_f(
                    model = .
                    )
                
                # format the genotype regression results table
                GenotypeRegressionResults_RSE <-
                  # create a tidy summary of the linear robust standard errors regression
                  TidyLinearRobustStandardErrorsRegressionResults(
                    DataForRegression = DataForRegression,
                    RegressionModel = RegressionModel
                    ) %>%
                   # add a label for the current variant, phenotype, and gene
                  mutate(
                    .data = .,
                    Variant = currentVariantDataFrame$SNP,
                    Gene = currentVariantDataFrame$Gene.refGeneWithVer,
                    Phenotype = currentPhenotype
                    ) %>%
                  # add the P-value from the test for Heteroscedastiticy with astericks for levels of significance
                  mutate(
                    .data = .,
                    Heteroscedastiticy_p.value = VarianceHomogeneityTestResult$p %>%
                                                 signif(x = .,digits = 3) %>%
                                                 AddAstericksToPvalues(columnVector = .)
                    ) %>%
                  # filter out the intercept term
                  filter(
                    .data = .,
                    term != "(Intercept)"
                    )

                # append the suffix "_RSE" to each column name to indicate robust standard errors
                names(x = GenotypeRegressionResults_RSE) <- GenotypeRegressionResults_RSE %>%
                                                    names(x = .) %>%
                                                    sapply(
                                                      X = .,
                                                      FUN = function(currentColumnName)
                                                      {
                                                         ValueToReturn <-
                                                           currentColumnName %>%
                                                           paste0(.,"_RSE")
                                                      },
                                                      simplify = TRUE
                                                        )

                # create a table of genotype codes (0,1,2) and their corresponding nucleotide genotypes
                GenotypeLabelTable <-
                  data.frame(
                             "Genotype_Numeric" = c(0,1,2),
                             "Genotype_alleles" = c(
                                                    paste0(currentVariantDataFrame$Ref,"/",currentVariantDataFrame$Ref),
                                                    paste0(currentVariantDataFrame$Ref,"/",currentVariantDataFrame$Alt),
                                                    paste0(currentVariantDataFrame$Alt,"/",currentVariantDataFrame$Alt)
                                                    )
                            )

              # create a variant ID label for the plot
                # If the variant ID has an rs-number, use the rs-number
                # otherwise use the variant ID
              VariantLabel <-
                VariantIDrsNumberMap %>%
                filter(
                  .data = .,
                  VariantID_fromInputFile == currentVariantDataFrame$SNP
                  ) %>%
                mutate(
                  .data = .,
                  Label = if_else(
                                  condition = (rsNumber_VEP=="*rsNA") | (is.na(x = rsNumber_VEP)),
                                  true = VariantID_fromInputFile,
                                  false = rsNumber_VEP
                                    )
                 ) %>%
                pull(
                  .data = .,
                  Label
                  )
              
              # create a directory for storing the boxplots if it doesn't exist
              if(!dir.exists(paths = "./genotypeBoxplots/"))
              {
                dir.create(path = "./genotypeBoxplots/")
                dir.create(path = "./genotypeBoxplots/noCovariates/")
              }
              
              # set a scatter plot point jitter width
              PointJitterWidth <- 0.40
              PointJitterHeight <- 0
              
              # create the genotype boxplot without any covariates
              PlotToReturn <-
                  DataForGenotypeBoxplot %>%
                  select(
                    .data = .,
                    VariantGenotype,
                    PhenotypeValue
                    ) %>%
                  na.omit(object = .) %>%
                  ggplot(
                    data = .,
                    mapping = aes(
                                  x = factor(x = VariantGenotype),
                                  y = PhenotypeValue
                                  )
                      ) +
                    geom_boxplot(
                      lwd = 1.5,
                      fatten = 1,
                      outlier.shape = NA
                      ) +
                    scale_y_continuous(
                      n.breaks = 10
                      # limits = c(0.01,2000),
                      # breaks = c(1,10,50,100,250,500,1000,2000),
                      # labels = c(1,10,50,100,250,500,1000,2000)
                    ) +
                    geom_point(
                      shape = 21,
                      colour = "black",
                      fill = "darkgoldenrod1",
                      size = 3,
                      stroke = 2/3,
                      alpha = 1,
                      position = position_jitter(
                                                 seed = 777,
                                                 width = PointJitterWidth,
                                                 height = PointJitterHeight
                                                )
                      ) +
                    theme_classic() +
                    # convert the genotypes of 0/1/2 on the x-axis to genotypes of ref,ref/ref,alt/alt,alt
                    # using the genotype label table define above
                    scale_x_discrete(
                      labels =  GenotypeLabelTable %>%
                                filter(
                                  .data = .,
                                  Genotype_Numeric %in% sort(x = unique(x = DataForGenotypeBoxplot$VariantGenotype))
                                ) %>%
                                pull(.data = .,Genotype_alleles)
                    ) +
                    xlab(label = paste(currentVariantDataFrame$Gene.refGeneWithVer,VariantLabel)) +
                    ylab(label = currentPhenotype) +
                    theme(
                      axis.title = element_text(
                                                face = "bold",
                                                size = 11,
                                                family = "Arial"
                                                ),
                      axis.text.x = element_text(
                                                face = "bold",
                                                size = 11,
                                                family = "Arial",
                                                angle = 45,
                                                hjust = 1
                                                ),
                      axis.text.y = element_text(
                                                face = "bold",
                                                size = 11,
                                                family = "Arial"
                                                ),
                      legend.position = c(.20,.15),
                      legend.box.just = "right",
                      legend.margin = margin(6, 6, 6, 6),
                      legend.text = element_text(family = "Arial",face = "bold",size = 11),
                      legend.title = element_text(family = "Arial",face = "bold",size = 11,hjust = 0.5),
                      legend.box.background = element_rect(color = "black",size = 1)
                      ) #+
                # annotate(
                #   geom = "text",
                #   x = 2,
                #   y = 0.1,
                #   label = paste0(
                #                 # F-statistic P-value for overall regression significance
                #                 paste0("F-stat. P-value = ",GenotypeRegressionResults_RSE$F_stat_p_value_RSE),
                #                 "\n",
                #                 # Genotype term significance
                #                 paste0("P-value = ",GenotypeRegressionResults_RSE$p.value_RSE),
                #                 "\n",
                #                 # Bonferroni adjusted alpha value for significance
                #                 paste0("alpha = ",signif(x = currentAlphaValue,digits = 3)),
                #                 "\n",
                #                 # beta coefficient +/- standard error
                #                 paste0("beta = ",GenotypeRegressionResults_RSE$estimate_std.error_RSE),
                #                 "\n",
                #                 # unadjusted coefficient of determination (unadjusted R-squared)
                #                 paste0("R2 = ",GenotypeRegressionResults_RSE$R2_RSE),
                #                 "\n",
                #                 # number of observations included
                #                 paste0("N = ",GenotypeRegressionResults_RSE$N_RSE)
                #                 ),
                #   fontface = "bold",
                #   size = 4
                #     )

              # save the plot to the directory labeled with the current variant ID
              tiff(
                filename = glue("./genotypeBoxplots/noCovariates/{currentPhenotype}_{currentVariant}_plot.tiff"),
                units = "in",
                compression = "none",
                res = 300,
                width = 7,
                height = 5
                )
              # print the plot to return so it can be saved to the file system
              print(x = PlotToReturn)
              dev.off()
              
              return(GenotypeRegressionResults_RSE)
        }
        ) %>%
      # bind each dataframe of regression results together by row
      do.call(
        what = "rbind",
        args = .
        )

##############

# ############# Perform multiple linear regressions with robust standard errors
# of continuous outcomes versus all significant variant genotypes from the GWAS for the variants
# that had a significant association with a phenotype after correction for multiple comparisons.
# The only phenotypes that had a significant variant association in the GWAS were:
    # DaysOnRisperidone, DaysToReadmittance, LengthOfStay, and MaxDoseRisperidone_mg.

Multiple_GenotypeRegressionResults_RSE <-
        # loop through the significant variants from the GWAS for each phenotype
        SignificantVariantsFromGWAS %>%
          group_by(
            .data = .,
            Phenotype
            ) %>%
          group_split(.tbl = .) %>%
          # loop through each set of significant variants
          lapply(
            X = .,
            FUN = function(currentVariantDataFrame)
            {
              # select the current set of variants as an array
              currentVariantSet <- currentVariantDataFrame$SNP %>%
                                  # loop through each variant in the array and
                                  # reformat the variant ID to match the genotype matrix format
                                   sapply(
                                     X = .,
                                     FUN = function(currentVariant)
                                     {
                                       ValueToReturn <-
                                         currentVariant %>%
                                         # replace colons with periods
                                         gsub(
                                           pattern = ":",
                                           replacement = ".",
                                           x = .
                                           )

                                       return(ValueToReturn)
                                     },
                                     simplify = FALSE
                                       ) %>%
                                    unlist(x = .) %>%
                                    unname(obj = .)

              # select the current phenotype
              currentPhenotype <- currentVariantDataFrame$Phenotype %>%
                                  unique(x = .)
              
              DataForRegression <-
                  # select the columns from the genotype matrix with clinical data
                  # for the current set of variants and current phenotype
                  currentVariantSet %>%
                  # select the column corresponding to each variant from the GenotypeMatrixAndClinicalData
                  lapply(
                    X = .,
                    FUN = function(currentVariant)
                    {
                      DataToReturn <-
                        GenotypeMatrixAndClinicalData %>%
                          select(
                            .data = .,
                            contains(match = currentVariant)
                            )

                      return(DataToReturn)
                    }
                      ) %>%
                    do.call(
                      what = "cbind",
                      args = .
                      ) %>%
                    cbind(
                      .,
                      select(
                        .data = GenotypeMatrixAndClinicalData,
                        "PhenotypeValue" = all_of(x = currentPhenotype)
                        )
                      ) %>%
                      # remove observations with missing values
                      na.omit(object = .)
                  
                  # specify the regression model formula with the significant variant genotypes: 
                      # Phenotype ~ variant genotypes
                  RegressionFormula <-
                    DataForRegression %>%
                    select(
                      .data = .,
                      -PhenotypeValue
                      ) %>%
                    names(x = .) %>%
                    paste(
                      .,
                      collapse = " + "
                    ) %>%
                    paste(
                      "PhenotypeValue",
                      .,
                      sep = " ~ "
                    )
                  
                  # create a robust standard errors regression model of the phenotype value versus significant variant genotypes
                  RegressionModel <-
                    DataForRegression %>%
                    lm_robust(
                      formula = as.formula(object = RegressionFormula),
                      data = .
                      )

                  # perform an F-test for variance homogeneity to see if there is any heteroscedasticity with
                  # the lm() function; the lm_robust() function cannot be used in the ols_test_f() function
                  VarianceHomogeneityTestResult <-
                    DataForRegression %>%
                    lm(
                      formula = as.formula(object = RegressionFormula),
                      data = .
                    ) %>%
                    ols_test_f(
                      model = .
                      )
                  
                  GenotypeRegressionResults_RSE <-
                    # create a tidy summary of the linear robust standard errors regression
                    TidyLinearRobustStandardErrorsRegressionResults(
                      DataForRegression = DataForRegression,
                      RegressionModel = RegressionModel
                      ) %>%
                     # add a label for the current phenotype and regression model formula
                    mutate(
                      .data = .,
                      Phenotype = currentPhenotype,
                      Regression = RegressionFormula
                      ) %>%
                    # add the P-value from the test for Heteroscedastiticy with astericks for levels of significance
                    mutate(
                      .data = .,
                      Heteroscedastiticy_p.value = VarianceHomogeneityTestResult$p %>%
                                                   signif(x = .,digits = 3) %>%
                                                   AddAstericksToPvalues(columnVector = .)
                      ) %>%
                    # filter out the intercept term
                    filter(
                      .data = .,
                      term != "(Intercept)"
                      )

                  return(GenotypeRegressionResults_RSE)



            }
            ) %>%
        do.call(
          what = "rbind",
          args = .
          ) %>%
        # label the gene for each variant
        mutate(
          .data = .,
          Gene = term %>%
                 lapply(
                   X = .,
                   FUN = function(currentVariant)
                   {
                     # reformat the string of the current variant to match the
                     # format in the SignificantVariantsFromGWAS table
                     currentVariant <-
                       currentVariant %>%
                       gsub(pattern = "[.]",replacement = ":",x = .) %>%
                       str_remove(string = .,pattern = "_[:alpha:]*")

                     ValueToReturn <-
                        SignificantVariantsFromGWAS %>%
                        # extract the gene label corresponding to the current variant
                        filter(
                          .data = .,
                          str_detect(string = currentVariant,pattern = SNP)==TRUE
                          ) %>%
                        pull(
                          .data = .,
                          Gene.refGeneWithVer
                          )

                     return(ValueToReturn)
                   }
                   ) %>%
                  unlist(x = .)
          ) %>%
        # select useful columns
        # and update column names
        select(
          .data = .,
          Phenotype,
          "Term" = term,
          Gene,
          "Beta ± S.E." = estimate_std.error,
          "P-value" = p.value,
          "F-stat. P-value" = F_stat_p_value,
          R2,
          "R2-adj." = R2_adj,
          N,
          Regression
          ) 

# save the results to the file system
if(!dir.exists(paths = "./MultipleGenotypeRegressions/"))
{
        dir.create(path = "./MultipleGenotypeRegressions/")
}

Multiple_GenotypeRegressionResults_RSE %>%
  write_delim(
    x = .,
    file = "./MultipleGenotypeRegressions/Multiple_GenotypeRegressionResults_RSE.txt",
    delim = "\t"
      )

##############

# ############# Perform multiple linear regressions with robust standard errors with adjustment for demographic/clinical covariates
# of continuous outcomes versus all significant variant genotypes from the GWAS for the variants
# that had a significant association with a phenotype after correction for multiple comparisons.
# The only phenotypes that had a significant variant association in the GWAS were:
    # DaysOnRisperidone, DaysToReadmittance, LengthOfStay, and MaxDoseRisperidone_mg.

# The demographic/clinical covariates that are included in this model are:
    # Age_years, Gender_MorF, BMI_kgPerSqM, Residential_visit, Insurance, PsychiatricHistory_family, NeglectAndAbuseHistory, and SubstanceAbuseHistory; and cyp2d6, dopaminergic, and serotonergic modulators

Multiple_GenotypeRegressionResults_RSE_withDemographicCovariates <-
        # loop through the significant variants from the GWAS for each phenotype
        SignificantVariantsFromGWAS %>%
          group_by(
            .data = .,
            Phenotype
            ) %>%
          group_split(.tbl = .) %>%
          # loop through each set of significant variants
          lapply(
            X = .,
            FUN = function(currentVariantDataFrame)
            {
              # select the current set of variants as an array
              currentVariantSet <- currentVariantDataFrame$SNP %>%
                                  # loop through each variant in the array and
                                  # reformat the variant ID to match the genotype matrix format
                                   sapply(
                                     X = .,
                                     FUN = function(currentVariant)
                                     {
                                       ValueToReturn <-
                                         currentVariant %>%
                                         # replace colons with periods
                                         gsub(
                                           pattern = ":",
                                           replacement = ".",
                                           x = .
                                           )

                                       return(ValueToReturn)
                                     },
                                     simplify = FALSE
                                       ) %>%
                                    unlist(x = .) %>%
                                    unname(obj = .)

              # select the current phenotype
              currentPhenotype <- currentVariantDataFrame$Phenotype %>%
                                  unique(x = .)
              
              DataForRegression <-
                  # select the columns from the genotype matrix with clinical data
                  # for the current set of variants, current phenotype, and demographic covariates
                  currentVariantSet %>%
                  # select the column corresponding to each variant from the GenotypeMatrixAndClinicalData
                  lapply(
                    X = .,
                    FUN = function(currentVariant)
                    {
                      DataToReturn <-
                        GenotypeMatrixAndClinicalData %>%
                          select(
                            .data = .,
                            contains(match = currentVariant)
                            )

                      return(DataToReturn)
                    }
                      ) %>%
                    do.call(
                      what = "cbind",
                      args = .
                      ) %>%
                    cbind(
                      .,
                      select(
                        .data = GenotypeMatrixAndClinicalData,
                        "PhenotypeValue" = all_of(x = currentPhenotype),
                        # demographic information
                        Age_years,
                        Gender_MorF,
                        BMI_kgPerSqM,
                        GAF_admit,
                        Residential_visit,
                        Insurance,
                        PsychiatricHistory_family,
                        NeglectAndAbuseHistory,
                        SubstanceAbuseHistory,
                        cyp2d6_modulator,
                        #dopaminergic_modulator,
                        #serotonergic_modulator,
                        C1:C10
                        )
                      ) %>%
                      # remove observations with missing values
                      na.omit(object = .)
              
                  # specify the regression model formula with the significant variant genotypes: 
                    # Phenotype ~ variant genotypes + demographic covariates
                  RegressionFormula_noMDS <-
                    DataForRegression %>%
                    select(
                      .data = .,
                      -PhenotypeValue,
                      -C1,-C2,-C3,-C4,-C5,-C6,-C7,-C8,-C9,-C10
                      ) %>%
                    names(x = .) %>%
                    paste(
                      .,
                      collapse = " + "
                      ) %>%
                    paste(
                      "PhenotypeValue",
                      .,
                      sep = " ~ "
                    )
                  
                  # create a robust standard errors regression model of the phenotype value versus significant variant genotypes
                  # with adjustment for demographic/clinical covariates
                  RegressionModel_noMDS <-
                    DataForRegression %>%
                    lm_robust(
                      formula = as.formula(object = RegressionFormula_noMDS),
                      data = .
                      )

                  # perform an F-test for variance homogeneity to see if there is any heteroscedasticity with
                  # the lm() function; the lm_robust() function cannot be used in the ols_test_f() function
                  VarianceHomogeneityTestResult_noMDS <-
                    DataForRegression %>%
                    lm(
                      formula = as.formula(object = RegressionFormula_noMDS),
                      data = .
                    ) %>%
                    ols_test_f(
                      model = .
                      )

                  GenotypeRegressionResults_RSE_noMDS <-
                    # create a tidy summary of the linear robust standard errors regression
                    TidyLinearRobustStandardErrorsRegressionResults(
                      DataForRegression = DataForRegression,
                      RegressionModel = RegressionModel_noMDS
                      ) %>%
                     # add a label for the current phenotype and regression model formula
                    # and a label for whether or not MDS components were included
                    mutate(
                      .data = .,
                      Phenotype = currentPhenotype,
                      MDS_components = "MDS components NOT included",
                      Regression = RegressionFormula_noMDS
                      ) %>%
                    # add the P-value from the test for Heteroscedastiticy with astericks for levels of significance
                    mutate(
                      .data = .,
                      Heteroscedastiticy_p.value = VarianceHomogeneityTestResult_noMDS$p %>%
                                                   signif(x = .,digits = 3) %>%
                                                   AddAstericksToPvalues(columnVector = .)
                      ) %>%
                    # filter out the intercept term
                    filter(
                      .data = .,
                      term != "(Intercept)"
                      )
                  
                  # specify the same regression model formula with 10 MDS components included to adjust for population
                  # stratification
                  RegressionFormula_withMDScomponents <-
                    RegressionFormula_noMDS %>%
                    paste(
                      .,
                      "+ C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10"
                      )
                  
                  # create the same regression model with 10 MDS components included to adjust for population stratification
                  RegressionModel_withMDScomponents <-
                    DataForRegression %>%
                    lm_robust(
                      formula = as.formula(object = RegressionFormula_withMDScomponents),
                      data = .
                      )
                  
                  # perform an F-test for variance homogeneity to see if there is any heteroscedasticity with
                  # the lm() function; the lm_robust() function cannot be used in the ols_test_f() function
                  VarianceHomogeneityTestResult_withMDScomponents <-
                    DataForRegression %>%
                    lm(
                      formula = as.formula(object = RegressionFormula_withMDScomponents),
                      data = .
                    ) %>%
                    ols_test_f(
                      model = .
                      )
                  
                  GenotypeRegressionResults_RSE_withMDScomponents <-
                    # create a tidy summary of the linear robust standard errors regression
                    TidyLinearRobustStandardErrorsRegressionResults(
                      DataForRegression = DataForRegression,
                      RegressionModel = RegressionModel_withMDScomponents
                      ) %>%
                     # add a label for the current phenotype and regression model formula
                    # and a label for whether or not MDS components were included
                    mutate(
                      .data = .,
                      Phenotype = currentPhenotype,
                      MDS_components = "MDS components ARE INCLUDED",
                      Regression = RegressionFormula_withMDScomponents
                      ) %>%
                    # add the P-value from the test for Heteroscedastiticy with astericks for levels of significance
                    mutate(
                      .data = .,
                      Heteroscedastiticy_p.value = VarianceHomogeneityTestResult_withMDScomponents$p %>%
                                                   signif(x = .,digits = 3) %>%
                                                   AddAstericksToPvalues(columnVector = .)
                      ) %>%
                    # filter out the intercept term
                    filter(
                      .data = .,
                      term != "(Intercept)"
                      ) %>%
                    # filter out the MDS components terms
                    filter(
                      .data = .,
                      str_detect(
                        string = term,
                        pattern = "^C[:digit:]+$"
                        )==FALSE
                      )
                  
                  # bind the regression tables without and with adjustment for 10 MDS components together by row
                  # and return the result
                  RegressionResultsToReturn <-
                    rbind(
                      GenotypeRegressionResults_RSE_noMDS,
                      GenotypeRegressionResults_RSE_withMDScomponents
                    )
                  
                  return(RegressionResultsToReturn)

            }
            ) %>%
        do.call(
          what = "rbind",
          args = .
          ) %>%
        # label the gene for each variant
        mutate(
          .data = .,
          Gene = term %>%
                 lapply(
                   X = .,
                   FUN = function(currentTerm)
                   {
                     # if the current term is a variant term (the term starts with an "X"),
                     # identify the gene for the current variant, 
                     # otherwise return a missing value
                     if(startsWith(x = currentTerm,prefix = "X"))
                     {
                         # reformat the string of the current variant to match the
                         # format in the SignificantVariantsFromGWAS table
                         currentVariant <-
                           currentTerm %>%
                           gsub(pattern = "[.]",replacement = ":",x = .) %>%
                           str_remove(string = .,pattern = "_[:alpha:]*")
  
                         ValueToReturn <-
                            SignificantVariantsFromGWAS %>%
                            # extract the gene label corresponding to the current variant
                            filter(
                              .data = .,
                              str_detect(string = currentVariant,pattern = SNP)==TRUE
                              ) %>%
                            pull(
                              .data = .,
                              Gene.refGeneWithVer
                              )
                         
                     } else {
                       ValueToReturn <- as.character(x = NA)
                     }
  
                     return(ValueToReturn)
                   }
                   ) %>%
                  unlist(x = .)
          ) %>%
        # select useful columns
        # and update column names
        select(
          .data = .,
          Phenotype,
          MDS_components,
          "Term" = term,
          Gene,
          "Beta ± S.E." = estimate_std.error,
          "P-value" = p.value,
          "F-stat. P-value" = F_stat_p_value,
          R2,
          "R2-adj." = R2_adj,
          N,
          Regression
          ) 

# save the results to the file system
Multiple_GenotypeRegressionResults_RSE_withDemographicCovariates %>%
  write_delim(
    x = .,
    file = "./MultipleGenotypeRegressions/Multiple_GenotypeRegressionResults_RSE_withDemographicCovariates.txt",
    delim = "\t"
      )

##############

# add an "_OLS" suffix to the column names of the SignificantVariantsFromGWAS
# to indicate ordinary least squares regression
names(x = SignificantVariantsFromGWAS) <- SignificantVariantsFromGWAS %>%
                                          names(x = .) %>%
                                          sapply(
                                            X = .,
                                            FUN = function(currentColumnName)
                                            {
                                                ValueToReturn <-
                                                  currentColumnName %>%
                                                  paste0(.,"_OLS")

                                                return(ValueToReturn)
                                            },
                                            simplify = TRUE
                                              )

# create a table of ordinary least squares regression summary statistics with the
# SignificantVariantsFromGWAS data and combine it with the robust standard errors
# regression summary statistics
GenotypeRegressionResults_OLS_RSE <-
  SignificantVariantsFromGWAS %>%
  left_join(
    x = .,
    y = GenotypeRegressionResults_RSE,
    by = c("SNP_OLS" = "Variant_RSE")
    ) %>%
  # add astericks to the OLS p-values
  mutate(
    .data = .,
    P_OLS = P_OLS %>% AddAstericksToPvalues(columnVector = .)
  ) 

if(!dir.exists(paths = "./primaryGenotypeRegressions/"))
{
        dir.create(path = "./primaryGenotypeRegressions/")
}

# save the table of ordinary least squares and robust standard errors regression
# results to the file system
GenotypeRegressionResults_OLS_RSE %>%
  write_delim(
    x = .,
    file = "./primaryGenotypeRegressions/SignificantVariantsFromGWAS_RegressionResults_OLS_RSE.txt",
    delim = "\t"
    )

############## perform survival analysis of antipsychotic drug discontinuation
  # The antipsychotic drugs being considered here are:
      # (1) aripiprazole, (2) chlorpromazine, (3) haloperidol, (4) olanzapine
      # (5) paliperidone, (6) quetiapine, (7) risperidone, and (8) ziprasidone

  # create a directory for saving survival analysis results if it does not exist
  if(!dir.exists(paths = "./SurvivalAnalysisResults/"))
  {
    dir.create(path = "./SurvivalAnalysisResults/")
  }

  # create a survival analysis dataset
  SurvivalAnalysisData <-
    GenotypeMatrixAndClinicalData %>%
      select(
        .data = .,
        ParticipantID,
        # select columns with the days on each antipsychotic medication
        starts_with(match = "DaysOn"),
        # select the columns that indicate if the medication was taken after discharge (a censored observation) or
        # discontinued during the visit
        contains(match = "AfterDischarge"),
        # select the column for risperidone discontinuation during the visit or not (censored observation)
        Risperidone_discontinued,
        # demographic information
        Age_years,
        Gender_MorF,
        BMI_kgPerSqM,
        GAF_admit,
        Residential_visit,
        Insurance,
        PsychiatricHistory_family,
        NeglectAndAbuseHistory,
        SubstanceAbuseHistory,
        cyp2d6_modulator,
        #dopaminergic_modulator,
        #serotonergic_modulator,
        ReasonForRisperidoneDiscontinuation
        ) %>%
    # make sure the covariate columns are the appropriate data type
    mutate(
      .data = .,
      Age_years = Age_years %>% as.numeric(x = .),
      Gender_MorF = Gender_MorF %>% factor(x = .),
      BMI_kgPerSqM = BMI_kgPerSqM %>% as.numeric(x = .),
      GAF_admit = GAF_admit %>% as.numeric(x = .),
      Residential_visit = Residential_visit %>% factor(x = .),
      Insurance = Insurance %>% factor(x = .),
      PsychiatricHistory_family = PsychiatricHistory_family %>% factor(x = .),
      NeglectAndAbuseHistory = NeglectAndAbuseHistory %>% factor(x = .),
      SubstanceAbuseHistory = SubstanceAbuseHistory %>% factor(x = .),
      cyp2d6_modulator = cyp2d6_modulator %>% factor(x = .),
      #dopaminergic_modulator = dopaminergic_modulator %>% factor(x = .),
      #serotonergic_modulator = serotonergic_modulator %>% factor(x = .),
      ReasonForRisperidoneDiscontinuation = ReasonForRisperidoneDiscontinuation %>% factor(x = .)
      ) %>%
    # set days on the antipsychotic columns to numerics
    mutate_at(
      .tbl = .,
      .vars = vars(starts_with(match = "DaysOn")),
      .funs = function(currentColumn)
      {
        ValuesToReturn <-
          currentColumn %>%
          as.numeric(x = .)

        return(ValuesToReturn)
      }
        ) %>%
    # convert medication after discharge columns to medication discontinued columns
    # by changing yes to no and no to yes
    mutate_at(
      .tbl = .,
       vars(contains(match = "AfterDischarge")),
      .funs = function(currentColumn)
      {

          ValuesToReturn <-
            currentColumn %>%
            # loop through each value in the column
            lapply(
              X = .,
              FUN = function(currentValue)
              {
                  # if the value is not missing
                  if(!is.na(x = currentValue))
                  {
                     # if the current value is "y" for yes,
                     # convert it to a numeric 0
                     if(currentValue=="y")
                     {
                       ValueToReturn <- as.numeric(x = 0)
                     }

                    # if the current value is "n" for no,
                    # convert it to a numeric 1
                    if (currentValue=="n")
                    {
                       ValueToReturn <- as.numeric(x = 1)
                    }

                  } else {

                    # otherwise, leave the value as NA
                    ValueToReturn <- as.numeric(x = NA)

                  }

                  return(ValueToReturn)
              }
                ) %>%
             unlist(x = .) %>%
            # convert the array to a logical vector
            as.logical(x = .)

          return(ValuesToReturn)
      }
      ) %>%
    # change the names of the medication after discharge columns to medication_discontinued columns
    select(
      .data = .,
      "Aripiprazole_discontinued" = AripiprazoleAfterDischarge,
      "Chlorpromazine_discontinued" = ChlorpromazineAfterDischarge,
      "Haloperidol_discontinued" = HaloperidolAfterDischarge,
      "Olanzapine_discontinued" = OlanzapineAfterDischarge,
      "Paliperidone_discontinued" = PaliperidoneAfterDischarge,
      "Quetiapine_discontinued" = QuetiapineAfterDischarge,
      "Ziprasidone_discontinued" = ZiprasidoneAfterDischarge,
      everything()
      ) %>%
    # convert risperidone discontinued to a logical vector
    mutate(
      .data = .,
      Risperidone_discontinued = Risperidone_discontinued %>%
                                 as.character(x = .) %>%
                                 as.numeric(x = .) %>%
                                 as.logical(x = .)
      )
  
  # define an array of antipsychotic drugs
  AntipsychoticDrugs <-
    c(
      "Aripiprazole",
      "Chlorpromazine",
      "Haloperidol",
      "Olanzapine",
      "Paliperidone",
      "Quetiapine",
      "Risperidone",
      "Ziprasidone"
      )
  
  # loop through each medication and prepare a dataframe with medication 
  # discontinuation data as well as clinical/demographic covariates
  KaplanMeierSurvivalModelData <-
    AntipsychoticDrugs %>%
      lapply(
        X = .,
        FUN = function(currentMedication)
        {
            # define the columns to select based on the medication
            ColumnsToSelect <-
              c(
                "ParticipantID",
                paste0("DaysOn",currentMedication),
                paste0(currentMedication,"_discontinued")
              ) %>%
              c(
                .,
                "Age_years",
                "Gender_MorF",
                "BMI_kgPerSqM",
                "GAF_admit",
                "Residential_visit",
                "Insurance",
                "PsychiatricHistory_family",
                "NeglectAndAbuseHistory",
                "SubstanceAbuseHistory",
                "cyp2d6_modulator",
                #"dopaminergic_modulator",
                #"serotonergic_modulator",
                "ReasonForRisperidoneDiscontinuation"
              )
            
            DataToReturn <-
              ColumnsToSelect %>%
              # select the columns for the current medication and demographic information
              SurvivalAnalysisData[.] %>%
              # create a column with a label for the current medication
              mutate(
                .data = .,
                antipsychotic = currentMedication
                ) %>%
              # remove the specific medication from the column name
              select(
                .data = .,
                "DaysOnMedication" = paste0("DaysOn",currentMedication),
                "Medication_discontinued" = paste0(currentMedication,"_discontinued"),
                antipsychotic,
                everything()
                )
          
          return(DataToReturn)
        }
        ) %>%
      do.call(
        what = "rbind",
        args = .
        )
  
  # create a Kaplan-Meier survival model with a curve for each antipsychotic drug without any adjustment for covariates
  ggsurvObjectAllAntipsychotics <-
     ggsurvplot(
       fit = survfit(
                     formula = Surv(DaysOnMedication,Medication_discontinued) ~ antipsychotic,
                     data = KaplanMeierSurvivalModelData %>%
                            select(
                              .data = .,
                              DaysOnMedication,
                              Medication_discontinued,
                              antipsychotic
                              ) %>%
                            na.omit(object = .)
                     ),
       xlab = "Days on medication",
       ylab = "Predicted patients on medication (%)",
       legend.title = "",
       legend.labs = AntipsychoticDrugs,
       font.legend = c(11, "bold"),
       font.main = c(11, "bold"),
       title = "",
       surv.scale = "percent",
       legend = c(.8,.725),
       risk.table = "absolute"
       #cumevents = TRUE,
       #cumcensor = TRUE
       )

  # plot the Kaplan-Meier survival curve for each antipsychotic drug
  ggsurvObjectAllAntipsychotics$plot <-
    ggsurvObjectAllAntipsychotics$plot +
      theme(
        plot.title = element_blank(),
        # axis.text = element_text(family = "Arial",face = "bold",size = 11),
        # axis.text.x = element_text(family = "Arial",face = "bold",size = 11),
        # axis.text.y = element_text(family = "Arial",face = "bold",size = 11),
        # axis.title = element_text(family = "Arial",face = "bold",size = 11),
        # axis.title.x = element_text(family = "Arial",face = "bold",size = 11),
        # axis.title.y = element_text(family = "Arial",face = "bold",size = 11),
        legend.title = element_blank(),
        # legend.text = element_text(family = "Arial",face = "bold",size = 11),
        # text = element_text(family = "Arial",face = "bold",size = 11),
        legend.box.background = element_rect(color="black", size=1)
           )

  # save the plot to the file system
  tiff(
    filename = "./SurvivalAnalysisResults/KaplanMeierSurvivalPlot_allAntipsychotics_noCovariates.tiff",
    width = 9,
    height = 9,
    units = "in",
    compression = "none",
    res = 300
    )

  print(x = ggsurvObjectAllAntipsychotics)
  
  dev.off()
  
  # plot the Kaplan-Meier survival curve for risperidone only stratified by 
  # reason for risperidone discontinuation (lack of efficacy, side effect, or other)
  ggsurvObjectReasonForRisperidoneDiscontinuation <-
    ggsurvplot(
      fit = survfit(
                    formula = Surv(DaysOnMedication,Medication_discontinued) ~ ReasonForRisperidoneDiscontinuation,
                    data = KaplanMeierSurvivalModelData %>%
                           filter(
                             .data = .,
                             antipsychotic == "Risperidone"
                             ) %>%
                           select(
                             .data = .,
                             DaysOnMedication,
                             Medication_discontinued,
                             ReasonForRisperidoneDiscontinuation
                             ) %>%
                           na.omit(object = .)
                    ),
      xlab = "Days on medication",
      ylab = "Predicted patients on medication (%)",
      legend.title = "",
      legend.labs = c("lack of efficacy","other reason","side effects"),
      font.legend = c(11, "bold"),
      font.main = c(11, "bold"),
      title = "",
      surv.scale = "percent",
      legend = c(.8,.725),
      risk.table = "absolute"
      )
  
  ggsurvObjectReasonForRisperidoneDiscontinuation$plot <-
    ggsurvObjectReasonForRisperidoneDiscontinuation$plot +
    theme(
      plot.title = element_blank(),
      # axis.text = element_text(family = "Arial",face = "bold",size = 11),
      # axis.text.x = element_text(family = "Arial",face = "bold",size = 11),
      # axis.text.y = element_text(family = "Arial",face = "bold",size = 11),
      # axis.title = element_text(family = "Arial",face = "bold",size = 11),
      # axis.title.x = element_text(family = "Arial",face = "bold",size = 11),
      # axis.title.y = element_text(family = "Arial",face = "bold",size = 11),
      legend.title = element_blank(),
      # legend.text = element_text(family = "Arial",face = "bold",size = 11),
      # text = element_text(family = "Arial",face = "bold",size = 11),
      legend.box.background = element_rect(color="black", size=1)
         )
  
  # save the plot to the file system
  tiff(
    filename = "./SurvivalAnalysisResults/KaplanMeierSurvivalPlot_ReasonForRisperidoneDiscontinuation.tiff",
    width = 9,
    height = 9,
    units = "in",
    compression = "none",
    res = 300
    )
  
  print(x = ggsurvObjectReasonForRisperidoneDiscontinuation)
  
  dev.off()

  # save the Kaplan Meier model output for all antipsychotic drugs with no covariates to the file system
  KaplanMeierSurvivalModelData %>%
  select(
    .data = .,
    DaysOnMedication,
    Medication_discontinued,
    antipsychotic
    ) %>%
  na.omit(object = .) %>%
  survfit(
    formula = Surv(DaysOnMedication,Medication_discontinued) ~ antipsychotic,
    data = .
  ) %>%
  # capture the survival model output
  capture.output(.) %>%
  # convert to a text file and read the text into R
  textConnection(object = .) %>%
  read.table(
    file = .,
    skip = 2,
    header = TRUE
    ) %>%
  # create a column for the antipsychotic
  mutate(
    .data = .,
    antipsychotic = rownames(x = .)
    ) %>%
  # remove the rownames
  data.frame(.,row.names = NULL) %>%
  # save to the file system
  write_delim(
    x = .,
    file = "./SurvivalAnalysisResults/KaplanMeierSurvivalTable_allAntipsychotics_noCovariates.txt",
    delim = "\t"
      )
  
  # save the Kaplan Meier model output the reason for risperidone discontinuation to the file system
  KaplanMeierSurvivalModelData %>%
    filter(
      .data = .,
      antipsychotic == "Risperidone"
    ) %>%
    select(
      .data = .,
      DaysOnMedication,
      Medication_discontinued,
      ReasonForRisperidoneDiscontinuation
    ) %>%
    na.omit(object = .) %>%
    survfit(
      formula = Surv(DaysOnMedication,Medication_discontinued) ~ ReasonForRisperidoneDiscontinuation,
      data = .
    ) %>%
    # capture the survival model output
    capture.output(.) %>%
    # convert to a text file and read the text into R
    textConnection(object = .) %>%
    read.table(
      file = .,
      skip = 2,
      header = TRUE
    ) %>%
    # create a column for the antipsychotic
    mutate(
      .data = .,
      ReasonForDiscontinuation = rownames(x = .)
    ) %>%
    # remove the rownames
    data.frame(.,row.names = NULL) %>%
    # save to the file system
    write_delim(
      x = .,
      file = "./SurvivalAnalysisResults/KaplanMeierSurvivalTable_ReasonForRispDiscontinuation.txt",
      delim = "\t"
    )
    
  
  # create a table for labeling the plot legends of kaplan-meier log rank tests
  KaplanMeierLogRankTestLegendLableTable <-
    tibble(
      "BinaryCovariate" = c(
                            "Gender_MorF",
                            "Residential_visit",
                            "Insurance",
                            "PsychiatricHistory_family",
                            "NeglectAndAbuseHistory",
                            "SubstanceAbuseHistory",
                            "cyp2d6_modulator",
                            #"dopaminergic_modulator",
                            #"serotonergic_modulator",
                            SignificantVariantsFromGWAS %>%
                            filter(
                              .data = .,
                              Phenotype_OLS=="DaysOnRisperidone"
                              ) %>%
                            pull(
                              .data = .,
                              SNP_OLS
                              ) %>%
                            paste0("VariantGenotype_",.)
                            ),
      "LegendTitle" = c(
                        "Gender",
                        "Visit",
                        "Insurance",
                        "Family psych. hist.",
                        "Neglect and abuse hist.",
                        "Substance abuse hist.",
                        "CYP2D6 modulator",
                        #"dopaminergic modulator",
                        #"serotonergic modulator",
                        "7:158034920:T"
                        # SignificantVariantsFromGWAS %>%
                        #   filter(
                        #     .data = .,
                        #     Phenotype_OLS=="DaysOnRisperidone"
                        #   ) %>%
                        #   pull(
                        #     .data = .,
                        #     SNP_OLS
                        #   )
                       ),
      "LegendLabels" = list(
                            c("Female","Male"),
                            c("acute","residential"),
                            c("medicaid","private"),
                            c("no","yes"),
                            c("no","yes"),
                            c("no","yes"),
                            c("no","yes"),
                            #c("no","yes"),
                            #c("no","yes"),
                            c("0 copy","1 copy")
                            )
    )
  
    # perform the Kaplan Meier log-rank test to see if there are differences in antipsychotic
    # medication discontinuation for binary covariates 
    KaplanMeierLogRankTestResults <-
        KaplanMeierSurvivalModelData %>%
        filter(
          .data = .,
          antipsychotic=="Risperidone"
          ) %>%
          group_by(
            .data = .,
            antipsychotic
            ) %>%
          group_split(
            .tbl = .
            ) %>%
          # loop through the data for each antipsychotic medication
          lapply(
            X = .,
            FUN = function(currentMedicationData)
            {
  
              # If the current medication is risperidone, include the variant that was
              # significantly associated with days on risperidone in the analysis
              if(unique(x = currentMedicationData$antipsychotic)=="Risperidone")
              {
  
                  VariantAssociatedWithDaysOnRisperidone <-
                      SignificantVariantsFromGWAS %>%
                      filter(
                        .data = .,
                        Phenotype_OLS=="DaysOnRisperidone"
                        ) %>%
                      pull(
                        .data = .,
                        SNP_OLS
                        )
                  
                  # define the binary covariates to loop through:
                    #gender, residential_visit, insurance, family psychiatric history, neglect and abuse history, substance abuse history,
                    #cyp2d6 modulator, dopaminergic modulator, serotonergic modulator, and the variant associated with days on risperidone
                  BinaryCovariatesToLoopThrough <-
                    c(
                    "Gender_MorF",
                    "Residential_visit",
                    "Insurance",
                    "PsychiatricHistory_family",
                    "NeglectAndAbuseHistory",
                    "SubstanceAbuseHistory",
                    "cyp2d6_modulator"#,
                    #"dopaminergic_modulator",
                    #"serotonergic_modulator"
                    ) %>%
                    c(
                      .,
                      VariantAssociatedWithDaysOnRisperidone
                     )
  
              } else {
  
                # define the binary covariates to loop through
                # (cyp2d6, dopaminergic, and serotonergic modulators are not included in this
                # because the majority of antipsychotics fit these categories)
                BinaryCovariatesToLoopThrough <-
                       c(
                       "Gender_MorF",
                       "Residential_visit",
                       "Insurance",
                       "PsychiatricHistory_family",
                       "NeglectAndAbuseHistory",
                       "SubstanceAbuseHistory"
                       )

                VariantAssociatedWithDaysOnRisperidone <- FALSE
  
              }
  
               ResultsToReturn <-
                   # loop through each binary covariate
                   BinaryCovariatesToLoopThrough %>%
                     lapply(
                       X = .,
                       FUN = function(currentBinaryCovariate)
                       {
                             # if the current binary covariate is the variant associated with days
                             # on risperidone, add a column with the variant genotype to the data for the
                             # log rank test
                             if(currentBinaryCovariate==VariantAssociatedWithDaysOnRisperidone)
                             {
  
                               # update the currentBinaryCovariate to match the VariantGenotype column name
                               # in the code below
                               currentBinaryCovariate <- paste0("VariantGenotype_",VariantAssociatedWithDaysOnRisperidone)
  
                               DataForLogRankTest <-
                                 currentMedicationData %>%
                                 # select days on medication, medication discontinued, and
                                 # variant genotype columns
                                 select(
                                   .data = .,
                                   ParticipantID,
                                   DaysOnMedication,
                                   Medication_discontinued
                                   ) %>%
                                 left_join(
                                   x = .,
                                   y = GenotypeMatrixAndClinicalData %>%
                                       select(
                                         .data = .,
                                         ParticipantID,
                                         "VariantGenotype" = contains(
                                                                      match = VariantAssociatedWithDaysOnRisperidone %>%
                                                                              gsub(
                                                                                pattern = ":",
                                                                                replacement = ".",
                                                                                x = .
                                                                                )
                                                                       )
                                        ),
                                   by = "ParticipantID"
                                   ) %>%
                                 # drop the participant ID column
                                 select(
                                   .data = .,
                                   -ParticipantID
                                   ) %>%
                                 # remove missing values
                                 na.omit(object = .) %>%
                                 # make sure the genotype column is heterozygous dominant
                                 mutate(
                                   .data = .,
                                   VariantGenotype = if_else(
                                                             condition = VariantGenotype>=1,
                                                             true = as.integer(x = 1),
                                                             false = as.integer(x = 0)
                                                               )
                                   ) %>%
                                 # rename the variant genotype column to the currentBinaryCovariate
                                 select(
                                   .data = .,
                                   "Covariate" = VariantGenotype,
                                   everything()
                                   )
                               
                             } else {
  
                               DataForLogRankTest <-
                                 currentMedicationData %>%
                                   # select days on medication, medication discontinued, and binary covariate columns
                                   select(
                                     .data = .,
                                     DaysOnMedication,
                                     Medication_discontinued,
                                     "Covariate" = all_of(x = currentBinaryCovariate)
                                   ) %>%
                                   # remove missing values
                                   na.omit(object = .)
                             }
  
                         # identify the number of covariate categories
                         NumberOfCovariateCategories <-
                           DataForLogRankTest$Covariate %>%
                           unique(x = .) %>%
                           length(x = .)
                         
                         # identify the number of outcome categories
                         NumberOfOutcomeCategories <-
                           DataForLogRankTest$Medication_discontinued %>%
                           unique(x = .) %>%
                           length(x = .)
                         
                         # if there is more than one category of the covariate, 
                         # and there is more than one category for the outcome,
                         # perform the log-rank test
                         if(NumberOfCovariateCategories>1 & NumberOfOutcomeCategories>1)
                         {
                              LogRankTestResult <-
                                 DataForLogRankTest %>%
                                 # define the survival model with the current medication data and covariate
                                 survdiff(
                                   formula = Surv(DaysOnMedication,Medication_discontinued) ~ Covariate,
                                   data = .
                                     )
  
                               # extract the chi-square P-value from the survdiff object
                              LogRankPvalue <- 1 - pchisq(
                                                          LogRankTestResult$chisq,
                                                          length(x = LogRankTestResult$n) - 1
                                                          )
                              
                              # create a results table to return
                              ResultsToReturn <-
                                  DataForLogRankTest %>%
                                  survfit(
                                    formula = Surv(DaysOnMedication,Medication_discontinued) ~ Covariate,
                                    data = .
                                      ) %>%
                                  # capture the survival model output
                                  capture.output(.) %>%
                                  # convert to a text file and read the text into R
                                  textConnection(object = .) %>%
                                  read.table(
                                    file = .,
                                    skip = 2,
                                    header = TRUE
                                  ) %>%
                                  # create a column for the antipsychotic
                                  mutate(
                                    .data = .,
                                    covariateCategory = rownames(x = .)
                                  ) %>%
                                  # remove the rownames
                                  data.frame(.,row.names = NULL) %>%
                                  mutate(
                                    .data = .,
                                    Pvalue_logRankTest = LogRankPvalue %>%
                                                         signif(x = .,digits = 3) %>%
                                                         AddAstericksToPvalues(columnVector = .),
                                    covariate = currentBinaryCovariate,
                                    antipsychotic = currentMedicationData$antipsychotic %>%
                                                    unique(x = .)
                                    )
                              
                                  # create a Kaplan-Meier survival plot of the current antipsychotic stratified 
                                  # by the current binary covariate
                          
                                  # identify the current antipsychotic medication
                                  currentMedication <-
                                    currentMedicationData$antipsychotic %>%
                                    unique(x = .)
                                  
                                  # plot the Kaplan-Meier curve of antipsychotic discontinuation stratified
                                  # by the binary covariate
                                  StratifiedKaplanMeierPlotObject <-
                                    ggsurvplot(
                                      fit = survfit(
                                                    formula = Surv(DaysOnMedication,Medication_discontinued) ~ Covariate,
                                                    data = DataForLogRankTest
                                                    ),
                                      data = DataForLogRankTest,
                                      xlab = glue("Days on {currentMedication}"),
                                      ylab = glue("Predicted patients on {currentMedication} (%)"),
                                      font.legend = c(11,"bold"),
                                      font.main = c(11,"bold"),
                                      legend.title = KaplanMeierLogRankTestLegendLableTable %>%
                                                     filter(
                                                       .data = .,
                                                       BinaryCovariate == currentBinaryCovariate
                                                       ) %>%
                                                      pull(
                                                        .data = .,
                                                        LegendTitle
                                                        ),
                                      legend.labs = KaplanMeierLogRankTestLegendLableTable %>%
                                                    filter(
                                                      .data = .,
                                                      BinaryCovariate == currentBinaryCovariate
                                                      ) %>%
                                                    pull(
                                                      .data = .,
                                                      LegendLabels
                                                      ) %>%
                                                    unlist(x = .),
                                      title = "",
                                      surv.scale = "percent",
                                      legend = c(0.8,0.85)#,
                                      #pval = TRUE
                                      )
                                    
                                  StratifiedKaplanMeierPlot <-
                                    StratifiedKaplanMeierPlotObject$plot +
                                    theme(
                                      plot.title = element_blank(),
                                      axis.text = element_text(family = "Arial",face = "bold",size = 11),
                                      axis.text.x = element_text(family = "Arial",face = "bold",size = 11),
                                      axis.text.y = element_text(family = "Arial",face = "bold",size = 11),
                                      axis.title = element_text(family = "Arial",face = "bold",size = 11),
                                      axis.title.x = element_text(family = "Arial",face = "bold",size = 11),
                                      axis.title.y = element_text(family = "Arial",face = "bold",size = 11),
                                      legend.title = element_text(family = "Arial",face = "bold",size = 11,hjust = 0.5),
                                      legend.text = element_text(family = "Arial",face = "bold",size = 11),
                                      text = element_text(family = "Arial",face = "bold",size = 11),
                                      legend.box.background = element_rect(color="black", size=1)
                                    )
                                  
                                    tiff(
                                      filename = glue("./SurvivalAnalysisResults/{currentMedication}_{currentBinaryCovariate}_KaplanMeierLogRankTest.tiff"),
                                      width = 7,
                                      height = 5,
                                      units = "in",
                                      compression = "none",
                                      res = 300
                                    )
                                  
                                    print(x = StratifiedKaplanMeierPlot)
                                    
                                    dev.off()
  
                                  return(ResultsToReturn)
                         }
  
                           }) %>%
                           do.call(
                             what = "rbind",
                             args = .
                             )
  
  
  
  
               return(ResultsToReturn)
            }
              ) %>%
          do.call(
            what = "rbind",
            args = .
            )
    
    # save the Kaplan Meier log-rank test results to the file system
    KaplanMeierLogRankTestResults %>%
      write_delim(
        x = .,
        file = "./SurvivalAnalysisResults/KaplanMeierLogRankTestResults_allAntipsychotics.txt",
        delim = "\t"
      )
    
    # temporarily suppress non-problematic warning messages in the cox proportional hazard regressions
    options(warn = -1)
    
    #### perform Cox proportional hazards regression analysis of discontinuation of each
    # antipsychotic versus individual and multiple demographic covariates
  CoxRegressionResults <-
      KaplanMeierSurvivalModelData %>%
        filter(
          .data = .,
          antipsychotic == "Risperidone"
          ) %>%
        group_by(
          .data = .,
          antipsychotic
          ) %>%
        group_split(
          .tbl = .
          ) %>%
        # loop through the data for each antipsychotic medication
        lapply(
          X = .,
          FUN = function(currentMedicationData)
          {
            
            # define the variant that was significantly associated
            # with days on risperidone from the GWAS analysis
            VariantSignificantlyAssociatedWithDaysOnRisperidone <-
               SignificantVariantsFromGWAS %>%
               filter(
                 .data = .,
                 Phenotype_OLS == "DaysOnRisperidone"
                 ) %>%
               pull(
                 .data = .,
                 SNP_OLS
                 )
            
            # if the current medication is risperidone, add the variant that was significantly associated with
            # days on risperidone to the covariates to loop through
            if(unique(x = currentMedicationData$antipsychotic)=="Risperidone")
            {
              CovariatesToLoopThrough <-
              # define covariate terms to loop through
                c(
                  "Gender_MorF",
                  "Age_years",
                  "BMI_kgPerSqM",
                  "GAF_admit",
                  "Residential_visit",
                  "Insurance",
                  "PsychiatricHistory_family",
                  "NeglectAndAbuseHistory",
                  "SubstanceAbuseHistory",
                  "cyp2d6_modulator"#,
                  #"dopaminergic_modulator",
                  #"serotonergic_modulator"
                 ) %>%
                c(
                  .,
                  VariantSignificantlyAssociatedWithDaysOnRisperidone
                  )
              
            } else {
              
              CovariatesToLoopThrough <-
                # define covariate terms to loop through (cyp2d6, dopamine, and serotonin modulators are not
                # included in this since nearly all antipsychotic drugs fit these categories)
                c(
                  "Gender_MorF",
                  "Age_years",
                  "BMI_kgPerSqM",
                  "GAF_admit",
                  "Residential_visit",
                  "Insurance",
                  "PsychiatricHistory_family",
                  "NeglectAndAbuseHistory",
                  "SubstanceAbuseHistory"
                )
              
            }
                
            primaryCoxRegressionModelResults <-
                # loop through each covariate and include in a primary 
                # cox-proportional hazard regression model
                CovariatesToLoopThrough %>%
                  lapply(
                    X = .,
                    FUN = function(currentCovariate)
                    {
                        
                        # if the current covariate is the variant that was significantly associated
                        # with days on risperidone from the GWAS analysis,
                        # add the variant genotype to the data for the regression
                        if(currentCovariate==VariantSignificantlyAssociatedWithDaysOnRisperidone)
                        {
                            DataForRegression <-
                              currentMedicationData %>%
                              select(
                                .data = .,
                                ParticipantID,
                                DaysOnMedication,
                                Medication_discontinued
                                ) %>%
                              # join the survival analysis data
                              # with the variant genotype that was associated with days on risperidone
                              left_join(
                                x = .,
                                y = GenotypeMatrixAndClinicalData %>%
                                    select(
                                      .data = .,
                                      ParticipantID,
                                      "Covariate" = contains(
                                                             match = VariantSignificantlyAssociatedWithDaysOnRisperidone %>%
                                                              gsub(
                                                                pattern = ":",
                                                                replacement = ".",
                                                                x = .
                                                                )
                                                              )
                                    ),
                                by = "ParticipantID"
                                ) %>%
                              # make sure the genotype column is a numeric
                              mutate(
                                .data = .,
                                Covariate = Covariate %>% as.numeric(x = .)
                                ) %>%
                              # make sure the genotype column is heterozygous dominant
                              mutate(
                                .data = .,
                                Covariate = if_else(
                                                    condition = Covariate>=1,
                                                    true = 1,
                                                    false = 0
                                                      )
                              ) %>%
                              # drop the participant ID column
                              select(
                                .data = .,
                                -ParticipantID
                                 ) %>%
                              na.omit(object = .)
                                  
                        } else {
                          
                            DataForRegression <-
                              currentMedicationData %>%
                              # select medication discontinued, days on medication, and
                              # covariate columns
                              select(
                                .data = .,
                                DaysOnMedication,
                                Medication_discontinued,
                                "Covariate" = all_of(x = currentCovariate)
                                ) %>%
                              # remove missing values
                              na.omit(object = .)
                          
                        }
                        
                        CoxRegressionModelResults <-
                            # define the Cox proportional hazard model
                            coxph(
                              formula = Surv(DaysOnMedication,Medication_discontinued) ~ Covariate,
                              data = DataForRegression
                                ) %>%
                            # tidy the regression output and exponentiate the regression estimate
                              # to obtain a hazard ratio
                            broom::tidy(
                              x = .,
                              exp = TRUE
                              ) %>%
                            # add a column for the antipsychotic and the covariate
                           mutate(
                             .data = .,
                             medication = currentMedicationData$antipsychotic %>%
                                          unique(x = .),
                             covariate = currentCovariate,
                             # add astericks to p-values
                             p.value = p.value %>%
                                       signif(x = .,digits = 3) %>%
                                       AddAstericksToPvalues(columnVector = .),
                             # paste estimate and standard error together
                             estimate_std.error = paste0(
                                                         round(x = estimate,digits = 3),
                                                         "±",
                                                         round(x = std.error,digits = 3)
                                                         ),
                             # add the interpretation of the hazard ratio:
                             # if the ratio is > 1, increased hazard
                             # if the ratio is < 1, decreased hazard
                             hazard = if_else(
                                              condition = estimate>1,
                                              true = "increased",
                                              false = if_else(
                                                              condition = estimate<1,
                                                              true = "decreased",
                                                              false = "neutral"
                                                                )
                                                )
                             )
                        
                        return(CoxRegressionModelResults)
                    }
                      ) %>%
                  do.call(
                    what = "rbind",
                    args = .
                  ) %>%
              # add a column for the regression type 
              mutate(
                .data = .,
                RegressionType = "primary"
                )
            
            # if the current antipsychotic medication is risperidone,
            # join the outcome and covariate data for the regression with 
            # the variant genotype that was significantly associated with days on risperidone
            # and perform the Cox proportional hazards regression with all demographic covariates
            # and the variant genotype included,
            # otherwise, perform the regression with demographic covariates only
            if(unique(x = currentMedicationData$antipsychotic)=="Risperidone")
            {
                  DataForMultipleRegression <-
                    c(
                      "ParticipantID",
                      "DaysOnMedication",
                      "Medication_discontinued",
                      "Gender_MorF",
                      "Age_years",
                      "BMI_kgPerSqM",
                      "GAF_admit",
                      "Residential_visit",
                      "Insurance",
                      "PsychiatricHistory_family",
                      "NeglectAndAbuseHistory",
                      "SubstanceAbuseHistory",
                      "cyp2d6_modulator"#,
                      #"dopaminergic_modulator",
                      #"serotonergic_modulator"
                    ) %>%
                    # select columns with demographic covariates, the outcome, and days on the medication 
                    currentMedicationData[.] %>%
                    # join with the variant genotype data
                    left_join(
                      x = .,
                      y = GenotypeMatrixAndClinicalData %>% 
                          select(
                            .data = .,
                            ParticipantID,
                            "VariantGenotype" = contains(
                                                         match = VariantSignificantlyAssociatedWithDaysOnRisperidone %>% 
                                                         gsub(
                                                           pattern = ":",
                                                           replacement = ".",
                                                           x = .
                                                           )
                                     )
                            ),
                      by = "ParticipantID"
                      ) %>%
                    # remove missing values
                    na.omit(object = .) %>%
                    # make sure the variant genotype column is heterozygous dominant
                    mutate(
                      .data = .,
                      VariantGenotype = if_else(
                                                condition = VariantGenotype>=1,
                                                true = as.integer(x = 1),
                                                false = as.integer(x = 0)
                                                  )
                      ) %>%
                    # remove the participant ID column
                    select(
                      .data = .,
                      -ParticipantID
                      )
                  
                  # define the multiple model with demographic covariates and the variant genotype
                  MultipleRegressionResults <-
                      coxph(
                        formula = Surv(DaysOnMedication,Medication_discontinued)~Gender_MorF+Age_years+BMI_kgPerSqM+GAF_admit+Residential_visit+Insurance+PsychiatricHistory_family+NeglectAndAbuseHistory+SubstanceAbuseHistory+cyp2d6_modulator+VariantGenotype,
                        data = DataForMultipleRegression
                          )
              
            } else {
              
                  # select columns with the outcome, days on the medication, and demographic covariates
                  DataForMultipleRegression <-
                      CovariatesToLoopThrough %>%
                      c(
                        .,
                        "DaysOnMedication",
                        "Medication_discontinued"
                        ) %>%
                      # remove missing values
                      currentMedicationData[.] %>%
                      na.omit(object = .)
                  
                  # define the cox regression model with demographic covariates
                  MultipleRegressionResults <-
                      coxph(
                        formula = Surv(DaysOnMedication,Medication_discontinued)~Gender_MorF+Age_years+BMI_kgPerSqM+GAF_admit+Residential_visit+Insurance+PsychiatricHistory_family+NeglectAndAbuseHistory+SubstanceAbuseHistory,
                        data = DataForMultipleRegression
                          )
            }
            
              
                MultipleRegressionResults <-
                  MultipleRegressionResults %>%
                   # tidy the regression output and exponentiate the regression estimate
                     # to obtain a hazard ratio
                   broom::tidy(
                     x = .,
                     exp = TRUE
                     ) %>%
                     # add a column for the antipsychotic and the covariate
                    mutate(
                      .data = .,
                      medication = currentMedicationData$antipsychotic %>%
                                   unique(x = .),
                      covariate = "all covariates",
                      # add astericks to p-values
                      p.value = p.value %>%
                                signif(x = .,digits = 3) %>%
                                AddAstericksToPvalues(columnVector = .),
                      # paste estimate and standard error together
                      estimate_std.error = paste0(
                                                  round(x = estimate,digits = 3),
                                                  "±",
                                                  round(x = std.error,digits = 3)
                                                  ),
                      # add the interpretation of the hazard ratio:
                      # if the ratio is > 1, increased hazard
                      # if the ratio is < 1, decreased hazard
                      hazard = if_else(
                                       condition = estimate>1,
                                       true = "increased",
                                       false = if_else(
                                                       condition = estimate<1,
                                                       true = "decreased",
                                                       false = "neutral"
                                                         )
                                         )
                      ) %>%
                     # add a column for the regression type
                     mutate(
                       .data = .,
                       RegressionType = "Multiple"
                       )
                
                 # combine the primary and multiple regression results tables
                 # and return the results
                 ResultsToReturn <-
                   rbind(
                     primaryCoxRegressionModelResults,
                     MultipleRegressionResults
                   )
                 
                 return(ResultsToReturn)
          }
          ) %>%
        do.call(
          what = "rbind",
          args = .
          ) %>%
        # select useful columns
        # and update column names
        select(
          .data = .,
          medication,
          covariate,
          "Term" = term,
          "Hazard ratio ± S.E." = estimate_std.error,
          hazard,
          "P-value" = p.value,
          RegressionType
          )
  
   # turn warning messages back on
   options(warn = 0)
  
   # save the results to the file system
   CoxRegressionResults %>%
     write_delim(
       x = .,
       file = "./SurvivalAnalysisResults/CoxProportionalHazardRegressions.txt",
       delim = "\t"
         )
  
