# Risperidone response genome-wide association study (GWAS) pipeline

**Author:** Jack Staples

**Overview:**

This is a repository for the code for a genome-wide association study (GWAS) of how genomic variation and other clinical/demographic factors [age, body mass index (BMI), sex, insurance, family psychiatric history, history of neglect and abuse, substance abuse history, CYP2D6 modulator use, and treatment program (acute or residential)] influence risperidone treatment outcomes (duration of risperidone treatment, maximum risperidone dose prescribed, risperidone discontinuation, frequency of patient readmittance after the initial hospital stay, time to readmittance after discharge, frequency of readmittance while on risperidone, and duration of hospital stay). Although GAF score is a measure of the severity of psychiatric illness, and not a measure of risperidone response, we included GAF score at admittance and discharge and change in GAF between the two as endpoints of this study.

This pipeline consists of a combination of unix bioinformatics command line tools and R scripts wrapped in bash scripts so the pipeline can be run on a bash shell terminal.

**Required software:**

-   ANNOVAR version date 06/08/2020 (<https://annovar.openbioinformatics.org/en/latest/>)

-   BCFtools version 1.11 from htslib version 1.11 (<http://www.htslib.org/download/>)

-   bgzip version 1.11 (included with htslib)

-   Ensembl variant effect predictor (VEP) web interface (<http://grch37.ensembl.org/Homo_sapiens/Tools/VEP>)

-   PGx-POP version 1.0 (<https://github.com/PharmGKB/PGxPOP>)

-   PLINK version 1.90 (<https://www.cog-genomics.org/plink/>)

-   Pyenv version 1.2.23 for switching between python version 2 and 3 (<https://github.com/pyenv/pyenv>)

-   pypgx version 0.16.0 (<https://github.com/sbslee/pypgx>)

-   Python version 2.7.18 and 3.8.2

-   R version 4.0.4

-   stargazer version 1.0.8 (<https://stargazer.gs.washington.edu/stargazerweb/index.html>)

-   tabix (included with htslib)

-   Trans-omics for precision medicine (TOPMed) imputation server (<https://imputation.biodatacatalyst.nhlbi.nih.gov/#!>): performs haplotype phasing with Eagle version 2.4 and Minimac4 version 1.6.6

-   VCFtools version 0.1.16 (<https://vcftools.github.io/index.html>)

**Required R packages:**

-   aod version 1.3.1

-   estimatr version 0.30.2

-   factoextra version 1.0.7

-   GGally version 2.1.2

-   ggpubr version 0.4.0

-   glue version 1.4.2

-   olsrr version 0.5.3

-   parallel version 4.0.4

-   qqman version 0.1.8

-   readxl version 1.3.1

-   Rtsne version 0.16

-   survival version 3.2.7

-   survminer version 0.4.9

-   tidyverse version 1.3.0

**The input files for this pipeline are:**

-   PLINK format .ped and .map files of genotype data (GSA2016_125_025.ped) and participant IDs (GSA2016_125_025.map) from microarray genotyping

-   Variant call format (VCF) (ALL.2of4intersection.20100804.genotypes.vcf.gz) and participant information file (20100804.ALL.panel) from phase 3 of the 1000 Genomes project downloaded from the 1000 genomes ftp (<ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz>)

-   Genomic regions with high linkage disequilibrium to exclude from linkage disequilibrium pruning calculations in quality control of microarray genotyping data (inversion.txt) (downloaded from Andries Marees GWAS tutorial github page: <https://github.com/MareesAT/GWA_tutorial>)

-   Tables of GRCh37 (NC_000007.13) reference assembly variant rs-identification numbers and star alleles for all pharmacogenes downloaded from the Pharmacogene Variation Consortium website (<https://www.pharmvar.org/genes>)

-   Clinical outcomes (phenotypes), additional demographic/clinical information, medication data, and HCPCS diagnostic code data from study participant electronic medical records (EHR_phenotypes.xlsx)

**Pipeline scripts and descriptions:**

-   **shodair_genetic_analysis_pipeline_master.R** lists the order that pipeline scripts need to be run in order to complete the analysis

    This runs the following scripts:

    -   **shodair_quality_control.sh** performs quality control specific for GWASs with microarray genotyping data with code adapted from Andries Marees GWAS tutorial github page: <https://github.com/MareesAT/GWA_tutorial>).

        **Corresponding paper citation:** Marees, A. T., *et al*., *Int J Methods Psychiatr Res.*, *27*(2), e1608, 2018. Available from: <https://doi.org/10.1002/mpr.1608>

    -   **prepare_VCF_for_imputation.sh** prepares VCF files (creates a separate vcf.gz file for each chromosome, sorts variants by GRCh37 genomic position, and identifies variants with strand flip issues) for phasing and imputation with the TOPMed imputation server.

    -   **convert_imputed_files_to_VCF.sh** combines all files that are output from the TOPMed imputation server into a single VCF file and filters variants based on commonly used post imputation quality control thresholds for GWAS analysis.

    -   **ANNOVAR_annotation.sh** annotates variant calling data with ANNOVAR, the Ensembl variant effect predictor, clinical variant annotation data from the pharmGKB, and star allele haplotypes from pharmVAR. This also performs *in-silico* predictions of deleterious variants using the ADME optimized prediction framework.

    -   **call_star_alleles.sh** assigns star allele haplotypes with starGazer, PGx-POP, and pypgx software.

    -   **association_analysis.sh** performs all the genotype-phenotype statistical association analysis, which includes:

        -   Computes basic summary statistics of risperidone treatment outcomes and clinical/demographic data.

        -   Linear regressions to test individual genomic variants and other clinical/demographic factors for associations with continuous risperidone treatment outcomes.

        -   Logistic regressions to test individual genomic variants and other clinical/demographic factors for associations with dichotomous risperidone treatment outcomes.

        -   Kaplan-Meier survival and Cox proportional hazard regression analysis of risperidone discontinuation.

        -   Linkage disequilibrium pruning of association results.

        -   Boxplots of continous phenotypes stratified by genotype for significant variants.

        -   Manhattan and QQ plots of association results for each variant tested.

    -   **load_necessary_packages.R** is run at the beginning of every R script for loading necessary R packages.
