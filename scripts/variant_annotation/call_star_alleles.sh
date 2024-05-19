# this is a script for calling star alleles on the imputed variant call data

# make a working directory for using stargazer if it does not exist
if [[ -d ~/Psychiatric_GWAS/Psychiatric_GWAS/starGazer ]]
then
  echo "stargazer haplotype assignment directory exists"
else
  echo "Created working directory for star allele haplotype assignment with stargazer"
  mkdir ~/Psychiatric_GWAS/Psychiatric_GWAS/starGazer/
fi

# working directory for this is:
cd ~/Psychiatric_GWAS/Psychiatric_GWAS/starGazer/

# make sure that python version 3 is running (this is required for starGazer)
# add the pyenv executable to the PATH
export PATH="$HOME/.pyenv/bin:$PATH"
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"

# set the current version of python as python 3
pyenv global 3.8.2

#### running stargazer with the imputed VCF file 
# run stagazer genotyping tool for all possible genes with the imputed vcf with the SNP-array generated data option
for currentGene in {cacna1s,cftr,cyp1a1,cyp1a2,cyp1b1,cyp2a6,cyp2a13,cyp2b6,cyp2c8,cyp2c9,cyp2c19,cyp2d6,cyp2e1,cyp2f1,cyp2j2,cyp2r1,cyp2s1,cyp2w1,cyp3a4,cyp3a5,cyp3a7,cyp3a43,cyp4b1,cyp26a1,cyp4f2,cyp19a1,dpyd,g6pd,gstm1,gstp1,gstt1,ifnl3,nat1,nat2,nudt15,por,ryr1,slc15a2,slc22a2,slco1b1,slco1b3,slco2b1,sult1a1,tbxas1,tpmt,ugt1a1,ugt1a4,ugt2b7,ugt2b15,ugt2b17,vkorc1}
do
  python ~/bioinformatics_software/Stargazer_v1.0.8/stargazer.py genotype -o "${currentGene}_results" -d chip -t "${currentGene}" --vcf "../TOPMEDimputationResults/imputed_variants_allCHROMs_INFO_0.3_MAF_0.05_maxMissing20.dose.recode.vcf"
done

# set the current version of python as the masOS system version (python 2)
pyenv global system

# plot (make haplotype and genotype-predicted phenotype bar charts and tables for each gene) with the stargazer star allele calling results
Rscript --no-save ../scripts/variant_annotation/summarize_stargazer_results.R

# make a working directory for using pgx-pop if it does not exist
if [[ -d ~/Psychiatric_GWAS/Psychiatric_GWAS/pgx_pop_annotation ]]
then
  echo "PGx-POP directory exists"
else
  echo "Created working directory for star allele haplotype assignment with PGx-POP"
  mkdir ~/Psychiatric_GWAS/Psychiatric_GWAS/pgx_pop_annotation/
fi

# working directory for this is:
cd ~/Psychiatric_GWAS/Psychiatric_GWAS/pgx_pop_annotation/
# add the directory with pgx-pop python scripts to the path
PATH=$PATH:~/bioinformatics_software/PGxPOP/bin/

# add the pyenv executable to the PATH for ensuring that the python version 3 is running for PGx-POP
export PATH="$HOME/.pyenv/bin:$PATH"
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"

# set the current version of python as python 3
pyenv global 3.8.2

# compress and zip the VCF file with all variants to be used for PGx-POP
bgzip -c ../TOPMEDimputationResults/imputed_variants_allCHROMs_INFO_0.3_MAF_0.05_maxMissing20.dose.recode.vcf > imputed_variants_allCHROMs_INFO_0.3_MAF_0.05_maxMissing20.dose.recode.vcf.gz
tabix -p vcf imputed_variants_allCHROMs_INFO_0.3_MAF_0.05_maxMissing20.dose.recode.vcf.gz
# delete the uncompressed vcf file now that it is no longer needed
rm ../TOPMEDimputationResults/imputed_variants_allCHROMs_INFO_0.3_MAF_0.05_maxMissing20.dose.recode.vcf

# run PGxPOP with python 3 (python 3.6 or higher required for PGx-POP) for all possible genes
for currentGene in {CFTR,CYP2B6,CYP2C19,CYP2C9,CYP2D6,CYP3A5,CYP4F2,DPYD,IFNL3,NUDT15,SLCO1B1,TPMT,UGT1A1,VKORC1}
do
  python ~/bioinformatics_software/PGxPOP/bin/PGxPOP.py --vcf imputed_variants_allCHROMs_INFO_0.3_MAF_0.05_maxMissing20.dose.recode.vcf.gz -g "${currentGene}" --phased --build hg19 -o "./PGx_POP_results_${currentGene}.txt"
done

# set the current version of python as the masOS system version (python 2)
pyenv global system

# summarize the PGx-POP results
Rscript --no-save ../scripts/variant_annotation/summarize_PGxPOP_results.R

# make a working directory for using pypgx if it does not exist
if [[ -d ~/Psychiatric_GWAS/Psychiatric_GWAS/pypgx ]]
then
  echo "pypgx haplotype assignment directory exists"
else
  echo "Created working directory for star allele haplotype assignment with pypgx"
  mkdir ~/Psychiatric_GWAS/Psychiatric_GWAS/pypgx/
fi

# working directory for this is:
cd ~/Psychiatric_GWAS/Psychiatric_GWAS/pypgx/

# add the pyenv executable to the PATH for ensuring that the python version 3 is running for pypgx
export PATH="$HOME/.pyenv/bin:$PATH"
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"

# set the current version of python as python 3
pyenv global 3.8.2

# move vcf files to the pypgx directory
mv ../pgx_pop_annotation/imputed_variants_allCHROMs_INFO_0.3_MAF_0.05_maxMissing20.dose.recode.vcf.gz ./imputed_variants_allCHROMs_INFO_0.3_MAF_0.05_maxMissing20.dose.recode.vcf.gz
mv ../pgx_pop_annotation/imputed_variants_allCHROMs_INFO_0.3_MAF_0.05_maxMissing20.dose.recode.vcf.gz.tbi ./imputed_variants_allCHROMs_INFO_0.3_MAF_0.05_maxMissing20.dose.recode.vcf.gz.tbi

# call star alleles with pypgx for all possible pharmacogenes
# with the chip (array) sequencing pipeline
for currentGene in {ABCB1,ABCG2,CACNA1S,CFTR,CYP1A1,CYP1A2,CYP1B1,CYP2A6,CYP2A13,CYP2B6,CYP2C8,CYP2C9,CYP2C19,CYP2D6,CYP2F1,CYP2J2,CYP2R1,CYP2S1,CYP2W1,CYP3A4,CYP3A5,CYP3A7,CYP3A43,CYP4A11,CYP4A22,CYP4B1,CYP4F2,CYP17A1,CYP19A1,CYP26A1,DPYD,F5,G6PD,GSTM1,GSTP1,GSTT1,IFNL3,NAT1,NAT2,NUDT15,POR,PTGIS,RYR1,SLC15A2,SLC22A2,SLCO1B1,SLCO1B3,SLCO2B1,SULT1A1,TBXAS1,TPMT,UGT1A1,UGT1A4,UGT2B7,UGT2B15,UGT2B17,VKORC1,XPC}
do
  pypgx run-chip-pipeline "${currentGene}" "${currentGene}-pipeline" imputed_variants_allCHROMs_INFO_0.3_MAF_0.05_maxMissing20.dose.recode.vcf.gz
done

# set the current version of python as the masOS system version (python 2)
pyenv global system

# uncompress the pypgx results file for each gene
for currentGene in {ABCB1,ABCG2,CACNA1S,CFTR,CYP1A1,CYP1A2,CYP1B1,CYP2A6,CYP2A13,CYP2B6,CYP2C8,CYP2C9,CYP2C19,CYP2D6,CYP2F1,CYP2J2,CYP2R1,CYP2S1,CYP2W1,CYP3A4,CYP3A5,CYP3A7,CYP3A43,CYP4A11,CYP4A22,CYP4B1,CYP4F2,CYP17A1,CYP19A1,CYP26A1,DPYD,F5,G6PD,GSTM1,GSTP1,GSTT1,IFNL3,NAT1,NAT2,NUDT15,POR,PTGIS,RYR1,SLC15A2,SLC22A2,SLCO1B1,SLCO1B3,SLCO2B1,SULT1A1,TBXAS1,TPMT,UGT1A1,UGT1A4,UGT2B7,UGT2B15,UGT2B17,VKORC1,XPC}
do
  cd "${currentGene}-pipeline/"
  unzip results.zip 
  rm *zip
  mv */data.tsv ./data.tsv
  cd ../
done

# summarize the pypgx results
Rscript --no-save ../scripts/variant_annotation/summarize_pypgx_results.R

# move vcf files back to the ../TOPMEDimputationResults/ directory
mv ./imputed_variants_allCHROMs_INFO_0.3_MAF_0.05_maxMissing20.dose.recode.vcf.gz ../TOPMEDimputationResults/imputed_variants_allCHROMs_INFO_0.3_MAF_0.05_maxMissing20.dose.recode.vcf.gz
mv ./imputed_variants_allCHROMs_INFO_0.3_MAF_0.05_maxMissing20.dose.recode.vcf.gz.tbi ../TOPMEDimputationResults/imputed_variants_allCHROMs_INFO_0.3_MAF_0.05_maxMissing20.dose.recode.vcf.gz.tbi


