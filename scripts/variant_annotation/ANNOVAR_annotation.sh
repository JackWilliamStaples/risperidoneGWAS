# add the bioinformatics software directory to the path
PATH=$PATH:~/bioinformatics_software/

# if a working directory for annotating genetic data does not exist, create it
if [[ -d ~/Psychiatric_GWAS/Psychiatric_GWAS/annovar_annotation ]]
then
  echo "Variant annotation working directory exists"
else
  echo "Created variant annnotation working directory"
  mkdir ~/Psychiatric_GWAS/Psychiatric_GWAS/annovar_annotation/
fi

# working directory for this is:
cd ~/Psychiatric_GWAS/Psychiatric_GWAS/annovar_annotation/

##### code for downloading ANNOVAR data for annotation, only run this code if necessary
#	#Downloaded refGeneWithVer hg19 annotation information
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --hgvs --downdb --webfrom annovar refGeneWithVer ~/bioinformatics_software/annovar/humandb/
#	#Downloaded hg19 cytogenetic band identifiers
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --downdb --buildver hg19 cytoBand ~/bioinformatics_software/annovar/humandb/
#	#Downloaded current 1000 genomes annotation for whole genome data
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar 1000g2015aug ~/bioinformatics_software/annovar/humandb/
#	#Downloaded current exac annotation for whole exome data
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar exac03 ~/bioinformatics_software/annovar/humandb/
#	#Downloaded dbnsfp30a: includes SIFT, PolyPhen2 HDIV, PolyPhen2 HVAR, LRT, MutationTaster, MutationAssessor, FATHMM, MetaSVM, MetaLR, VEST, CADD, GERP++, DANN, fitCons, PhyloP and SiPhy scores, but ONLY on coding variants
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar dbnsfp35a ~/bioinformatics_software/annovar/humandb/
#	#Downloaded dbscsnv11: dbscSNV version 1.1 for splice site prediction by AdaBoost and Random Forest, which score how likely that the variant may affect splicing
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar dbscsnv11 ~/bioinformatics_software/annovar/humandb/
#	#Downloaded dbsnp 138 database and avsnp150 to see if there were any differences in the version/to see if one was better than the other for recognizing indels (the avsnp data bases are dbsnp but with left normalization)
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar snp138 ~/bioinformatics_software/annovar/humandb/
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar avsnp138 ~/bioinformatics_software/annovar/humandb/
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar avsnp142 ~/bioinformatics_software/annovar/humandb/
#	# Downloaded clinVar registry
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar clinvar_20190305 ~/bioinformatics_software/annovar/humandb/
#	# Downloaded catalog of variants previously reported in GWAS studies
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl -build hg19 -downdb gwasCatalog ~/bioinformatics_software/annovar/humandb/
#	# Download transcription factor binding site motif annotation
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl -buildver hg19 -downdb tfbsConsSites ~/bioinformatics_software/annovar/humandb/
#	#Downloaded GWAVA scores (tentative) [too many GB]
#		#perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar gwava ~/bioinformatics_software/annovar/humandb/
#	# Downloaded CADD scores (tentative) [350GB]
#		#perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar cadd13 ~/bioinformatics_software/annovar/humandb/
#	# Downloaded DANN scores (tentative) [200GB]
#		#perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar dann ~/bioinformatics_software/annovar/humandb/
#	# Downloaded COSMIC database
#		perl ~/bioinformatics_software/annovar/annotate_variation.pl --buildver hg19 --downdb --webfrom annovar cosmic70 ~/bioinformatics_software/annovar/humandb/
	
perl ~/bioinformatics_software/annovar/table_annovar.pl "../TOPMEDimputationResults/imputed_variants_allCHROMs_INFO_0.3_MAF_0.05_maxMissing20.dose.recode.vcf" ~/bioinformatics_software/annovar/humandb/ --outfile "ANNOVAR_annotation" --buildver hg19 --onetranscript --nastring '.' --protocol refGeneWithVer,cytoBand,gwasCatalog,tfbsConsSites,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,exac03,dbnsfp35a,dbscsnv11,snp138,avsnp138,avsnp142,clinvar_20190305,cosmic70 --operation g,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f --remove --vcfinput
# delete unneccesary/temporary files produced by ANNOVAR
rm *avinput
rm *vcf

# create an annotation output directory if it does not exist and
# move all of the final annotated files to the annotation_output directory
if [[ -d ../annotation_output ]]
then
  echo "Annotation output directory exists"
else
  echo "Created annotation output directory"
  mkdir ../annotation_output/
fi
mv *multianno* ../annotation_output/
  
# make a directory for storing a pharmGKB annotation file from the pharmGKB website (https://www.pharmgkb.org/downloads) if it does not exist
if [[ -d ./pharmGKB_annotation_file ]]
then
  echo "Directory for pharmGKB file exists"
else
  echo "Created directory for pharmGKB annotation file"
  mkdir ./pharmGKB_annotation_file/
fi

# change directories to the pharmGKB_annotation_file directory
cd ./pharmGKB_annotation_file/
# download the pharmGKB clinicalAnnotations.zip  file to this directory
wget https://api.pharmgkb.org/v1/download/file/data/clinicalAnnotations.zip 
# unzip the clinicalAnnotations.zip file
unzip clinicalAnnotations.zip
# change directories back to the annovar_annotation directory
cd ../

# Annotate the annovar text files additionally with ADME optimized prediction framework scrores as 
# described in Zhou, et al., 2018, The Pharmacogenomics Journal, 19:115-126.
# the tidyverse package must be installed for this script to run
Rscript --no-save ../scripts/variant_annotation/ADME_optimized_framework_prediction.R
# delete files that are no longer needed
rm ../annotation_output/ANNOVAR_annotation.hg19_multianno.txt

### prepare a compressed VCF file for each chromosome for annotation with the variant effect predictor (VEP)
# for each chromosome in the VCF file
for currentChromosome in {1..23}
do
  # split the vcf file by chromosome
  plink --vcf ../TOPMEDimputationResults/imputed_variants_allCHROMs_INFO_0.3_MAF_0.05_maxMissing20.dose.recode.vcf --keep-allele-order --chr ${currentChromosome} --recode vcf --out chr${currentChromosome}
  # compress the resulting file with gzip
  gzip chr${currentChromosome}.vcf
  # delete temporary files
  rm *bed
  rm *bim
  rm *fam
  rm *log
  rm *nosex
done

# create a directory for storing variant effect predictor (VEP) annotation output if it doesn't exist
if [[ -d ../VEP_annotation_output ]]
then
  echo "VEP annotation directory exists"
else
  echo "Created VEP annotation directory"
  mkdir ../VEP_annotation_output/
fi
    
# manually perform an ANNOTATION of the VCF file for each chromosome using the VEP web interface
# see VEP_settings.pdf screenshot in ../VEP_annotation_output for the settings used for the annotation
  # the input file used were:
      # ./chr1.vcf.gz
      # ./chr2.vcf.gz
      # ./chr3.vcf.gz
      # ./chr4.vcf.gz
      # ./chr5.vcf.gz
      # ./chr6.vcf.gz
      # ./chr7.vcf.gz
      # ./chr8.vcf.gz
      # ./chr9.vcf.gz
      # ./chr10.vcf.gz
      # ./chr11.vcf.gz
      # ./chr12.vcf.gz
      # ./chr13.vcf.gz
      # ./chr14.vcf.gz
      # ./chr15.vcf.gz
      # ./chr16.vcf.gz
      # ./chr17.vcf.gz
      # ./chr18.vcf.gz
      # ./chr19.vcf.gz
      # ./chr20.vcf.gz
      # ./chr21.vcf.gz
      # ./chr22.vcf.gz
      # ./chr23.vcf.gz
# the results were manually downloaded from the VEP website interface and saved to ../VEP_annotation_output/ with the file names:
      # ../VEP_annotation_output/chr1.txt
      # ../VEP_annotation_output/chr2.txt
      # ../VEP_annotation_output/chr3.txt
      # ../VEP_annotation_output/chr4.txt
      # ../VEP_annotation_output/chr5.txt
      # ../VEP_annotation_output/chr6.txt
      # ../VEP_annotation_output/chr7.txt
      # ../VEP_annotation_output/chr8.txt
      # ../VEP_annotation_output/chr9.txt
      # ../VEP_annotation_output/chr10.txt
      # ../VEP_annotation_output/chr11.txt
      # ../VEP_annotation_output/chr12.txt
      # ../VEP_annotation_output/chr13.txt
      # ../VEP_annotation_output/chr14.txt
      # ../VEP_annotation_output/chr15.txt
      # ../VEP_annotation_output/chr16.txt
      # ../VEP_annotation_output/chr17.txt
      # ../VEP_annotation_output/chr18.txt
      # ../VEP_annotation_output/chr19.txt
      # ../VEP_annotation_output/chr20.txt
      # ../VEP_annotation_output/chr21.txt
      # ../VEP_annotation_output/chr22.txt
      # ../VEP_annotation_output/chr23.txt

# delete the vcf files for each chromosome after the annotation is complete
rm *vcf.gz

# Parse the VEP results with an R script for any interesting annotations that are not included as part of the ANNOVAR annotation
# and join the useful columns from the VEP annotation with ANNOVAR and pharmGKB annotated data.
# Additionally, add pharmGKB clinical annotations and star alleles (based on rs-number matching) to the annotation file
Rscript --no-save ../scripts/variant_annotation/combine_ANNOVAR_and_pharmGKBdata_with_VEP_data.R

# filter the variant annotation data to variants in cyp2d6, cyp2c9, and cyp2c19 loci
Rscript --no-save ../scripts/variant_annotation/identify_cyp2d6_2c9_2c19_alleles.R



          