#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

#export DATE=`date +"%F"`
script_version=4.13_xcat

in_file=$1 # IN
pedigree_file=$2


# GLOBAL VARIABLES
in_file_base=$(basename $in_file)
DIR=$(dirname "${in_file}")


# Tool variables  
export intermediate_files=$DIR"/intermediate_files"
export externe_festplatte='/gpfs01/homedirs/aschmidt'
export anno_sources=$externe_festplatte'/annotation_sources'

# Parameter Variables

# REMOVE old files / folders

rm -rf $intermediate_files # make sure no old files will corrupt the output
mkdir $intermediate_files

######################################################################################################################################################################################################################################
# 0 Download recent ClinVar file, for parallel processing: Do this first and skip this step.

#(cd $intermediate_files && cp /gpfs01/homedirs/aschmidt/M14-1693/intermediate_files/clinvar.vcf.gz .) 
#(cd $intermediate_files && cp /gpfs01/homedirs/aschmidt/M14-1693/intermediate_files/clinvar.vcf.gz.tbi .) 

(cd $intermediate_files && cp /ceph01/projects/bioinformatics_resources/Genomes/Human/GATK/b37/clinvar.vcf.gz .) 
(cd $intermediate_files && cp /ceph01/projects/bioinformatics_resources/Genomes/Human/GATK/b37/clinvar.vcf.gz.tbi .) 

#(cd $intermediate_files && curl -O ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz)
#(cd $intermediate_files && curl -O ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi)


######################################################################################################################################################################################################################################
# in_file=$1 # IN
trimmed_vcf=$intermediate_files'/'$in_file_base'.tri.vcf' #OUT

# 1. Trim Alternative Alleles
#cat $in_file | $anno_sources/bcftools/bcftools view --trim-alt-alleles > $trimmed_vcf
cat $in_file | bcftools view --trim-alt-alleles > $trimmed_vcf


######################################################################################################################################################################################################################################
#IN: trimmed_vcf=$intermediate_files'/'$in_file_base'_tri.vcf' 
trimmed_annotated_vcf=$in_file'.tri_ann.VCF' #OUT
#Parameter: anno_sources

# 2. VEP Annotation
vep -v  \
--cache \
--offline \
--dir_cache $anno_sources \
--force_overwrite \
--fork 6 \
--hgvs \
--format vcf \
--numbers \
--vcf \
--gene_phenotype \
--pick_allele_gene \
--af \
--af_esp \
--check_existing \
-a GRCh37 \
--custom $intermediate_files'/clinvar.vcf.gz',ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
--custom $anno_sources'/genes_annotation.bed.gz',genes,bed,overlap,0 \
--custom $anno_sources'/counts_in_house2019-11-08_fix.vcf.gz',in_house,vcf,exact,0,num_parent,num_index,num_parent_not_index \
--plugin CADD,$anno_sources'/whole_genome_SNVs.tsv.gz',$anno_sources'/InDels.tsv.gz' \
--fasta $anno_sources'/hg19_ref_fasta/human_g1k_v37.fasta' \
--custom $anno_sources'/gnomad.exomes.r2.1.1.sites.vcf.bgz',gnomAD_ex,vcf,exact,0,AF,AN,AC,AF_afr,AF_amr,AF_asj,AF_eas,AF_fin,AF_nfe,AF_oth,AF_sas,nhomalt \
--custom $anno_sources'/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz',gnomAD_ge,vcf,exact,0,AF,AN,AC,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_SAS \
--plugin dbNSFP,$anno_sources'/dbNSFP41/dbNSFP4.1a.gz',Ensembl_proteinid,MutationTaster_pred,MutationTaster_score,PROVEAN_pred,PROVEAN_score,PrimateAI_pred,PrimateAI_score,Polyphen2_HVAR_pred,SIFT_pred,UK10K_AF,REVEL_score,REVEL_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,BayesDel_noAF_score,BayesDel_noAF_pred,gnomAD_genomes_AC,gnomAD_genomes_nhomalt \
--custom $anno_sources'/splice_ai_CrossMap_v38tov37_lift_sorted_stripped.vcf.gz',SpliceAI,vcf,exact,0,ALLELE,SYMBOL,DS_AG,DS_AL,DS_DG,DS_DL  \
-i $trimmed_vcf \
-o $trimmed_annotated_vcf 
