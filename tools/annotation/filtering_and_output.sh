#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

#export DATE=`date +"%F"`
script_version=4.13_xcat

function add_filter_tag
{
#filter_string should filter for the variants that should be kept; the variants NOT kept will get a filter-value
in_vcf=$1
out=$2
filter_stringx=$3
tag_name=$4


#TEMP-Files
temp_filtered=$in_vcf"temp_filtered.vcf"
temp_filtered_complement=$in_vcf"temp_filtered_complement.vcf"
temp_for_anno=$in_vcf"temp_for_anno.vcf"

echo "filter_vep -i  "$in_vcf" --format vcf -o $temp_filtered --force_overwrite --filter \""$filter_stringx"\"" | bash

bgzip $in_vcf -fc > $in_vcf".gz"
tabix $in_vcf".gz"

bgzip $temp_filtered -f
tabix $temp_filtered".gz" -f

# calculate complement of filtered variants
#$anno_sources/bcftools/bcftools isec -C -o $temp_filtered_complement -O v $in_vcf".gz" $temp_filtered".gz" -w1
bcftools isec -C -o $temp_filtered_complement -O v $in_vcf".gz" $temp_filtered".gz" -w1

# add filter tags
#$anno_sources/bcftools/bcftools filter $temp_filtered_complement -o $temp_for_anno -e CHROM=CHROM -s "$tag_name"
bcftools filter $temp_filtered_complement -o $temp_for_anno -e CHROM=CHROM -s "$tag_name"

bgzip $temp_for_anno -f
tabix $temp_for_anno".gz" -f

# use this file for annotating
#$anno_sources/bcftools/bcftools annotate -c =FILTER -a $temp_for_anno".gz" -o $out $in_vcf".gz" 
bcftools annotate -c =FILTER -a $temp_for_anno".gz" -o $out $in_vcf".gz" 

}


function add_vase_annotation
{
# starts Vase with a pedigree-file and adds INFO tag_name to the lines that were output by VASE
#USAGE: add_vase_annotation in_vcf_vase_function ped_file_vase_function out filter_stringx_vase_function tag_name
in_vcf_vase_function=$1
ped_file_vase_function=$2
out_vase_function=$3
filter_stringx_vase_function=$4
tag_name=$5
 	
vase_output=$in_vcf_vase_function"_VASE_filtered"

export vase="/ceph01/homedirs/sivalingam/vase-0.2.4/bin/vase"

if echo "python $vase -i  "$in_vcf_vase_function" -ped $ped_file_vase_function -o $vase_output $filter_stringx_vase_function" | bash; then

bgzip $vase_output -f
tabix $vase_output".gz" -f


temp_for_anno=$in_vcf_vase_function"_temp_VASE"
bcftools annotate -m +$tag_name -a $vase_output".gz" $vase_output".gz" -o $temp_for_anno # adds TAGGG

bgzip $temp_for_anno -f
tabix $temp_for_anno".gz" -f
bgzip $in_vcf_vase_function -f
tabix $in_vcf_vase_function".gz" -f

bcftools annotate -c INFO/$tag_name -a $temp_for_anno".gz" -o $out_vase_function $in_vcf_vase_function".gz"
return 0

else

return 1
fi
}


function write_output_w_header
{

anno_tab_in=$1
anno_tab_out=$2
split_string1=$3


#VARIABLES define which FIELDS should be PRINTED

split_string2='§?\t§?\t§?\t§?\t§?\t§?\t§?\t§?\t§?\t§?\t§?\t'
split_string2_header='PEDIA\tID_Liste\tGen-Klinik\tGen-Klinik-Anmerkung\tVarianten-Beurteilung\tVariante_einschätzen/raus?\tACMG\t2.Auswerter\t§?\t§?\tX=HYPERLINK("http://localhost:60151//goto?locus=chr"&A2&":"&B2&"-"&B2)\t'
split_string3='%FILTER\t%in_house_num_parent\t%gnomAD_ex_AF\t%gnomAD_ex_AC\t%gnomAD_ex_nhomalt\t%gnomAD_ge_AF\t%gnomAD_ge_AC\t%AF_1000G\t%AA_AF\t%EA_AF\t%in_house_num_parent_not_index\t%in_house_num_index\t%gnomAD_ex_AF_afr\t%gnomAD_ex_AF_amr\t%gnomAD_ex_AF_asj\t%gnomAD_ex_AF_eas\t%gnomAD_ex_AF_fin\t%gnomAD_ex_AF_nfe\t%gnomAD_ex_AF_oth\t%gnomAD_ex_AF_sas\t%gnomAD_ge_AF_AFR\t%gnomAD_ge_AF_AMR\t%gnomAD_ge_AF_ASJ\t%gnomAD_ge_AF_EAS\t%gnomAD_ge_AF_FIN\t%gnomAD_ge_AF_NFE\t%gnomAD_ge_AF_OTH\t%gnomAD_ge_AF_SAS\t%UK10K_AF\t%gnomAD_genomes_AC\t%gnomAD_genomes_nhomalt\t%Existing_variation\t%Ensembl_proteinid\t%REVEL_score\t%REVEL_rankscore\t%MutationTaster_pred\t%PROVEAN_pred\t%PROVEAN_score\t%PrimateAI_pred\t%PrimateAI_score\t%Polyphen2_HVAR_pred\t%SIFT_pred\t%BayesDel_noAF_score\t%BayesDel_noAF_pred\t%phyloP100way_vertebrate\t%phyloP100way_vertebrate_rankscore\t%phastCons100way_vertebrate\t%phastCons100way_vertebrate_rankscore\t%SpliceAI_ALLELE\t%SpliceAI_SYMBOL\t%SpliceAI_DS_AG\t%SpliceAI_DS_AL\t%SpliceAI_DS_DG\t%SpliceAI_DS_DL\t%QD\t%FS\t%MQ\t%MQRankSum\t%ReadPosRankSum\t%END\%tSVLEN\t%ClinVar_CLNDN\t%ClinVar_CLNREVSTAT\t%ClinVar_CLNSIG\t%ClinVar\t%GENE_PHENO\t%CLIN_SIG\t%PHENO\t%Gene'


###########Create Output
# save COLUMN-HEADINGs to the output file
#header1=$($anno_sources/bcftools/bcftools query -u -H -f $split_string1'\n' $anno_tab_in | head -n1) ||  echo ""
header1=$(bcftools query -u -H -f $split_string1'\n' $anno_tab_in | head -n1) ||  echo ""
header2=$(echo $split_string2_header$split_string3 | sed 's:\\t:        :g')
echo $header1\t$header2 | awk -v OFS="\t" '$1=$1' | cut -f2- > $anno_tab_out

# save the DATA to the output file
cat $anno_tab_in | sed 's:||:|0|:g'| sed 's:||:|0|:g' | bcftools +$split_vep -f $split_string1$split_string2$split_string3'\n' -d -A tab | sed 's:\§?: :g' >> $anno_tab_out


########### QC-Check: Count the variants in the output files
num_variants=0
variants_concat="" #first row in output file - the file directory
for file_name in $in_file $trimmed_vcf $trimmed_annotated_vcf $Fvep_in $bcfF_out $FILTER_OUT $anno_tab_in $anno_tab_out;do
	#echo $file_name
	num_variants=$(cat $file_name | grep -v "^#" | cut -f1,2,4,5 | sort| uniq | wc -l) #Zählen..
	variants_concat="${variants_concat} ${file_name} ${num_variants} "
done

#In die erste Zeile Infos zur Annotation speichern.
in_file_info=$(stat --printf="Name: %n size: %s bytes, created: %w last modified %y" $in_file)
ped_file_info=$(stat --printf="Name: %n size: %s bytes, created: %w last modified %y" $pedigree_file)

echo "Erstellt mit annotate_and_filter-Skript, Version $script_version, Input: $in_file_info und $ped_file_info, gefiltert: Filter-Stings: bcftools: $bcf_filter_str, filter_vep: $filter_vep_str, das gefilterte vcf heißt: $anno_tab_in ;Anzahl-Varianten/Dateien: $variants_concat" | cat - $anno_tab_out > temp && mv temp $anno_tab_out
}


function return_AF_string_w_genomes
{
echo "(AF_1000G < "$1" or not AF_1000G) and (gnomAD_ex_AF_afr < "$1" or not gnomAD_ex_AF_afr) and (gnomAD_ex_AF_amr < "$1" or not gnomAD_ex_AF_amr) and (gnomAD_ex_AF_eas < "$1" or not gnomAD_ex_AF_eas) and (gnomAD_ex_AF_fin < "$1" or not gnomAD_ex_AF_fin) and (gnomAD_ex_AF_nfe < "$1" or not gnomAD_ex_AF_nfe) and (gnomAD_ex_AF_oth < "$1" or not gnomAD_ex_AF_oth) and (gnomAD_ex_AF_sas < "$1" or not gnomAD_ex_AF_sas) and (gnomAD_ge_AF_AFR < "$1" or not gnomAD_ge_AF_AFR) and (gnomAD_ge_AF_AMR < "$1" or not gnomAD_ge_AF_AMR) and (gnomAD_ge_AF_EAS < "$1" or not gnomAD_ge_AF_EAS) and (gnomAD_ge_AF_FIN < "$1" or not gnomAD_ge_AF_FIN) and (gnomAD_ge_AF_NFE < "$1" or not gnomAD_ge_AF_NFE) and (gnomAD_ge_AF_SAS < "$1" or not gnomAD_ge_AF_SAS) and (AA_AF < "$1" or not AA_AF) and (EA_AF < "$1" or not EA_AF)"
}

function return_AF_string_wo_genomes
{
echo "(AF_1000G < "$1" or not AF_1000G) and (gnomAD_ex_AF_afr < "$1" or not gnomAD_ex_AF_afr) and (gnomAD_ex_AF_amr < "$1" or not gnomAD_ex_AF_amr) and (gnomAD_ex_AF_eas < "$1" or not gnomAD_ex_AF_eas) and (gnomAD_ex_AF_fin < "$1" or not gnomAD_ex_AF_fin) and (gnomAD_ex_AF_nfe < "$1" or not gnomAD_ex_AF_nfe) and (gnomAD_ex_AF_oth < "$1" or not gnomAD_ex_AF_oth) and (gnomAD_ex_AF_sas < "$1" or not gnomAD_ex_AF_sas) and (AA_AF < "$1" or not AA_AF) and (EA_AF < "$1" or not EA_AF)"
}


in_file=$1 # IN
pedigree_file=$2


# GLOBAL VARIABLES
in_file_base=$(basename $in_file)
DIR=$(dirname "${in_file}")


# Tool variables  
export intermediate_files=$DIR"/intermediate_files"
export externe_festplatte='/gpfs01/homedirs/aschmidt'
export anno_sources=$externe_festplatte'/annotation_sources'
export tsv_annotator2='/ceph01/projects/bioinformatics_resources/Genomes/Human/VEP/b37/tsv_annotator2.py'
#export BCFTOOLS_PLUGINS='/home/sivalingam/scratch/conda/envs/anno/libexec/bcftools'
#export split_vep='/home/sivalingam/scratch/conda/envs/anno/libexec/bcftools/split-vep.so'
export BCFTOOLS_PLUGINS='/home/sivalingam/scratch/conda/envs/bcftools/libexec/bcftools/'
export split_vep='/home/sivalingam/scratch/conda/envs/bcftools/libexec/bcftools/split-vep.so'

# Parameter Variables
# Pindel Adjusted
bcf_filter_str="QUAL>30 || QUAL='.'"
filter_vep_str="$(return_AF_string_w_genomes 0.02)"



######################################################################################################################################################################################################################################
trimmed_vcf=$intermediate_files'/'$in_file_base'.tri.vcf' #OUT
trimmed_annotated_vcf=$in_file'.tri_ann.VCF' #in
Fvep_in=$intermediate_files'/'$in_file_base"_Fvep_IN.vcf" #OUT

# 3. change AF field to not have arbitrary INFO tags
cat $trimmed_annotated_vcf | sed -e 's/|AF|/|AF_1000G|/g' > $Fvep_in 


######################################################################################################################################################################################################################################
#CREATE TAB-delimited FILE for use in EXCEL

# trimmed_annotated_vcf #IN
out_file_pre=$in_file".post_annot.anno_tab" #OUT

# REMOVE old files
rm -f $out_file_pre -f
string='%CHROM\t%POS\t%ID\t%REF\t%ALT\t%Allele\t%QUAL\t[%TGT\t%GT\t%AD\t%DP\t%GQ\t]%IMPACT\t%CADD_PHRED\t%Consequence\t%HGVSc\t%HGVSp\t%SYMBOL\t%genes\t'
bcfF_out="" 
FILTER_OUT=""

write_output_w_header $Fvep_in $out_file_pre $string


######################################################################################################################################################################################################################################
# Fvep_in=$intermediate_files'/'$in_file_base"_Fvep_IN.vcf" #IN
bcfF_out=$intermediate_files'/'$in_file_base"_bcfF.vcf" #OUT

# 4. DELETE variants with low QUAL
#cat $Fvep_in | $anno_sources/bcftools/bcftools filter -i $bcf_filter_str -o $bcfF_out
cat $Fvep_in | bcftools filter -i $bcf_filter_str -o $bcfF_out


######################################################################################################################################################################################################################################
#IN bcfF_out=$intermediate_files'/'$in_file_base"_bcfF.vcf" #IN
Fvep_out=$intermediate_files'/'$in_file_base"_Fvep_OUT.vcf" #OUT

# 5. DELETE variants with high POPULATION FREQUENCY
echo "filter_vep -i  "$bcfF_out" --format vcf -o "$Fvep_out" --force_overwrite --filter \""$filter_vep_str"\"" | bash

######################################################################################################################################################################################################################################
#Fvep_out=$intermediate_files'/'$in_file_base"_Fvep_OUT.vcf" #IN
Refseq_ann_interim=$intermediate_files'/'$in_file_base"_Refseq_ann.vcf" #OUT
Refseq_ann=$intermediate_files'/'$in_file_base"_Refseq_ann_head.vcf" #OUT
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
--vcf_info_field vep_refseq \
-a GRCh37 \
--refseq \
--use_given_ref \
--fasta $anno_sources'/hg19_ref_fasta/human_g1k_v37.fasta' \
-i $Fvep_out \
-o $Refseq_ann_interim

# Remove the VEP header lines to not confuse bcftools by two VEP annotations
cat $Fvep_out | grep "^##" > $intermediate_files'/'$in_file_base"_header.vcf" # take header from before annotation
echo "##INFO=<ID=vep_refseq,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP.">"  >> $intermediate_files'/'$in_file_base"_header.vcf" # add line for the second VEP annotation which won't be misunderstood by bcftools
cat $Refseq_ann_interim | grep -v "^##" | cat $intermediate_files'/'$in_file_base"_header.vcf" - > $Refseq_ann # replace header 


######################################################################################################################################################################################################################################
# Fvep_out=$intermediate_files'/'$in_file_base"_Fvep_OUT.vcf" #IN
Fvep_soft1=$intermediate_files'/'$in_file_base"_Fvep_soft1.vcf" #OUT

# 6. Add FILTER-columns to variants with LOW IMPACT
add_filter_tag $Refseq_ann $Fvep_soft1 "((IMPACT is HIGH) or (IMPACT is MODERATE) or (IMPACT is LOW) or (SpliceAI_SYMBOL))" "MOD"


######################################################################################################################################################################################################################################
# Fvep_soft1=$intermediate_files'/'$in_file_base"_Fvep_soft1.vcf" #IN
Fvep_soft2=$intermediate_files'/'$in_file_base"_Fvep_soft2.vcf" #OUT

# 7. Add FILTER-columns to variants with HIGH POPULATION FREQUENCY
add_filter_tag $Fvep_soft1 $Fvep_soft2 "(gnomAD_ex_nhomalt < 4 or not gnomAD_ex_nhomalt)" "HOM4" 


######################################################################################################################################################################################################################################
# Fvep_soft1=$intermediate_files'/'$in_file_base"_Fvep_soft2.vcf" #IN
Fvep_soft3=$intermediate_files'/'$in_file_base"_Fvep_soft3.vcf" #OUT

# 8. Add FILTER-columns to variants with HIGH POPULATION FREQUENCY
add_filter_tag $Fvep_soft2 $Fvep_soft3 "(not gnomAD_ex_AC or gnomAD_ex_AC < 9)" "AC8" 


######################################################################################################################################################################################################################################
# Fvep_soft4=$intermediate_files'/'$in_file_base"_Fvep_soft3.vcf" #IN
Fvep_soft4=$intermediate_files'/'$in_file_base"_Fvep_soft4.vcf" #OUT

# 8. Add FILTER-columns to variants with HIGH POPULATION FREQUENCY
add_filter_tag $Fvep_soft3 $Fvep_soft4 "(in_house_num_parent < 7 or not in_house_num_parent)" "InH6" 

######################################################################################################################################################################################################################################
# Fvep_soft4=$intermediate_files'/'$in_file_base"_Fvep_soft4.vcf" #IN
Fvep_soft5=$intermediate_files'/'$in_file_base"_Fvep_soft5.vcf" #OUT

# 8. Add FILTER-columns to variants with HIGH POPULATION FREQUENCY
add_filter_tag $Fvep_soft4 $Fvep_soft5 "(in_house_num_parent < 2 or not in_house_num_parent)" "InH1" 


######################################################################################################################################################################################################################################
# Fvep_soft5=$intermediate_files'/'$in_file_base"_Fvep_soft5.vcf" #IN
Fvep_soft6=$intermediate_files'/'$in_file_base"_Fvep_soft6.vcf" #OUT

# 9. Add FILTER-columns to variants with BAD QUALITY
bcftools filter -m+ -o $Fvep_soft6 -e "INFO/QD < 2.0 || INFO/FS > 60.0 || INFO/MQ < 40.0 || INFO/MQRankSum < -12.5 || INFO/ReadPosRankSum < -8.0" -s "Q" $Fvep_soft5


######################################################################################################################################################################################################################################
#10. NEW VASE Inheritance filtering https://github.com/david-a-parry/vase
# Fvep_soft6=$intermediate_files'/'$in_file_base"_Fvep_soft6.vcf" #IN
#Intermediate files
VASE_TEMP1=$intermediate_files'/'$in_file_base"_VASE_TEMP1.vcf"
VASE_TEMP2=$intermediate_files'/'$in_file_base"_VASE_TEMP2.vcf"
VASE_TEMP3=$intermediate_files'/'$in_file_base"_VASE_TEMP3.vcf"

FILTER_OUT=$intermediate_files'/'$in_file_base"_FILTER_OUT.vcf" #OUT
#USAGE: add_vase_annotation in_vcf ped_file out filter_stringx tag_name

#test if vase annotation works out:
if add_vase_annotation $Fvep_soft6 $pedigree_file $VASE_TEMP1 "--recessive --gq -1" "VASE_AR"; then 
add_vase_annotation $VASE_TEMP1 $pedigree_file $VASE_TEMP2 "--recessive --gq -1 --keep_filters PASS AC8 Q InH1" "VASE_STRICT_AR" # do not use variants with all other than AC8 and Q

# Try to add autosomal dominant / de novo filter
	if add_vase_annotation $VASE_TEMP2 $pedigree_file $VASE_TEMP3 "--dominant --de_novo --gq -1" "VASE_AD"; then
		add_vase_annotation $VASE_TEMP3 $pedigree_file $FILTER_OUT "--dominant --de_novo --gq 50 --pass_filters" "VASE_STRICT_AD"; 
		string_VASE='%CHROM\t%POS\t%ID\t%REF\t%ALT\t%Allele\t%QUAL\t[%TGT\t%GT\t%AD\t%DP\t%GQ\t]%vep_refseq\t%SVLEN\t%VASE_AR\t%VASE_AD\t%VASE_STRICT_AR\t%VASE_STRICT_AD\t%IMPACT\t%CADD_PHRED\t%Consequence\t%HGVSc\t%HGVSp\t%SYMBOL\t%genes\t' # add inheritance column
		string_janno='%CHROM\t%POS\t%ID\t%REF\t%ALT\t%Allele\t%QUAL\t[%TGT\t%GT\t%AD\t%DP\t%GQ\t]%INHERITANCE\t%vep_refseq\t%SVLEN\t%VASE_AR\t%VASE_AD\t%VASE_STRICT_AR\t%VASE_STRICT_AD\t%IMPACT\t%CADD_PHRED\t%Consequence\t%HGVSc\t%HGVSp\t%SYMBOL\t%genes\t' # add inheritance column
	else
		cp -f $VASE_TEMP2 $FILTER_OUT	
		string_VASE='%CHROM\t%POS\t%ID\t%REF\t%ALT\t%Allele\t%QUAL\t[%TGT\t%GT\t%AD\t%DP\t%GQ\t]%vep_refseq\t%SVLEN\t%VASE_AR\t%VASE_STRICT_AR\t%IMPACT\t%CADD_PHRED\t%Consequence\t%HGVSc\t%HGVSp\t%SYMBOL\t%genes\t' # add inheritance column
        	string_janno='%CHROM\t%POS\t%ID\t%REF\t%ALT\t%Allele\t%QUAL\t[%TGT\t%GT\t%AD\t%DP\t%GQ\t]%INHERITANCE\t%vep_refseq\t%SVLEN\t%VASE_AR\t%VASE_STRICT_AR\t%IMPACT\t%CADD_PHRED\t%Consequence\t%HGVSc\t%HGVSp\t%SYMBOL\t%genes\t' # add inheritance column
	fi

else # is executed if VASE completely failed:

string_VASE='%CHROM\t%POS\t%ID\t%REF\t%ALT\t%Allele\t%QUAL\t[%TGT\t%GT\t%AD\t%DP\t%GQ\t]%vep_refseq\t%SVLEN\t%IMPACT\t%CADD_PHRED\t%Consequence\t%HGVSc\t%HGVSp\t%SYMBOL\t%genes\t' #no inheritance column
string_janno='%CHROM\t%POS\t%ID\t%REF\t%ALT\t%Allele\t%QUAL\t[%TGT\t%GT\t%AD\t%DP\t%GQ\t]%INHERITANCE\t%vep_refseq\t%SVLEN\t%IMPACT\t%CADD_PHRED\t%Consequence\t%HGVSc\t%HGVSp\t%SYMBOL\t%genes\t' 
cp -f $Fvep_soft6 $FILTER_OUT

fi

######################################################################################################################################################################################################################################
#11. CREATE TAB-delimited FILE for use in EXCEL

# janno_out=$in_file"_janno.vcf" #IN
out_file_FILTERED_VASE=$in_file".VASE_FILTERED.anno_tab" #OUT
out_file_VASE_TSV=$in_file".VASE_FILTERED_TSV.anno_tab"
out_file_VASE_TSV_TEMP=$intermediate_files"/VASE_FILTERED_TSV_TEMP"

# REMOVE old files
rm -f $out_file_FILTERED_VASE -f

write_output_w_header $FILTER_OUT $out_file_FILTERED_VASE $string_VASE

######################################################################################################################################################################################################################################
# EXPERIMENTAL: add custom annotations to TSV file
#out_file_VASE_TSV=$in_file".VASE_FILTERED_TSV.anno_tab"
#out_file_VASE_TSV_TEMP=$intermediate_files"/VASE_FILTERED_TSV_TEMP"

#The file tsv_annotator2.py is a python script that simply loads two text files as pandas data frames and applies a join opperation
# INPUT file: -i
# Output file name: -o
# source of annotations (tab delimited): -a
# Column name that is used for joining in the input-file: -k
# Column name that is used for joining in the source of annotations: -z
# columns to append to the input file to form the output file: -n
# -l is depracted - it was used to indicate where in the TSV file the new columns should be inserted.

#Add the gnomAD constraint metrices
#python $anno_sources'/tsv_annotator2.py' \
python $tsv_annotator2 \
-i $out_file_FILTERED_VASE \
-o $out_file_VASE_TSV \
-a $anno_sources'/gnomad.v2.1.1.lof_metrics.by_gene.txt' \
-k %Gene \
-z gene_id \
-n "pLI","mis_z","oe_mis","oe_lof" \
-l "%Gene"

#Add some annotations from the dbNSFP4.1_gene file to the output
#python $anno_sources'/tsv_annotator2.py' \
python $tsv_annotator2 \
-i $out_file_VASE_TSV \
-o $out_file_VASE_TSV \
-a $anno_sources'/dbNSFP41/dbNSFP4.1_gene' \
-k %Gene \
-z Ensembl_gene \
-n "Tissue_specificity(Uniprot)","Expression(egenetics)","Expression(GNF/Atlas)","MGI_mouse_phenotype","GO_biological_process","GO_cellular_component","GO_molecular_function" \
-l "%Gene"

cut -f 3- $out_file_VASE_TSV > $out_file_VASE_TSV_TEMP
cat $out_file_VASE_TSV_TEMP > $out_file_VASE_TSV


######################################################################################################################################################################################################################################
# 12. pre-filter Jannovar files to possibly get rid of problematic coordinates
# FILTER_OUT=$intermediate_files'/'$in_file_base"_FILTER_OUT.vcf" #IN
janno_in=$intermediate_files'/'$in_file_base"janno_in_prefilter.vcf" #OUT

#$anno_sources/bcftools/bcftools filter $FILTER_OUT -o $janno_in -e 'FILTER="InH6"'
bcftools filter $FILTER_OUT -o $janno_in -e 'FILTER="InH6"'


######################################################################################################################################################################################################################################
# 13. !!!!TRY!!!! Add INHERITANCE-columns with JANNOVAR using the ped-file
#janno_in=$intermediate_files'/'$in_file_base"janno_in_prefilter.vcf" #OUT
janno_out=$intermediate_files'/'$in_file_base".janno.vcf"
out_file_FILTERED_janno="" #set the variables that they are not unbound
out_file_janno_TSV=""
out_file_janno_TSV_TEMP=""

if jannovar annotate-vcf -i $janno_in -o $janno_out --pedigree-file $pedigree_file -d $anno_sources/referenceGenome/hg19_refseq.ser ; then


######################################################################################################################################################################################################################################

# 14. CREATE TAB-delimited FILE for use in EXCEL
# janno_out=$in_file"_janno.vcf" #IN
out_file_FILTERED_janno=$in_file".janno_FILTERED.anno_tab" #OUT

# REMOVE old files
rm -f $out_file_FILTERED_janno

write_output_w_header $janno_out $out_file_FILTERED_janno $string_janno

#####################################################################################################################################################################################################################################
# EXPERIMENTAL: add custom annotations to TSV-file
out_file_janno_TSV=$in_file".janno_FILTERED_TSV.anno_tab"
out_file_janno_TSV_TEMP=$intermediate_files"/janno_FILTERED_TSV_TEMP"

#see above for more informations

#Add the gnomAD constraint metrices
#python $anno_sources'/tsv_annotator2.py' \
python $tsv_annotator2 \
-i $out_file_FILTERED_janno \
-o $out_file_janno_TSV \
-a $anno_sources'/gnomad.v2.1.1.lof_metrics.by_gene.txt' \
-k %Gene \
-z gene_id \
-n "pLI","mis_z","oe_mis","oe_lof" \
-l "%Gene"

#Add some annotations from the dbNSFP4.1_gene file to the output
#python $anno_sources'/tsv_annotator2.py' \
python $tsv_annotator2 \
-i $out_file_janno_TSV \
-o $out_file_janno_TSV \
-a $anno_sources'/dbNSFP41/dbNSFP4.1_gene' \
-k %Gene \
-z Ensembl_gene \
-n "Tissue_specificity(Uniprot)","Expression(egenetics)","Expression(GNF/Atlas)","MGI_mouse_phenotype","GO_biological_process","GO_cellular_component","GO_molecular_function" \
-l "%Gene"

cut -f 3- $out_file_janno_TSV > $out_file_janno_TSV_TEMP
cat $out_file_janno_TSV_TEMP > $out_file_janno_TSV 
fi

######################################################################################################################################################################################################################################
# Clean UP

zip -m $in_file".annotated.zip" $out_file_VASE_TSV $out_file_FILTERED_VASE $out_file_janno_TSV $out_file_FILTERED_janno $out_file_pre $trimmed_annotated_vcf $trimmed_annotated_vcf"_summary.html"

rm -rf $intermediate_files

######################################################################################################################################################################################################################################

