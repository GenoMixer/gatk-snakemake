Variant calling workflow
------------------------

This workflow can be used as end-to-end pipeline for variant calling from
fastq input. If follows the `GATK best practices guidelines`.

This workflow is using GATK 4, so it may be considered incompatible to
earlier versions using prior versions (i.e. 3.8)

Paired end reads are aligned using bwa mem and filtered for duplicates
utilizing the picard toolchain. Also base recalibration is applied using
the GATK toolchain.

For variant calling the Haplotype Caller is used with GVCF mode, prompting
inclusion of non-variant sites. The resulting sites are genotyped using
a joint genotyping approach over all samples included in the sample list.

This workflow does not filter the resulting vcf files.

Quality control is performed at various stages throughout the pipeline. The
input is evaluated using fastqc. Bams are evaluated using the results from
picard's CollectHsMetrics and MarkDuplicates. Stats for vcf files are
collected from bcftools stats. Finally, the peddy tool can be applied to
check for the validity of the supplied pedigree.

Quality Control Tools:

- fastqc
- verifybamid
- CollectHsMetrics
- MarkDuplicates
- bcftools stats
- peddy

Overview of the workflow:

- [Optional] Download input from remote
- Alignment:

  - bwa
  - mark duplicates
  - bqsr

- Variant Calling:

  - Haplotype Caller
  - GenotypeGVCF

Results of the workflow include:

- per sample vcf, bam
- per trio vcf (families are specified in the sample config)
- quality control reports

.. _GATK best practices workflow: https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145
