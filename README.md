## General Description

This snakemake pipeline for variant calling of whole-exome (WES) or whole-genome (WGS) sequencing data follows basic gatk best practice recommendations. In principal the data is mapped to the decoy reference genome (GRCh37) and variant calling is performed on the core refseq targets from Twist (mainly used in Bonn) using a interval padding of +/-100bp in case of WES or against the complete genome in case of WGS data. Different configurations regarding the target regions and interval padding are possible by modifying the config.yml. Sophisticated QC is carried out on multiple levels of the data such as fastq, bam and vcfs. Variant annotation is carried out using vep. For variant filtering custom scripts, population-wide AF and inhouse variant databases are applied. For the Excel-based analysis, tsv files are generated, and for the web-based analysis, all project-specific vcfs are uploaded into Varfish.



## Standard operating procedure (SOP)

1. Clone the repository:

    ```bash
    git clone https://git.meb.uni-bonn.de/workflows/gatk_snakemake.git
    ```

2. Create a directory e.g. "input" to save the raw fastq files:

    ```bash
    mkdir input
    ```

3. Generate a raw data sheet with the corresponding fastq files in a format provided in the "files.tsv"

4. Generate a sample sheet with the corresponding phenotypic data in a format provided in the "samples.tsv"


5. Create symlinks to the reference data: 
    
    ```bash
    ln -s /ceph01/projects/bioinformatics_resources/Genomes/Human/GATK/b37/* library/
    ```

6. Edit the config.yml by setting the corresponding "Flowcell_ID"


7. Start the session in a screen:

    ```bash
    screen -S Flowcell_ID
    ```

8. Activate the snakemake environemnt:

    ```bash
    conda activate snakemake
    ```

9. Start a dry run (-n) and estimate the numbers of jobs provided by snakemake:

    ```bash
    snakemake --use-conda --profile slurm -k --rerun-incomplete --cluster-config cluster.yml --configfile config.yml -n
    ```

10. Start real execution with estimated numbers of jobs that can be run in parallel (-j XXX):

    ```bash
    snakemake --use-conda --profile slurm -j XXX -k --rerun-incomplete --cluster-config cluster.yml --configfile config.yml
    ```

## Software Version

| Tool      | Version |
| :---      | :--- |
| BWA       | 3.7  / 0.7.15  |
| GATK      | 3.7  / 4.1.0.0 |
| BCFtools  | 1.10.2    |
| FastQC    | 0.11.9    |
| GLnexus   | 1.3.2     |
| Mosdepth  | 0.2.6     |
| MultiQC   | 1.9       |
| Peddy     | 0.4.6     |
| Picard    | 2.20.1    |
| Pindel    | 0.2.5b9   |
| Samtools  | 1.9       |
| VCFtools  | 0.1.16    |
| VEP       | 101.0       |
| VerifyBamID  | 1.1.3    |


## Database Version

| Database      | Version |
| :---      | :--- |
| Genome       | GRCh37_decoy   |
| GnomAD_WES       | 2.1.1    |
| GnomAD_WGS       | 2.0.1    |
| CADD      | 1.4   |
| ClinVar  |  20.03.22   |
| SpliceAI    |  1.4   |
| ENSEMBL  | 101     |
| dbNSFP  | 4.1a     |
| Inhouse  | 20.04.22     |
