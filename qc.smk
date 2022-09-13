rule fastqc:
    input: unpack(get_fastq)
    output:
        html=[
            temp("output/{sample}/{sample}_1.html"),
            temp("output/{sample}/{sample}_2.html")
        ],
        zip=[
            temp("output/{sample}/{sample}_1_fastqc.zip"),
            temp("output/{sample}/{sample}_2_fastqc.zip")
        ]
    params:
        extra=""
    log: "logs/{sample}/{sample}.fastqc.log"
    group: "alignment"
    resources:
        runtime=120,
        mem_mb=16192
    wrapper: "file:wrapper/fastqc"

rule peddy:
    input:
        vcf=rules.merge_vcfs.output.vcf,
        ped=rules.create_all_pedigree.output
    output:
        ped=temp("output/{project}/{project}.peddy.ped"),
        het_check=temp("output/{project}/{project}.het_check.csv"),
        ped_check=temp("output/{project}/{project}.ped_check.csv"),
        sex_check=temp("output/{project}/{project}.sex_check.csv"),
        background_pca=temp("output/{project}/{project}.background_pca.json")
    params:
        prefix="{project}"
    log: "logs/{project}/peddy.log"
    threads: 4
    resources:
        runtime=120,
        mem_mb=16192
    wrapper: "file:wrapper/peddy"

rule collect_hs_metrics:
    input:
        bam=rules.apply_bqsr.output.bam,
        bai=rules.apply_bqsr.output.bai,
        target_intervals=config["ref"]["target-region"],
        bait_intervals=config["ref"]["bait-region"],
        ref=config["ref"]["genome"]
    output:
        metrics=temp("output/{sample}/{sample}.hsmetrics.txt")
    params:
        extra="",
        jvm_args=""
    log: "logs/{sample}/collect_hs_metrics.log"
    resources:
        runtime=300,
        mem_mb=32192
    wrapper: "file:wrapper/picard/CollectHsMetrics"

rule mosdepth:
    input:
        exons=config["ref"]["exons"],    
        bam=rules.apply_bqsr.output.bam,
        bai=rules.apply_bqsr.output.bai
    output:
        dist="output/{sample}/{sample}.mosdepth.dist.txt",
        thresholds="output/{sample}/{sample}.mosdepth.thresholds.bed.gz",
        summary="output/{sample}/{sample}.mosdepth.summary.txt"
    log: "logs/{sample}/{sample}.mosdepth.log"
    params:
        extra="--no-per-base --fast-mode --mapq 20 --thresholds 0,10,15,20,30"
    threads: 4
    resources:
        runtime=300,
        mem_mb=32192
    wrapper: "file:wrapper/mosdepth"

rule variant_stats_all:
    input: rules.merge_vcfs.output.vcf
    output: temp("output/{project}/{project}.stats.txt")
    log: "logs/{project}/variant_stats.log"
    resources:
        runtime=120,
        mem_mb=16096
    wrapper: "file:wrapper/bcftools/stats"

rule variant_stats:
    input: rules.select_variants.output.vcf
    output: temp("output/{sample}/{sample}.stats.txt")
    log: "logs/{sample}/variant_stats.log"
    resources:
        runtime=120,
        mem_mb=16096
    wrapper: "file:wrapper/bcftools/stats"

rule verifybamid:
    input: 
        vcf=config["ref"]["1kg-variants"],
        bam=rules.apply_bqsr.output.bam,
        bai=rules.apply_bqsr.output.bai
    output:
        self_sm=temp("output/{sample}/{sample}.selfSM")
    params:
        extra="--ignoreRG --verbose"
    log: "logs/{sample}/verifybamid.log"
    resources:
        runtime=300,
        mem_mb=32384
    wrapper: "file:wrapper/verifybamid"

rule multiqc:
    input: get_qc_input
    output: report("output/{project}/{project}.html", caption="report/multiqc.rst", category="Quality control")
    params:
        extra=""
    log: "logs/{project}/multiqc.log"
    resources:
        runtime=60,
        mem_mb=16096
    wrapper: "file:wrapper/multiqc"
    
rule variant_eval:
    input: 
        vcf=rules.select_variants.output.vcf,
        ref=config["ref"]["genome"]
    output: eval="output/{sample}/{sample}.eval.txt"
    params:
        dbsnp=config["ref"]["dbsnp"]
    log: "logs/{sample}/variant_eval.log"
    resources:
        runtime=300,
        mem_mb=16096
    wrapper: "file:wrapper/gatk/VariantEvaluation"
    
rule happy:
    input:
        vcf=rules.select_variants.output.vcf,
        giab=config["ref"]["giab"],
        ref=config["ref"]["genome2"]
    output:
        summary="output/{sample}/{sample}.summary.csv"
    params:
        prefix="{sample}"
    log: "logs/{sample}/happy.log"
    threads: 4
    resources:
        runtime=120,
        mem_mb=16192
    wrapper: "file:wrapper/happy"

rule collect_insertsize_metrics:
    input:
        bam=rules.apply_bqsr.output.bam,
        bai=rules.apply_bqsr.output.bai,
    output:
        metrics=temp("output/{sample}/{sample}.insert.size.metrics.txt"),
        histogram=temp("output/{sample}/{sample}.insert.size.histogram.pdf")
    params:
        extra="M=0.5",
        jvm_args=""
    log: "logs/{sample}/collect_insertsize_metrics.log"
    resources:
        runtime=300,
        mem_mb=16192
    wrapper: "file:wrapper/picard/CollectInsertSizeMetrics"
    
rule relatedness2:
    input:
        vcf=rules.merge_vcfs.output.vcf
    output:
        relatedness2="output/{project}/{project}.relatedness2"
    params:
        prefix="{project}"
    log: "logs/{project}/relatedness2.log"
    resources:
        runtime=120,
        mem_mb=16192
    wrapper: "file:wrapper/vcftools/relatedness2"

rule mosdepth_xlsx:
    input:
        thres="output/{sample}/{sample}.mosdepth.thresholds.bed.gz",
        mos=config["ref"]["mos"]
    output:
        xlsx="output/{sample}/{sample}.mosdepth.thresholds.xlsx"
    log:
        "logs/{sample}/mosdepth.xlsx.log"
    resources:
        runtime=120,
        mem_mb=16192
    script:
        "./tools/mosdepth_to_xlsx/mosdepth.to.xlsx.R"
        
rule roh:
    input:
        vcf=rules.select_variants.output.vcf,
        gnomad=config["ref"]["gnomad-af"]
    output:
        roh="output/{sample}/{sample}.roh.txt"
    log:
        "logs/{sample}/bcftools.roh.log"
    params:
        extra="-M 1e-8"
    threads: 1
    resources:
        runtime=60,
        mem_mb=4096
    wrapper:
        "file:wrapper/bcftools/roh"

rule fastqscreen:
    input:
        unpack(get_fastq),
        config=config["fastqscreen"]
    output:
        html=[
            temp("output/{sample}/{sample}_1_screen.html"),
            temp("output/{sample}/{sample}_2_screen.html")
        ],
        txt=[
            temp("output/{sample}/{sample}_1_screen.txt"),
            temp("output/{sample}/{sample}_2_screen.txt")
        ],
        png=[
            temp("output/{sample}/{sample}_1_screen.png"),
            temp("output/{sample}/{sample}_2_screen.png")
        ]
    params:
        extra="--subset 100000 --aligner bowtie2"
    log: "logs/{sample}/{sample}.fastqscreen.log"
    threads: 4
    wrapper:
        "file:wrapper/fastqscreen"