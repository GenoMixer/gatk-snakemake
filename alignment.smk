rule ngs_cli_download:
    input: unpack(get_ngs_cli_files)
    output:
        os.path.join(work_dir, "input/{sample}_1.fastq.gz"),
        os.path.join(work_dir, "input/{sample}_2.fastq.gz")
    params:
        sample="{sample}",
        project=config["project"]
    log: "logs/{sample}/ngs_cli_download.log"
    resources:
        runtime=120
    group: "alignment"
    wrapper: "file:wrapper/ngs-cli/download"

rule bwa_map:
    input: unpack(get_fastq)
    output:
        bam=temp(os.path.join(work_dir, "output/{sample}/{sample}.sorted.bam")),
        bai=temp(os.path.join(work_dir, "output/{sample}/{sample}.sorted.bai")),
        # md5=temp("output/{sample}/{sample}.sorted.bam.md5")
    params:
        index=config["ref"]["genome"],
        extra="-M -R \"@RG\\tID:{sample}\\tSM:{sample}\\tLB:{sample}\\tPL:Illumina\"",
        sort_jvm_args="-Xmx8g"
    group: "alignment"
    threads: 8
    resources:
        mem_mb=32384, # 16GB
        runtime=1440 # 3 hours
    log: "logs/{sample}/{sample}.bwa_map.log"
    wrapper: "file:wrapper/bwa/mem"

rule mark_duplicates:
    input: 
        bam=rules.bwa_map.output.bam,
        bai=rules.bwa_map.output.bai
    output:
        bam=temp(os.path.join(work_dir, "output/{sample}/{sample}.marked.sorted.bam")),
        bai=temp(os.path.join(work_dir, "output/{sample}/{sample}.marked.sorted.bai")),
        # md5=temp("output/{sample}/{sample}.marked.sorted.bam.md5"),
        metrics=temp("output/{sample}/{sample}.marked.metrics.txt")
    params:
        jvm_args="",
        extra="VALIDATION_STRINGENCY=SILENT OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500"
    group: "alignment"
    threads: 1
    resources:
        mem_mb=32196,
        runtime=1440
    log: "logs/{sample}/{sample}.mark_duplicates.log"
    wrapper: "file:wrapper/picard/MarkDuplicates"

rule apply_bqsr:
    input:
        bam=rules.mark_duplicates.output.bam,
        bai=rules.mark_duplicates.output.bai,
        ref=config["ref"]["genome"],
        known_variants=[*config["ref"]["known-variants"], config["ref"]["dbsnp"]]
    output:
        bam="output/{sample}/{sample}.bam",
        bai="output/{sample}/{sample}.bai",
        md5="output/{sample}/{sample}.bam.md5"
    params:
        jvm_args="",
        extra=""
    group: "alignment"
    threads: 1
    resources:
        mem_mb=32384,
        runtime=1440
    wrapper: "file:wrapper/gatk/BaseRecalibrator"
