include: "main.smk"
include: "alignment.smk"
include: "upload.smk"
include: "variant-calling.smk"
include: "trio-calling.smk"
include: "qc.smk"
include: "vep.smk"
include: "pindel.smk"

families = samples.groupby('family_id').groups.keys()

rule all:
    input:
        expand("output/{sample}/{sample}.bam", sample=samples.index),
        expand("output/{sample}/{sample}.bai", sample=samples.index),
        expand("output/{sample}/{sample}.vcf.gz", sample=samples.index),
        expand("output/{sample}/{sample}.vcf.gz.tbi", sample=samples.index),
        expand("output/{sample}/{sample}.eval.txt",sample=samples.index),
        expand("output/{sample}/{sample}.roh.txt", sample=samples.index),
        expand("output/{trio}/{trio}.vcf.gz", trio=families),
        expand("output/{trio}/{trio}.vcf.gz.tbi", trio=families),
        expand("output/{project}/{project}.html", project=config["project"]),
        expand("output/{project}/{project}.relatedness2", project=config["project"]),
        expand("output/{sample}/{sample}.insert.size.metrics.txt", sample=samples.index),
        expand("output/{sample}/{sample}.insert.size.histogram.pdf", sample=samples.index),
        expand("output/{sample}/{sample}.mosdepth.thresholds.xlsx", sample=samples.index),
        expand("output/{trio}/{trio}_pindel.vcf.gz", trio=families),
        expand("output/{trio}/{trio}_pindel.vcf.gz.tbi", trio=families),
        expand("output/{trio}/{trio}_HC_Pindel.vcf.gz", trio=families),
        expand("output/{trio}/{trio}_HC_Pindel.vcf.gz.tbi", trio=families),
        expand("output/{trio}/{trio}_HC_Pindel.vcf.gz.annotated.zip", trio=families)
