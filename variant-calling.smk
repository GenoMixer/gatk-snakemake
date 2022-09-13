rule partition_bam:
    input: unpack(get_input_bam),
    output:
        bins=temp(expand(
            "output/{{sample}}/bins/{i}.list",
            i=map(lambda num: "{:04d}".format(num), range(0, config['bins']))
        ))
    log: "logs/{sample}/partition_bam.log"
    threads: config['bins']
    params:
        extra="--length"
    resources:
        runtime=5
    wrapper: "file:tools/split-bam"

rule haplotype_caller_gvcf:
    input: unpack(get_hc_in)
    output:
        gvcf=temp("output/{sample}/bins/{i}.gvcf"),
        idx=temp("output/{sample}/bins/{i}.gvcf.idx")
    params:
        jvm_args="",
        extra=""
    log: "logs/{sample}/{sample}.{i}.haplotype_caller_gvcf.log"
    resources:
        mem_mb=32768,
        runtime=300
    wrapper: "file:wrapper/gatk/HaplotypeCaller"

rule combine_gvcf_parts:
    input:
        gvcfs=expand("output/{sample}/bins/{{i}}.gvcf", sample=samples.index),
        idx=expand("output/{sample}/bins/{{i}}.gvcf.idx", sample=samples.index),
        ref=config["ref"]["genome"]
    output:
        gvcf=temp(os.path.join(work_dir, "output/{project}/{i}.gvcf")),
        idx=temp(os.path.join(work_dir, "output/{project}/{i}.gvcf.idx")),
        # md5=temp(os.path.join(work_dir, "output/{project}/{i}.gvcf.md5"))
    log: "logs/{project}/combine_gvcfs.{i}.log"
    group: "joint_genotyping"
    resources:
        mem_mb=16384,
        runtime=120
    wrapper: "file:wrapper/gatk/CombineGVCFs"

rule genotype_gvcf_parts:
    input:
        gvcf=rules.combine_gvcf_parts.output.gvcf,
        idx=rules.combine_gvcf_parts.output.idx,
        ref=config["ref"]["genome"]
    output:
        vcf=temp("output/{project}/{i}.vcf.gz"),
        idx=temp("output/{project}/{i}.vcf.gz.tbi"),
        md5=temp("output/{project}/{i}.vcf.gz.md5")
    params:
        extra=(" --dbsnp {dbsnp}".format(dbsnp=config["ref"]["dbsnp"])
            + (" -L {}".format(config["ref"]["bait-region"]) if "bait-region" in config["ref"] else ""))
    log: os.path.join("logs/{project}/genotype_gvcfs.{i}.log")
    group: "joint_genotyping"
    resources:
        mem_mb=65536,
        runtime=300
    wrapper: "file:wrapper/gatk/GenotypeGVCFs"

def aggregate_vcfs(wildcards):
    return expand(
        "output/{project}/{i}.vcf.gz",
        i=map(lambda num: "{:04d}".format(num), range(0, config['bins'])),
        project=wildcards.project
    )

rule merge_vcfs:
    input: 
        vcfs=aggregate_vcfs,
        ref=config["ref"]["genome"]
    output: 
        vcf=report("output/{project}/{project}.vcf.gz", caption="report/vcf.rst", category="Variant calling"),
        idx="output/{project}/{project}.vcf.gz.tbi",
        # md5="output/{project}/{project}.vcf.gz.md5"
    params:
        jvm_args="",
        extra=""
    log: "logs/{project}/merge_vcfs.log"
    resources:
        mem_mb=8192,
        runtime=60
    wrapper: "file:wrapper/picard/MergeVcfs"

rule select_variants:
    input:
        vcf=expand("output/{project}/{project}.vcf.gz", project=config["project"]),
        ref=config["ref"]["genome"]
    output:
        vcf="output/{sample}/{sample}.vcf.gz",
        idx="output/{sample}/{sample}.vcf.gz.tbi"
    params:
        jvm_args="",
        extra=("--exclude-non-variants" 
            + (" -L {}".format(config["ref"]["bait-region"]) if "bait-region" in config["ref"] else "")
            + " --sample-name {sample}"
        )
    log: "logs/{sample}/select_variants.log"
    resources:
        mem_mb=8192,
        runtime=60
    wrapper: "file:wrapper/gatk/SelectVariants"
