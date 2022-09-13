rule vep_pindel:
    input:
        vcf="output/{trio}/{trio}_HC_Pindel.vcf.gz",
        ped="output/{trio}/{trio}.ped"
    output:
        "output/{trio}/{trio}_HC_Pindel.vcf.gz.tri_ann.VCF",
        "output/{trio}/{trio}_HC_Pindel.vcf.gz.tri_ann.VCF_summary.html"
    params:
        clinvar=config["ref"]["clinvar"]
    conda:
        "envs/vep_101.yml"
    log:
        log="logs/{trio}/{trio}.vep.log"
    resources:
        mem_mb=16192,
        runtime=300
    shell:
        "sh 'tools/annotation/vep_annotation.sh' {input.vcf} {input.ped} {params.clinvar} {output} &> {log.log}"


rule filtering_pindel:
    input:
        vcf="output/{trio}/{trio}_HC_Pindel.vcf.gz",
        anno="output/{trio}/{trio}_HC_Pindel.vcf.gz.tri_ann.VCF",
        ped="output/{trio}/{trio}.ped"
    output:
        "output/{trio}/{trio}_HC_Pindel.vcf.gz.annotated.zip"
    params:
        clinvar=config["ref"]["clinvar"]
    conda:
        "envs/vep_101.yml"
    log:
        log="logs/{trio}/{trio}.filtering.log"
    resources:
        mem_mb=16192,
        runtime=300
    shell:
        "sh 'tools/annotation/filtering_and_output.sh' {input.vcf} {input.ped} {params.clinvar} {output} &> {log.log}"
