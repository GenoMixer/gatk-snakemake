rule create_pedigree:
    input: expand("output/{project}/{project}.vcf.gz", project=config["project"])
    output: "output/{trio}/{trio}.ped"
    params:
        trio="{trio}"
    log: "logs/{trio}/create_pedigree.log"
    resources:
        runtime=5
    run:
        import os
        family = samples[samples['family_id'] == params.trio]
        if not os.path.exists(os.path.dirname(output[0])):
            os.path.makedirs(os.path.dirname(output[0]))
        with open(output[0], 'w+') as out_f:
            for idx, row in family.iterrows():
                out_f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    row.family_id,
                    row.sample_id,
                    row.father_id,
                    row.mother_id,
                    '1' if row.sex == 'M' else '2' if row.sex == 'F' else '.',
                    row.is_affected    
                ))

rule create_all_pedigree:
    input: "output/{project}/{project}.vcf.gz"
    output: "output/{project}/{project}.ped"
    log: "logs/{project}/create_all_pedigree.log"
    resources:
        runtime=5
    run:
        import os
        if not os.path.exists(os.path.dirname(output[0])):
            os.path.makedirs(os.path.dirname(output[0]))
        with open(output[0], 'w+') as out_f:
            for idx, row in samples.iterrows():
                out_f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    row.family_id,
                    row.sample_id,
                    row.father_id if not pandas.isna(row.father_id) else '',
                    row.mother_id if not pandas.isna(row.mother_id) else '',
                    '1' if row.sex == 'M' else '2' if row.sex == 'F' else '.',
                    row.is_affected    
                ))

rule gather_family:
    input:
        vcf=expand("output/{project}/{project}.vcf.gz", project=config["project"]),
        ped=rules.create_pedigree.output,
        ref=config["ref"]["genome"]
    output:
        vcf="output/{trio}/{trio}.vcf.gz",
        idx="output/{trio}/{trio}.vcf.gz.tbi",
        md5="output/{trio}/{trio}.vcf.gz.md5"
    params:
        jvm_args="",
        extra=lambda wildcards, input: "-ped " + str(input.ped) + " " + " ".join(
            map(lambda x: "-sn " + x[1].sample_id, samples[samples['family_id'] == wildcards.trio].iterrows())
        )
    log: "logs/{trio}/gather_family.log"
    resources:
        mem_mb=8192,
        runtime=60
    wrapper: "file:wrapper/gatk/SelectVariants"