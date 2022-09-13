rule create_pindel_config:
    input: expand("output/{project}/{project}.vcf.gz", project=config["project"])
    output: temp("output/{trio}/{trio}_pindel.config")
    params: trio="{trio}"
    log: "logs/{trio}/create_pindel_config.log"
    resources:
        runtime=5
    run:
        import os
        family = samples[samples['family_id'] == params.trio]
        if not os.path.exists(os.path.dirname(output[0])):
            os.path.makedirs(os.path.dirname(output[0]))
        with open(output[0], 'w+') as out_f:
            for idx, row in family.iterrows():
                out_f.write("{}\t{}\t{}\n".format(
                    os.path.join("output/", row.sample_id, row.sample_id + ".bam"),
                    "250",
                    row.sample_id
                ))

rule pindel:
    input:
        ref=config["ref"]["genome2"],
        config="output/{trio}/{trio}_pindel.config"
    output:
        bp=temp("output/{trio}/{trio}_BP"),
        close=temp("output/{trio}/{trio}_CloseEndMapped"),
        d=temp("output/{trio}/{trio}_D"),
        int=temp("output/{trio}/{trio}_INT_final"),
        inv=temp("output/{trio}/{trio}_INV"),
        li=temp("output/{trio}/{trio}_LI"),
        rp=temp("output/{trio}/{trio}_RP"),
        si=temp("output/{trio}/{trio}_SI"),
        td=temp("output/{trio}/{trio}_TD"),
        vcf=temp("output/{trio}/{trio}_pindel_initial.vcf"),
        vcf2="output/{trio}/{trio}_pindel.vcf.gz",
        tbi="output/{trio}/{trio}_pindel.vcf.gz.tbi"
    params:
        prefix=os.path.join("output/", "{trio}/{trio}"),
        extra="-e 20 -is 20 -b -G -co 10 -he 0.2 -ho 0.8"
    log:
        "logs/{trio}/pindel.log"
    threads: 6
    resources:
        runtime=1440,
        mem_mb=64384
    wrapper:
        "file:wrapper/pindel"

rule pindel_GQ:
    input:
        vcf="output/{trio}/{trio}_pindel.vcf.gz"
    output:
        vcf="output/{trio}/{trio}_pindel_w_GQ.vcf.gz",
        idx="output/{trio}/{trio}_pindel_w_GQ.vcf.gz.tbi"
    params:
        jvm_args="",
        extra=""
    log: "logs/{trio}/pindel_format.log"
    resources:
        mem_mb=8192,
        runtime=60
    run:
        import os
        import pysam
        vcf_in = pysam.VariantFile(input[0])
        vcf_in.header.formats.add('GQ','1','Integer','Genotype Quality')
        vcf_out = pysam.VariantFile(output[0], 'w', header=vcf_in.header)
        for record in vcf_in:
            for sample in record.samples:
                record.samples[sample]['GQ'] = 99
            vcf_out.write(record)
        vcf_out.close()
        pysam.tabix_index(filename=output[0], preset="vcf")

rule pindel_merge:
    input:
        vcf1="output/{trio}/{trio}.vcf.gz",
        vcf2="output/{trio}/{trio}_pindel_w_GQ.vcf.gz"
    output:
        vcf="output/{trio}/{trio}_HC_Pindel.vcf.gz",
        idx="output/{trio}/{trio}_HC_Pindel.vcf.gz.tbi"
    params:
        jvm_args="",
        extra="-a"
    log: "logs/{trio}/pindel_concat.log"
    resources:
        mem_mb=8192,
        runtime=60
    wrapper: "file:wrapper/bcftools/concat"
