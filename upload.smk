rule ngs_cli_upload:
    input: get_output_files
    output: temp("output/{sample}/{sample}.upload-successful")
    params:
        sample="{sample}",
        project=config["project"]
    log: "logs/{sample}/ngs_cli_upload.log"
    group: "upload"
    wrapper: "file:wrapper/ngs-cli/upload"

rule ngs_cli_upload_family:
    input: get_family_output
    output: temp("output/{trio}/{trio}.upload-successful")
    params:
        sample=get_family_index,
        project=config["project"]
    log: "logs/{trio}/ngs_cli_upload_family.log"
    wrapper: "file:wrapper/ngs-cli/upload"
