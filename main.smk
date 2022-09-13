import pandas
from snakemake.utils import validate
from snakemake.exceptions import IncompleteCheckpointException
from snakemake.io import checkpoint_target
import os
import sys

##
# Report
report: "report/workflow.rst"

##
# Load configuration and sample sheet
#
configfile: 'config.yaml'
validate(config, schema='schemas/config.schema.yaml')

samples = pandas.read_csv(
    config['sample-sheet'], sep='\t'
).set_index('sample_id', drop=False)
validate(samples, schema='schemas/samples.schema.yaml')

##
# Setup scratch dir handling
#
work_dir = os.getenv("SCRATCH_DIR", os.getcwd()) # workdir is a reserved keyword
print('Work_dir:\t{}'.format(work_dir), file=sys.stderr)

##
# Set wildcard_constraints identifying samples and trios
#
wildcard_constraints:
    sample="|".join(samples.index),
    trio="|".join(samples.groupby('family_id').groups.keys()),
    project=config["project"],
    i="|".join(map(lambda num: "{:04d}".format(num), range(0, config['bins'])))

##
# Helper functions to find input files
#
def get_ngs_cli_files(wildcards):
    data = samples[
        (samples['sample_id'] == wildcards.sample)
        & (samples['remote'] == 'ngs-cli')
    ].loc[:, ['fastq1', 'fastq2']].dropna()
    if not os.path.exists(os.path.join(work_dir, 'input', '.mock')):
        os.makedirs(os.path.join(work_dir, 'input', '.mock'))
    if len(data) == 2:
        path = os.path.join(work_dir, 'input', '.mock', data.fastq1.values[0])
        if not os.path.exists(path):
            open(path, 'a').close()
        return {
            'fq1': path
        }
    else:
        path1 = os.path.join(work_dir, 'input', '.mock', data.fastq1.values[0])
        if not os.path.exists(path1):
            open(path1, 'a').close()
        path2 = os.path.join(work_dir, 'input', '.mock', data.fastq2.values[0])
        if not os.path.exists(path2):
            open(path2, 'a').close()
        return {
            'fq1': path1,
            'fq2': path2
        }

def get_fastq(wildcards):
    fastqs = samples[
        samples['sample_id'] == wildcards.sample
    ].loc[:, ['remote', 'fastq1', 'fastq2']].dropna(axis=1)
    reads = list()
    for _, row in fastqs.iterrows():
        if row.remote == 'ngs-cli':
            reads.append(os.path.join(work_dir, 'input', 
                                      wildcards.sample + '_1.fastq.gz'))
            if len(*fastqs.values) > 2:
                reads.append(os.path.join(work_dir, 'input',
                                          wildcards.sample + '_2.fastq.gz'))
        else:
            reads.append(os.path.join(work_dir, 'input',
                                      fastqs.fastq1.values[0]))
            if len(*fastqs.values) > 2:
                reads.append(os.path.join(work_dir, 'input',
                                          fastqs.fastq2.values[0]))
    return {
        'reads': reads
    }


def get_output_files(wildcards):
    return list(map(
        lambda ext: '{}/{}/{}.{}'.format(
            "output",
            wildcards.sample,
            wildcards.sample,
            ext
        ),
        config['output-extensions']['single']
    ))

def get_family_output(wildcards):
    return list(map(
        lambda ext: '{}/{}/{}.{}'.format(
            "output",
            wildcards.trio,
            wildcards.trio,
            ext
        ),
        config['output-extensions']['family']
    ))

def get_family_index(wildcards):
    family = samples[samples['family_id'] == wildcards.trio]
    if len(family) == 1:
        print(family, family[0].sample_id)
        return family[0].sample_id
    children = family[(family['mother_id'] | family['father_id'])]
    print(children)
    if len(children) > 0:
        return children[0].sample_id
    else:
        raise Exception(
            "No children found for family {}".format(wildcards.trio)
        )


def get_input_bam(wildcards):
    return {
        'bam': 'output/{0}/{0}.bam'.format(wildcards.sample),
        'bai': 'output/{0}/{0}.bai'.format(wildcards.sample)
    }

def get_hc_in(wildcards):
    interval = "output/{0}/bins/{1}.list".format(wildcards.sample, wildcards.i)
    if "bait-region" in config["ref"] and config["ref"]["bait-region"]:
        interval = [interval, config["ref"]["bait-region"]]
    return {
        **get_input_bam(wildcards),
        "interval": interval,
        "ref": config["ref"]["genome"],
        "dbsnp": config["ref"]["dbsnp"]
    }

def gvcf_files(wildcards):
    return expand("output/{sample}/{sample}.gvcf", sample=samples.index)

def get_qc_input(wildcards):
    ##
    # Expected QC:
    #   fastqc          [sample]            _{1,2}.fastqc.zip
    #   peddy           [project]           {het,ped,sex}_check.csv, peddy.ped
    #   hsMetrics       [sample]            hsmetrics.txt
    #   markDuplicates  [sample]            marked.metrics.txt
    #   stats           [project,sample]    stats.txt
    #   verifybamid     [sample]            selfSM
    #   mosdepth        [sample]            mosdepth.dist.txt, mosdepth.summary.txt
    #   insertsize      [sample]            insert.size.metrics.txt
    #   fastqscreen     [sample]            screen.txt
    return [
        *expand(
            "output/{sample}/{sample}.{ending}",
            sample=samples.index,
            ending=[
                "hsmetrics.txt",
                "selfSM",
                "stats.txt",
                "marked.metrics.txt",
                "mosdepth.dist.txt",
                "mosdepth.summary.txt",
                "eval.txt",
                "insert.size.metrics.txt"
            ]
        ),
        *expand(
            "output/{project}/{project}.{ending}",
            project=config["project"],
            ending=[
                "peddy.ped",
                "het_check.csv",
                "ped_check.csv",
                "sex_check.csv",
                "background_pca.json",
                "stats.txt",
                "relatedness2"
            ]
        ),
        *expand(
            "output/{sample}/{sample}_{ending}",
            sample=samples.index,
            ending=[
                "1_fastqc.zip",
                "2_fastqc.zip"
            ]
        )
    ]
