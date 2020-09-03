from snakemake.utils import validate
import pandas as pd
import os.path
import glob

configfile: "config.yaml"

samples_df = pd.read_table(config["samples"],sep = '\t')
## will want to validate
samples_df = samples_df.set_index("Sample")
samples = list(samples_df.index.unique())
runs = list(samples_df.loc[:,"Project_dir"])
project_dirs = {}
for s in samples:
    project_dirs[s] = list(samples_df.loc[s,"Project_dir"])


wildcard_constraints:
    sample = "|".join(samples)
         
def get_fast5(wildcards):
      
    f5 = glob.glob(os.path.join(config["raw_data"],wildcards.sample,"2*","fast5_pass"))
    return(f5)

rule all:
    input: 
        expand("results/Methylation/{sample}_frequency.tsv",sample=samples)


rule combine_tech_reps:
    input:
        fqs = lambda wildcards: glob.glob(os.path.join(config["raw_data"],"{sample}","2*","{sample}_fastq_pass.gz").format(sample=wildcards.sample))

    output:
        fq = os.path.join(config["raw_data"],"{sample}","{sample}_fastq_pass.gz")

    shell: """
        zcat {input} > {output}
    """

rule build_index:
    input:
        genome = config["genome"]

    output:
        index = "resources/index/genome.mmi"

    threads: config["resources"]["index_genome"]["threads"]
         

    shell:'''
        minimap2 -t {threads} -d {output.index} {input.genome}
        '''

rule map_reads:
    input:
        fq = os.path.join(config["raw_data"],"{sample}","{sample}_fastq_pass.gz"),
        index = "resources/index/genome.mmi",

    output:
        bam = "results/alignments/{sample}.bam"
     
    params:
        map = config["resources"]["align"]["threads"], 
        sam =  config["resources"]["samtools_sort"]["threads"]

    log:
        "results/logs/minimap2/{sample}.log"

    shell:"""
        minimap2 -t {params.map} -ax map-ont {input.index} {input.fq} \
        2> {log} | samtools view -bh - | samtools sort -@ {params.sam} -T {wildcards.sample}.tmp -o {output.bam} - 
        samtools index {output.bam}
        """
rule flagstat:
    input:
        bam = 'results/alignments/{sample}.bam'

    output:
        stat = "results/alignments/{sample}_flagstat.txt"        

    shell: """
        samtools flagstat {input.bam} > {output.stat}
        """
        
rule QC:
    input:
        summary = lambda wildcards: glob.glob(os.path.join(config["raw_data"],"{sample}", "2*", "sequencing_summary.txt").format(sample=wildcards.sample)),
        bam = "results/alignments/{sample}.bam"
    output:
        html = "results/QC/{sample}.html",
        json = "resources/QC/{sample}.json"

    shell: '''
        pycoQC -f {input.summary} -a {input.bam} \
        -o {output.html} -j {output.json}
        '''

rule collect_summaries:
    input:
        summary = lambda wildcards: glob.glob(os.path.join(config["raw_data"],"{sample}", "2*", "sequencing_summary.txt").format(sample=wildcards.sample))

    output:
        out = os.path.join(config["raw_data"],"{sample}","{sample}_sequencing_summaries.txt")

    run:
        with open(output.out,'w') as out:
            for item in input:
                 out.write("%s\n" % os.path.abspath(item))


rule index_fastq:
    input:
        summary = os.path.join(config["raw_data"],"{sample}","{sample}_sequencing_summaries.txt"),
        fq = os.path.join(config["raw_data"],"{sample}","{sample}_fastq_pass.gz"),
        f5 = get_fast5 

    output:
        fastq_index = os.path.join(config["raw_data"],"{sample}","{sample}_fastq_pass.gz.index"),
        fastq_index_fai = os.path.join(config["raw_data"],"{sample}","{sample}_fastq_pass.gz.index.fai"),
        fastq_index_gzi = os.path.join(config["raw_data"],"{sample}","{sample}_fastq_pass.gz.index.gzi"),
        fastq_index_readdb = os.path.join(config["raw_data"],"{sample}","{sample}_fastq_pass.gz.index.readdb")

    run:
         fast5_dirs=""
         for d in input.f5:
             fast5_dirs += "-d " + d
         
         c = "nanopolish index -f {input.summary} " + "fast5_dirs" + " {input.fq}"
         shell(c)
    

rule call_meth:
    input:
        genome=config["genome"],
        fq = os.path.join(config["raw_data"],"{sample}","{sample}_fastq_pass.gz"),
        bam = "results/alignments/{sample}.bam"

    output:
        meth = "results/Methylation/{sample}_methylation_calls.tsv"    

    threads:
        config["resources"]["call_meth"]["threads"]

    shell: """
        nanopolish call-methylation -t {threads} -g {input.genome} -r {input.fq} -b {input.bam} >\
        {output.meth}
        """

rule meth_freq:
    input:
        meth = "results/Methylation/{sample}_methylation_calls.tsv"    

    output:
        freq = "results/Methylation/{sample}_frequency.tsv"


    shell:"""
    scripts/calculate_methylation_frequency.py -s {input.meth} > {output.freq}
    """    





            