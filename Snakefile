from snakemake.utils import validate
import pandas as pd
import os.path
import glob

configfile: "config.yaml"

samples_df = pd.read_table(config["samples"],sep = '\t')
samples_df = samples_df.set_index("Sample")
samples = list(samples_df.index.unique())


wildcard_constraints:
    sample = "|".join(samples)
         
def get_fast5(wildcards):
      
    f5 = glob.glob(os.path.join(config["raw_data"],wildcards.sample,"2*","fast5_pass"))
    return(f5)

localrules: all,build_index

rule all:
    input: 
        expand("results/Methylation/{sample}_frequency.tsv",sample=samples),
        expand("results/alignments/{sample}_flagstat.txt",sample=samples),
        expand("results/QC/{sample}_pycoQC.html",sample=samples),
        "results/QC/multiqc_report.html"     



rule combine_tech_reps:
    input:
        fqs = lambda wildcards: glob.glob(os.path.join(config["raw_data"],"{sample}","2*","{sample}_fastq_pass.gz").format(sample=wildcards.sample))

    output:
        fq = os.path.join(config["raw_data"],"{sample}","{sample}_fastq_pass.gz")

    shell: """
        cat {input} > {output}
    """

rule build_index:
    input:
        genome = config["genome"]

    output:
        index = "resources/index/genome.mmi"

    shell:'''
        minimap2 -d {output.index} {input.genome}
        '''

rule map_reads:
    input:
        fq = os.path.join(config["raw_data"],"{sample}","{sample}_fastq_pass.gz"),
        index = "resources/index/genome.mmi",

    output:
        bam = "results/alignments/{sample}.bam"
     
    resources:
        cpus = config["resources"]["align"]["threads"], 
        time =  config["resources"]["align"]["time"],
        mem = config["resources"]["align"]["mem"]

    log:
        "results/logs/minimap2/{sample}.log"

    shell:"""
        minimap2 -t {resources.cpus} -ax map-ont {input.index} {input.fq} \
        2> {log} | samtools view -bh - | samtools sort -@ {resources.cpus} -T {wildcards.sample}.tmp -o {output.bam} - 
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
        summary = lambda wildcards: glob.glob(os.path.join(config["raw_data"],"{sample}", "2*", "sequencing_summary*.txt").format(sample=wildcards.sample)),
        bam = "results/alignments/{sample}.bam"
    output:
        html = "results/QC/{sample}_pycoQC.html",
        json = "resources/QC/{sample}_pycoQC.json"

    shell: '''
        pycoQC -f {input.summary} -a {input.bam} -o {output.html} -j {output.json}
        '''

rule collect_summaries:
    input:
        summary = lambda wildcards: glob.glob(os.path.join(config["raw_data"],"{sample}", "2*", "sequencing_summary*.txt").format(sample=wildcards.sample))

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

    log: "results/logs/{sample}_nanopolish_index.log"

    params:
        config["nanopolish_path"] + "/nanopolish"       

    resources:
        time = "4:00:00"
    run:
         fast5_dirs=""
         for d in input.f5:
             fast5_dirs += " -d " + d
         
         c = "{params} index -f {input.summary} " + fast5_dirs + " {input.fq}"
         shell(c)
    
rule call_meth:
    input:
        genome=config["genome"],
        fastq_index_readdb = os.path.join(config["raw_data"],"{sample}","{sample}_fastq_pass.gz.index.readdb"),
        fq = os.path.join(config["raw_data"],"{sample}","{sample}_fastq_pass.gz"),
        bam = "results/alignments/{sample}.bam"

    output:
        meth = "results/Methylation/{sample}_methylation_calls.tsv" 

    params:
        config["nanopolish_path"] + "/nanopolish"       

    resources:
        cpus = config["resources"]["call_meth"]["threads"],
        time = config["resources"]["call_meth"]["time"],
        mem = config["resources"]["call_meth"]["mem"]

    shell: """
        LD_LIBRARY_PATH="$CONDA_PREFIX/lib"
        {params} call-methylation -t {resources.cpus} -g {input.genome} -r {input.fq} -b {input.bam} >\
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

rule multiqc:
    input:
        expand("resources/QC/{sample}_pycoQC.json",sample=samples)

    output:
        "results/QC/multiqc_report.html"     

    shell:
        "multiqc -n {output} {input}"



            