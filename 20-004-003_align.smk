from snakemake.utils import validate
import pandas as pd
import os.path
import glob

samples = "20-004-003"

configfile: "config.yaml"


rule all:
    input:
        "results/alignments/20-004-003.bam",
        "results/alignments/20-004-003_flagstat.txt"        


rule map_reads:
    input:
        fq = os.path.join(config["raw_data"],"{sample}","{sample}_fastq_pass_fix"),
        index = "resources/index/genome.mmi",

    output:
        sam = "results/alignments/{sample}.sam"
     
    resources:
        cpus = config["resources"]["align"]["threads"], 
        time =  config["resources"]["align"]["time"],
        mem = config["resources"]["align"]["mem"]

    log:
        "results/logs/minimap2/{sample}.log"

    shell:"""
        minimap2 -t {resources.cpus} -ax map-ont {input.index} {input.fq} > {output.sam} \
        2> {log} 
        """

rule samtools_sort:
    input:
        sam = "results/alignments/{sample}.sam"

    output:
        bam =  "results/alignments/{sample}.bam"   

    resources:
        cpus = config["resources"]["align"]["threads"],
        time =  config["resources"]["align"]["time"],
        mem = config["resources"]["align"]["mem"]


    shell:"""
     samtools view -bh {input.sam} | samtools sort -@ {resources.cpus} -T {wildcards.sample}.tmp -o {output.bam} - 
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
