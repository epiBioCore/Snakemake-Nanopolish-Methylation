from snakemake.utils import validate
import pandas as pd

configfile: "config.yaml"

samples = pd.read_table(configfile["sample"],sep = '\t',header=None)

rule build_index:
    input:
        genome = config["genome"]

    output:
        index = "resources/index/genome.mmi"

    params:
        opts = config["minimap_index_opts"]  

    threads:
        config["threads"]     

    shell:'''
        minimap2 -t {threads} {params.opt} -d {output.index} \
        {input.genome}
        '''

rule map_reads:
    input:
        fq = lamba ## figure this out
        index = rules.build_index.output.index

    output:
        bam = "results/alignments/{sample}.bam"

    params:
        opts = config["minimap2_opts"]     

    threads:
        config["threads"]   

    log:
        "results/logs/minimap2/{sample}.log"

    shell:"""
        minimap2 -t {threads} -ax map-ont {params.opt} {input.index} {input.fq} \
        2> {log} | samtools view -bh - | samtools sort -@ {threads} -T {sample}.tmp -o out.bam - 
        samtools index output.bam
        """
rule flagstat:
    input:
        bam = rules.map_reads.output.bam

    output:
        stat = "results/alignments/{sample}_flagstat.txt"        

    shell: """
        samtools flagstat input.bam > output.stat
        """
        
rule QC:
    input:
        summary = "[samples]/sequencing_summary.txt"
        ## handle multiple directories
        bam = "[samples].bam"
    output:
        html = "results/QC/[sample].html"
        json = "resources/QC/[sample].json"

    shell: '''
        pycoQC -f input.summary -a input.bam \
        -o output.html -j output.json
        '''
