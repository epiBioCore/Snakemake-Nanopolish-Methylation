from snakemake.utils import validate
import pandas as pd

configfile: "config.yaml"

samples_df = pd.read_table(configfile["sample"],sep = '\t')
## will want to validate
samples_df = samples_df.set_index("Sample")
samples = list(samples_df.index.unique)
project_dirs = {}
for s in samples:
    project_dirs[s] = list(samples_df.loc[s,["Project_dir"]])


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

rule index_fastq:
    input:
        summary = rules.QC.input.summary
        fq = fq ## check
        f5 = {sample}.fast5 # check

    output:
        index = "results/fastq_index/{sample}"    ## check

    shell:"""
        nanopolish index -f input.summary input.f5 input.fq
     """

rule call_meth:
    input:
        genome=config["genome"]
        fq= fq ## check
        bam = rules.map_reads.output.bam ## check

    output:
        meth = "results/Methylation/{sample}_methylation_calls.tsv"    

    threads:
        config["threads"]

    shell: """
        nanopolish call-methylation -t {threads} -g input.genome -r input.fq -b input.bam >\
        output.meth
        """

rule meth_freq:
    input:
        meth = {sample}_meth ## check

    output:
        freq = results/Methylation/{sample}_frequency.tsv

    script: ##check    

    shell:"""
    calculate_methylation_frequency.py -s input.meth > output.freq
    """    





            