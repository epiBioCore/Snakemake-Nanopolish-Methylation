from snakemake.utils import validate
import pandas as pd
import os.path
import glob
import csv
from pathlib import Path

include: "Common.smk"

## for testing


configfile: "config.yaml"

samples = []
with open(config["samples"]) as f:
    reader = csv.reader(f,delimiter = "\t")
    next(reader)
    for row in reader:
        samples.append(row[0])


contrasts = {}
with open(config["comparisons"]) as f:
    reader = csv.reader(f,delimiter = "\t")
    next(reader)
    for row in reader:
        c = "%s_%s_vs_%s" % (row[2],row[0],row[1])
        contrasts[c] = row


def get_contrast(wildcards):
    return contrasts[wildcards.contrast]

def get_bioc_species_pkg(wildcards):
    """Get the package bioconductor package name for the the species in config.yaml"""
    species_letters = config["txdb_annotation"]["species"][0:2].capitalize()
    return "org.{species}.eg.db".format(species=species_letters)

def get_bioc_TxDb_pkg(wildcards):
    """Get the package bioconductor package name for the the species in config.yaml"""
    species = config["txdb_annotation"]["species"].capitalize()
    annotation = config["txdb_annotation"]["annotation"]
    build = config["txdb_annotation"]["build"]
    version = config["txdb_annotation"]["build"]

    if config["txdb_annotation"]["annotation"] == "UCSC":
        pkg = "TxDb.{species}.{annotation}.{build}.knownGene".format(species = species,annotation = annotation,build=build)
    elif config["txdb_annotation"]["annotation"] == "Ensembl":
        pkg = "EnsDb.{species}.v{version}".format(species=species,version=version) 
    return pkg       

def get_bioc_pkg_path(wildcards):
    return "resources/bioconductor_anno/lib/R/library/{pkg}".format(pkg=get_bioc_species_pkg(wildcards))

def get_bioc_txdb_pkg_path(wildcards):
    return "resources/bioconductor_txdb/lib/R/library/{pkg}".format(pkg=get_bioc_TxDb_pkg(wildcards))

localrules: all, download_bioconductor_annotation_packages,download_bioconductor_txdb_packages

rule all:
    input:
        expand("beta_distribution_colored_by_{var}.pdf",var=config["beta_distribution"]["col_by"]),
        expand("PCA_PC1vsPC2_colored_by_{var}.pdf",var=config["pca"]["col"]),
        expand("results/Differential_Methylation/{contrast}_SIG_Meth_CpGs_annotated.txt",contrast=dict.keys(contrasts)),
        expand("results/Differential_Methylation/{contrast}_SIG_Meth_regions_annotated.txt",contrast=dict.keys(contrasts)),
        expand("results/Differential_Methylation/{contrast}_ALL_CpgGs_annotated.txt",contrast=dict.keys(contrasts)),


rule bsseq_init:
    input:
        meth = expand("results/Methylation/{sample}_frequency.tsv",sample=samples),
        meta = config["metadata"]

    output:
        full = "results/R_Objects/bsseq.rda",
        filt = "results/R_Objects/bsseq_filtered.rda"

    params:
        partition="talon"
    

    script:
        "scripts/bsseq_init.R"   

rule DML_test:
    input:
        filt = "results/R_Objects/bsseq_filtered.rda"

    output:
        protected("results/R_Objects/{contrast}_DMLtest.rda"),
        
    params:
        get_contrast, 
        partition="talon"

    resources:
        cpus = config["resources"]["dml"]["threads"],
        time = config["resources"]["dml"]["time"], 
        mem = config["resources"]["dml"]["mem"]


    log:
        "results/logs/R/DML_test_{contrast}.txt"
    script:
        "scripts/DML.R"


rule callDML:
    input:
        "results/R_Objects/{contrast}_DMLtest.rda"

    output:
        dml = "results/Differential_Methylation/{contrast}_SIG_Meth_CpGs.txt",
        dmr = "results/Differential_Methylation/{contrast}_SIG_Meth_regions.txt",
        all = "results/Differential_Methylation/{contrast}_ALL_CpgGs.txt"

    params:
        partition = "talon",
        thres=config["delta"]

    
    log:
        "results/logs/R/callDML_{contrast}.txt"
    
    script:
        "scripts/callDML.R"

rule download_bioconductor_annotation_packages:
    output:
        directory("resources/bioconductor_anno/lib/R/library/{package}")

    params:
        path=lambda wc, output: Path(output[0]).parents[3],
    
    log:
        "results/logs/R/install_{package}.txt"
    shell: """
     conda create --quiet --yes -p {params.path} --channel bioconda bioconductor-{wildcards.package}
     """
rule download_bioconductor_txdb_packages:
    output:
        directory("resources/bioconductor_txdb/lib/R/library/{package}")

    params:
        path=lambda wc, output: Path(output[0]).parents[3],
    
    log:
        "results/logs/R/install_{package}.txt"
    shell: """
     conda create --quiet --yes -p {params.path} --channel bioconda bioconductor-{wildcards.package}
     """

rule annotate_DM:
    input:
        list="results/Differential_Methylation/{list}.txt",
        species_anno = get_bioc_pkg_path,
        txdb = get_bioc_txdb_pkg_path

    output:
        gene_annot_pie = "results/Differential_Methylation/{list}_gene_annotation_pie.pdf",
        list = "results/Differential_Methylation/{list}_annotated.txt",
        cpg_annot_pie = "results/Differential_Methylation/{list}_CpG_island_annotation_barchart.pdf"


    params:
        txdb = get_bioc_TxDb_pkg,
        species = get_bioc_species_pkg,
        build = config["txdb_annotation"]["build"],
        partition="talon"

    log:
        "results/logs/R/annotate_DM_{list}_out.txt"    
    conda:
        "R_env.yaml"
    script:
        "scripts/annotate_lists.R"

rule bssmooth:
    input:  
        filt = "results/R_Objects/bsseq_filtered.rda"

    output:
        s = "results/R_Objects/bsseq_smoothed.rda"  

    log:
        "results/logs/bssmooth/bssmooth.txt"      

    resources:
        cpus = config["resources"]["bssmooth"]["threads"],
        time = config["resources"]["bssmooth"]["time"],
        mem = config["resources"]["bssmooth"]["mem"],
         
    params:
        partition="talon-fat"

    script:
        "scripts/BSSmooth.R"     

rule beta:
    input:
        s = "results/R_Objects/bsseq_smoothed.rda"  

    output:
        beta= "beta_distribution_colored_by_{var}.pdf"

    params:
        vars=config["beta_distribution"]["col_by"],
        partition="talon"

    log:
        "results/logs/R/beta_distribution_by_{var}.txt"
    script:
        "scripts/beta.R"    

rule pca:
    input:
        s = "results/R_Objects/bsseq_smoothed.rda"  

    output:
        pca = "PCA_{dim}_colored_by_{var}.pdf"

    params:
        labs = config["pca"]["sample_labels"],
        col = config["pca"]["col"],
        partition="talon"

    log:
        "results/logs/R/PCA_{dim}_colored_by_{var}.txt"

    script:
        "scripts/PCA.R"

