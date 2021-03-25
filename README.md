# Snakemake methylation calling pipeline using Nanopolish
This pipeline will align nanopore reads to reference genome using minimap2 and call CpG methylation using Nanopolish.

## How to install pipeline:  
1. Download environment.yaml file and create conda environment.
```
wget https://raw.githubusercontent.com/epiBioCore/Snakemake-Nanopolish-Methylation/master/environment.yaml
conda env create --file environment.yaml
conda activate nanopolish
```
2. install nanopolish from github.
```
git clone --recursive https://github.com/jts/nanopolish.git
cd nanopolish
make
```

## How to run pipeline
Modify the config.yaml file to include path to reference genome, sample file, and path to raw data. This pipeline expects the raw data to be in a particular directory structure.  
- `Raw_data` - put this path in config file  
    - `Sample1`  
        - `20200106_1957_X2_FAL31326_48b7243e` - include all run folders for a sample  
            - `Sample1_fastq_pass.gz` - gzip all fastq pass files into one file
            - `fast5_pass` - directory of fast5 files
            - `sequencing_summary_FAL31236_2aaf377f.txt` - sequencing summary file   
        - `20200107_1735_X2_FAL31326_6bd8a10e`  
    - `Sample 2`
        - `20200106_1957_X1_FAL31236_9553b13e`

Pipeline will concatenate all the fastq files of technical replicates, then use combined fastq file for alignment step and methylation calling. Currently, pipeline requires that all fastq files for a run be concatenated and gzip'ed before running. 

Pipeline can be run with:  
``` 
Snakemake
 ```

If running on cluster, you must a directory for the cluster logs ex. "logs_slurm" and create a cluster config.yaml (an example of this file is the "slurm_config.yaml") file and place it in your home directory:  
```
 ~/.config/snakemake/slurm/config.yaml
 ```  

Then pipeline can be run with:  
```
 snakemake --profile slurm  
 ```
