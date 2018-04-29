# proteinClusteringPipeline
a pipeline to cluster homologous proteins into families

## How to install the pipeline
You need to install Mmseqs2, Hhblits, Snakemake and the biopython module

## how to run the pipeline
It's a two steps pipeline. A first clustering is done using the MMseqs2 software, this step defines the subfamilies. From the subfamilies, an HMM-HMM comparison is performed using the Hhblits software.

### How to create subfamilies using subfamilies.py
    ./subfamilies.py FASTA_FILENAME OUTPUT_DIRECTORY

### How to create families using hhblits.smk
    snakemake -p --snakefile hhblits.smk --configfile CONFIG_FILENAME 

## how to assess the quality of the clustering
