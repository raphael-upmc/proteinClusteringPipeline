# proteinClusteringPipeline
a pipeline to cluster homologous proteins into families

## How to install the pipeline
You need to install Mmseqs2, MCL, HH-pred suite and the biopython module

## how to run the pipeline
It's a two steps pipeline. A first clustering is done using the MMseqs2 software, this step defines the subfamilies. From the subfamilies, an HMM-HMM comparison is performed using the Hhblits software from the HHpred suite.

### How to create subfamilies using subfamilies.py
    ./subfamilies.py FASTA_FILENAME --output-directory OUTPUT_DIRECTORY

### How to create families 
    ./hhblits.py CONFIG_FILENAME  # creating the HMMs and perfoming the HMM-HMM comparison
    ./runningMclClustering.py CONFIG_FILENAME # running MCL to create the families

## how to assess the quality of the clustering
