# ADR-2_localization
This repository contains scripts to process FASTQ files as described in Eliad et al., 2023. 
To run the scripts, activating our “py3” environment is essential.
First, create the environment from the py3_environment.yml file:

`conda env create -f py3_environment.yml`

Second, activate the new environment:

`conda activate py3`

While `/Processing/fastq_trimming.py` enables trimming low-quality edges of reads in a FASTQ file, `/Processing/fastq_collapse.py` enables merging identical reads originated in PCR duplications. 

To understand how to run the scripts, please run the following:

`python /Processing/fastq_trimming.py –-help`

`python /Processing/fastq_trimming.py example`

`python /Processing/fastq_trimming.py param`

`python /Processing/fastq_collapse.py –-help`

`python /Processing/fastq_collapse.py example`

`python /Processing/fastq_collapse.py param`
