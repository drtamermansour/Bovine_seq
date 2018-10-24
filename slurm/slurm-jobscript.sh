#!/bin/sh
# properties = {properties}
export PATH=$HOME/miniconda3/bin:$PATH 
source activate snakemake
{exec_job}
