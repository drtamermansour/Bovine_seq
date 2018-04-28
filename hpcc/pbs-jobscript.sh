#!/bin/sh
# properties = {properties}
source ~/miniconda3/bin/activate snakemake 
export PATH=$HOME/miniconda3/bin:$PATH
{exec_job}
