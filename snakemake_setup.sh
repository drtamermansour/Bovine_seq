############################################################
## SnakeMake with conda and generic snakemake profiles
######################################################
## Install Snakemake using miniconda    ## http://snakemake.readthedocs.io/en/stable/getting_started/installation.html
# install miniconda  ## https://github.com/taylorreiter/olive_genome
cd ~
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh  # Note: accept adding to .bashrc but you need to source .bashrc or restart or you can export the path #version is conda 4.4.10
export PATH=$HOME/miniconda3/bin:$PATH
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda create -n snakemake python==3.6 snakemake==4.8.1  ## I am not sure why using "=="

## Deploy pbs-torque profile  ## https://github.com/Snakemake-Profiles/pbs-torque
source ~/miniconda3/bin/activate snakemake
pip install cookiecutter        ## https://cookiecutter.readthedocs.io/en/latest/installation.html
## This is to make a global profile
# mkdir -p ~/.config/snakemake
# cd ~/.config/snakemake
# cookiecutter https://github.com/Snakemake-Profiles/pbs-torque.git
## This is to make project specific pbs-torque profile
work_dir=$(pwd)  ## This is the root directory of your project where you have the Snakefile
cd $work_dir
cookiecutter https://github.com/Snakemake-Profiles/pbs-torque.git ## use "hpcc" for profile_name
cd hpcc
## update the jobscript
# loading module in jobscript: https://bitbucket.org/snakemake/snakemake/issues/155/feature-request-integration-with-the
#sed -i 's|#!/bin/sh|#!/bin/bash|' pbs-jobscript.sh
sed -i '3i export PATH=$HOME/miniconda3/bin:$PATH \
source activate snakemake' pbs-jobscript.sh
## update the config file ## this config file is saved as config_backup.yaml. It will submit the jobs, keep track and resubmit on failure
## currently I am using the default config (just change the true to True
sed -i 's/cluster:.*/cluster: "pbs-submit.py "/' config.yaml 
sed -i 's/true/True/' config.yaml 
echo -e 'rerun-incomplete: True
keep-going: True
latency-wait: 10
max-jobs-per-second: 1
restart-times: 2' >> config.yaml
## create bash script to run your project
echo -e '
snakemake                               \
    --use-conda                         \
    --profile ./hpcc' > submit.sh

cd ../ && chmod u+x hpcc/*


############
cookiecutter https://github.com/Snakemake-Profiles/slurm.git ## keep it slurm
cd slurm
## update the jobscript
# loading module in jobscript: https://bitbucket.org/snakemake/snakemake/issues/155/feature-request-integration-with-the
#sed -i 's|#!/bin/sh|#!/bin/bash|' slurm-jobscript.sh
sed -i '3i export PATH=$HOME/miniconda3/bin:$PATH \
source activate snakemake' slurm-jobscript.sh
## update the config file ## this config file is saved as config_backup.yaml. It will submit the jobs, keep track and resubmit on failure
## currently I am using the default config (just change the true to True
mv config.yaml config.yaml_backup
cp ../hpcc/config.yaml .
sed -i 's/pbs/slurm/g' config.yaml
sed -i 's/--depend/--dependency/' config.yaml
## create bash script to run your project
echo -e '
snakemake                               \
    --use-conda                         \
    --profile ./slurm' > submit-slurm.sh

cd ../ && chmod u+x slurm/*

echo -e '
snakemake -np                           \
    --use-conda                         \
    --profile ./slurm' > submit-slurm_dry.sh    

cd ../ && chmod u+x slurm/*

