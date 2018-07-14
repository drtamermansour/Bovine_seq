## install miniconda for Spacecatgraph & snakemake
cd ~
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh  # Note: accept adding to .bashrc but you still need to source or restart or you can export the path
export PATH=$HOME/miniconda3/bin:$PATH
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda create -n catsenv python==3.6 snakemake==4.8.1

## startup your system (re-run each time you want to run the tools)
module unload Python
module swap GNU/4.4.5 GNU/6.2
module load CMake/3.1.0
source activate catsenv

## Installation of BCALM
cd ~/src/
git clone --recursive https://github.com/GATB/bcalm
cd bcalm
mkdir build;  cd build;  cmake ..;  make -j 8

## install Spacecatgraph
git clone https://github.com/spacegraphcats/spacegraphcats/ -b mphf
cd spacegraphcats
pip install Cython
pip install -r requirements.txt

## Deploy pbs-torque profile  ## https://github.com/Snakemake-Profiles/pbs-torque
pip install cookiecutter      ## https://cookiecutter.readthedocs.io/en/latest/installation.html

work_dir=$(pwd)  ## This is the root directory of your project where you have the Snakefile
cd $work_dir
cookiecutter https://github.com/Snakemake-Profiles/pbs-torque.git ## use "hpcc_cats" for profile_name
```
you might be asked one or more of the questions:
Is it okay to delete and re-download it? [yes]: no
Do you want to re-use the existing version? [yes]: yes
profile_name [pbs-torque]: hpcc_cats
immediate_submit [False]: 
```
cd hpcc_cats
## update the jobscript
sed -i 's|#!/bin/sh|#!/bin/bash|' pbs-jobscript.sh
sed -i '3i module unload Python \
module load GNU/6.2 \
module load CMake/3.1.0 \
export PATH=$HOME/miniconda3/bin:$HOME/spacegraphcats:$PATH \
source activate catsenv' pbs-jobscript.sh
## update the config file
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
    --snakefile Snakefile_cats          \
    --use-conda                         \
    --profile ./hpcc_cats' > submit.sh

cd ../ && chmod u+x hpcc_cats/*



