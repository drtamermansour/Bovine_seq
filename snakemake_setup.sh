############################################################
## SnakeMake with conda and generic snakemake profiles
######################################################

work_dir=$1  ## This is the root directory of your project where you have the Snakefile

## Install Snakemake using miniconda    ## http://snakemake.readthedocs.io/en/stable/getting_started/installation.html
# install miniconda  ## https://github.com/taylorreiter/olive_genome
cd ~
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh  # Note: accept adding to .bashrc but you still need to source or restart or you can export the path
export PATH=$HOME/miniconda3/bin:$PATH
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda create -n snakemake python==3.6 snakemake=4.3.0


## Deploy pbs-torque profile  ## https://github.com/Snakemake-Profiles/pbs-torque
source ~/miniconda3/bin/activate snakemake
pip install cookiecutter        ## https://cookiecutter.readthedocs.io/en/latest/installation.html
## This is to make a global profile
# mkdir -p ~/.config/snakemake
# cd ~/.config/snakemake
# cookiecutter https://github.com/Snakemake-Profiles/pbs-torque.git
## This is to make project specific pbs-torque profile
cd $work_dir
cookiecutter https://github.com/Snakemake-Profiles/pbs-torque.git ## use "hpcc" for profile_name
cd hpcc
## update the jobscript
sed -i '3i source ~/miniconda3/bin/activate snakemake \
export PATH=$HOME/miniconda3/bin:$PATH' pbs-jobscript.sh
## create bash script to run your project
echo -e '
logdir=hpcc/log
mkdir -p $logdir \n
QSUB="qsub -l nodes=1:ppn={threads} "
QSUB="$QSUB -o $logdir -e $logdir -A ged" \n
snakemake                               \
    -j 100                              \
    --profile ./hpcc                    \
    --rerun-incomplete                  \
    --keep-going                        \
    --latency-wait 10                   \
    --max-jobs-per-second 1             \
    --use-conda                         \
    --cluster "$QSUB" $@' > submit.sh

cd ../ && chmod u+x hpcc/*

