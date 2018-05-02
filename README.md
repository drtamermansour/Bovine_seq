## start up the project 
```
## Clone the repo to your hard disk
git clone https://github.com/drtamermansour/Bovine_seq.git
cd Bovine_seq
## delete the hpcc (Do not worry, the setup script will make it again)
rm -r hpcc
## get the data folder 
mv /path/to/your/data . 
```

## Install Snakemake using conda & deploy PBS profile
```
bash snakemake_setup.sh 
```

## Run Snakemake
```
# 1. add conda and the hpcc folder to you path
export PATH=$HOME/miniconda3/bin:$(pwd)/hpcc:$PATH
# 2. turn on the environment
source activate snakemake
# 3. run by the submoit script
. hpcc/submit.sh
```
