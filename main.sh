## download the data of the 1st 3 patches from the hard disk
mkdir Bovine_seq/data/fastq
cd Bovine_seq/data/fastq

## download the data from the 4th patch (ftp server
cd Bovine_seq/data/fastq/5_2_18/
for f in *.gz;do
newf=$(echo $f | sed 's/_/0/') 
echo $newf;
mv $f $newf;
done


## download data from SRP072240 (Carlson et al. 2016 paper)
mkdir -p Bovine_seq2/data/fastq/SRP072240
cd Bovine_seq2/data/fastq/SRP072240
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR329/002/SRR3290632/SRR3290632_1.fastq.gz -o SRR3290632_cellline2122_L000_R1_001.fastq.gz
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR329/002/SRR3290632/SRR3290632_2.fastq.gz -o SRR3290632_cellline2122_L000_R2_001.fastq.gz

curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR329/001/SRR3290631/SRR3290631_1.fastq.gz -o SRR3290631_cellline2120_L000_R1_001.fastq.gz
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR329/001/SRR3290631/SRR3290631_2.fastq.gz -o SRR3290631_cellline2120_L000_R2_001.fastq.gz

curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR329/005/SRR3290615/SRR3290615_1.fastq.gz -o SRR3290615_RCI002_L000_R1_001.fastq.gz
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR329/005/SRR3290615/SRR3290615_2.fastq.gz -o SRR3290615_RCI002_L000_R2_001.fastq.gz

curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR329/005/SRR3290535/SRR3290535_1.fastq.gz -o SRR3290535_RCI001_L000_R1_001.fastq.gz
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR329/005/SRR3290535/SRR3290535_2.fastq.gz -o SRR3290535_RCI001_L000_R2_001.fastq.gz

## setup and run variant calling pipeline 
bash snakemake_setup.sh
. hpcc/submit.sh

## try to do the target assembly analysis
bash targetRead_Mapsembler.sh

bash snakemake_cats_setup.sh 
. hpcc_cats/submit.sh

bash targetRead_grep.sh


## generate stats
bash summary_stats.sh

