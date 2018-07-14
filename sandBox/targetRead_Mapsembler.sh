cd ~
wget http://www.irisa.fr/symbiose/people/ppeterlongo/mapsembler2_2.2.4.zip
unzip mapsembler2_2.2.4.zip
cd mapsembler2_pipeline
./compile_all_tools.sh
##A
#1/ copy executables "tools/mapsembler2_extremities", "tools/mapsembler2_extend", "tools/kissreads_graph", "tools/kissreads" and "tools/GSVDesktop in a directory member of the PATH environment variable (e.g. /usb/local/bin)
#2/ replace PATH_RS="./tools" by PATH_RS="" in the "run_mapsembler_and_phaser.sh" configuration file
##B
#Leave it as is. In this case, if working outside this current directory ("/mnt/home/mansourt/mapsembler2_pipeline/visu"), you will have to indicate in the "run_mapsembler_and_phaser.sh" where executables "mapsembler", "phaser" and "GSV" are located by changing the value of the PATH_RS variable
sed -i 's|PATH_TOOLS="./tools/"|PATH_TOOLS="$HOME/mapsembler2_pipeline/tools/"|' run_mapsembler2_pipeline.sh

## define the program path
mapsembler=$HOME/mapsembler2_pipeline

## Create your project
work_dir=$SCRATCH/testProject  ## the folder that has your input data
mkdir $work_dir && cd $work_dir
wget wget https://de.cyverse.org/dl/d/3CE425D7-ECDE-46B8-AB7F-FAF07048AD42/samples.tar.gz
tar xvzf samples.tar.gz
rm samples.tar.gz

wget https://de.cyverse.org/dl/d/A9330898-FC54-42A5-B205-B1B2DC0E91AE/dog_chr5.fa.gz
gunzip dog_chr5.fa.gz

echo ">starter1" > starters.fa
head -n1000004 dog_chr5.fa | tail -n4 >> starters.fa
echo ">starter2" >> starters.fa
head -n1020000 dog_chr5.fa | tail -n4 >> starters.fa
echo ">starter3" >> starters.fa
head -n1040000 dog_chr5.fa | tail -n4 >> starters.fa

# Run full pipeline
$mapsembler/run_mapsembler2_pipeline.sh -s starters.fa -r "BD143_TGACCA_L005_R1_001.pe.fq.gz" -t 2 -p full_pipe -k 25 -c 2 -d 3 -g 1000000 &> log1

```
/mnt/home/mansourt/mapsembler2_pipeline/tools/\mapsembler2_extremities --k 25 --min-solid-subkmer 2 --starters starters.fa --reads "BD143_TGACCA_L005_R1_001.pe.fq.gz" --output starter_extremities.fa
/mnt/home/mansourt/mapsembler2_pipeline/tools/mapsembler2_extend starter_extremities.fa BD143_TGACCA_L005_R1_001.pe.fq.gz -t 2 -k 25 -c 2 -g 1000000 -x 40 -y 10000 -f 1 -i full_pipe -o full_pipe
/mnt/home/mansourt/mapsembler2_pipeline/tools/kissreads full_pipe\_k_25\_c_2\_t_2.fasta BD143_TGACCA_L005_R1_001.pe.fq.gz -f -o full_pipe\_coherent_k_25\_c_2\_t_2.fasta -u full_pipe\_uncoherent_k_25\_c_2\_t_2.fasta
```

## cut the input sequence into starters sequences with k interval, each starter = 2k  
cd $SCRATCH/Tamer/Bovine_seq/targetReads_Mapsembler
cat bait.fa > starters.fa
k=25
while read header; read line;do
  window=$(($k * 2));
  end=$((${#line} +1 - $window));
  for i in $(seq $k $k $end);do
    read=`expr substr $line $i $window`
    echo -e "$header.position.$i\n$read";done
done < bait.fa >> starters.fa
# use http://www.bioinformatics.org/sms2/rev_comp.html to get starters_reverse.fa
# edit the starters.fa to generate starters_kmers.fa


ln -s ../data/trimmed/12_19_17/AVEAY001A_S20_L003_R1_001.fastq.gz .
$mapsembler/run_mapsembler2_pipeline.sh -s starters.fa -r "AVEAY001A_S20_L003_R1_001.fastq.gz" -t 2 -p full_pipe -k 25 -c 2 -d 3 -g 1000000 &> log1

rm AVEAY001A_S20_L003_R1_001.fastq.gz 
cp ../data/trimmed/12_19_17/AVEAY001A_S20_L003_R1_001.fastq.gz .
gunzip AVEAY001A_S20_L003_R1_001.fastq.gz
module load jellyfish2/2.1.1
jellyfish count -m 25 -s 100M -t 10 -C AVEAY001A_S20_L003_R1_001.fastq
jellyfish query mer_counts.jf -s bait.fa > bait.count
jellyfish query mer_counts.jf -s starters_kmers.fa > starters_kmers.count


$mapsembler/tools/mapsembler2_extremities --k 25 --min-solid-subkmer 2 --starters starters.fa --reads "AVEAY001A_S20_L003_R1_001.fastq" --output starters_sub2_extremities.fa
$mapsembler/tools/mapsembler2_extremities --k 25 --min-solid-subkmer 1 --starters starters.fa --reads "AVEAY001A_S20_L003_R1_001.fastq" --output starters_sub1_extremities.fa

$mapsembler/tools/mapsembler2_extremities --k 25 --min-solid-subkmer 2 --starters starters_reverse.fa --reads "AVEAY001A_S20_L003_R1_001.fastq" --output startersR_sub2_extremities.fa
$mapsembler/tools/mapsembler2_extremities --k 25 --min-solid-subkmer 1 --starters starters_reverse.fa --reads "AVEAY001A_S20_L003_R1_001.fastq" --output startersR_sub1_extremities.fa

$mapsembler/tools/mapsembler2_extremities --k 17 --min-solid-subkmer 1 --starters starters.fa --reads "AVEAY001A_S20_L003_R1_001.fastq" --output starters17_sub1_extremities.fa

$mapsembler/tools/mapsembler2_extremities --k 25 --min-solid-subkmer 2 --starters bait.fa --reads "AVEAY001A_S20_L003_R1_001.fastq" --output bait_sub2_extremities.fa

```
## My experince with mapsembler
- mapsembler2 is in some sense comparable to spacegraphcats so maybe it is useful to share my recent experience of using it with you 
- Mapsembler should be able of targeted NGS assembly using starter sequences. It was first published in 2012 by Pierre Peterlongo and Rayan Chikhi (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3514201/). This is the group who made bcalm, Mania and GATB and also made Kissplice for reference free variant calling in RNAseq. The original implementation does not index the reads but does mapping of the starters to the read. It is subject to coverage constraints and can accommodate a bounded number of substitutions (They use Hamming distance).
- The new version (Mapsembler2; https://gatb.inria.fr/software/mapsembler/   &&  https://colibread.inria.fr/software/mapsembler2/) uses  Minia data structure  for indexing reads and the starters are only extended instead of searching for all (possibly numerous) sub-starters. 
- The README file coming with Mapsembler2 has couple errors that indicate lack of recent documentation.
- The example coming with Mapsembler2 works fine. I also tried the test data made for GWAS study in DIBSI 2017 and it seems to work. I tried to do a run using the Bovine sequence but all my trials failed to find starters extremities (kmers at the edge of starters that have enough abundance in the input reads). I did independent kmer count using jellyfish and did a query of the kmers at the starters ends and found many of them with good abundance.     
- I might try to go around this step to see if the Software will succeed to finish the target assembly task. If anybody interested, just let me know and I will keep you updated.
```


