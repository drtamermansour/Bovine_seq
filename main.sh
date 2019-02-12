## Data downlowd 
## Raw data was organizied by batches in Bovine_seq/data/fastq directory. The directory structure of raw data in Bovine_seq/sandBox/rawData_dir.txt
## Raw data can be downloaded from SRA. SRA info can be found in Bovine_seq/sandBox/metadata-4898878-processed-ok.tsv

## download data from SRP072240 (Carlson et al. 2016 paper)
## Whole genome sequencing of two edited calves along with those of their progenitor cell lines, 2122 and 2120
mkdir -p Bovine_seq/external/SRP072240/data/fastq/SRP072240
cd Bovine_seq/external/SRP072240/data/fastq/SRP072240
## parental cellline_2122, isolate: Holstein bull
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR329/002/SRR3290632/SRR3290632_1.fastq.gz -o SRR3290632_cellline2122_L000_R1_001.fastq.gz
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR329/002/SRR3290632/SRR3290632_2.fastq.gz -o SRR3290632_cellline2122_L000_R2_001.fastq.gz
## parental cellline_2120, isolate: Holstein bull
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR329/001/SRR3290631/SRR3290631_1.fastq.gz -o SRR3290631_cellline2120_L000_R1_001.fastq.gz
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR329/001/SRR3290631/SRR3290631_2.fastq.gz -o SRR3290631_cellline2120_L000_R2_001.fastq.gz
## RCI_002, isolate: Cloned bull (Buri)
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR329/005/SRR3290615/SRR3290615_1.fastq.gz -o SRR3290615_RCI002_L000_R1_001.fastq.gz
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR329/005/SRR3290615/SRR3290615_2.fastq.gz -o SRR3290615_RCI002_L000_R2_001.fastq.gz
## RCI_001, isolate: Cloned bull
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR329/005/SRR3290535/SRR3290535_1.fastq.gz -o SRR3290535_RCI001_L000_R1_001.fastq.gz
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR329/005/SRR3290535/SRR3290535_2.fastq.gz -o SRR3290535_RCI001_L000_R2_001.fastq.gz

## downlaod the software you need on the server
## multiQC (for Slurm configuration)
conda create -n multiQC multiqc==1.6
## GATK
# https://software.broadinstitute.org/gatk/documentation/quickstart
# 1. Requirements
module load Java/jdk1.8.0 ## java -version ## java version "1.8.0_102"
# 2. Get GATK
# 2A. download: https://software.broadinstitute.org/gatk/download/
wget https://github.com/broadinstitute/gatk/releases/download/4.0.9.0/gatk-4.0.9.0.zip
unzip gatk-4.0.9.0.zip
cd gatk-4.0.9.0
# 2B. Add the path to the .bashrc file 
echo "export PATH=$PATH:$(pwd)/gatk-4.0.9.0/gatk" >> $HOME/.bashrc
source $HOME/.bashrc
# 2B. run the conda env (I have installed miniconda already & I have miniconda3/bin in the PATH)
# https://software.broadinstitute.org/gatk/documentation/article?id=12836
conda update conda  ## conda --version ## conda 4.5.11
conda env create -n gatk -f gatkcondaenv.yml

###############
## Install Snakemake using conda & deploy PBS or slurm profile
bash snakemake_setup.sh

### Run Snakemake on hpc
## 1. add conda and the hpcc folder to you path
#export PATH=$HOME/miniconda3/bin:$(pwd)/hpcc:$PATH
## 2. turn on the environment
#source activate snakemake
## 3. run by the submoit script
#. hpcc/submit.sh

## Run Snakemake on slurm
# 1. add conda and the hpcc folder to you path
export PATH=$HOME/miniconda3/bin:$(pwd)/slurm:$PATH
# 2. turn on the environment
source activate snakemake
# 3. run by the submoit script
. slurm/submit-slurm.sh 

###############
## continue after Snakemake analysis
###############

## statistics per sample
#module load tabix/0.2.6
#bgzip -c vc/hapCaller_raw_withExternal.pass_SNPs.vcf > vc/hapCaller_raw_withExternal.pass_SNPs.vcf.gz
#tabix -p vcf vc/hapCaller_raw_withExternal.pass_SNPs.vcf.gz
#module load vcftools/0.1.14
#vcf-stats vc/hapCaller_raw_withExternal.pass_SNPs.vcf.gz

## Calc the SNP and INDEL rates
mkdir -p var_stat
#i=10
for var in SNP INDEL;do #echo $vcf;done
  vcf=vc/hapCaller_pass.${var}s.vcf
  for i in {10..44};do
    id=$(head -n10000 $vcf | grep "^#CHROM" | awk -F'\t' -v i=$i '{print $i}');
    echo $id.$var.stat
    echo "genotypes $id" > var_stat/$id.$var.stat
    cat $vcf | awk -F '\t' -v i=$i '/^#/{next;}{print $i;}' | awk -F ':' '{A[$1]++} END {for(x in A) print x,A[x]}' >> var_stat/$id.$var.stat
done;done;

first=0;
for var in SNP INDEL;do #echo $vcf;done
  for f in var_stat/*.$var.stat;do
    echo $f
    if [ "$first" -eq 0 ];then 
      inputA=$f;header=$(head -n1 $f);first=1;echo "done if;"
    else
      echo "start else"
      inputB=$f;
      header=$header" "$(head -n1 $f | cut -d" " -f2)
      tail -n+2 $inputA | sort -k 1b,1 | join - <(tail -n+2 $inputB | sort -k 1b,1) > temp
      echo $header > temp2
      cat temp >> temp2
      inputA="temp2";fi
  done
  first=0;
  mv temp2 var_stat/$var.stat
done &> log

## merge SNPs and psudo indels
module load VCFtools/0.1.15-Perl-5.26.1 #module load vcftools/0.1.14
vcf-concat vc/hapCaller_pass.SNPs.vcf vc/hapCaller_pass_INDELs.monoAllel_edit.vcf | vcf-sort > vc/allSnp.vcf
## genome annotation map
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/263/795/GCA_002263795.2_ARS-UCD1.2/GCA_002263795.2_ARS-UCD1.2_assembly_report.txt -O refGenome/assembly_report.txt
grep -v "^#" refGenome/assembly_report.txt | awk '{if($2=="assembled-molecule")print $5,$1}' > refGenome/assembly_map
module load bcftools/1.9.64 #module load bcftools/1.2
bcftools annotate --rename-chrs refGenome/assembly_map vc/allSnp.vcf > vc/allSnp_chr.vcf
## prepare Plink input
module load VCFtools/0.1.15-Perl-5.26.1 #module load vcftools/0.1.14
mkdir -p plink
vcftools --vcf vc/allSnp_chr.vcf --plink --out plink/allSnp 


## create binary inputs
module load PLINK/1.9b_4.1-x86_64  #module load plink/1.9
plink --file plink/allSnp --allow-no-sex --cow --make-bed --out plink/allSnp.binary #14084653 variants,35 cattle//genotyping rate 0.996337.
## create file of alternative alleles 
#cat allSnp.vcf | awk 'BEGIN{FS="\t";OFS="\t";}/#/{next;}{{if($3==".")$3=$1":"$2;}print $3,$5;}'  > alt_alleles

## pruning
plink --bfile plink/allSnp.binary --cow --maf 0.05 --geno 0.05 --indep 50 5 2 --out plink/allSnp_filter # 258044 variants removed due to missing genotype data // 2419167 variants removed due to minor allele threshold // 11407442 variants and 35 cattle pass filters and QC. // Pruning complete.  10686280 of 11407442 variants removed
plink --bfile plink/allSnp.binary --cow --extract plink/allSnp_filter.prune.in --out plink/pruned_allSnp --recode #721162 variants remaining // genotyping rate is 0.996072
plink --file plink/pruned_allSnp --allow-no-sex --cow --make-bed --out plink/pruned_allSnp.binary

## IBS similarity matrix (typically used for population stratification): distances are expressed in "shared allele counts", "proportions of alleles IBS", "1 minus the identity-by-state value" 
plink --bfile plink/pruned_allSnp.binary --distance square allele-ct ibs 1-ibs --cow --allow-no-sex --out plink/distance ## Excluding 25414 variants on non-autosomes
awk 'FNR==NR{a[$1]=$0;next}{ print a[$1]}' myIDs plink/distance.mdist.id > plink/distance.myIDs
## R (using the 1 minus the identity-by-state value distances)
## distance.myIDs is a version of "distance.mdist.id" with more columns for alt ids
Rscript -e 'args=(commandArgs(TRUE));'\
'library(ape)'\
'setwd("plink")'\
'mdist_id=read.table(args[1])'\
'mdist_table=read.table(args[2], fill=T, col.names=mdist_id$V3)'\
'mdist=as.dist(mdist_table)'\
'hc = hclust(mdist)'\
'mypal = c("red", "blue", "green", "gold", "black", "blueviolet","coral4","cadetblue", "darkred", "darkslategray", "darkorange3", "firebrick")'\
'clus = cutree(hc, 12)'\
'png(file="fan600.png",width=6500,height=6500,res=600, type="cairo")'\
'plot(as.phylo(hc), type = "fan", cex = 1, font=2, use.edge.length = TRUE, tip.color = mypal[clus], label.offset = 0.002)'\
'dev.off()' distance.myIDs distance.mdist

## IBS similarity matrix using all variants
plink --bfile plink/allSnp.binary --distance square allele-ct ibs 1-ibs --cow --allow-no-sex --out plink/distanceAll ## Excluding 387045 variants on non-autosomes
## R (using the 1 minus the identity-by-state value distances)
Rscript -e 'args=(commandArgs(TRUE));'\
'library(ape)'\
'setwd("plink")'\
'mdist_id=read.table(args[1])'\
'mdist_table=read.table(args[2], fill=T, col.names=mdist_id$V3)'\
'mdist=as.dist(mdist_table)'\
'hc = hclust(mdist)'\ 
'mypal = c("red", "blue", "green", "gold", "black", "blueviolet","coral4","cadetblue", "darkred", "darkslategray", "darkorange3", "firebrick")'\
'clus = cutree(hc, 12)'\ 
'png(file="fan600_All.png",width=6000,height=6000,res=600, type="cairo")'\
'plot(as.phylo(hc), type = "fan", cex = 1, font=2, use.edge.length = TRUE, tip.color = mypal[clus], label.offset = 0.002)'\
'dev.off()' distance.myIDs distanceAll.mdist


## statistics of Mendel errors in trios using pruned marker
cp plink/pruned_allSnp.binary.bed plink/pruned_allSnp_forMendel.binary.bed
cp plink/pruned_allSnp.binary.bim plink/pruned_allSnp_forMendel.binary.bim
cp plink/pruned_allSnp.binary.fam plink/tempFam
awk '{if($1=="AVEAY012A")print "Trio1",$2,"AVEAY0010B","AVEAY013A",$5,$6;else print $0;}' plink/tempFam > plink/tempFam1
awk '{if($1=="AVEAY012B")print "Trio1",$2,"AVEAY0010B","AVEAY013B",$5,$6;else print $0;}' plink/tempFam1 > plink/tempFam2
awk '{if($1=="AVEAY14A")print "Trio1",$2,"AVEAY0010B","AVEAY14B",$5,$6;else print $0;}' plink/tempFam2 > plink/tempFam3
awk '{if($1=="AVEAY0010B")print "Trio1",$2,$3,$4,"1",$6;else print $0;}' plink/tempFam3 > plink/tempFam4
awk '{if($1=="AVEAY013A")print "Trio1",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam4 > plink/tempFam5
awk '{if($1=="AVEAY013B")print "Trio1",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam5 > plink/tempFam6
awk '{if($1=="AVEAY14B")print "Trio1",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam6 > plink/tempFam

awk '{if($1=="AVEAY001B")print "Trio2",$2,"AVEAY001A","AVEAY004B",$5,$6;else print $0;}' plink/tempFam > plink/tempFam1
awk '{if($1=="AVEAY002A")print "Trio2",$2,"AVEAY001A","AVEAY006A",$5,$6;else print $0;}' plink/tempFam1 > plink/tempFam2
awk '{if($1=="AVEAY002B")print "Trio2",$2,"AVEAY001A","AVEAY005A",$5,$6;else print $0;}' plink/tempFam2 > plink/tempFam3
awk '{if($1=="AVEAY003A")print "Trio2",$2,"AVEAY001A","AVEAY007A",$5,$6;else print $0;}' plink/tempFam3 > plink/tempFam4
awk '{if($1=="AVEAY003B")print "Trio2",$2,"AVEAY001A","AVEAY005B",$5,$6;else print $0;}' plink/tempFam4 > plink/tempFam5
awk '{if($1=="AVEAY004A")print "Trio2",$2,"AVEAY001A","AVEAY006B",$5,$6;else print $0;}' plink/tempFam5 > plink/tempFam6
awk '{if($1=="AVEAY001A")print "Trio2",$2,$3,$4,"1",$6;else print $0;}' plink/tempFam6 > plink/tempFam7
awk '{if($1=="AVEAY004B")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam7 > plink/tempFam8
awk '{if($1=="AVEAY006A")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam8 > plink/tempFam9
awk '{if($1=="AVEAY005A")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam9 > plink/tempFam10
awk '{if($1=="AVEAY007A")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam10 > plink/tempFam11
awk '{if($1=="AVEAY005B")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam11 > plink/tempFam12
awk '{if($1=="AVEAY006B")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam12 > plink/tempFam

awk '{if($1=="AVEAY007B")print "Trio3",$2,"AVEAY011B","AVEAY0010A",$5,$6;else print $0;}' plink/tempFam > plink/tempFam1
awk '{if($1=="AVEAY011B")print "Trio3",$2,$3,$4,"1",$6;else print $0;}' plink/tempFam1 > plink/tempFam2
awk '{if($1=="AVEAY0010A")print "Trio3",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam2 > plink/tempFam

awk '{if($1=="AVEAY008A")print "Trio4",$2,"AVEAY011A","AVEAY009B",$5,$6;else print $0;}' plink/tempFam > plink/tempFam1
awk '{if($1=="AVEAY008B")print "Trio4",$2,"AVEAY011A","AVEAY009A",$5,$6;else print $0;}' plink/tempFam1 > plink/tempFam2
awk '{if($1=="AVEAY011A")print "Trio4",$2,$3,$4,"1",$6;else print $0;}' plink/tempFam2 > plink/tempFam3
awk '{if($1=="AVEAY009B")print "Trio4",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam3 > plink/tempFam4
awk '{if($1=="AVEAY009A")print "Trio4",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam4 > plink/tempFam
cp plink/tempFam plink/pruned_allSnp_forMendel.binary.fam
plink --bfile plink/pruned_allSnp_forMendel.binary --cow --allow-no-sex --mendel --out "plink/pruned_mendel_Buri"
tail -n+2 plink/pruned_mendel_Buri.mendel | awk '{A[$1" "$2]++}END{for(i in A)print i,A[i]}'

cp plink/pruned_allSnp_forMendel.binary.fam plink/pruned_allSnp_forMendel.binary.fam_backup
#awk '{if($2=="AVEAY001A")print "AVEAY001A",$2,$3,$4,"0",$6;else print $0;}' plink/tempFam > plink/tempFam1
#awk '{if($3=="AVEAY001A"){{gsub("AVEAY001A","SRR3290615",$3)} print;}else print $0;}' plink/tempFam1 > plink/tempFam2
#awk '{if($1=="SRR3290615")print "Trio2",$2,$3,$4,"1",$6;else print $0;}' plink/tempFam2 > plink/tempFam3
#cp plink/tempFam3 plink/pruned_allSnp_forMendel.binary.fam
#plink --bfile plink/pruned_allSnp_forMendel.binary --cow --allow-no-sex --mendel --out "plink/pruned_mendel_RCI002"
#tail -n+2 plink/pruned_mendel_RCI002.mendel | awk '{A[$1" "$2]++}END{for(i in A)print i,A[i]}'
#rm plink/tempFam*


################################
## create binary inputs but after exclusion of old bad samples
plink --file plink/allSnp --keep good_samples --allow-no-sex --cow --make-bed --out plink/goodSnp.binary #14084653 variants, 32 cattle // genotyping rate 0.997505

## pruning
plink --bfile plink/goodSnp.binary --cow --maf 0.05 --geno 0.05 --indep 50 5 2 --out plink/goodSnp_filter # 218070 variants removed due to missing genotype data // 2537388 variants removed due to minor allele threshold // 11329195 variants and 32 cattle pass filters and QC. // Pruning complete.  10631741 of 11329195 variants removed
plink --bfile plink/goodSnp.binary --cow --extract plink/goodSnp_filter.prune.in --out plink/pruned_goodSnp --recode #697454 variants remaining // genotyping rate is 0.997634
plink --file plink/pruned_goodSnp --allow-no-sex --cow --make-bed --out plink/pruned_goodSnp.binary

## IBS similarity matrix but after exclusion of old bad samples
plink --bfile plink/pruned_goodSnp.binary --distance square allele-ct ibs 1-ibs --cow --allow-no-sex --out plink/distanceGood ## Excluding 24683 variants on non-autosomes
awk 'FNR==NR{a[$1]=$0;next}{ print a[$1]}' myIDs plink/distanceGood.mdist.id > plink/distanceGood.myIDs
## R (using the 1 minus the identity-by-state value distances)
Rscript -e 'args=(commandArgs(TRUE));'\
'library(ape)'\
'setwd("plink")'\
'mdist_id=read.table(args[1])'\
'mdist_table=read.table(args[2], fill=T, col.names=mdist_id$V3)'\
'mdist=as.dist(mdist_table)'\
'hc = hclust(mdist)'\
'mypal = c("red", "blue", "green", "gold", "black", "blueviolet","coral4","cadetblue", "darkred", "darkslategray", "darkorange3", "firebrick")'\
'clus = cutree(hc, 12)'\
'png(file="fan600_Good.png",width=6500,height=6500,res=600, type="cairo")'\
'plot(as.phylo(hc), type = "fan", cex = 1, font=2, use.edge.length = TRUE, tip.color = mypal[clus], label.offset = 0.002)'\
'dev.off()' distanceGood.myIDs distanceGood.myIDs


## statistics of Mendel errors in trios using 645697 pruned marker
cp plink/pruned_goodSnp.binary.bed plink/pruned_goodSnp_forMendel.binary.bed
cp plink/pruned_goodSnp.binary.bim plink/pruned_goodSnp_forMendel.binary.bim
cp plink/pruned_goodSnp.binary.fam plink/tempFam
awk '{if($1=="AVEAY012A")print "Trio1",$2,"AVEAY0010B","AVEAY013A",$5,$6;else print $0;}' plink/tempFam > plink/tempFam1
awk '{if($1=="AVEAY012B")print "Trio1",$2,"AVEAY0010B","AVEAY013B",$5,$6;else print $0;}' plink/tempFam1 > plink/tempFam2
awk '{if($1=="AVEAY14A")print "Trio1",$2,"AVEAY0010B","AVEAY14B",$5,$6;else print $0;}' plink/tempFam2 > plink/tempFam3
awk '{if($1=="AVEAY0010B")print "Trio1",$2,$3,$4,"1",$6;else print $0;}' plink/tempFam3 > plink/tempFam4
awk '{if($1=="AVEAY013A")print "Trio1",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam4 > plink/tempFam5
awk '{if($1=="AVEAY013B")print "Trio1",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam5 > plink/tempFam6
awk '{if($1=="AVEAY14B")print "Trio1",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam6 > plink/tempFam

awk '{if($1=="AVEAY001B")print "Trio2",$2,"AVEAY15B","AVEAY004B",$5,$6;else print $0;}' plink/tempFam > plink/tempFam1
awk '{if($1=="AVEAY002A")print "Trio2",$2,"AVEAY15B","AVEAY006A",$5,$6;else print $0;}' plink/tempFam1 > plink/tempFam2
awk '{if($1=="AVEAY002B")print "Trio2",$2,"AVEAY15B","AVEAY15A",$5,$6;else print $0;}' plink/tempFam2 > plink/tempFam3
awk '{if($1=="AVEAY003A")print "Trio2",$2,"AVEAY15B","AVEAY007A",$5,$6;else print $0;}' plink/tempFam3 > plink/tempFam4
awk '{if($1=="AVEAY003B")print "Trio2",$2,"AVEAY15B","AVEAY16A",$5,$6;else print $0;}' plink/tempFam4 > plink/tempFam5
awk '{if($1=="AVEAY004A")print "Trio2",$2,"AVEAY15B","AVEAY006B",$5,$6;else print $0;}' plink/tempFam5 > plink/tempFam6
awk '{if($1=="AVEAY15B")print "Trio2",$2,$3,$4,"1",$6;else print $0;}' plink/tempFam6 > plink/tempFam7
awk '{if($1=="AVEAY004B")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam7 > plink/tempFam8
awk '{if($1=="AVEAY006A")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam8 > plink/tempFam9
awk '{if($1=="AVEAY15A")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam9 > plink/tempFam10
awk '{if($1=="AVEAY007A")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam10 > plink/tempFam11
awk '{if($1=="AVEAY16A")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam11 > plink/tempFam12
awk '{if($1=="AVEAY006B")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam12 > plink/tempFam

awk '{if($1=="AVEAY007B")print "Trio3",$2,"AVEAY011B","AVEAY0010A",$5,$6;else print $0;}' plink/tempFam > plink/tempFam1
awk '{if($1=="AVEAY011B")print "Trio3",$2,$3,$4,"1",$6;else print $0;}' plink/tempFam1 > plink/tempFam2
awk '{if($1=="AVEAY0010A")print "Trio3",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam2 > plink/tempFam

awk '{if($1=="AVEAY008A")print "Trio4",$2,"AVEAY011A","AVEAY009B",$5,$6;else print $0;}' plink/tempFam > plink/tempFam1
awk '{if($1=="AVEAY008B")print "Trio4",$2,"AVEAY011A","AVEAY009A",$5,$6;else print $0;}' plink/tempFam1 > plink/tempFam2
awk '{if($1=="AVEAY011A")print "Trio4",$2,$3,$4,"1",$6;else print $0;}' plink/tempFam2 > plink/tempFam3
awk '{if($1=="AVEAY009B")print "Trio4",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam3 > plink/tempFam4
awk '{if($1=="AVEAY009A")print "Trio4",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam4 > plink/tempFam
cp plink/tempFam plink/pruned_goodSnp_forMendel.binary.fam
plink --bfile plink/pruned_goodSnp_forMendel.binary --cow --allow-no-sex --mendel --out "plink/pruned_mendel_NewSeq"
tail -n+2 plink/pruned_mendel_NewSeq.mendel | awk '{A[$1" "$2]++}END{for(i in A)print i,A[i]}'
#rm plink/tempFam*


## create binary inputs but after filteration of bad variants (no pruning) [in the good samples only]
plink --bfile plink/goodSnp.binary --maf 0.05 --geno 0.05 --allow-no-sex --cow --make-bed --out plink/goodSnp_filter.binary # 218070 variants removed due to missing genotype data // 2537388 variants removed due to minor allele threshold // 11329195 variants and 32 cattle pass filters and QC // genotyping rate is 0.997505.

## statistics of Mendel errors in trios using filtered (but not pruned) markers
cp plink/goodSnp_filter.binary.bed plink/goodSnp_filter_forMendel.binary.bed
cp plink/goodSnp_filter.binary.bim plink/goodSnp_filter_forMendel.binary.bim
cp plink/goodSnp_filter.binary.fam plink/tempFam
awk '{if($1=="AVEAY012A")print "Trio1",$2,"AVEAY0010B","AVEAY013A",$5,$6;else print $0;}' plink/tempFam > plink/tempFam1
awk '{if($1=="AVEAY012B")print "Trio1",$2,"AVEAY0010B","AVEAY013B",$5,$6;else print $0;}' plink/tempFam1 > plink/tempFam2
awk '{if($1=="AVEAY14A")print "Trio1",$2,"AVEAY0010B","AVEAY14B",$5,$6;else print $0;}' plink/tempFam2 > plink/tempFam3
awk '{if($1=="AVEAY0010B")print "Trio1",$2,$3,$4,"1",$6;else print $0;}' plink/tempFam3 > plink/tempFam4
awk '{if($1=="AVEAY013A")print "Trio1",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam4 > plink/tempFam5
awk '{if($1=="AVEAY013B")print "Trio1",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam5 > plink/tempFam6
awk '{if($1=="AVEAY14B")print "Trio1",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam6 > plink/tempFam

awk '{if($1=="AVEAY001B")print "Trio2",$2,"AVEAY15B","AVEAY004B",$5,$6;else print $0;}' plink/tempFam > plink/tempFam1
awk '{if($1=="AVEAY002A")print "Trio2",$2,"AVEAY15B","AVEAY006A",$5,$6;else print $0;}' plink/tempFam1 > plink/tempFam2
awk '{if($1=="AVEAY002B")print "Trio2",$2,"AVEAY15B","AVEAY15A",$5,$6;else print $0;}' plink/tempFam2 > plink/tempFam3
awk '{if($1=="AVEAY003A")print "Trio2",$2,"AVEAY15B","AVEAY007A",$5,$6;else print $0;}' plink/tempFam3 > plink/tempFam4
awk '{if($1=="AVEAY003B")print "Trio2",$2,"AVEAY15B","AVEAY16A",$5,$6;else print $0;}' plink/tempFam4 > plink/tempFam5
awk '{if($1=="AVEAY004A")print "Trio2",$2,"AVEAY15B","AVEAY006B",$5,$6;else print $0;}' plink/tempFam5 > plink/tempFam6
awk '{if($1=="AVEAY15B")print "Trio2",$2,$3,$4,"1",$6;else print $0;}' plink/tempFam6 > plink/tempFam7
awk '{if($1=="AVEAY004B")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam7 > plink/tempFam8
awk '{if($1=="AVEAY006A")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam8 > plink/tempFam9
awk '{if($1=="AVEAY15A")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam9 > plink/tempFam10
awk '{if($1=="AVEAY007A")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam10 > plink/tempFam11
awk '{if($1=="AVEAY16A")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam11 > plink/tempFam12
awk '{if($1=="AVEAY006B")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam12 > plink/tempFam

awk '{if($1=="AVEAY007B")print "Trio3",$2,"AVEAY011B","AVEAY0010A",$5,$6;else print $0;}' plink/tempFam > plink/tempFam1
awk '{if($1=="AVEAY011B")print "Trio3",$2,$3,$4,"1",$6;else print $0;}' plink/tempFam1 > plink/tempFam2
awk '{if($1=="AVEAY0010A")print "Trio3",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam2 > plink/tempFam

awk '{if($1=="AVEAY008A")print "Trio4",$2,"AVEAY011A","AVEAY009B",$5,$6;else print $0;}' plink/tempFam > plink/tempFam1
awk '{if($1=="AVEAY008B")print "Trio4",$2,"AVEAY011A","AVEAY009A",$5,$6;else print $0;}' plink/tempFam1 > plink/tempFam2
awk '{if($1=="AVEAY011A")print "Trio4",$2,$3,$4,"1",$6;else print $0;}' plink/tempFam2 > plink/tempFam3
awk '{if($1=="AVEAY009B")print "Trio4",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam3 > plink/tempFam4
awk '{if($1=="AVEAY009A")print "Trio4",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam4 > plink/tempFam
cp plink/tempFam plink/goodSnp_filter_forMendel.binary.fam
plink --bfile plink/goodSnp_filter_forMendel.binary --cow --allow-no-sex --mendel --out "plink/filtered_mendel_NewSeq"
tail -n+2 plink/filtered_mendel_NewSeq.mendel | awk '{A[$1" "$2]++}END{for(i in A)print i,A[i]}'
#rm plink/tempFam*


## statistics of Mendel errors in trios using all (not even filtered) markers
cp plink/goodSnp.binary.bed plink/goodSnp_forMendel.binary.bed
cp plink/goodSnp.binary.bim plink/goodSnp_forMendel.binary.bim
cp plink/goodSnp.binary.fam plink/tempFam
awk '{if($1=="AVEAY012A")print "Trio1",$2,"AVEAY0010B","AVEAY013A",$5,$6;else print $0;}' plink/tempFam > plink/tempFam1
awk '{if($1=="AVEAY012B")print "Trio1",$2,"AVEAY0010B","AVEAY013B",$5,$6;else print $0;}' plink/tempFam1 > plink/tempFam2
awk '{if($1=="AVEAY14A")print "Trio1",$2,"AVEAY0010B","AVEAY14B",$5,$6;else print $0;}' plink/tempFam2 > plink/tempFam3
awk '{if($1=="AVEAY0010B")print "Trio1",$2,$3,$4,"1",$6;else print $0;}' plink/tempFam3 > plink/tempFam4
awk '{if($1=="AVEAY013A")print "Trio1",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam4 > plink/tempFam5
awk '{if($1=="AVEAY013B")print "Trio1",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam5 > plink/tempFam6
awk '{if($1=="AVEAY14B")print "Trio1",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam6 > plink/tempFam

awk '{if($1=="AVEAY001B")print "Trio2",$2,"AVEAY15B","AVEAY004B",$5,$6;else print $0;}' plink/tempFam > plink/tempFam1
awk '{if($1=="AVEAY002A")print "Trio2",$2,"AVEAY15B","AVEAY006A",$5,$6;else print $0;}' plink/tempFam1 > plink/tempFam2
awk '{if($1=="AVEAY002B")print "Trio2",$2,"AVEAY15B","AVEAY15A",$5,$6;else print $0;}' plink/tempFam2 > plink/tempFam3
awk '{if($1=="AVEAY003A")print "Trio2",$2,"AVEAY15B","AVEAY007A",$5,$6;else print $0;}' plink/tempFam3 > plink/tempFam4
awk '{if($1=="AVEAY003B")print "Trio2",$2,"AVEAY15B","AVEAY16A",$5,$6;else print $0;}' plink/tempFam4 > plink/tempFam5
awk '{if($1=="AVEAY004A")print "Trio2",$2,"AVEAY15B","AVEAY006B",$5,$6;else print $0;}' plink/tempFam5 > plink/tempFam6
awk '{if($1=="AVEAY15B")print "Trio2",$2,$3,$4,"1",$6;else print $0;}' plink/tempFam6 > plink/tempFam7
awk '{if($1=="AVEAY004B")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam7 > plink/tempFam8
awk '{if($1=="AVEAY006A")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam8 > plink/tempFam9
awk '{if($1=="AVEAY15A")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam9 > plink/tempFam10
awk '{if($1=="AVEAY007A")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam10 > plink/tempFam11
awk '{if($1=="AVEAY16A")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam11 > plink/tempFam12
awk '{if($1=="AVEAY006B")print "Trio2",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam12 > plink/tempFam

awk '{if($1=="AVEAY007B")print "Trio3",$2,"AVEAY011B","AVEAY0010A",$5,$6;else print $0;}' plink/tempFam > plink/tempFam1
awk '{if($1=="AVEAY011B")print "Trio3",$2,$3,$4,"1",$6;else print $0;}' plink/tempFam1 > plink/tempFam2
awk '{if($1=="AVEAY0010A")print "Trio3",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam2 > plink/tempFam

awk '{if($1=="AVEAY008A")print "Trio4",$2,"AVEAY011A","AVEAY009B",$5,$6;else print $0;}' plink/tempFam > plink/tempFam1
awk '{if($1=="AVEAY008B")print "Trio4",$2,"AVEAY011A","AVEAY009A",$5,$6;else print $0;}' plink/tempFam1 > plink/tempFam2
awk '{if($1=="AVEAY011A")print "Trio4",$2,$3,$4,"1",$6;else print $0;}' plink/tempFam2 > plink/tempFam3
awk '{if($1=="AVEAY009B")print "Trio4",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam3 > plink/tempFam4
awk '{if($1=="AVEAY009A")print "Trio4",$2,$3,$4,"2",$6;else print $0;}' plink/tempFam4 > plink/tempFam
cp plink/tempFam plink/goodSnp_forMendel.binary.fam
plink --bfile plink/goodSnp_forMendel.binary --cow --allow-no-sex --mendel --out "plink/mendel_NewSeq"
tail -n+2 plink/mendel_NewSeq.mendel | awk '{A[$1" "$2]++}END{for(i in A)print i,A[i]}'
#rm plink/tempFam*


###############
## generate stats
bash summary_stats.sh
###############

## Run the kmer bait analysis
bash sandBox/targetRead_grep.sh

