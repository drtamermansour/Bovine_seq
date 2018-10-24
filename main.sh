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

## ftp download
wget --no-verbose --no-parent --recursive --level=1 --no-directories --user=gslftp --password=gsl23ftp! ftp://gslserver.qb3.berkeley.edu/181004_150PE_HS4K2A/L235678/Young/AVEAY*

cd Bovine_seq/data/fastq/10_9_18/
for f in *.gz;do
newf=$(echo $f | sed 's/_/0/')
echo $newf;
mv $f $newf;
done
	

## download data from SRP072240 (Carlson et al. 2016 paper)
mkdir -p Bovine_seq2/data/fastq/SRP072240
cd Bovine_seq2/data/fastq/SRP072240 ## Whole genome sequencing of two edited calves along with those of their progenitor cell lines, 2122 and 2120
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

## Install Snakemake using conda & deploy PBS or slurm profile
bash snakemake_setup.sh

## Run Snakemake
# 1. add conda and the hpcc folder to you path
export PATH=$HOME/miniconda3/bin:$(pwd)/hpcc:$PATH
# 2. turn on the environment
source activate snakemake
# 3. run by the submoit script
. hpcc/submit.sh


## Run Snakemake on slurm
# 1. add conda and the hpcc folder to you path
export PATH=$HOME/miniconda3/bin:$(pwd)/slurm:$PATH
# 2. turn on the environment
source activate snakemake
# 3. run by the submoit script
. slurm/submit.sh


## continue after Snakemake analysis
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
  vcf=vc/hapCaller_raw_withExternal.pass_${var}s.vcf
  for i in {10..41};do
    id=$(grep "^#CHROM" $vcf | awk -F'\t' -v i=$i '{print $i}');
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
module load vcftools/0.1.14
vcf-concat vc/hapCaller_raw_withExternal.pass_SNPs.vcf vc/hapCaller_raw_withExternal.pass_INDELs.monoAllel_edit.vcf | vcf-sort > vc/allSnp.vcf
## genome annotation map
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/263/795/GCA_002263795.2_ARS-UCD1.2/GCA_002263795.2_ARS-UCD1.2_assembly_report.txt -O refGenome/assembly_report.txt
grep -v "^#" refGenome/assembly_report.txt | awk '{if($2=="assembled-molecule")print $5,$1}' > refGenome/assembly_map
module load bcftools/1.2
bcftools annotate --rename-chrs refGenome/assembly_map vc/allSnp.vcf > vc/allSnp_chr.vcf
## prepare Plink input
module load vcftools/0.1.14
mkdir -p plink
vcftools --vcf vc/allSnp_chr.vcf --plink --out plink/allSnp 
## create binary inputs
module load plink/1.9
plink --file plink/allSnp --allow-no-sex --cow --make-bed --out plink/allSnp.binary
## create file of alternative alleles 
#cat allSnp.vcf | awk 'BEGIN{FS="\t";OFS="\t";}/#/{next;}{{if($3==".")$3=$1":"$2;}print $3,$5;}'  > alt_alleles
## pruning
plink --bfile plink/allSnp.binary --cow --maf 0.05 --geno 0.05 --indep 50 5 2 --out plink/allSnp_filter
plink --bfile plink/allSnp.binary --cow --extract plink/allSnp_filter.prune.in --out plink/pruned_allSnp --recode
plink --file plink/pruned_allSnp --allow-no-sex --cow --make-bed --out plink/pruned_allSnp.binary

## IBS similarity matrix (typically used for population stratification): distances are expressed in "shared allele counts", "proportions of alleles IBS", "1 minus the identity-by-state value" 
plink --bfile plink/pruned_allSnp.binary --distance square allele-ct ibs 1-ibs --cow --allow-no-sex --out plink/distance
## R (using the 1 minus the identity-by-state value distances)
library(ape)
setwd("plink")
mdist_id=read.table("distance.myIDs")
mdist_table=read.table("distance.mdist", fill=T, col.names=mdist_id$V4)
mdist=as.dist(mdist_table)
hc = hclust(mdist) ## setwd("~/Desktop/FileZilla_client/breedSp/clustering_brachy"); save(hc, file="hclust.rda"); load("hclust.rda"); 
# vector of colors
#mypal = c("#556270", "#4ECDC4", "#1B676B", "#FF6B6B", "#C44D58")
mypal = c("red", "blue", "green", "gold", "black", "blueviolet","coral4","cadetblue", "darkred", "darkslategray", "darkorange3", "firebrick")
# cutting dendrogram in 5 clusters
clus = cutree(hc, 12)
# plot
png(file="fan600.png",width=6000,height=6000,res=600)
plot(as.phylo(hc), type = "fan", cex = 1, font=2, use.edge.length = TRUE, tip.color = mypal[clus], label.offset = 0.002)
dev.off()

## IBS similarity matrix using all variants
plink --bfile plink/allSnp.binary --distance square allele-ct ibs 1-ibs --cow --allow-no-sex --out plink/distanceAll
## R (using the 1 minus the identity-by-state value distances)
library(ape)
setwd("plink")
mdist_id=read.table("distance.myIDs")
mdist_table=read.table("distanceAll.mdist", fill=T, col.names=mdist_id$V4)
mdist=as.dist(mdist_table)
hc = hclust(mdist) 
# vector of colors
#mypal = c("#556270", "#4ECDC4", "#1B676B", "#FF6B6B", "#C44D58")
mypal = c("red", "blue", "green", "gold", "black", "blueviolet","coral4","cadetblue", "darkred", "darkslategray", "darkorange3", "firebrick")
# cutting dendrogram in 5 clusters
clus = cutree(hc, 12)
# plot
png(file="fan600_All.png",width=6000,height=6000,res=600)
plot(as.phylo(hc), type = "fan", cex = 1, font=2, use.edge.length = TRUE, tip.color = mypal[clus], label.offset = 0.002)
dev.off()

## statistics of Mendel errors in trios using 645697 pruned marker
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
awk '{if($2=="AVEAY001A")print "AVEAY001A",$2,$3,$4,"0",$6;else print $0;}' plink/tempFam > plink/tempFam1
awk '{if($3=="AVEAY001A"){{gsub("AVEAY001A","SRR3290615",$3)} print;}else print $0;}' plink/tempFam1 > plink/tempFam2
awk '{if($1=="SRR3290615")print "Trio2",$2,$3,$4,"1",$6;else print $0;}' plink/tempFam2 > plink/tempFam3
cp plink/tempFam3 plink/pruned_allSnp_forMendel.binary.fam
plink --bfile plink/pruned_allSnp_forMendel.binary --cow --allow-no-sex --mendel --out "plink/pruned_mendel_RCI002"
tail -n+2 plink/pruned_mendel_RCI002.mendel | awk '{A[$1" "$2]++}END{for(i in A)print i,A[i]}'
rm plink/tempFam*
#sed -i 's/^BD411/Trio1/' pruned_allSnp_KeepTrios_forMendel.binary.fam







## try to do the target assembly analysis
bash targetRead_Mapsembler.sh

bash snakemake_cats_setup.sh 
. hpcc_cats/submit.sh

bash targetRead_grep.sh


## generate stats
bash summary_stats.sh

