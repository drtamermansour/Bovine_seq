mkdir -p BovineSeq_LR/data 
work_dir=$(pwd)/BovineSeq_LR
# data download
cd BovineSeq_LR/data
wget -r --level=10 -nH -nc --cut-dirs=3 --no-parent --reject "wget_index.html" --no-check-certificate --header "Cookie: sessionid=ln5lvp001udo9vn9nn5wuhaj6v1n9wmj;" https://bioshare.bioinformatics.ucdavis.edu/bioshare/wget/2doso2j3dbwnzbh/wget_index.html
# uncompress the fastq data
tar xvzf fastq.tar.gz

#Install [assembly-stats](https://github.com/sanger-pathogens/assembly-stats):
source activate workEnv1
git clone https://github.com/sanger-pathogens/assembly-stats.git
cd assembly-stats/
mkdir build
cd build
cmake -DINSTALL_DIR:PATH=$HOME/bin/ ..
make
make test
make install  #sudo make install

cd $work_dir
assembly-stats <(cat data/fastq_pass/*.fastq)
#sum = 29927808118, n = 2427368, ave = 12329.32, largest = 197990
#N50 = 24635, n = 416831
#N60 = 20797, n = 548924
#N70 = 16868, n = 708242
#N80 = 12490, n = 912722
#N90 = 7034, n = 1224953
#N100 = 70, n = 2427368
#N_count = 0
#Gaps = 0

assembly-stats <(cat data/fastq_fail/*.fastq)
#sum = 7266887390, n = 1875249, ave = 3875.16, largest = 160817
#N50 = 15533, n = 141728
#N60 = 11869, n = 195146
#N70 = 8317, n = 267997
#N80 = 5099, n = 379186
#N90 = 2429, n = 584087
#N100 = 5, n = 1875249
#N_count = 0
#Gaps = 0


bash targetRead_grep_plasmid.sh  ## find reads aligned to the empty TOPO plasmid
assembly-stats  targetReads_grep_plasmid/mergedFiles/fastq_*
#stats for targetReads_grep_plasmid/mergedFiles/fastq_fail.fa
#sum = 454934, n = 26, ave = 17497.46, largest = 55562
#N50 = 27872, n = 6
#N60 = 25087, n = 8
#N70 = 19364, n = 10
#N80 = 14206, n = 12
#N90 = 9004, n = 16
#N100 = 1586, n = 26
#N_count = 0
#Gaps = 0
#stats for targetReads_grep_plasmid/mergedFiles/fastq_pass.fa
#sum = 1567576, n = 55, ave = 28501.38, largest = 60727
#N50 = 40779, n = 17
#N60 = 35029, n = 21
#N70 = 30810, n = 26
#N80 = 26975, n = 32
#N90 = 17861, n = 39
#N100 = 521, n = 55
#N_count = 0
#Gaps = 0


conda create --name kProcessV1 python=3.7
source activate kProcessV1
mkdir ~/kPro_v1 && cd ~/kPro_v1
wget https://github.com/drtamermansour/nu-ngs02/raw/master/Day5/kProcessor/dist/kProcessor-0.1-cp37-cp37m-linux_x86_64.whl
pip install kProcessor-0.1-cp37-cp37m-linux_x86_64.whl

cd $work_dir
cat targetReads_grep_plasmid/bait_unwrapped.fa | sed 's/>/@/' > targetReads_grep_plasmid/bait_unwrapped.fq
cat targetReads_grep_plasmid/bait_unwrapped.fa | sed 's/>/+/' | tr  "atcgATCG" "........" >> targetReads_grep_plasmid/bait_unwrapped.fq
mkdir -p targetReads_grep_plasmid/{plasmid_match_fq,plasmid_match_fa}
IFS=$'\t'; cat targetReads_grep_plasmid/singleFiles/fastq_*/*.fq | paste - - - - | \
 while read header seq empty qual;do readID=$(echo $header | awk -F"[ =]" '{print $5}');
  echo -e $header"\n"$seq"\n"$empty"\n"$qual > targetReads_grep_plasmid/plasmid_match_fq/$readID.fastq;
  newHeader=$(echo $header | sed 's/^@/>/');
  echo -e $newHeader"\n"$seq > targetReads_grep_plasmid/plasmid_match_fa/$readID.fasta;
done
for f in targetReads_grep_plasmid/plasmid_match_fq/*.fastq;do if [ -s $f ];then echo $f;python shared_with_plasmid.py -i $f;fi;done | paste - - > targetReads_grep_plasmid/plasmid_matchingKmers
cat targetReads_grep_plasmid/plasmid_matchingKmers | awk '{if($2>10)print $0}' > targetReads_grep_plasmid/plasmid_matchingKmers.sig
cat targetReads_grep_plasmid/plasmid_matchingKmers | awk '{if($2>10)print $1}' | xargs cat > targetReads_grep_plasmid/plasmid_matchingKmers.fastq

assembly-stats  targetReads_grep_plasmid/plasmid_matchingKmers.fastq
#sum = 175631, n = 8, ave = 21953.88, largest = 47455
#N50 = 42151, n = 2
#N60 = 27507, n = 3
#N70 = 18680, n = 4
#N80 = 17690, n = 5
#N90 = 14522, n = 6
#N100 = 521, n = 8

## Test mapping
source activate workEnv1
conda install -c bioconda minimap2
minimap2 -a targetReads_grep_plasmid/plasmid_match_fa/13566.fasta targetReads_grep_plasmid/bait_unwrapped.fa > 13566_plasmid.sam


##  Create a new version of the genome with the suspected insertion (the plasmid with 2 clones)
# 1. get a copy of plasmid.fa & clone_From_Carlson2016.fa
# 2. get a copy (or short cut) for the reference genome 
ln -s $HOME/Tamer/Bovine_seq/refGenome/ARS-UCD1.2_chr.fa .
# 3. Test extraction 
module load  BEDTools/2.27.1
echo -e "CM008168.2\t2429108\t2429118" > test.bed  #start of the 1st 10pb of the 1st 212bp Polled (591bp) CM008168.2:2429109-2429118
echo -e "CM008168.2\t2429329\t2429335" >> test.bed # the 6pb of deletion CM008168.2:2429330-2429335  (should be replaced by 208)
echo -e "CM008168.2\t2428499\t2429891" >> test.bed # the reference seq replaced by Carlson 1594bp clone CM008168.2:2428500-2429891 (start-830 TO end+556)
bedtools getfasta -fi ARS-UCD1.2_chr.fa -bed test.bed
# 4. create the suggestedInsert.fa (clone.forward__TOPO.seq:1-293.reverse__TOPO.seq:294-3931.reverse__clone.forward) 1594+3931+1594=7119
echo -e "pCR2.1_TOPO\t0\t293" > insertedPlasmid.bed # TOPO.seq:1-293
echo -e "pCR2.1_TOPO\t293\t3931" >> insertedPlasmid.bed # TOPO.seq:294-3931
bedtools getfasta -fi plasmid.fa -bed insertedPlasmid.bed | paste - - | awk '{print $2}' | tr "atcg" "ATCG" | rev | tr "ATGC" "TACG" | awk '{ printf("%s",$0);}' > insertedPlasmid.seq
cat clone_From_Carlson2016.fa  | awk '/^>/ {next; } { printf("%s",$0); }' > clone_unwrapped.seq
echo ">suggestedInsert clone.forward__TOPO.seq:1-293.reverse__TOPO.seq:294-3931.reverse__clone.forward" > suggestedInsert.fa
cat clone_unwrapped.seq insertedPlasmid.seq clone_unwrapped.seq | awk '{ printf("%s",$0);}' >> suggestedInsert.fa
# Test mapping
minimap2 -a targetReads_grep_plasmid/plasmid_match_fa/13566.fasta suggestedInsert.fa > 13566_suggestedInsert.sam
# 5. Edit the genome
Bef="$(bedtools getfasta -fi ARS-UCD1.2_chr.fa -bed <(echo -e "CM008168.2\t2428499\t2429891") -tab | awk '{print $2}')"
#Aft="$(cat clone_unwrapped.seq)"
Aft="$(cat suggestedInsert.fa | tail -n1)"
conda install -c bioconda seqkit
cat ARS-UCD1.2_chr.fa | seqkit replace -p $Bef -r $Aft -s -i > ARS-UCD1.2_chr_suggestedInsert.fa


## Alignment
minimap2 -a ARS-UCD1.2_chr_suggestedInsert.fa targetReads_grep_plasmid/plasmid_matchingKmers.fastq > matchReads_sugRef.sam
cat matchReads_sugRef.sam | grep -v "^@" | awk '{print $3, $4, $5, length($10)}'

module load SAMtools
samtools view -bo matchReads_sugRef.bam matchReads_sugRef.sam 
samtools sort matchReads_sugRef.bam -o matchReads_sugRef.sorted.bam
samtools index matchReads_sugRef.sorted.bam

## The expected co-oridinats of the new features
## The started of the suggested insert should be CM008168.2:2428500
bedtools getfasta -fi ARS-UCD1.2_chr_suggestedInsert.fa -bed <(echo -e "CM008168.2\t2428499\t2428509") ## 10bp
## the end of the suggested insert should be CM008168.2:2435618   (2428500+7119-1=2435618)
bedtools getfasta -fi ARS-UCD1.2_chr_suggestedInsert.fa -bed <(echo -e "CM008168.2\t2435608\t2435618") ## 10bp
## The position of the 208 bp replacing the 6bp in the 1st clone CM008168.2:2429330-2429537
bedtools getfasta -fi ARS-UCD1.2_chr_suggestedInsert.fa -bed <(echo -e "CM008168.2\t2429329\t2429537") ## 208bp ## The start of the change 
## The plasmid (start: 2428500+1594=2430094) & (End: 2430094+3931-1=2434024)
bedtools getfasta -fi ARS-UCD1.2_chr_suggestedInsert.fa -bed <(echo -e "CM008168.2\t2430093\t2430103") ## 10bp
bedtools getfasta -fi ARS-UCD1.2_chr_suggestedInsert.fa -bed <(echo -e "CM008168.2\t2434014\t2434024") ## 10bp
## The position of the 208 bp replacing the 6bp in the 2st clone (start 2434024+1+(2429330-2428500)) CM008168.2:2434855-2435062
bedtools getfasta -fi ARS-UCD1.2_chr_suggestedInsert.fa -bed <(echo -e "CM008168.2\t2434854\t2435062") ## 208bp ## The end  of the change
## The position of the actual suggested mutation CM008168.2:2429330-2435062 (which should have replaced the 6 bps) 

###############################################
bash targetRead_grep_clone.sh

source activate kProcessV1
cd $work_dir
cat targetReads_grep_clone/bait_unwrapped.fa | sed 's/>/@/' > targetReads_grep_clone/bait_unwrapped.fq
cat targetReads_grep_clone/bait_unwrapped.fa | sed 's/>/+/' | tr  "atcgATCG" "........" >> targetReads_grep_clone/bait_unwrapped.fq

> targetReads_grep_clone/log
> targetReads_grep_clone/clone_matchingKmers.fastq
let i=0
IFS=$'\t'; cat targetReads_grep_clone/mergedFiles/fastq_*.fq | paste - - - - | \
 while read header seq empty qual;do 
  let i++
  x=$(python shared_with_reference.py -k 31 -i <(echo -e $header"\n"$seq"\n"$empty"\n"$qual) -r targetReads_grep_clone/bait_unwrapped.fq)
  echo $i $header $x
  if [ "$x" -gt 10 ];then
    echo -e $header"\n"$seq"\n"$empty"\n"$qual >> targetReads_grep_clone/clone_matchingKmers.fastq;fi
done &> targetReads_grep_clone/log


##  Create a new version of the genome with the clone inserted. The co-ordinated of this insertion in this new version of the genome should be CM008168.2:2429330-2429537).  
source activate workEnv1
module load  BEDTools/2.27.1
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < clone_From_Carlson2016.fa | tail -n +2 > clone_unwrapped.fa
Bef="$(bedtools getfasta -fi ARS-UCD1.2_chr.fa -bed <(echo -e "CM008168.2\t2428499\t2429891") -tab | awk '{print $2}')"
Aft="$(cat clone_unwrapped.fa | tail -n1)"
cat ARS-UCD1.2_chr.fa | seqkit replace -p $Bef -r $Aft -s -i > ARS-UCD1.2_chr_cloneInsert.fa


## Alignment
minimap2 -a ARS-UCD1.2_chr_cloneInsert.fa targetReads_grep_clone/clone_matchingKmers.fastq > matchReads_insRef.sam
cat matchReads_insRef.sam | grep -v "^@" | awk '{print $3, $4, $5, length($10)}' | sort -k1,1 -k2,2n > matchReads_insRef.sam.summary

#module load SAMtools
#samtools view -bo matchReads_insRef.bam matchReads_insRef.sam
#samtools sort matchReads_insRef.bam -o matchReads_insRef.sorted.bam
#samtools index matchReads_insRef.sorted.bam

# reads previously aligned to the plasmid
cat matchReads_sugRef.sam | grep -v "^@" | awk '{print $1}' | grep -Fwf - matchReads_insRef.sam > plasmid_reads.sam
cat matchReads_sugRef.sam | grep -v "^@" | awk '{print $1,$3, $4, $5, length($10)}'
cat plasmid_reads.sam | grep -v "^@" | awk '{print $1, $3, $4, $5, length($10)}' ## most of the reads showed supplementary alignment 
# reads other than those previously aligned to the plasmid
cat matchReads_sugRef.sam | grep -v "^@" | awk '{print $1}' | grep -vFwf - matchReads_insRef.sam > matchReads_insRef_noplasmid.sam
cat matchReads_insRef_noplasmid.sam | grep -v "^@" | awk '{print $1, $3, $4, $5, length($10)}' | grep CM008168.2 | awk '$3<2429330' | awk '$3+$5>2429537'
module load SAMtools
samtools view -bo matchReads_insRef_noplasmid.bam matchReads_insRef_noplasmid.sam
samtools sort matchReads_insRef_noplasmid.bam -o matchReads_insRef_noplasmid.sorted.bam
samtools index matchReads_insRef_noplasmid.sorted.bam


#minimap2 -a ARS-UCD1.2_chr_suggestedInsert.fa targetReads_grep_clone/clone_matchingKmers.fastq > matchReads_sugRef2.sam
#cat matchReads_sugRef2.sam | grep -v "^@" | awk '{print $1,$3, $4, $5, length($10)}'

#samtools view -bo matchReads_sugRef2.bam matchReads_sugRef2.sam
#samtools sort matchReads_sugRef2.bam -o matchReads_sugRef2.sorted.bam
#samtools index matchReads_sugRef2.sorted.bam

