module load picard/2.18.1-Java-1.8.0_152
mkdir -p data/temp_sorted_reads/
for r1 in data/trimmed/10_*/*_R1_001.fastq.gz;do #echo $r1;#done
name=$(basename $r1)
SM=$(echo $name | cut -d "_" -f1)
LB=$(echo $name | cut -d"_" -f1,2)  
batch=$(basename "$(dirname $r1)")
if [ "$batch" != "trimmed" ];then LB=$batch.$LB;fi
PL="Illumina"
header=$(head -n1 <(zcat $r1) | grep ':*:*:*:*:*:*')
if [ "$header" != "" ]; then
RGID=$(echo "$header" | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)
else # "make unique ID and PU using checksum"
checksum=$(shasum {input.r1} | awk '{ print $1 }')
RGID="UnChrPU_"$checksum
fi
PU=$RGID.$LB
echo RGID $RGID LB $LB PL $PL PU $PU SM $SM
dir=$(echo $r1 | cut -d"/" -f3)
id=$(basename $r1 _R1_001.fastq.gz)
input=data/sorted_reads/$dir/$id.bam
output=data/temp_sorted_reads/$dir/$id.bam
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=$input O=$output RGID=$RGID RGLB=$LB RGPL=$PL RGPU=$PU RGSM=$SM
done

