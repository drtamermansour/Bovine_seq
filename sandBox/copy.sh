current_dir="Tamer/Bovine_seq"
new_dir="Tamer2/Bovine_seq" ## "Tamer/Bovine_seq"
mkdir -p ../../$new_dir/{refGenome,data,qc,vc}
cp * ../../$new_dir/.
cp -R hpcc* ../../$new_dir/.
cp -R .git ../../$new_dir/.
cp refGenome/*.fa ../../$new_dir/refGenome/.
cp -R refGenome/*Index ../../$new_dir/refGenome/.
cp -R data/fastq ../../$new_dir/data/.
cp -R data/trimmed ../../$new_dir/data/.
cp -R qc/fastqc ../../$new_dir/qc/.
cp -R qc/multiqc ../../$new_dir/qc/.
cp -R data/mapped_reads ../../$new_dir/data/.
cp -R data/sorted_reads ../../$new_dir/data/.
cp -R data/merged_reps ../../$new_dir/data/.  ## replaced by the "select merged_reps" module
cp -R qc/mappingQC ../../$new_dir/qc/.        ## replaced by the "select mappingQC" module
cp -R data/dedup ../../$new_dir/data/.        ## replaced by the "select dedup" module
cp -R vc/hapCaller_single ../../$new_dir/vc/.        ## replaced by the "select gvcf" module
cp -R vc/combineGVCF ../../$new_dir/vc/.        ## replaced by the "select combined gvcf" module
cp -R reports ../../$new_dir/.
cp -R logs* ../../$new_dir/.

cp -R targetRead* ../../$new_dir/.
cp -R singleFiles ../../$new_dir/.

cp -R .snakemake ../../$new_dir/.

#####################
ext="external/SRP072240"
mkdir -p ../../$new_dir/$ext/{refGenome,data,qc,vc}
cp $ext/* ../../$new_dir/$ext/.
cp -R $ext/hpcc* ../../$new_dir/$ext/.
#cp -R $ext/.git ../../$new_dir/$ext/.
cp $ext/refGenome/*.fa ../../$new_dir/$ext/refGenome/.
cp -R $ext/refGenome/*Index ../../$new_dir/$ext/refGenome/.
cp -R $ext/data/fastq ../../$new_dir/$ext/data/.
cp -R $ext/data/trimmed ../../$new_dir/$ext/data/.
cp -R $ext/qc/fastqc ../../$new_dir/$ext/qc/.
cp -R $ext/qc/multiqc ../../$new_dir/$ext/qc/.
cp -R $ext/data/mapped_reads ../../$new_dir/$ext/data/.
cp -R $ext/data/sorted_reads ../../$new_dir/$ext/data/.
#cp -R $ext/data/merged_reps ../../$new_dir/$ext/data/. 
for singleSample in data/sorted_reads/SRP072240/*.bam;do 
  newSample=$(basename $singleSample); 
  ln -s $(pwd)/$singleSample data/merged_reps/$newSample;done
cp -R $ext/qc/mappingQC ../../$new_dir/$ext/qc/.
cp -R $ext/data/dedup ../../$new_dir/$ext/data/.
cp -R $ext/vc/hapCaller_single ../../$new_dir/$ext/vc/. 
cp -R $ext/vc/combineGVCF ../../$new_dir/$ext/vc/.
cp -R $ext/logs* ../../$new_dir/$ext/.

#cp -R $ext/.snakemake ../../$new_dir/$ext/.

####################
## select merged_reps
mkdir -p ../../$new_dir/data/merged_reps
cp data/merged_reps/AVEAY001A.bam ../../$new_dir/data/merged_reps/.
cd ../../$new_dir/
for sample in AVEAY011A AVEAY011B AVEAY012A AVEAY012B AVEAY013A AVEAY013B;do
  ln -s $(pwd)/data/sorted_reads/5_2_18/${sample}_*.bam data/merged_reps/${sample}.bam
done
cd ../../$current_dir

## select mappingQC
mkdir -p ../../$new_dir/qc/mappingQC
for sample in AVEAY011A AVEAY011B AVEAY012A AVEAY012B AVEAY013A AVEAY013B AVEAY001A;do
  cp qc/mappingQC/$sample.cov ../../$new_dir/qc/mappingQC/.
  cp qc/mappingQC/$sample.stat ../../$new_dir/qc/mappingQC/.
done

## select dedup
mkdir -p ../../$new_dir/data/dedup
for sample in AVEAY011A AVEAY011B AVEAY012A AVEAY012B AVEAY013A AVEAY013B AVEAY001A;do
  cp data/dedup/$sample.bam ../../$new_dir/data/dedup/.
  cp data/dedup/$sample.bai ../../$new_dir/data/dedup/.
  cp data/dedup/$sample.txt ../../$new_dir/data/dedup/.
done

## select gvcf
mkdir -p ../../$new_dir/vc/hapCaller_single
for sample in AVEAY011A AVEAY011B AVEAY012A AVEAY012B AVEAY013A AVEAY013B AVEAY001A;do
  cp vc/hapCaller_single/$sample.g.vcf ../../$new_dir/vc/hapCaller_single/.
  cp vc/hapCaller_single/$sample.g.vcf.idx ../../$new_dir/vc/hapCaller_single/.
done

## select combined gvcf
mkdir -p ../../$new_dir/vc/combineGVCF
cp vc/combineGVCF/combined_patch_[67].vcf ../../$new_dir/vc/combineGVCF/.
cp vc/combineGVCF/combined_patch_[67].vcf.idx ../../$new_dir/vc/combineGVCF/.

