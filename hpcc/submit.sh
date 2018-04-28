
logdir=hpcc/log
mkdir -p $logdir 

QSUB="qsub -l nodes=1:ppn={threads} "
QSUB="$QSUB -o $logdir -e $logdir -A ged" 

snakemake                               \
    -j 100                              \
    --profile ./hpcc                    \
    --rerun-incomplete                  \
    --keep-going                        \
    --latency-wait 10                   \
    --max-jobs-per-second 1             \
    --use-conda                         \
    --cluster "$QSUB" $@
