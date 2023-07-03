#!/bin/bash

## Record the start time
start=`date +%s`

## Record the host being run on
echo "Hostname: $(eval hostname)"

THREADS=4
MEM=10

echo "Allocated threads: " $THREADS
echo "Allocated memory: " $MEM

## Where cellranger executable is located
export PATH=/share/workshop/scRNA_workshop/Software/cellranger-7.1.0/bin:$PATH
## Set the parameters for the run
basedir="/share/workshop/scRNA_workshop/${USER}"
transcriptome="/share/workshop/scRNA_workshop/Software/refdata-gex-GRCh38-2020-A"
fastqs="${basedir}/scrnaseq_example/00-RawData"
outdir="${basedir}/01-Cellranger"
mkdir -p $outdir
cd $outdir
## loop over samples in sample sheet, running cellranger on each
for sample in `cat ${basedir}/samples.txt`
do
  ## https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome
  ## Create the call
  call="cellranger count \
    --id=${sample} \
    --sample=${sample} \
    --transcriptome=${transcriptome} \
    --fastqs=${fastqs} \
    --nosecondary \
    --localcores=${THREADS} \
    --localmem=${MEM}"

  ## Some other parameters that may be usefull/needed
  ## --expect-cells=NUM, number of cells expected
  ## --include-introns         Include intronic reads in count
  ## --nosecondary, skip the unnecessary secondary analysis
  ## --r2-length=NUM, if your R2 qualities are really poor
  ## --chemistry=CHEM, should it fail chemistry detection

  ## Echo the call
  echo $call
  ## Evaluate the call
  #eval $call
done

## Record the start time, and output runtime
end=`date +%s`
runtime=$((end-start))
echo $runtime
