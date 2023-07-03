#!/bin/bash
#SBATCH --time=2-12 # days-hours
#SBATCH --job-name=cellrngr  # Job name
#SBATCH --array=1-3
#SBATCH --nodes=1
#SBATCH --ntasks=4 # Number of cores
#SBATCH --mem=16000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # Partition to submit to
#SBATCH --reservation=rnaworkshop
#SBATCH --account=workshop
#SBATCH --output=slurmout/counts-cellrngr-%A-%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/counts-cellrngr-%A-%a.err # File to which STDERR will be written
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=youremail@whatever.edu # Email to which notifications will be sent

## Record the start time
start=`date +%s`

## Record the host being run on
echo "Hostname: $(eval hostname)"


THREADS=${SLURM_NTASKS}
MEM=$(expr ${SLURM_MEM_PER_NODE} / 1024)

echo "Allocated threads: " $THREADS
echo "Allocated memory: " $MEM

## Where cellranger executable is located
## a) by loading a module
#module load cellranger/7.0.0

## b) or, by placing the location of the executables on the path (edit to your location)
export PATH=/share/workshop/scRNA_workshop/Software/cellranger-7.1.0/bin:$PATH

## c) or if they are already on the path, do nothing

## Set the parameters for the run
basedir="/share/workshop/scRNA_workshop/$USER/scrnaseq_example"
transcriptome=/share/workshop/scRNA_workshop/Software/refdata-gex-GRCh38-2020-A
fastqs="${basedir}/00-RawData"
outdir="${basedir}/01-Cellranger"

## provide the script the row # of the sample to be run
sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" ${basedir}/samples.txt`

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
## --expect-cells=NUM, override auto-estimate of expected number of recovered cells
## --force-cells=NUM, force pipeline to use this number of cells, bypassing cell detection algorithm
## --include-introns=false         exclude intronic reads in count (new in cellranger v7.0)
## --nosecondary, skip the unnecessary secondary analysis
## --r2-length=NUM, if your R2 qualities are really poor
## --chemistry=CHEM, should it fail chemistry detection

## Create output directory and change directory
mkdir -p $outdir
cd $outdir
## Echo the call
echo $call
## Evaluate the call
eval $call

## Record the start time, and output runtime
end=`date +%s`
runtime=$((end-start))
echo $runtime
