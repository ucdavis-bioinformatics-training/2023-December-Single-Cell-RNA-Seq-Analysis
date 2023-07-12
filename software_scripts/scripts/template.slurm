#!/bin/bash
# options for sbatch
#SBATCH --job-name=name # Job name
#SBATCH --nodes=1 # should never be anything other than 1
#SBATCH --ntasks=1 # number of cpus to use
#SBATCH --time=10 # Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds".
#SBATCH --mem=500 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # cluster partition
#SBATCH --account=workshop # cluster account to use for the job
#SBATCH --reservation=workshop # cluster account reservation
#SBATCH --array=1 # Task array indexing, see https://slurm.schedmd.com/job_array.html
#SBATCH --output=template-%A-%a.out # File to which STDOUT will be written
#SBATCH --error=template-%A-%a.err # File to which STDERR will be written
#SBATCH --mail-type=END # when to email user: END, FAIL, ALL, and others. See https://slurm.schedmd.com/sbatch.html
#SBATCH --mail-user=myemail@email.com # edit to receive email from the scheduler

# for calculating the amount of time the job takes and echo the hostname
begin=`date +%s`
echo $HOSTNAME

# Sleep for 5 minutes
sleep 300

# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed
